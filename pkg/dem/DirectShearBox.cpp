/*************************************************************************
*  Copyright (C) 2026 by Chang GUO                                       *
*  chang.guo@connect.polyu.hk                                            *
*                                                                        *
*  DirectShearBox / DirectShearBox6D: direct-shear boundary controller,  *
*  added in V2.0.0.                                                      *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "DirectShearBox.hpp"
#include <lib/high-precision/Constants.hpp>
#include <core/Clump.hpp>
#include <core/Scene.hpp>
#include <core/State.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/CohesiveFrictionalContactLaw.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <assert.h>
#include <preprocessing/dem/Shop.hpp>
#include <core/Omega.hpp>
#include <pkg/dem/DirectShearBoxRecorder.hpp>
#include <pkg/dem/NewtonIntegrator.hpp>
#include <preprocessing/dem/SpherePack.hpp>
#include <unordered_set>

namespace yade { // Cannot have #include directive inside.

CREATE_CPP_LOCAL_LOGGER("DirectShearBox.cpp");
YADE_PLUGIN((DirectShearBox)(DirectShearBox6D));

Vector3r DirectShearBox::getStress(int boundId) const 
{
	/*
	return stress. index: wallsX0 = 0, wallsX1=1, wallsY0=2, wallsY1=3, wallsZ0=4, wallsZ1=5. stress[0][0], [1][0], [2][1], [3][1], [4][2] and [5][2] may be normally needed and always positive.
	*/
	assert(boundId >= 0 && boundId <= 5);
	return stress[boundId];
}

void DirectShearBox::updateStiffness() 
{
	/*
	Sum the normal stiffness kn of all load-bearing (fn != 0) contacts touching the pressor
	(wall_id_Z1). Result is stored in 'stiffness' and used by the stress servo control.
	Called every stiffnessUpdateInterval iterations (modified from the triaxial controller).
	*/
	stiffness = 0;
	Real local_stiffness = 0;
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided)
#endif
	for (long idx = 0; idx < (long)(scene->interactions->size()); ++idx) {
		const shared_ptr<Interaction>& contact = (*scene->interactions)[idx];
		if (!contact || !contact->isReal()) continue;
		
		Real fn = (static_cast<FrictPhys*>(contact->phys.get()))->normalForce.norm();
		if (fn == 0) continue;
		
		int id1 = contact->getId1(), id2 = contact->getId2();
		if (wall_id_Z1 == id1 || wall_id_Z1 == id2) {
			FrictPhys* currentContactPhysics = static_cast<FrictPhys*>(contact->phys.get());
#ifdef YADE_OPENMP
			#pragma omp critical
#endif
			local_stiffness += currentContactPhysics->kn;
		}
	}
	stiffness = local_stiffness;
}

void DirectShearBox::controlExternalStress(Real resultantForceZ) // controls walls such that Sum Forces from Sample on Wall = resultantForce
{
	/*
	Servo-control the pressor wall (Z1) so that the total normal force exerted by the
	sample on the wall equals resultantForceZ (= pressorGoal * wall area). Uses the
	steady-flow law with stressDamping / steadyflowFactor and clamps the resulting
	velocity to max_vel. Key routine for external compaction (modified from triaxial controller).
	*/
	// Ensure cached pointers are valid
	if (!cachedWallZ1 || !cachedStateZ1) {
		cachedWallZ1 = Body::byId(wall_id_Z1, scene);
		cachedStateZ1 = cachedWallZ1->state.get();
	}
	
	Real translation = getForce(scene, wall_id_Z1).z() - resultantForceZ;
	
	if (translation != 0) {
		if (stiffness!= 0) 	translation /= stiffness;
		else translation = max_vel * math::sign(translation) * scene->dt;
	}
	previousTranslation = (1 - stressDamping) * translation + steadyflowFactor * previousTranslation; // formula for "steady-flow" evolution with fluctuations	
	previousTranslation	 = math::min(math::abs(previousTranslation), max_vel * scene->dt) * math::sign(previousTranslation);
	externalWork += previousTranslation*getForce(scene, wall_id_Z1)[2];
	cachedStateZ1->vel[2] = previousTranslation / scene->dt;
}

void DirectShearBox::action()
{
	/*
	Main entry point, called once per iteration by the simulation loop.
	It (1) initializes volumes on the first step, (2) periodically refreshes
	porosity/height/stiffness/stress, then (3) runs the compaction or shear
	controller depending on the current mode.
	*/
	// sync thread storage of ForceContainer
	scene->forces.sync();

	if (first) {
		// On the very first iteration: compute the solid (particle) volume and
		// record the reference free height; then start internal compaction from initiald50.
		particlesVolume = 0;
		const size_t nbodies = scene->bodies->size();
		Real local_volume = 0;
#ifdef YADE_OPENMP
#pragma omp parallel for reduction(+:local_volume) schedule(guided)
#endif
		for (size_t i = 0; i < nbodies; ++i) {
			const auto& b = (*scene->bodies)[i];
			if (b && b->shape->getClassName() == "Sphere" && b->isDynamic() && !b->isClumpMember()) {
				const shared_ptr<Sphere>& sphere = YADE_PTR_CAST<Sphere>(b->shape);
				local_volume += 1.3333333 * Mathr::PI * pow(sphere->radius, 3);
			}
		}
		particlesVolume = local_volume;
		// Initialize cached pointers for performance
		cachedWallZ1 = Body::byId(wall_id_Z1, scene);
		cachedWallZ0 = Body::byId(wall_id_Z0, scene);
		cachedStateZ1 = cachedWallZ1->state.get();
		cachedStateZ0 = cachedWallZ0->state.get();
		heightZ0 = cachedStateZ1->se3.position.z() - cachedStateZ0->se3.position.z() - thicknessZ;
		currentd50 = initiald50;
		first = false;
	}
	// Refresh porosity and free height. Done only at a coarse interval because it
	// is only needed for data output, not for the control loop.
	if (scene->iter % porosityHeightZUpdateInterval == 0) 
	{
		if (!cachedWallZ1 || !cachedStateZ1 || !cachedWallZ0 || !cachedStateZ0) {
			cachedWallZ1 = Body::byId(wall_id_Z1, scene);
			cachedWallZ0 = Body::byId(wall_id_Z0, scene);
			cachedStateZ1 = cachedWallZ1->state.get();
			cachedStateZ0 = cachedWallZ0->state.get();
		}
		heightZ = cachedStateZ1->se3.position.z() - cachedStateZ0->se3.position.z() - thicknessZ;
		boxVolume = depthY0 * widthX0 * heightZ;
		porosity = (boxVolume - particlesVolume - pipeVolume) / (boxVolume - pipeVolume);
	}

	//Update stiffness only if it has been computed by StiffnessCounter (see "stiffnessUpdateInterval")
	if (scene->iter % stiffnessUpdateInterval == 0 || scene->iter < 100) updateStiffness();
	//Update stress
	if (scene->iter % computeStressInterval == 0 ) computeStress();

	if (!internalCompaction)   // ---- EXTERNAL compaction (walls move, particles fixed) ----
	{
		// Servo-control the pressor so its normal force matches pressorGoal*area.
		controlExternalStress(pressorGoal * widthX0 * depthY0); //move pressor. always do this in external condition.
		if (sampling) // step through targetPressures once the current target is reached
		{
			pressorGoal=targetPressures[samplingNextIndex];
			if ((ComputeUnbalancedForce(false)<=samplingUnbThreshold) and (abs((stress[wallsZ1][2]-pressorGoal)/pressorGoal)<=samplingPressureThreshold))
				if (abs(cachedStateZ1->vel[2])<=samplingVelThreshold)
				{
					std::cout<<"Target pressure "<<pressorGoal<<" has been reached."<< endl;
					sampling =false;
					Omega::instance().saveSimulation(savePath + "externaledP" + boost::lexical_cast<string>(pressorGoal) + ".yade.bz2");
					sampling = true;
					samplingNextIndex += 1;
					if (samplingNextIndex==(int)(targetPressures.size()))
					{
						std::cout<<"Sampling finished."<< endl;
						sampling=false;

					}
					else
					{
						pressorGoal=targetPressures[samplingNextIndex];
						std::cout<<"pressorGoal is set as "<<pressorGoal<<"."<< endl;
					}
				}
		}
		else if (!axialloading && autoSimulation)
		{
			// No sampling: once the single pressorGoal is reached, auto-start shearing.
			if ((ComputeUnbalancedForce(false)<=samplingUnbThreshold) and (abs((stress[wallsZ1][2]-pressorGoal)/pressorGoal)<=samplingPressureThreshold))
				if (abs(cachedStateZ1->vel[2])<=samplingVelThreshold)
				{
					std::cout<<"Target pressure "<<pressorGoal<<" has been reached. Axial loading begin."<< endl;
					setExternal();
					setAxial();
					axialAction();
				}
		}
		else if (scene->iter % controllerInterval == 0 && axialloading) controllerAction();
	}
	else // ---- INTERNAL compaction (particles grow, walls fixed) ----
	{
		if (!cachedStateZ1 || !cachedWallZ1) {
			cachedWallZ1 = Body::byId(wall_id_Z1, scene);
			cachedStateZ1 = cachedWallZ1->state.get();
		}
		cachedStateZ1->vel  = Vector3r::Zero(); // don't move pressor
		if (scene->iter % radiusControlInterval == 0)
		{
			if ((targetd50>0) and (targetd50>=currentd50))
			{
				// Grow particles only once the assembly has re-settled (low unbalanced force).
				Real unb=ComputeUnbalancedForce(false);
				if (isnan(unb) or (unb<internalUnbThreshold)){
					// ratio = progress from initiald50 toward targetd50 in [0,1).
					// Pick a progressively smaller enlarging factor (f0 > f1 > f2 > f3) as
					// the target is approached, so the final packing is not over-stressed.
					Real ratio=(currentd50-initiald50)/(targetd50-initiald50);
					if (ratio<1) {
						if (ratio<r1) previousMultiplier=1+f0*(targetd50-initiald50)/currentd50;
						else if (ratio<r2) previousMultiplier=1+f1*(targetd50-initiald50)/currentd50;
						else if (ratio<r3) previousMultiplier=1+f2*(targetd50-initiald50)/currentd50;
						else previousMultiplier=1+f3*(targetd50-initiald50)/currentd50;
					}
					controlInternalStress(previousMultiplier);
				}
			}
			else if (!internalCompactionSaved) 
			{
				// Internal compaction target reached (or explicitly skipped): save the state.
				if (targetd50<0)
				{
					std::cout<<"targetd50<0. no need to internal compact. "<< endl << "You may want to set setAllFrictions(), setAllEta() if 6D, sampling=True, targetPressures=[17e3,34e3,50e3,100e3] or pressorGoal=17e3, internalCompaction = false."<< endl;
					Omega::instance().saveSimulation(savePath + "nointernal.yade.bz2");
				}
				else
				{
					std::cout<<"Internal Compaction finished. d50 = "<<currentd50<<", porosity = "<< porosity << endl << "You may want to set setAllFrictions(), setAllEta() if 6D, sampling=True, targetPressures=[17e3,34e3,50e3,100e3] or pressorGoal=17e3, internalCompaction = false."<< endl;
					Omega::instance().saveSimulation(savePath + "internaledD50" + boost::lexical_cast<string>(currentd50) + ".yade.bz2");
					sp->fromSimulation();
					sp->toFile(savePath + "Par_internaledD50" + boost::lexical_cast<string>(currentd50) + ".txt");
				}
				internalCompactionSaved=true;
			}
			else if (autoSimulation) 
			{
				// After internal compaction, automatically switch to external (wall) compaction.
				std::cout<<"External Compaction begins."<< endl;
				setExternal();
				internalCompaction=false;
			}
		}
	}
	
}

void DirectShearBox::computeStress()
{
	/*
	Recompute resultant stresses/forces on each active boundary from the current
	contact forces, update meanStress (averaged over present boundaries), and
	accumulate externalWork. Called every computeStressInterval iterations.
	*/
	// Ensure cached pointers are valid
	if (!cachedStateZ1 || !cachedStateZ0 || !cachedWallZ1 || !cachedWallZ0) {
		cachedWallZ1 = Body::byId(wall_id_Z1, scene);
		cachedWallZ0 = Body::byId(wall_id_Z0, scene);
		cachedStateZ1 = cachedWallZ1->state.get();
		cachedStateZ0 = cachedWallZ0->state.get();
	}
	
	Real newHeightZ = cachedStateZ1->se3.position.z() - cachedStateZ0->se3.position.z() - thicknessZ;
	
	// Only recalculate surface areas if height changed (performance optimization)
	if (!cachedHeightZ || newHeightZ != cachedHeightZ || !cachedInvXSurface || !cachedInvYSurface) {
		cachedHeightZ = newHeightZ;
		cachedInvXSurface = 1.f / (depthY0 * newHeightZ);
		cachedInvYSurface = 1.f / (widthX0 * newHeightZ);
	}
	cachedInvZSurface = 1.f / (depthY0 * widthX0);
	heightZ = cachedHeightZ;

	if (wall_ids_Y0.size()>0) {
		force[wallsY0]=Vector3r::Zero();
		for (const int& id : wall_ids_Y0) force[wallsY0]  += getForce(scene, id);
		stress[wallsY0] = -force[wallsY0] * cachedInvYSurface;
	}

	if (wall_ids_Y1.size()>0) {
		force[wallsY1]=Vector3r::Zero();
		for (const int& id : wall_ids_Y1) force[wallsY1]  += getForce(scene, id);
		stress[wallsY1] = force[wallsY1] * cachedInvYSurface;
	}

	if (wall_ids_X0.size()>0) {
		force[wallsX0]=Vector3r::Zero();
		for (const int& id : wall_ids_X0) force[wallsX0]  += getForce(scene, id);
		stress[wallsX0] = -force[wallsX0] * cachedInvXSurface;
	}

	if (wall_ids_X1.size()>0) {
		force[wallsX1]=Vector3r::Zero();
		for (const int& id : wall_ids_X1) force[wallsX1]  += getForce(scene, id);
		stress[wallsX1] = force[wallsX1] * cachedInvXSurface;
	}

	force[wallsZ1]   = getForce(scene, wall_id_Z1);
	stress[wallsZ1]  = force[wallsZ1] * cachedInvZSurface;
	force[wallsZ0]    = getForce(scene, wall_id_Z0);
	stress[wallsZ0]   = -force[wallsZ0] * cachedInvZSurface;

	meanStress = stress[wallsX0][0] + stress[wallsX1][0] + stress[wallsY0][1] + stress[wallsY1][1] + stress[wallsZ0][2] + stress[wallsZ1][2];
	meanStress /= (2.0+(wall_ids_X0.size()!=0)+(wall_ids_X1.size()!=0)+(wall_ids_Y0.size()!=0)+(wall_ids_Y1.size()!=0));// ( sXX(xLeft) + sXX(xRight) + sYY(yBottom) + sYY(yTop) + sZZ(zBack) + sZZ(zFront) ) / 6
}

void DirectShearBox::controlInternalStress(Real multiplier)
{
	/*
	Apply the internal-compaction growth step: scale every particle's radius by
	'multiplier' (via growParticles) so the packing densifies without moving walls.
	Also scales the tracked particle volume and current d50 consistently.
	*/
	particlesVolume *= pow(multiplier, 3);
	growParticles(multiplier, true, true);
	currentd50 *= multiplier;
}

Real DirectShearBox::ComputeUnbalancedForce(bool maxUnbalanced) {
	/* Wrapper around Shop::unbalancedForce(): returns the current unbalanced-force ratio of
	   the assembly (maxUnbalanced=true -> max over contacts, otherwise the mean). */
    return Shop::unbalancedForce(maxUnbalanced, scene);
}

void DirectShearBox::growParticles(Real multiplier, bool updateMass, bool dynamicOnly)
{
	/*
	Grow (or shrink) every particle radius by 'multiplier'. When updateMass is true,
	also scale mass and inertia consistently. When dynamicOnly is true, clumps and
	their members are skipped. The contact reference radii (refR1/refR2) are updated
	so existing interactions keep their relative overlap.
	*/
	for (const auto& b : *scene->bodies) {
		if (dynamicOnly && !b->isDynamic()) continue;
		if (b->isClump() or b->isClumpMember()) continue;
		//We grow only spheres and clumps
		if (b->shape->getClassName() == "Sphere") {
			if (updateMass) {
				b->state->mass *= pow(multiplier, 3);
				b->state->inertia *= pow(multiplier, 5);
			}
			// for spheres, we update radius
			(YADE_CAST<Sphere*>(b->shape.get()))->radius *= multiplier;
		}
	}
	FOREACH(const shared_ptr<Interaction>& ii, *scene->interactions)
	{
		if (ii->isReal()) {
			shared_ptr<Body>& b1 = (*(scene->bodies))[ii->getId1()];
			shared_ptr<Body>& b2 = (*(scene->bodies))[ii->getId2()];
			const shared_ptr<FrictMat>& sdec1 = YADE_PTR_CAST<FrictMat>(b1->material);
			const shared_ptr<FrictMat>& sdec2 = YADE_PTR_CAST<FrictMat>(b2->material);
			
			FrictPhys* contactPhysics = YADE_CAST<FrictPhys*>(ii->phys.get());
			ScGeom*     geom = YADE_CAST<ScGeom*>(ii->geom.get());
			if ((!dynamicOnly || b1->isDynamic()) && b1->shape->getClassName() == "Sphere")
				geom->refR1 = YADE_CAST<Sphere*>(b1->shape.get())->radius;
			if ((!dynamicOnly || b2->isDynamic()) && b2->shape->getClassName() == "Sphere")
				geom->refR2 = YADE_CAST<Sphere*>(b2->shape.get())->radius;
			Real          Ea             = sdec1->young;
			Real          Eb             = sdec2->young;
			Real          Va             = sdec1->poisson;
			Real          Vb             = sdec2->poisson;
			Real          Da             = geom->radius1;
			Real          Db             = geom->radius2;
			Real          Kn             = 2.0 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); 
			Real Ks;
			if (Va && Vb)
				Ks = 2.0 * Ea * Da * Va * Eb * Db * Vb / (Ea * Da * Va + Eb * Db * Vb); 
			else
				Ks = 0;

			contactPhysics->kn = Kn;
			contactPhysics->ks = Ks;
		}
	}
}


void DirectShearBox6D::growParticles(Real multiplier, bool updateMass, bool dynamicOnly)
{
	/* 
		Grow all non-clumpmember spheres by multiplier. always update radius in body's shape. and kn ks kr ktw maxRollPl maxTwistPl in contactPhysics. 
		Update mass and inertia if updateMass. Only grow dynamic sphere if dynamicOnly.
	*/

	const size_t nbodies = scene->bodies->size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided)
#endif
	for (size_t i = 0; i < nbodies; ++i) {
		const auto& b = (*scene->bodies)[i];
		if (!b) continue;  // Skip empty slots
		if (dynamicOnly && !b->isDynamic()) continue;
		if (b->isClump() or b->isClumpMember()) continue;
		//We grow only spheres and clumps
		if (b->shape->getClassName() == "Sphere") {
			if (updateMass) {
				b->state->mass *= pow(multiplier, 3);
				b->state->inertia *= pow(multiplier, 5);
			}
			// for spheres, we update radius
			(YADE_CAST<Sphere*>(b->shape.get()))->radius *= multiplier;
		}
	}
	FOREACH(const shared_ptr<Interaction>& ii, *scene->interactions)
	{
		if (ii->isReal()) {
			shared_ptr<Body>& b1 = (*(scene->bodies))[ii->getId1()];
			shared_ptr<Body>& b2 = (*(scene->bodies))[ii->getId2()];
			const shared_ptr<CohFrictMat>& sdec1 = YADE_PTR_CAST<CohFrictMat>(b1->material);
			const shared_ptr<CohFrictMat>& sdec2 = YADE_PTR_CAST<CohFrictMat>(b2->material);
			
			CohFrictPhys* contactPhysics = YADE_CAST<CohFrictPhys*>(ii->phys.get());
			ScGeom6D*     geom = YADE_CAST<ScGeom6D*>(ii->geom.get());
			if ((!dynamicOnly || b1->isDynamic()) && b1->shape->getClassName() == "Sphere")
				geom->refR1 = YADE_CAST<Sphere*>(b1->shape.get())->radius;
			if ((!dynamicOnly || b2->isDynamic()) && b2->shape->getClassName() == "Sphere")
				geom->refR2 = YADE_CAST<Sphere*>(b2->shape.get())->radius;
			Real          Ea             = sdec1->young;
			Real          Eb             = sdec2->young;
			Real          Va             = sdec1->poisson;
			Real          Vb             = sdec2->poisson;
			Real          Da             = geom->radius1;
			Real          Db             = geom->radius2;
			Real          Kn             = 2.0 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); //harmonic average of two stiffnesses
			// harmonic average of alphas parameters
			Real AlphaKr, AlphaKtw;
			if (sdec1->alphaKr && sdec2->alphaKr) AlphaKr = 2.0 * sdec1->alphaKr * sdec2->alphaKr / (sdec1->alphaKr + sdec2->alphaKr);
			else
				AlphaKr = 0;
			if (sdec1->alphaKtw && sdec2->alphaKtw) AlphaKtw = 2.0 * sdec1->alphaKtw * sdec2->alphaKtw / (sdec1->alphaKtw + sdec2->alphaKtw);
			else
				AlphaKtw = 0;

			Real Ks;
			if (Va && Vb)
				Ks = 2.0 * Ea * Da * Va * Eb * Db * Vb
				        / (Ea * Da * Va + Eb * Db * Vb); //harmonic average of two stiffnesses with ks=V*kn for each sphere
			else
				Ks = 0;

			contactPhysics->kr                     = Da * Db * Ks * AlphaKr;
			contactPhysics->ktw                    = Da * Db * Ks * AlphaKtw;

			contactPhysics->kn = Kn;
			contactPhysics->ks = Ks;

			contactPhysics->maxRollPl         = math::min(sdec1->etaRoll * Da, sdec2->etaRoll * Db);
			contactPhysics->maxTwistPl        = math::min(sdec1->etaTwist * Da, sdec2->etaTwist * Db);
		}
	}
}

void DirectShearBox::controllerAction()
{
	/*
	Periodic shear controller, called every controllerInterval while axialloading is true.
	Optionally re-applies pipe friction at a slower cadence, then checks the shear-box
	displacement against the current target (currentTargetY) to decide saving / direction switch.
	*/
	if (targetPipeFriction>0 && pipefrictionUpdateNum>0) {if (scene->iter % (controllerInterval * pipefrictionUpdateNum) == 0) setInteractionsFrictionOnClump(targetPipeFriction, shearbox_id);}
	State* s_pipe = Body::byId(shearbox_id, Omega::instance().getScene().get())->state.get();
	// relativePos in [0,1]: fraction of the current stroke (from the start target to currentTargetY) already travelled.
	Real relativePos = currentTargetY-s_pipe->displ()[1];
	relativePos = 1-relativePos/(targetYs[(int)(loading)]-targetYs[(int)(!loading)]);
	if ((saverNum>0) and (relativePos>=currentSaverNum/saverNum))
	{
		Omega::instance().saveSimulation(savePath+ "P" + boost::lexical_cast<string>(pressorGoal)+ "N" + boost::lexical_cast<string>(currentN) + ((loading) ? "Loading" : "Unloading")+ "n"+boost::lexical_cast<string>(currentSaverNum)+ "Y"+boost::lexical_cast<string>(s_pipe->displ()[1])+ ".yade.bz2");
		currentSaverNum+=1;
	}

	if (relativePos>=1) // reach target Y
	{
		if (loading==0) currentN+=1; //if loading==0, this cycle finish

		if (targetN!=0)
		{
			loading = not loading; // switch loading condition
			currentSaverNum=0;
		}
		if (targetN<=currentN) //cyclic loading finish?
		{

			axialloading=false;
			s_pipe->vel[1]=0;
			std::cout << "axial loading finish." <<endl;
			Omega::instance().pause();
			Omega::instance().saveSimulation(savePath + "N" + boost::lexical_cast<string>(currentN) + "finished_iter" + boost::lexical_cast<string>(scene->iter) + ".yade.bz2");
		}
		else axialAction();
	}

}

void DirectShearBox::axialAction() 
{
	/*
	Start (or continue) the shear phase. Called by setAxial() and by the
	autoSimulation branch in action() once compaction is finished.
	*/
	if (targetPipeFriction>0 && pipefrictionUpdateNum>0)setInteractionsFrictionOnClump(targetPipeFriction, shearbox_id);
	// Target Y position for the shear box: targetYs[1] while loading, targetYs[0] while unloading.
    currentTargetY=targetYs[(int) loading];
	std::cout << "current target Y is set to " << currentTargetY << ".";
	State* s_pipe = Body::byId(shearbox_id, Omega::instance().getScene().get())->state.get();
	// Drive the shear box at +/- targetShearRate along Y depending on the loading direction.
	s_pipe->vel[1]  = ((loading) ? 1.0 : -1.0)*targetShearRate;
	std::cout << "current shear rate is set to " << s_pipe->vel[1] <<endl;
	axialloading=true;
	// Forcefully end compaction/sampling: once shearing starts these modes no longer apply.
	if (internalCompaction)
	{
		internalCompaction=false;
		std::cout << "internal compaction are ended by force." << endl;
	}
	if (sampling)
	{
		sampling=false;
		std::cout << "sampling are ended by force." << endl;
	}
}

void DirectShearBox::setBodyFriction(Real angleRad, Body::id_t bodyID)
{
	/*
	Set the material friction angle [rad] of a single body. If the body is a clump,
	the angle is applied to every clump member (so the clump surface uses one friction).
	*/
	const shared_ptr<Body>& b   = Body::byId(bodyID,scene);
	if (b->isClump()) {
		FOREACH(const auto& mm, YADE_PTR_CAST<Clump>(b->shape)->members)
		{
			// mm.first is Body::id_t, mm.second is local Se3r of that body in the clump
			const Body::id_t&       memberId = mm.first;
			const shared_ptr<Body>& member   = Body::byId(memberId,scene);
			YADE_PTR_CAST<FrictMat>(member->material)->frictionAngle = angleRad;
		}
	}
	else YADE_PTR_CAST<FrictMat>(b->material)->frictionAngle = angleRad;
}

void DirectShearBox::setBodiesFriction(Real angleRad, vector<Body::id_t> bodyIDs){
	/*
	Apply setBodyFriction to every body in the list (convenience loop over setBodyFriction).
	*/
    for (Body::id_t bodyID : bodyIDs) setBodyFriction(angleRad,bodyID);}

void DirectShearBox::updateInteractionsFriction()
{
	/*
	Re-derive every contact's tangensOfFrictionAngle from the current materials of the
	two bodies (min of the two friction angles). Call this AFTER changing any material
	friction, otherwise existing contacts keep their old friction.
	*/
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	FOREACH(const shared_ptr<Interaction>& ii, *scene->interactions)
	{
		if (!ii->isReal()) continue;
		const shared_ptr<FrictMat>& sdec1 = YADE_PTR_CAST<FrictMat>((*bodies)[(Body::id_t)((ii)->getId1())]->material);
		const shared_ptr<FrictMat>& sdec2 = YADE_PTR_CAST<FrictMat>((*bodies)[(Body::id_t)((ii)->getId2())]->material);
		//FIXME - why dynamic_cast fails here?
		FrictPhys*  contactPhysics             = YADE_CAST<FrictPhys*>((ii)->phys.get());
		const Real& fa                         = sdec1->frictionAngle;
		const Real& fb                         = sdec2->frictionAngle;
		contactPhysics->tangensOfFrictionAngle = math::tan(math::min(fa, fb));
	}
}

void DirectShearBox::setWallFriction(Real angleRad, Body::id_t bodyID)
{
	/*
	Set one wall's material friction, then refresh its contacts.
	*/
	setBodyFriction(angleRad, bodyID);
	updateInteractionsFriction();
}

void DirectShearBox::setWallsFriction(Real angleRad, vector<Body::id_t> bodyIDs)
{
	/*
	Set a list of walls' material friction, then refresh their contacts.
	*/
	setBodiesFriction(angleRad, bodyIDs);
	updateInteractionsFriction();
}

void DirectShearBox::setSpheresFrictionOnly(Real angleRad)
{
	/*
	Set the material friction angle (radians) of every free (non-clump-member) dynamic
	sphere only (does not touch walls or the pipe).
	*/
	shared_ptr<BodyContainer>& bodies = Omega::instance().getScene().get()->bodies;
	const size_t nbodies = bodies->size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided)
#endif
	for (size_t i = 0; i < nbodies; ++i) {
		const auto& b = (*bodies)[i];
		if (b && b->isDynamic() && !b->isClumpMember() && b->shape && b->shape->getClassName() == "Sphere")
			YADE_PTR_CAST<FrictMat>(b->material)->frictionAngle = angleRad;
	}
}

void DirectShearBox::setSpheresFriction(Real angleRad)
{
	/*
	Set the material friction of all spheres, then refresh every contact so they pick up the new friction.
	*/
	setSpheresFrictionOnly(angleRad);
	updateInteractionsFriction();
}

void DirectShearBox::setAllFriction()
{
	/*
	Apply every target* friction attribute that is >= 0 (negative target means "leave as is"),
	then refresh all contacts once. This is the main entry used by setExternal().
	*/
	bool hasChanges = false;
	if (targetWallsFriction_X0>=0) 
	{
		setBodiesFriction(targetWallsFriction_X0,wall_ids_X0);
		std::cout << "Friction of walls x0 is set as "<< targetWallsFriction_X0 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallsFriction_X1>=0) 
	{
		setBodiesFriction(targetWallsFriction_X1,wall_ids_X1);
		std::cout << "Friction of walls x1 is set as "<< targetWallsFriction_X1 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallsFriction_Y0>=0)
	{
		setBodiesFriction(targetWallsFriction_Y0,wall_ids_Y0);
		std::cout << "Friction of walls y0 is set as "<< targetWallsFriction_Y0 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallsFriction_Y1>=0) 
	{
		setBodiesFriction(targetWallsFriction_Y1,wall_ids_Y1);
		std::cout << "Friction of walls y1 is set as "<< targetWallsFriction_Y1 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallFriction_Z0>=0) 
	{
		setBodyFriction(targetWallFriction_Z0,wall_id_Z0);
		std::cout << "Friction of wall z0 is set as "<< targetWallFriction_Z0 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallFriction_Z1>=0)
	{
		setBodyFriction(targetWallFriction_Z1,wall_id_Z1);
		std::cout << "Friction of wall z1 is set as "<< targetWallFriction_Z1 <<"." <<endl;
		hasChanges = true;
	}
	if (targetPipeFriction>=0)
	{
		setBodyFriction(targetPipeFriction,shearbox_id);
		std::cout << "Friction of shearbox(pipe) is set as "<< targetPipeFriction <<"." <<endl;
		hasChanges = true;
	}
	if (targetParticlesFriction>=0) 
	{
		setSpheresFrictionOnly(targetParticlesFriction);
		std::cout << "Friction of particles is set as "<< targetParticlesFriction <<"." <<endl;
		hasChanges = true;
	}
	if (hasChanges) updateInteractionsFriction();
}

void DirectShearBox::resetRefSe3()
{
	/*
	Reset every body's reference pose (refPos/refOri) to its current pose, and the
	periodic cell's refHSize. After this, reported displacements start from zero.
	Called by setExternal() just before shearing begins.
	*/
	const shared_ptr<Scene>& scene_ptr = Omega::instance().getScene();
	const size_t nbodies = scene_ptr->bodies->size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided)
#endif
	for (size_t i = 0; i < nbodies; ++i) {
		const auto& b = (*scene_ptr->bodies)[i];
		if (b) {
			b->state->refPos = b->state->pos;
			b->state->refOri = b->state->ori;
		}
	}
	if (scene_ptr->isPeriodic) { scene_ptr->cell->refHSize = scene_ptr->cell->hSize; }
	std::cout << "All refSe3 are reset." <<endl;
}

void DirectShearBox::setPipeDOFs(std::string DOFs)
{
	/*
	Lock the blocked degrees of freedom of the shear box / pipe (e.g. "xz" allows only Y motion)
	so the pipe can only travel along the shear direction.
	*/
	State* s_pipe = Body::byId(shearbox_id, Omega::instance().getScene().get())->state.get();
	s_pipe->blockedDOFs_vec_set(DOFs);
	std::cout << "Pipe DOFs have been set as " <<s_pipe->blockedDOFs_vec_get() <<endl;
}

void DirectShearBox::setCycle (int N)
{
	/*
	Set the target cyclic number targetN (0 => monotonic loading, >0 => periodic loading/unloading).
	*/
	targetN = N;
	std::cout << "target N is set as " << targetN << "." << endl;
}

void DirectShearBox::setLoading (bool L) 
{
	/*
	Set the loading direction flag: true drives the shear box toward targetYs[1] (+rate),
	false drives it back toward targetYs[0] (-rate, unloading).
	*/
	loading = L;
	std::cout << ((L) ? "Loading condition is set as loading." : "Loading condition is set as unloading.") << endl;
}

void DirectShearBox::setY (vector<Real> Ys)
{
	/*
	Set the cyclic shear targets targetYs. The two values are auto-ordered so that
	targetYs[0] (unloading) < targetYs[1] (loading), regardless of input order.
	*/
	if (Ys[0]>Ys[1])
	{
		Real Y = Ys[1];
		Ys[1] = Ys[0];
		Ys[0] = Y;
	}
	targetYs[0]= Ys[0] ;
	targetYs[1]= Ys[1] ;
	std::cout << "target Y is set as : " << targetYs[0] << "," << targetYs[1] << endl;
}

void DirectShearBox::setV (Real V) 
{
	/*
	Set the absolute target shear rate targetShearRate (stored as abs(V); the sign is
	applied later during shearing depending on the loading flag).
	*/
	targetShearRate = abs(V);
	std::cout << "target V is set as : " << abs(V) << endl;
}



void DirectShearBox::setCurrentVY (Real V)
{
	/*
	Override the current shear-box Y velocity for THIS loading/unloading stroke only
	(does not change targetShearRate). The sign follows the loading flag.
	*/
	State* s_pipe = Body::byId(shearbox_id, Omega::instance().getScene().get())->state.get();
	s_pipe->vel[1]  = ((loading) ? 1.0 : -1.0)*abs(V);
	std::cout << "current shear rate is set to " << s_pipe->vel[1] <<". Only for this loading/unloading. use setV() to change. "<<endl;
}

Real DirectShearBox::getCurrentVY () 
{
	/*
	Read back the current shear-box Y velocity (used to monitor the shear rate).
	*/
	return Body::byId(shearbox_id, Omega::instance().getScene().get())->state.get()->vel[1];
}

void DirectShearBox::setGravityZ (Real gravityZ)
{
	/*
	Set the gravity vector to [0,0,gravityZ] on the scene's NewtonIntegrator.
	The NewtonIntegrator is located lazily on first call if not yet cached.
	*/
	if (!newton) 
	{
		vector<shared_ptr<Engine>>::iterator itFirst = Omega::instance().getScene().get()->engines.begin();
		vector<shared_ptr<Engine>>::iterator itLast  = Omega::instance().getScene().get()->engines.end();
		for (; itFirst != itLast; ++itFirst) {
			if ((*itFirst)->getClassName() == "NewtonIntegrator") {
				LOG_DEBUG("NewtonIntegrator engine found");
				newton = YADE_PTR_CAST<NewtonIntegrator>(*itFirst);
			}
		}
	}
	if (!newton) std::cout << "NewtonIntegrator engine NOT found" <<endl;
	else
	{
		newton->gravity=Vector3r(0,0,gravityZ);
		std::cout << "Gravity is "<< newton->gravity <<"." <<endl;
	}
}

void DirectShearBox::setDamping (Real D)
{
	/*
	Set the global damping of the scene's NewtonIntegrator (located lazily on first call).
	*/
	if (!newton) 
	{
		vector<shared_ptr<Engine>>::iterator itFirst = Omega::instance().getScene().get()->engines.begin();
		vector<shared_ptr<Engine>>::iterator itLast  = Omega::instance().getScene().get()->engines.end();
		for (; itFirst != itLast; ++itFirst) {
			if ((*itFirst)->getClassName() == "NewtonIntegrator") {
				LOG_DEBUG("NewtonIntegrator engine found");
				newton = YADE_PTR_CAST<NewtonIntegrator>(*itFirst);
			}
		}
	}
	if (!newton) std::cout << "NewtonIntegrator engine NOT found" <<endl;
	else
	{
		newton->damping=D;
		std::cout << "Damping is set as "<< newton->damping <<"." <<endl;
	}
}

void DirectShearBox::setInteractionsFrictionOnClump(Real angleRad, Body::id_t bodyID)
{
	/*
	Set the contact friction of every interaction touching a clump body (e.g. the pipe),
	by directly overwriting tangensOfFrictionAngle on load-bearing contacts (fn != 0).
	Unlike setBodyFriction, this changes the contact physics, not the material.
	*/
	const shared_ptr<Body>& b   = Body::byId(bodyID,scene);
	if (b->isClump()) {
		InteractionContainer::iterator ii    = scene->interactions->begin();
		InteractionContainer::iterator iiEnd = scene->interactions->end();
		for (; ii != iiEnd; ++ii)
			if ((*ii)->isReal()) {
				const shared_ptr<Interaction>& contact = *ii;
				Real                           fn      = (static_cast<FrictPhys*>(contact->phys.get()))->normalForce.norm(); // TODO:capillary?
				if (fn != 0) {
					Body::id_t id1 = contact->getId1(), id2 = contact->getId2();
					FOREACH(const auto& mm, YADE_PTR_CAST<Clump>(b->shape)->members)
					{
						// mm.first is Body::id_t, mm.second is local Se3r of that body in the clump
						const Body::id_t&       memberId = mm.first;
						if (memberId == id1 || memberId == id2) {
							FrictPhys* currentContactPhysics = static_cast<FrictPhys*>(contact->phys.get());
							currentContactPhysics->tangensOfFrictionAngle = math::tan(angleRad);
						}
					}
				}
			}

	}

}

void DirectShearBox::setInteractionsFrictionOnBody(Real angleRad, Body::id_t bodyID)
{
	/*
	Set the contact friction of every interaction touching one body (fn != 0), by directly
	overwriting tangensOfFrictionAngle. Used to change friction on existing contacts without
	touching materials (e.g. for the pipe during shearing).
	*/
	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();
	for (; ii != iiEnd; ++ii)
		if ((*ii)->isReal()) {
			const shared_ptr<Interaction>& contact = *ii;
			Real                           fn      = (static_cast<FrictPhys*>(contact->phys.get()))->normalForce.norm(); // TODO:capillary?
			if (fn != 0) {
				Body::id_t id1 = contact->getId1(), id2 = contact->getId2();
				if (bodyID == id1 || bodyID == id2) {
					FrictPhys* currentContactPhysics = static_cast<FrictPhys*>(contact->phys.get());
					currentContactPhysics->tangensOfFrictionAngle = math::tan(angleRad);
				}
			}
		}
}
void DirectShearBox::setInteractionsFrictionOnBodies(Real angleRad, vector<Body::id_t> bodyIDs)
{
	// Performance optimization: avoid nested loops by using a set and single interaction pass
	std::unordered_set<Body::id_t> bodyIdSet(bodyIDs.begin(), bodyIDs.end());
	Real tangensAngle = math::tan(angleRad);
	
	for (auto& ii : *scene->interactions) {
		if (!ii->isReal()) continue;
		Body::id_t id1 = ii->getId1(), id2 = ii->getId2();
		if (bodyIdSet.count(id1) || bodyIdSet.count(id2)) {
			FrictPhys* phys = static_cast<FrictPhys*>(ii->phys.get());
			phys->tangensOfFrictionAngle = tangensAngle;
		}
	}
}

void DirectShearBox::setAllInteractionsFrictionOnWalls()
{
	/*
	Apply the target* wall frictions directly to the existing wall contacts (fn != 0),
	using setInteractionsFrictionOnBodies/Body. Unlike setAllFriction, this does NOT change
	wall materials, only the contact tangensOfFrictionAngle. Called by setAxial().
	*/
	if (targetWallsFriction_X0>=0) 
	{
		setInteractionsFrictionOnBodies(targetWallsFriction_X0,wall_ids_X0);
		std::cout << "Interactions friction with walls x0 is set as "<< targetWallsFriction_X0 <<"." <<endl;
	}
	if (targetWallsFriction_X1>=0) 
	{
		setInteractionsFrictionOnBodies(targetWallsFriction_X1,wall_ids_X1);
		std::cout << "Interactions friction with walls x1 is set as "<< targetWallsFriction_X1 <<"." <<endl;
	}
	if (targetWallsFriction_Y0>=0)
	{
		setInteractionsFrictionOnBodies(targetWallsFriction_Y0,wall_ids_Y0);
		std::cout << "Interactions friction with walls y0 is set as "<< targetWallsFriction_Y0 <<"." <<endl;
	}
	if (targetWallsFriction_Y1>=0) 
	{
		setInteractionsFrictionOnBodies(targetWallsFriction_Y1,wall_ids_Y1);
		std::cout << "Interactions friction with walls y1 is set as "<< targetWallsFriction_Y1 <<"." <<endl;
	}
	if (targetWallFriction_Z0>=0) 
	{
		setInteractionsFrictionOnBody(targetWallFriction_Z0,wall_id_Z0);
		std::cout << "Interactions friction with wall z0 is set as "<< targetWallFriction_Z0 <<"." <<endl;
	}
	if (targetWallFriction_Z1>=0)
	{
		setInteractionsFrictionOnBody(targetWallFriction_Z1,wall_id_Z1);
		std::cout << "Interactions friction with wall z1 is set as "<< targetWallFriction_Z1 <<"." <<endl;
	}
}

void DirectShearBox::setAxial()
{
	/*
	Switch the engine into shear mode: push the cycle/y/velocity/loading targets to
	the internal controllers, re-apply wall frictions, reset the reference height,
	and (optionally) re-point the DirectShearBoxRecorder to a "_shearing.csv" file so
	per-step data is logged. Called once before shearing starts.
	*/
	if (targetDamping>0) setDamping(targetDamping);
	setCycle(targetN);
	setY(targetYs);
	setV(targetShearRate);
	setLoading(loading);
	setAllInteractionsFrictionOnWalls();
	if (!cachedWallZ1 || !cachedWallZ0 || !cachedStateZ1 || !cachedStateZ0) {
		cachedWallZ1 = Body::byId(wall_id_Z1, scene);
		cachedWallZ0 = Body::byId(wall_id_Z0, scene);
		cachedStateZ1 = cachedWallZ1->state.get();
		cachedStateZ0 = cachedWallZ0->state.get();
	}
	heightZ0 = cachedStateZ1->se3.position.z() - cachedStateZ0->se3.position.z() - thicknessZ;
	std::cout << "heightZ0 is reset." << endl;
	if (!dsbr) {
		vector<shared_ptr<Engine>>::iterator itFirst = Omega::instance().getScene().get()->engines.begin();
		vector<shared_ptr<Engine>>::iterator itLast  = Omega::instance().getScene().get()->engines.end();
		for (; itFirst != itLast; ++itFirst) {
			if ((*itFirst)->getClassName() == "DirectShearBoxRecorder") {
				LOG_DEBUG("DirectShearBoxRecorder engine found");
				dsbr = YADE_PTR_CAST<DirectShearBoxRecorder>(*itFirst);
				break;  // Found recorder, exit loop early
			}
		}
		if (!dsbr) {
			LOG_ERROR("DirectShearBoxRecorder engine NOT found");
			return;  // Early exit if recorder not found
		}
	}
	// Modify filename and reopen stream
	std::string fileTem=dsbr->file;
	int pos = fileTem.find(".csv");
	if (pos >-1) fileTem.replace(pos,4,"_shearing.csv");
	dsbr->file=fileTem;
	std::cout << "Output file set to: " << dsbr->file << endl;
	
	dsbr->out.close();
	std::string fileTemp = dsbr->file;
	if (dsbr->addIterNum) fileTemp += "-" + boost::lexical_cast<string>(scene->iter);

	if (fileTemp.empty()) throw ios_base::failure(__FILE__ ": Empty filename.");
	dsbr->out.open(fileTemp.c_str(), dsbr->truncate ? std::fstream::trunc : std::fstream::app);
	if (!dsbr->out.good()) throw ios_base::failure(__FILE__ ": I/O error opening file `" + fileTemp + "'.");
	dsbr->action();
}

void DirectShearBox::setcolors()
{
	/*
	Color every free dynamic sphere by its Y layer: even layers -> colorOne, odd -> colorTwo,
	using colorlayerNum equal-thickness layers across the box depth (depthY0). Purely for
	visualisation; only runs when setColor is enabled and colorlayerNum > 0.
	*/
	const Real gap=depthY0/colorlayerNum;
	shared_ptr<BodyContainer>& bodies = Omega::instance().getScene().get()->bodies;
	const size_t nbodies = bodies->size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided)
#endif
	for (size_t i = 0; i < nbodies; ++i) 
	{
		const auto& b = (*bodies)[i];
		if (b && b->isDynamic() && !b->isClumpMember() && b->shape && b->shape->getClassName() == "Sphere")
		{
			Real temp=b->state->pos.y()/gap;
			int tem=(int) temp;
			if (tem%2) b->shape->color=colorTwo;
			else b->shape->color=colorOne;
		}
	}
	std::cout << "Colors are set in "<< colorlayerNum <<" layers " << endl;
}

void DirectShearBox::setExternal()
{
	/*
	Apply all boundary conditions needed before the shear phase: optionally lock the pipe
	DOFs, apply gravity, set all friction angles, reset the reference pose used for displacement
	output, and color particles by layer. This is the main entry that switches the box from
	compaction to the ready-to-shear state.
	*/
	if (changePipeDOFornot) setPipeDOFs(targetBlockedDOFs);
	if (gravityAction) setGravityZ(targetGravityZ);
	setAllFriction();
	resetRefSe3();
	if ((setColor) and (colorlayerNum>0)) setcolors();
}

void DirectShearBox::showState()
{
	/*
	Print a human-readable summary of the current phase (internal compaction / sampling /
	external compaction / axial loading) together with the relevant targets and the live
	stress / unbalanced-force. Called automatically each step when logornot is true.
	*/
	scene->forces.sync();

	if (internalCompaction) {
		std::cout << "Internal compacting."<< endl;
		std::cout << "--d50: "<< currentd50 << "->"<< targetd50 << endl;
		std::cout << "--unb: "<< Shop::unbalancedForce(false, scene) << "->"<< samplingUnbThreshold << endl;
	}
	else {
		if (sampling) 
		{
			std::cout << " Sampling. Target pressures: ";
			for (auto tp = targetPressures.begin(); tp != targetPressures.end();++tp){std::cout << *tp << " ";}
			std::cout << endl;
		}
		else {std::cout << "External compacting. " << endl;}
		std::cout << "--pressor: "<< stress[wallsZ1][2] << "->"<< pressorGoal << endl;
		std::cout << "--unb: "<< Shop::unbalancedForce(false, scene) << endl;
	}
	if (axialloading) {
		State* pipeState = Body::byId(shearbox_id, scene)->state.get();
		std::cout << "Axial loading. " << endl;
		std::cout << "Y: " << targetYs[(int) (not loading)] << ((loading) ? "->" : "<-") << pipeState->displ()[1] <<"->" << targetYs[(int) loading] << endl;
		std::cout << "V: " << pipeState->vel[1] << "->" << ((loading) ? 1.0 : -1.0)*targetShearRate << endl;
		std::cout << "N: " << currentN << "->" << targetN << ((loading) ? ": loading" : ": unloading") << endl;
	}
	else{
		if (changePipeDOFornot) std::cout << "targetBlockedDOFs: " << targetBlockedDOFs << ". " << endl;
		if (gravityAction) std::cout << "targetGravityZ: " << targetGravityZ << "." << endl;
		if (targetParticlesFriction>=0) std::cout << "targetParticlesFriction: "<< targetParticlesFriction <<"." <<endl; 
		if (targetWallFriction_Z1>=0) std::cout << "targetWallFriction_Z1: "<< targetWallFriction_Z1 <<"." <<endl;
		if (targetWallsFriction_X0>=0) std::cout << "targetWallsFriction_X0: "<< targetWallsFriction_X0 <<"." <<endl; 
		if (targetWallsFriction_X1>=0) std::cout << "targetWallsFriction_X1: "<< targetWallsFriction_X1 <<"." <<endl; 
		if (targetWallsFriction_Y0>=0) std::cout << "targetWallsFriction_Y0: "<< targetWallsFriction_Y0 <<"." <<endl; 
		if (targetWallsFriction_Y1>=0) std::cout << "targetWallsFriction_Y1: "<< targetWallsFriction_Y1 <<"." <<endl; 
		if (targetWallFriction_Z0>=0) std::cout << "targetWallFriction_Z0: "<< targetWallFriction_Z0 <<"." <<endl; 
		if (targetPipeFriction>=0) std::cout << "targetPipeFriction: "<< targetPipeFriction <<"." <<endl; 
		std::cout << "target Y: " << targetYs[0] << "," << targetYs[1] << endl;
		std::cout << "target V: " << targetShearRate << endl;
		std::cout << "target N: " << targetN << endl;
	}
}


void DirectShearBox6D::setBodyEtaRoll(Real er, Body::id_t bodyID)
{
	/*
	Set the rolling resistance coefficient (etaRoll) of one body (or every member of a clump)
	by writing into its CohFrictMat. etaRoll scales the contact radius that caps the rolling
	plastic moment.
	*/
	const shared_ptr<Body>& b   = Body::byId(bodyID,scene);
	if (b->isClump()) {
		FOREACH(const auto& mm, YADE_PTR_CAST<Clump>(b->shape)->members)
		{
			// mm.first is Body::id_t, mm.second is local Se3r of that body in the clump
			const Body::id_t&       memberId = mm.first;
			const shared_ptr<Body>& member   = Body::byId(memberId,scene);
			YADE_PTR_CAST<CohFrictMat>(member->material)->etaRoll = er;
		}
	}
	else YADE_PTR_CAST<CohFrictMat>(b->material)->etaRoll = er;
}

void DirectShearBox6D::setBodyEtaTwist(Real et, Body::id_t bodyID)
{
	/*
	Set the twisting resistance coefficient (etaTwist) of one body (or every member of a clump)
	by writing into its CohFrictMat. etaTwist scales the contact radius that caps the twisting
	plastic moment.
	*/
	const shared_ptr<Body>& b   = Body::byId(bodyID,scene);
	if (b->isClump()) {
		FOREACH(const auto& mm, YADE_PTR_CAST<Clump>(b->shape)->members)
		{
			// mm.first is Body::id_t, mm.second is local Se3r of that body in the clump
			const Body::id_t&       memberId = mm.first;
			const shared_ptr<Body>& member   = Body::byId(memberId,scene);
			YADE_PTR_CAST<CohFrictMat>(member->material)->etaTwist = et;
		}
	}
	else YADE_PTR_CAST<CohFrictMat>(b->material)->etaTwist = et;
}

void DirectShearBox6D::setBodiesEtaRoll(Real er, vector<Body::id_t> bodyIDs){
	/*
	Apply setBodyEtaRoll to every body in the list (convenience loop over setBodyEtaRoll).
	*/
    for (Body::id_t bodyID : bodyIDs) setBodyEtaRoll(er,bodyID);}
void DirectShearBox6D::setBodiesEtaTwist(Real et, vector<Body::id_t> bodyIDs){
	/*
	Apply setBodyEtaTwist to every body in the list (convenience loop over setBodyEtaTwist).
	*/
    for (Body::id_t bodyID : bodyIDs) setBodyEtaTwist(et,bodyID);}

void DirectShearBox6D::setInteractionsEta()
{
	/*
	Recompute every contact's maxRollPl / maxTwistPl from the current etaRoll / etaTwist of
	the two adjacent CohFrictMat materials (min of the two, scaled by the contact radius), and
	enable the moment-rotation law only if both materials request it. Call this AFTER changing
	any eta* so existing contacts pick up the new resistances.
	*/
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	FOREACH(const shared_ptr<Interaction>& ii, *scene->interactions)
	{
		if (!ii->isReal()) continue;
		ScGeom6D*    geom  = YADE_CAST<ScGeom6D*>(ii->geom.get());
		const shared_ptr<CohFrictMat>& sdec1 = YADE_PTR_CAST<CohFrictMat>((*bodies)[(Body::id_t)((ii)->getId1())]->material);
		const shared_ptr<CohFrictMat>& sdec2 = YADE_PTR_CAST<CohFrictMat>((*bodies)[(Body::id_t)((ii)->getId2())]->material);
		//FIXME - why dynamic_cast fails here?
		CohFrictPhys*  contactPhysics = YADE_CAST<CohFrictPhys*>((ii)->phys.get());
		const Real& Da = geom->radius1;
		const Real& Db = geom->radius2;
		contactPhysics->maxRollPl = math::min(sdec1->etaRoll * Da, sdec2->etaRoll * Db);
		contactPhysics->maxTwistPl = math::min(sdec1->etaTwist * Da, sdec2->etaTwist * Db);
		contactPhysics->momentRotationLaw = (sdec1->momentRotationLaw && sdec2->momentRotationLaw);
	}
}

void DirectShearBox6D::setWallEtaRoll(Real er, Body::id_t bodyID)
{
	/*
	Set the rolling resistance of one wall, then refresh all contacts (setInteractionsEta).
	*/
	setBodyEtaRoll(er, bodyID);
	setInteractionsEta();
}

void DirectShearBox6D::setWallEtaTwist(Real et, Body::id_t bodyID)
{
	/*
	Set the twisting resistance of one wall, then refresh all contacts (setInteractionsEta).
	*/
	setBodyEtaTwist(et, bodyID);
	setInteractionsEta();
}

void DirectShearBox6D::setWallsEtaRoll(Real er, vector<Body::id_t> bodyIDs)
{
	/*
	Set the rolling resistance of a list of walls, then refresh all contacts.
	*/
	setBodiesEtaRoll(er, bodyIDs);
	setInteractionsEta();
}

void DirectShearBox6D::setWallsEtaTwist(Real et, vector<Body::id_t> bodyIDs)
{
	/*
	Set the twisting resistance of a list of walls, then refresh all contacts.
	*/
	setBodiesEtaTwist(et, bodyIDs);
	setInteractionsEta();
}

void DirectShearBox6D::setSpheresEtaRollOnly(Real er)
{
	/*
	Set the rolling resistance (etaRoll) of every free dynamic sphere only (does not refresh
	contacts; use setSpheresEtaRoll for that).
	*/
	shared_ptr<BodyContainer>& bodies = Omega::instance().getScene().get()->bodies;
	const size_t nbodies = bodies->size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided)
#endif
	for (size_t i = 0; i < nbodies; ++i) {
		const auto& b = (*bodies)[i];
		if (b && b->isDynamic() && !b->isClumpMember() && b->shape && b->shape->getClassName() == "Sphere")
			YADE_PTR_CAST<CohFrictMat>(b->material)->etaRoll = er;
	}
}

void DirectShearBox6D::setSpheresEtaTwistOnly(Real et)
{
	/*
	Set the twisting resistance (etaTwist) of every free dynamic sphere only (does not refresh
	contacts; use setSpheresEtaTwist for that).
	*/
	shared_ptr<BodyContainer>& bodies = Omega::instance().getScene().get()->bodies;
	const size_t nbodies = bodies->size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided)
#endif
	for (size_t i = 0; i < nbodies; ++i) {
		const auto& b = (*bodies)[i];
		if (b && b->isDynamic() && !b->isClumpMember() && b->shape && b->shape->getClassName() == "Sphere")
			YADE_PTR_CAST<CohFrictMat>(b->material)->etaTwist = et;
	}
}

void DirectShearBox6D::setSpheresEtaRoll(Real et)
{
	/*
	Set the rolling resistance of all spheres, then refresh every contact (setInteractionsEta).
	*/
	setSpheresEtaRollOnly(et);
	setInteractionsEta();
}

void DirectShearBox6D::setSpheresEtaTwist(Real et)
{
	/*
	Set the twisting resistance of all spheres, then refresh every contact (setInteractionsEta).
	*/
	setSpheresEtaTwistOnly(et);
	setInteractionsEta();
}


void DirectShearBox6D::setAllEta()
{
	/*
	Apply the whole rolling/twisting resistance plan at once: for every target*EtaRoll /
	target*EtaTwist attribute that is >= 0, push it into the corresponding wall / particle /
	pipe material, then call setInteractionsEta() once to recompute maxRollPl / maxTwistPl
	for all contacts. This is the main entry used by setExternal() in the 6D variant.
	*/
	bool hasChanges = false;
	if (targetWallsEtaRoll_X0>=0) 
	{
		setBodiesEtaRoll(targetWallsEtaRoll_X0,wall_ids_X0);
		std::cout << "EtaRoll of walls x0 is set as "<< targetWallsEtaRoll_X0 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallsEtaTwist_X0>=0) 
	{
		setBodiesEtaTwist(targetWallsEtaTwist_X0,wall_ids_X0);
		std::cout << "EtaTwist of walls x0 is set as "<< targetWallsEtaTwist_X0 <<"." <<endl;
		hasChanges = true;
	}

	if (targetWallsEtaRoll_X1>=0) 
	{
		setBodiesEtaRoll(targetWallsEtaRoll_X1,wall_ids_X1);
		std::cout << "EtaRoll of walls x1 is set as "<< targetWallsEtaRoll_X1 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallsEtaTwist_X1>=0) 
	{
		setBodiesEtaTwist(targetWallsEtaTwist_X1,wall_ids_X1);
		std::cout << "EtaTwist of walls x1 is set as "<< targetWallsEtaTwist_X1 <<"." <<endl;
		hasChanges = true;
	}

	if (targetWallsEtaRoll_Y0>=0)
	{
		setBodiesEtaRoll(targetWallsEtaRoll_Y0,wall_ids_Y0);
		std::cout << "EtaRoll of walls y0 is set as "<< targetWallsEtaRoll_Y0 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallsEtaTwist_Y0>=0)
	{
		setBodiesEtaTwist(targetWallsEtaTwist_Y0,wall_ids_Y0);
		std::cout << "EtaTwist of walls y0 is set as "<< targetWallsEtaTwist_Y0 <<"." <<endl;
		hasChanges = true;
	}

	if (targetWallsEtaRoll_Y1>=0) 
	{
		setBodiesEtaRoll(targetWallsEtaRoll_Y1,wall_ids_Y1);
		std::cout << "EtaRoll of walls y1 is set as "<< targetWallsEtaRoll_Y1 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallsEtaTwist_Y1>=0) 
	{
		setBodiesEtaTwist(targetWallsEtaTwist_Y1,wall_ids_Y1);
		std::cout << "EtaTwist of walls y1 is set as "<< targetWallsEtaTwist_Y1 <<"." <<endl;
		hasChanges = true;
	}

	if (targetWallEtaRoll_Z0>=0) 
	{
		setBodyEtaRoll(targetWallEtaRoll_Z0,wall_id_Z0);
		std::cout << "EtaRoll of wall z0 is set as "<< targetWallEtaRoll_Z0 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallEtaTwist_Z0>=0) 
	{
		setBodyEtaTwist(targetWallEtaTwist_Z0,wall_id_Z0);
		std::cout << "EtaTwist of wall z0 is set as "<< targetWallEtaTwist_Z0 <<"." <<endl;
		hasChanges = true;
	}

	if (targetWallEtaRoll_Z1>=0)
	{
		setBodyEtaRoll(targetWallEtaRoll_Z1,wall_id_Z1);
		std::cout << "EtaRoll of wall z1 is set as "<< targetWallEtaRoll_Z1 <<"." <<endl;
		hasChanges = true;
	}
	if (targetWallEtaTwist_Z1>=0)
	{
		setBodyEtaTwist(targetWallEtaTwist_Z1,wall_id_Z1);
		std::cout << "EtaTwist of wall z1 is set as "<< targetWallEtaTwist_Z1 <<"." <<endl;
		hasChanges = true;
	}

	if (targetPipeEtaRoll>=0)
	{
		setBodyEtaRoll(targetPipeEtaRoll,shearbox_id);
		std::cout << "EtaRoll of shearbox(pipe) is set as "<< targetPipeEtaRoll <<"." <<endl;
		hasChanges = true;
	}
	if (targetPipeEtaTwist>=0)
	{
		setBodyEtaTwist(targetPipeEtaTwist,shearbox_id);
		std::cout << "EtaTwist of shearbox(pipe) is set as "<< targetPipeEtaTwist <<"." <<endl;
		hasChanges = true;
	}

	if (targetParticlesEtaRoll>=0) 
	{
		setSpheresEtaRollOnly(targetParticlesEtaRoll);
		std::cout << "EtaRoll of particles is set as "<< targetParticlesEtaRoll <<"." <<endl;
		hasChanges = true;
	}
	if (targetParticlesEtaTwist>=0) 
	{
		setSpheresEtaTwistOnly(targetParticlesEtaTwist);
		std::cout << "EtaTwist of particles is set as "<< targetParticlesEtaTwist <<"." <<endl;
		hasChanges = true;
	}

	if (hasChanges) setInteractionsEta();
}

void DirectShearBox6D::setExternal()
{
	/*
	Same as DirectShearBox::setExternal(), plus it also applies the planned rolling/twisting
	resistances (setAllEta) so the 6D moment law is active before the external stress servo runs.
	*/
	if (changePipeDOFornot) setPipeDOFs(targetBlockedDOFs);
	if (gravityAction) setGravityZ(targetGravityZ);
	setAllFriction();
	setAllEta();
	resetRefSe3();
	if ((setColor) and (colorlayerNum>0)) setcolors();
}

void DirectShearBox6D::showState() 
{
	/*
	Print the same phase/target summary as DirectShearBox::showState(), including the current
	6D etaRoll / etaTwist settings, to stdout. Called automatically each step when logornot is true.
	*/
	scene->forces.sync();
	
	if (internalCompaction) {
		std::cout << "Internal compacting."<< endl;
		std::cout << "--d50: "<< currentd50 << "->"<< targetd50 << endl;
		std::cout << "--unb: "<< Shop::unbalancedForce(false, scene) << "->"<< samplingUnbThreshold << endl;
	}
	else {
		if (sampling) 
		{
			std::cout << " Sampling. Target pressures: ";
			for (auto tp = targetPressures.begin(); tp != targetPressures.end();++tp){std::cout << *tp << " ";}
			std::cout << endl;
		}
		else {std::cout << "External compacting. " << endl;}
		std::cout << "--pressor: "<< stress[wallsZ1][2] << "->"<< pressorGoal << endl;
		std::cout << "--unb: "<< Shop::unbalancedForce(false, scene) << endl;
	}
	if (axialloading) {
		State* pipeState = Body::byId(shearbox_id, scene)->state.get();
		std::cout << "Axial loading. " << endl;
		std::cout << "Y: " << targetYs[(int) (not loading)] << ((loading) ? "->" : "<-") << pipeState->displ()[1] <<"->" << targetYs[(int) loading] << endl;
		std::cout << "V: " << pipeState->vel[1] << "->" << ((loading) ? 1.0 : -1.0)*targetShearRate << endl;
		std::cout << "N: " << currentN << "->" << targetN << ((loading) ? ": loading" : ": unloading") << endl;
	}
	else{
		if (changePipeDOFornot) std::cout << "targetBlockedDOFs: " << targetBlockedDOFs << ". " << endl;
		if (gravityAction) std::cout << "targetGravityZ: " << targetGravityZ << "." << endl;
		if (targetParticlesFriction>=0) std::cout << "targetParticlesFriction: "<< targetParticlesFriction <<"." <<endl; 
		if (targetParticlesEtaRoll>=0) std::cout << "targetParticlesEtaRoll: "<< targetParticlesEtaRoll <<"." <<endl; 
		if (targetParticlesEtaTwist>=0) std::cout << "targetParticlesEtaTwist: "<< targetParticlesEtaTwist <<"." <<endl; 

		if (targetWallFriction_Z1>=0) std::cout << "targetWallFriction_Z1: "<< targetWallFriction_Z1 <<"." <<endl;
		if (targetWallEtaRoll_Z1>=0) std::cout << "targetWallEtaRoll_Z1: "<< targetWallEtaRoll_Z1 <<"." <<endl;
		if (targetWallEtaTwist_Z1>=0) std::cout << "targetWallEtaTwist_Z1: "<< targetWallEtaTwist_Z1 <<"." <<endl;

		if (targetWallsFriction_X0>=0) std::cout << "targetWallsFriction_X0: "<< targetWallsFriction_X0 <<"." <<endl; 
		if (targetWallsEtaRoll_X0>=0) std::cout << "targetWallsEtaRoll_X0: "<< targetWallsEtaRoll_X0 <<"." <<endl; 
		if (targetWallsEtaTwist_X0>=0) std::cout << "targetWallsEtaTwist_X0: "<< targetWallsEtaTwist_X0 <<"." <<endl; 

		if (targetWallsFriction_X1>=0) std::cout << "targetWallsFriction_X1: "<< targetWallsFriction_X1 <<"." <<endl; 
		if (targetWallsEtaRoll_X1>=0) std::cout << "targetWallsEtaRoll_X1: "<< targetWallsEtaRoll_X1 <<"." <<endl; 
		if (targetWallsEtaTwist_X1>=0) std::cout << "targetWallsEtaTwist_X1: "<< targetWallsEtaTwist_X1 <<"." <<endl; 

		if (targetWallsFriction_Y0>=0) std::cout << "targetWallsFriction_Y0: "<< targetWallsFriction_Y0 <<"." <<endl; 
		if (targetWallsEtaRoll_Y0>=0) std::cout << "targetWallsEtaRoll_Y0: "<< targetWallsEtaRoll_Y0 <<"." <<endl; 
		if (targetWallsEtaTwist_Y0>=0) std::cout << "targetWallsEtaTwist_Y0: "<< targetWallsEtaTwist_Y0 <<"." <<endl; 

		if (targetWallsFriction_Y1>=0) std::cout << "targetWallsFriction_Y1: "<< targetWallsFriction_Y1 <<"." <<endl; 
		if (targetWallsEtaRoll_Y1>=0) std::cout << "targetWallsEtaRoll_Y1: "<< targetWallsEtaRoll_Y1 <<"." <<endl; 
		if (targetWallsEtaTwist_Y1>=0) std::cout << "targetWallsEtaTwist_Y1: "<< targetWallsEtaTwist_Y1 <<"." <<endl; 

		if (targetWallFriction_Z0>=0) std::cout << "targetWallFriction_Z0: "<< targetWallFriction_Z0 <<"." <<endl; 
		if (targetWallEtaRoll_Z0>=0) std::cout << "targetWallEtaRoll_Z0: "<< targetWallEtaRoll_Z0 <<"." <<endl; 
		if (targetWallEtaTwist_Z0>=0) std::cout << "targetWallEtaTwist_Z0: "<< targetWallEtaTwist_Z0 <<"." <<endl; 

		if (targetPipeFriction>=0) std::cout << "targetPipeFriction: "<< targetPipeFriction <<"." <<endl; 
		if (targetPipeEtaRoll>=0) std::cout << "targetPipeEtaRoll: "<< targetPipeEtaRoll <<"." <<endl; 
		if (targetPipeEtaTwist>=0) std::cout << "targetPipeEtaTwist: "<< targetPipeEtaTwist <<"." <<endl; 
		
		std::cout << "target Y: " << targetYs[0] << "," << targetYs[1] << endl;
		std::cout << "target V: " << targetShearRate << endl;
		std::cout << "target N: " << targetN << endl;
	}
}


} // namespace yade


