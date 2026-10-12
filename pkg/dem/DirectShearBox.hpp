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

#pragma once

#include <lib/base/Math.hpp>
#include <core/Scene.hpp>
#include <pkg/common/BoundaryController.hpp>
#include <pkg/dem/NewtonIntegrator.hpp>
#include <boost/array.hpp>
#include <preprocessing/dem/SpherePack.hpp>
#include <pkg/dem/DirectShearBoxRecorder.hpp>
namespace yade { // Cannot have #include directive inside.

class Scene;
class State;
class SpherePack;

// class DirectShearBoxRecorder;
class NewtonIntegrator;
/*! \brief Controls the stress on the boundaries of a box and compute strain-like and stress-like quantities for the packing. The algorithms used have been developed initialy for simulations reported in [Chareyre2002a] and [Chareyre2005]. They have been ported to Yade in a second step and used in e.g. [Kozicki2008],[Scholtes2009b],[Jerier2010b].
*/

class DirectShearBox : public BoundaryController {

public:
	inline const Vector3r getForce(Scene* rb, Body::id_t id) const
	{
		return rb->forces.getForce(id); /* needs sync, which is done at the beginning of action */
	}

	//! Internal indices for retrieving walls (each wall may be made of several bodies, see wall_ids_*).
	enum { wallsX0 = 0, wallsX1, wallsY0, wallsY1, wallsZ0, wallsZ1 };
	//! Translation of the stress-control at the previous time step (used by the steady-flow servo algorithm).
	Real previousTranslation=0;
	//! Current boundary-packing interaction stiffness, updated every stiffnessUpdateInterval iterations.
	Real stiffness=0;
	//! Resultant stresses on the six boundary walls (always positive, see computeStress()).
	Vector3r stress[6];
	//! Resultant forces on the six boundary walls.
	Vector3r force[6];
	//! Current volume of the sample box (depthY0 * widthX0 * heightZ).
	Real boxVolume;

	// Performance optimization: cached values and pointers
	Real cachedHeightZ;
	Real cachedInvXSurface, cachedInvYSurface, cachedInvZSurface;
	shared_ptr<Body> cachedWallZ1, cachedWallZ0;
	State* cachedStateZ1;
	State* cachedStateZ0;

	virtual ~DirectShearBox() { };

	void action() override;
	// ! Regulate the stress applied on walls with flag wall_XXX_activated = true
	void controlExternalStress(Real resultantForceZ);
	//! Regulate the mean stress by changing spheres size, WARNING : this function assumes that all dynamic bodies in the problem are spheres
	void controlInternalStress(Real multiplier);
	//! update the stiffness of boundary-packing interaction (sum of contacts stiffness on the boundary)
	void updateStiffness();
	//! Compute stresses on walls as "Vector3r stress[6]" and compute meanStress
	void computeStress();
	//! Compute the mean/max unbalanced force in the assembly (normalized by mean contact force)
	Real ComputeUnbalancedForce(bool maxUnbalanced = false);
	///! Getter for stress and rates in python
	Vector3r getStress(int boundId) const;
	virtual void growParticles(Real multiplier, bool updateMass, bool dynamicOnly);


	// ---- shear / controller helpers ----
	shared_ptr<DirectShearBoxRecorder> dsbr;   //!< Recorder engine inserted by setAxial() to log the run.
	shared_ptr<NewtonIntegrator> newton;       //!< Shortcut to the scene's NewtonIntegrator (used to set damping).
	SpherePack* sp = new SpherePack();         //!< Used to dump particle geometry to a text file after internal compaction.
	//! Periodic shear controller: advances currentTargetY and switches loading direction at stroke ends.
	void controllerAction();
	//! Drive the shear box at +/- targetShearRate along Y; also force-ends compaction/sampling when called.
	void axialAction();
	//! Set the material friction angle [rad] of a single body (wall, pipe, ...).
	void setBodyFriction(Real angleRad, Body::id_t bodyID);
	//! Set the material friction angle [rad] of a list of bodies.
	void setBodiesFriction(Real angleRad, vector<Body::id_t> bodyIDs);
	//! Refresh contact-level friction on every interaction using the current materials (call after changing friction).
	void updateInteractionsFriction();
	//! Set material friction angle [rad] of one wall body.
	void setWallFriction(Real angleRad, Body::id_t bodyID);
	//! Set material friction angle [rad] of a list of wall bodies.
	void setWallsFriction(Real angleRad, vector<Body::id_t> bodyIDs);
	//! Set material friction angle [rad] of all spheres only (not clumps / walls).
	void setSpheresFrictionOnly(Real angleRad);
	//! Set material friction angle [rad] of all spheres and clump members.
	void setSpheresFriction(Real angleRad);
	//! Convenience: apply all target* frictions (walls / pipe / particles) in one call.
	void setAllFriction();
	//! Reset the reference pose (se3) of all bodies, so displacement output starts from zero at shear begin.
	void resetRefSe3();
	//! Lock the blocked degrees of freedom of the shear box / pipe (string like "xz" or "").
	void setPipeDOFs(std::string DOFs);
	//! Set the target cyclic number targetN (0 = monotonic loading).
	void setCycle (int N);
	//! Set the loading direction flag (true = loading toward targetYs[1]).
	void setLoading (bool L);
	//! Set the shear target positions targetYs (targetYs[0] < targetYs[1]).
	void setY (vector<Real> Ys);
	//! Set the absolute target shear rate targetShearRate.
	void setV (Real V);
	//! Override the current shear velocity (sets currentTargetY and the pipe Y-velocity directly).
	void setCurrentVY (Real V);
	//! Read back the current shear-box Y velocity.
	Real getCurrentVY ();
	//! Apply gravity [0,0,gravityZ] to all dynamic bodies.
	void setGravityZ (Real gravityZ);
	//! Set the material friction angle [rad] of the interactions touching one body.
	void setInteractionsFrictionOnBody(Real angleRad, Body::id_t bodyID);
	//! Set the material friction angle [rad] of the interactions touching a list of bodies.
	void setInteractionsFrictionOnBodies(Real angleRad, vector<Body::id_t> bodyIDs);
	//! Set the material friction angle [rad] of all wall-wall / wall-particle interactions.
	void setAllInteractionsFrictionOnWalls();
	//! Switch the engine into shear mode: push targets, apply frictions, insert the recorder.
	void setAxial();
	//! Apply all boundary conditions before shearing (gravity, frictions, ref pose, colors).
	virtual void setExternal();
	//! Print a one-line human-readable state summary (used when logornot is true).
	virtual void showState();
	//! Set the global damping of the NewtonIntegrator.
	void setDamping (Real D);
	//! Set the material friction angle [rad] of the interactions touching a clump body (e.g. the pipe).
	void setInteractionsFrictionOnClump(Real angleRad, Body::id_t bodyID);
	//! Assign layer-based render colors (colorOne / colorTwo over colorlayerNum layers) to particles.
	void setcolors();

	// ---- end of controller helpers ----
	// clang-format off
 		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(DirectShearBox,BoundaryController,
		"Engine to run a direct-shear (or similar) test: interface shearing, axial soil-pipe interaction, etc.\n"
		"Two compaction modes are available, switched by :yref:`internalCompaction<DirectShearBox::internalCompaction>`:\n"
		"  - INTERNAL compaction (internalCompaction==1): particle diameters grow from :yref:`initiald50<DirectShearBox::initiald50>` toward :yref:`targetd50<DirectShearBox::targetd50>` via a staged schedule controlled by the critical-size factors (r1,r2,r3) and enlarging factors (f0..f3). Particle size, mass and inertia are scaled accordingly, and the contact geometry (refR1/refR2) is updated. This increases the solid fraction without moving the walls.\n"
		"  - EXTERNAL compaction (internalCompaction==0): the pressor (wall_id_Z1) is servo-controlled so that its normal stress tends to :yref:`pressorGoal<DirectShearBox::pressorGoal>`. During sampling, pressorGoal is stepped through :yref:`targetPressures<DirectShearBox::targetPressures>` once the pressure / unbalanced-force / velocity thresholds are met. The shear sequence can then be driven by axial loading or cyclic (Y,V,N) control.\n"
		"Typical usage: (1) construct the box walls and a shear body, (2) set internalCompaction + targetd50 (or pressorGoal + targetPressures), (3) call setExternal() to apply frictions/gravity/DOFs, (4) let action() run until the target state is reached, (5) call setAxial() to start shearing, (6) read results from DirectShearBoxRecorder."
		,
		((unsigned int,stiffnessUpdateInterval,5,,"Iteration period for measuring the resultant packing-boundary stiffness, used by the stress servo-control."))
   		((unsigned int,radiusControlInterval,5,,"Iteration period for the internal-radius growth step during internal compaction."))
		((unsigned int,computeStressInterval,5,,"Iteration period for recomputing boundary stresses/forces via computeStress()."))
		((unsigned int,porosityHeightZUpdateInterval,100,,"Iteration period for recomputing porosity and the free height heightZ (only when pressorGoal != 0)."))
		((Real,stressDamping,0.25,,"Steady-flow servo coefficient. Wall damping for the stress controller: 0 implies a (theoretical) perfect control, 1 means no boundary movement."))
		((Real,steadyflowFactor,0.8,,"Steady-flow servo coefficient. previousTranslation = (1 - stressDamping) * translation * normal[4] + steadyflowFactor * previousTranslation."))
		((Real,thicknessZ,0,,"Thickness of the boundary walls in the 2-axis (normal) direction; subtracted when computing the free height heightZ."))
		((vector<Body::id_t>,wall_ids_Y0,vector<Body::id_t>({4,5}),,"Body ids of the 1- boundary. Set to [] if this boundary does not exist."))
		((vector<Body::id_t>,wall_ids_Y1,vector<Body::id_t>({6,7}),,"Body ids of the 1+ boundary. Set to [] if this boundary does not exist."))
		((vector<Body::id_t>,wall_ids_X0,vector<Body::id_t>({0,1}),,"Body ids of the 0- boundary. Set to [] if this boundary does not exist."))
		((vector<Body::id_t>,wall_ids_X1,vector<Body::id_t>({2,3}),,"Body ids of the 0+ boundary. Set to [] if this boundary does not exist."))
		((int,wall_id_Z1,8,,"Body id of the pressor (2+ boundary). Default is valid if aabbWalls are appended BEFORE spheres."))
		((int,wall_id_Z0,9,,"Body id of the 2- boundary. Default is valid if aabbWalls are appended BEFORE spheres."))
		((int,shearbox_id,10,,"Body id of the shear box / pipe; MUST be given."))
		((Real,initiald50,1,,"Initial mean particle diameter d50; MUST be given. Reference for internal-compaction scaling."))
		((Real,currentd50,,Attr::readonly,"Current mean particle diameter d50 (updated during internal compaction). Read only |yupdate|"))
		((Real,targetd50,0,,"Target mean particle diameter d50. MUST be given for internal compaction. If targetd50<=0 or targetd50<currentd50, internal compaction does nothing."))
		((Real,pipeVolume,0,,"Volume of the pipe (if any); used together with particlesVolume to compute porosity."))
		((Real,depthY0,1,,"Box size along the 1-axis (shear direction reference); MUST be given."))
		((Real,widthX0,1,,"Box size along the 0-axis; MUST be given."))
		((Real,heightZ,,Attr::readonly,"Box size along the 2-axis (free height = wallZ1 - wallZ0 - thicknessZ). |yupdate|"))
		((Real,heightZ0,,,"Reference free height (2-axis) recorded at the start of shearing, used for vertical-displacement output."))
		((Real,pressorGoal,17000,,"Prescribed normal stress on the pressor (2-axis), active during external compaction. Set from targetPressures while sampling; equals targetPressures[0] if sampling is never used."))
		((Real,max_vel,0.5,,"Maximum allowed boundary velocity [m/s]. Overrides the servo-assigned velocity when the latter is higher. max_vel can be infinity in many cases, but sometimes helps stabilize packings. Per-axis maxima are derived from the sample dimensions so that moving every boundary at its max yields an isotropic strain rate (see e.g. :yref:`TriaxialStressController::max_vel1`)."))
		((Real,previousStress,0,Attr::readonly,"Previous boundary stress, used by the servo. |yupdate|"))
		((Real,previousMultiplier,1,Attr::readonly,"Previous growth multiplier, used by the servo. |yupdate|"))
		((bool,internalCompaction,true,,"Switch the compaction mode: true = INTERNAL (grow particles), false = EXTERNAL (move walls)."))
		((Real,meanStress,0,Attr::readonly,"Mean stress over the active boundaries (see computeStress()). |yupdate|"))
		((Real,particlesVolume,0,Attr::readonly,"Solid volume of all dynamic non-clump spheres (and clumps). |yupdate|"))
		((Real,porosity,,Attr::readonly,"Porosity of the packing, computed from particlesVolume, pipeVolume and boxVolume. |yupdate|"))
		((Real,externalWork,0,Attr::readonly,"Mechanical work of the boundary conditions: $\\int_{\\partial \\Omega} \\mat{T} \\cdot \\mat{u}\\, ds$ with $\\mat{T}$ the surface traction and $\\mat{u}$ the boundary displacement. |yupdate|"))
		((bool,logornot,false,,"Whether to print the state log via showState() during the run."))
		((string,savePath,"log/",,"Directory for saved simulation files. You may add a descriptive sub-path; please verify it exists before running."))
		((Real,internalUnbThreshold,0.05,,"Unbalanced-force threshold for internal compaction: once reached, particles are enlarged again."))
		((bool,first,true,,"Internal flag: true before the first internal-compaction growth step."))
		((Real,r1,0.7,,"Internal-compaction stage-1 critical-size factor (growth active when currentd50 > r1 * targetd50)."))
		((Real,r2,0.9,,"Internal-compaction stage-2 critical-size factor."))
		((Real,r3,0.95,,"Internal-compaction stage-3 critical-size factor."))
		((Real,f0,0.0002,,"Internal-compaction enlarging factor for stage 1 (largest growth)."))
		((Real,f1,0.0001,,"Internal-compaction enlarging factor for stage 2."))
		((Real,f2,0.00005,,"Internal-compaction enlarging factor for stage 3."))
		((Real,f3,0.000005,,"Internal-compaction enlarging factor for stage 4 (finest adjustment)."))
		((bool,internalCompactionSaved,false,,"Whether the internally-compacted state has been saved."))
		((unsigned int,controllerInterval,10,,"Iteration period for the axial-loading controller."))
		((unsigned int,pipefrictionUpdateNum,10,,"Update interval (controllerInterval * pipefrictionUpdateNum) for re-setting pipe friction during shearing; <=0 means never update."))
		((bool,axialloading,false,,"Whether axial (shear) loading control is active."))
		((bool,gravityAction,false,,"Whether to apply gravity. If true, gravity is set to [0,0,targetGravityZ]."))
		((Real,targetGravityZ,-9.81,,"Target gravity along Z applied when gravityAction is true."))
		((Real,targetWallsFriction_X0,-1,,"Target friction angle [rad] of the X0 walls. Negative means no change."))
		((Real,targetWallsFriction_X1,-1,,"Target friction angle [rad] of the X1 walls. Negative means no change."))
		((Real,targetWallsFriction_Y0,-1,,"Target friction angle [rad] of the Y0 walls. Negative means no change."))
		((Real,targetWallsFriction_Y1,-1,,"Target friction angle [rad] of the Y1 walls. Negative means no change."))
		((Real,targetWallFriction_Z0,-1,,"Target friction angle [rad] of the Z0 wall. Negative means no change."))
		((Real,targetWallFriction_Z1,-1,,"Target friction angle [rad] of the Z1 (pressor) wall. Negative means no change."))
		((Real,targetPipeFriction,-1,,"Target friction angle [rad] of the pipe. Negative means no change."))
		((Real,targetParticlesFriction,-1,,"Target friction angle [rad] of particles. Negative means no change."))
		((bool,changePipeDOFornot,false,,"Whether to change the shear box (pipe) blocked DOFs. If true, DOFs are changed once external compaction begins, then this is auto-set to false."))
		((std::string,targetBlockedDOFs,"",,"Blocked DOFs string for the shear box, applied once external compaction begins."))
		((int,saverNum,40,,"Number of yade-file saves per loading; <=0 means no periodic save."))
		((Real,currentSaverNum,0.0,,"Current count of saves performed."))
		((Real,currentTargetY,,,"Shear box's current target position in Y; used for cyclic loading."))
		((vector<Real>,targetYs,vector<Real>({0.0,0.01}),,"Shear box target positions in Y for cyclic loading. targetYs[0] (unloading) must be smaller than targetYs[1] (loading)."))
		((bool,loading,true,,"Loading (true) or unloading (false). targetYs[1] is the loading target, targetYs[0] the unloading target."))
		((int,currentN,0,,"Current cyclic number reached."))
		((int,targetN,0,,"Target cyclic number. 0 means monotonic (single-pass) loading."))
		((Real,targetShearRate,0.0005,,"Absolute target shear velocity in Y (must be positive). Sign is applied from loading: +targetShearRate when loading, -targetShearRate otherwise."))
		((bool,sampling,true,,"Whether to step through targetPressures (sampling). If true, pressorGoal is advanced through the list once thresholds are met."))
		((vector<Real>,targetPressures,vector<Real>({17000,34000,50000,100000}),,"Target pressures for sampling. If sampling is never used, pressorGoal is set to targetPressures[0]."))
		((Real,samplingPressureThreshold,0.01,,"Relative pressure threshold for sampling: once met, pressorGoal advances to the next targetPressures item."))
		((Real,samplingUnbThreshold,0.05,,"Unbalanced-force threshold for external-loading sampling: once met, pressorGoal advances to the next item."))
		((Real,samplingVelThreshold,0.001,,"Velocity threshold for sampling: once met, pressorGoal advances to the next item."))
		((int,samplingNextIndex,0,,"Index of the next pressure to apply in targetPressures."))
		((bool,autoSimulation,false,,"Whether the engine drives the whole simulation automatically (compaction -> shear) without external script stepping."))
		((Real,targetDamping,-1,,"Target damping ratio for the simulation. Negative means no change."))
		((bool,setColor,false,,"Whether to color particles in layers after internal compaction."))
		((Vector3r,colorOne,Vector3r(0.75,0.75,0.75),,"Render color for odd-numbered layers (used by setcolors())."))
		((Vector3r,colorTwo,Vector3r(0.3,0.3,0.3),,"Render color for even-numbered layers (used by setcolors())."))
		((int,colorlayerNum,6,,"Number of color layers used by setcolors()."))
		,
		/* extra initializers */
		,
   		/* constructor */
		// initRun=true;
		,
		.def("stress",&DirectShearBox::getStress,(boost::python::arg("id")),
			"Return the current normal stress (Vector3r) on the boundary wall identified by 'id'. "
			"Wall index map: wallsX0 = 0, wallsX1 = 1, wallsY0 = 2, wallsY1 = 3, wallsZ0 = 4, wallsZ1 = 5. "
			"The components that are normally needed are stress[0][0] (X0), [1][0] (X1), [2][1] (Y0), "
			"[3][1] (Y1), [4][2] (Z0) and [5][2] (Z1); these are always reported as positive (compressive).")
		.def("axialAction", &DirectShearBox::axialAction,
			"Begin the axial (shear) loading stage. Picks the target Y position (targetYs[1] while loading, "
			"targetYs[0] while unloading), sets the shear-box (pipe) velocity to +/- targetShearRate along Y, "
			"and forces internalCompaction / sampling to end since they no longer apply during shearing.")
		.def("_setWallFriction", &DirectShearBox::setWallFriction,(py::arg("angleRad"),py::arg("bodyID")),
			"Set the material friction angle (in radians) of the single wall/body with the given bodyID, "
			"and immediately refresh the friction of every contact touching that body.")
		.def("_setWallsFriction", &DirectShearBox::setWallsFriction,(py::arg("angleRad"),py::arg("bodyIDs")),
			"Set the material friction angle (in radians) of every body in the bodyIDs list, and refresh "
			"the contacts of all of them in one go.")
		.def("_setSpheresFriction", &DirectShearBox::setSpheresFriction,(py::arg("angleRad")),
			"Set the material friction angle (in radians) of ALL dynamic, non-clump-member spheres, then "
			"refresh every contact so existing interactions pick up the new friction.")
		.def("_setAllFriction", &DirectShearBox::setAllFriction,
			"Apply the whole friction plan at once: set wall, particle and pipe friction angles from their "
			"respective target* attributes (only those with target >= 0 are changed) and refresh all contacts.")
		.def("_resetRefSe3", &DirectShearBox::resetRefSe3,
			"Reset the reference Se3 (position+orientation) of every body to its current value, so that "
			"subsequent displacement/rotation measurements start from this instant.")
		.def("_setPipeDOFs", &DirectShearBox::setPipeDOFs,(py::arg("DOFs")),
			"Constrain the shear-box (pipe) wall by setting its blockedDOFs string (e.g. 'z' to lock the Z axis). "
			"Used to let the pipe move only along the shear direction.")
		.def("_setCycle", &DirectShearBox::setCycle,(py::arg("N")),
			"Set the number of loading/unloading cycles to perform: assigns targetN = N.")
		.def("_setLoading", &DirectShearBox::setLoading,(py::arg("L")),
			"Set the loading flag. L = True drives the shear box toward targetYs[1] at +targetShearRate; "
			"L = False drives it back toward targetYs[0] at -targetShearRate (unloading).")
		.def("_setY", &DirectShearBox::setY,(py::arg("Ys")),
			"Set the two target Y positions Ys = [unloadedY, loadedY] that bound the shear stroke.")
		.def("_setV", &DirectShearBox::setV,(py::arg("V")),
			"Set the shear-box velocity magnitude targetShearRate = V (m/s) used during axial loading.")
		.def("_setGravityZ", &DirectShearBox::setGravityZ,(py::arg("gravityZ")=-9.81),
			"Set the vertical gravity acceleration as [0, 0, gravityZ]. Defaults to -9.81 m/s^2. "
			"Only takes effect if gravityAction is enabled in setExternal().")
		.def("setExternal", &DirectShearBox::setExternal,
			"Prepare the box for the EXTERNAL (wall-stress servo) stage. In one call it applies the planned "
			"pipe DOFs, gravity, wall/particle/pipe friction, (for 6D) rolling/twisting resistances, resets the "
			"reference Se3, and optionally recolours the sample layers. After this, the stress servo drives the "
			"walls toward pressorGoal.")
		.def("setAxial", &DirectShearBox::setAxial,
			"Prepare the box for the AXIAL (shear) loading stage. In one call it sets the cycle count, the loading "
			"flag, the target Y stroke and the target shear velocity from the plan, then starts moving the pipe. "
			"Typical Python usage: engine.setAxial(loading=True, currentN=0).")
		.def("showState", &DirectShearBox::showState,
			"Print a human-readable summary of the current simulation state to stdout: whether internal compaction, "
			"sampling or axial loading is active, the relevant targets (Y, V, N, friction, gravity, DOFs) and the "
			"current stress / unbalanced force.")
		.def("setcolors", &DirectShearBox::setcolors,
			"Assign a colour to each particle according to the layer it belongs to (based on its initial Z position "
			"and colorlayerNum), for visualisation. Only runs if setColor is enabled and colorlayerNum > 0.")
		.add_property("currentVY",&DirectShearBox::getCurrentVY,&DirectShearBox::setCurrentVY,
			"Get/set the instantaneous shear velocity along Y (m/s) for the CURRENT loading/unloading half-cycle. "
			"For a permanent change of the target rate, use _setV() instead.")
		)
	// clang-format on
	// DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(DirectShearBox);

class DirectShearBox6D : public DirectShearBox {
public:
	virtual ~DirectShearBox6D() { };
	void growParticles(Real multiplier, bool updateMass, bool dynamicOnly) override;
	void setBodyEtaRoll(Real er, Body::id_t bodyID);
	void setBodyEtaTwist(Real et, Body::id_t bodyID);
	void setBodiesEtaRoll(Real er, vector<Body::id_t> bodyIDs);
	void setBodiesEtaTwist(Real et, vector<Body::id_t> bodyIDs);
	void setInteractionsEta();
	void setWallEtaRoll(Real er, Body::id_t bodyID);
	void setWallEtaTwist(Real et, Body::id_t bodyID);
	void setWallsEtaRoll(Real er, vector<Body::id_t> bodyIDs);
	void setWallsEtaTwist(Real et, vector<Body::id_t> bodyIDs);
	void setSpheresEtaRollOnly(Real er);
	void setSpheresEtaTwistOnly(Real et);
	void setSpheresEtaRoll(Real er);
	void setSpheresEtaTwist(Real et);
	void setAllEta();
	void setExternal() override;
	void showState() override;


	// clang-format off
 		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(DirectShearBox6D,DirectShearBox,
		"6D extension of :yref:`DirectShearBox`. In addition to the same compaction/shear workflow, it controls rolling (etaRoll) and twisting (etaTwist) resistances on walls, particles and pipe so that rotational degrees of freedom can be damped.\n"
		"Compaction modes are identical to DirectShearBox: if internalCompaction==1 the currentd50 grows toward targetd50 and then internalCompaction auto-resets to 0; if internalCompaction==0 the pressor stress servo-controls toward pressorGoal (stepped via targetPressures during sampling)."
		,
		/* attrs */
		((Real,targetWallsEtaRoll_X0,-1,,"target EtaRoll of walls of x0. negative means no need to change."))
		((Real,targetWallsEtaRoll_X1,-1,,"target EtaRoll of walls of x1. negative means no need to change."))
		((Real,targetWallsEtaRoll_Y0,-1,,"target EtaRoll of walls of y0. negative means no need to change."))
		((Real,targetWallsEtaRoll_Y1,-1,,"target EtaRoll of walls of y1. negative means no need to change."))
		((Real,targetWallEtaRoll_Z0,-1,,"target EtaRoll of wall of z0. negative means no need to change."))
		((Real,targetWallEtaRoll_Z1,-1,,"target EtaRoll of wall of z1. negative means no need to change."))
		((Real,targetPipeEtaRoll,-1,,"target EtaRoll of pipe. negative means no need to change."))
		((Real,targetParticlesEtaRoll,-1,,"target EtaRoll of particles. negative means no need to change."))
		((Real,targetWallsEtaTwist_X0,-1,,"target EtaTwist of walls of x0. negative means no need to change."))
		((Real,targetWallsEtaTwist_X1,-1,,"target EtaTwist of walls of x1. negative means no need to change."))
		((Real,targetWallsEtaTwist_Y0,-1,,"target EtaTwist of walls of y0. negative means no need to change."))
		((Real,targetWallsEtaTwist_Y1,-1,,"target EtaTwist of walls of y1. negative means no need to change."))
		((Real,targetWallEtaTwist_Z0,-1,,"target EtaTwist of wall of z0. negative means no need to change."))
		((Real,targetWallEtaTwist_Z1,-1,,"target EtaTwist of wall of z1. negative means no need to change."))
		((Real,targetPipeEtaTwist,-1,,"target EtaTwist of pipe. negative means no need to change."))
		((Real,targetParticlesEtaTwist,-1,,"target EtaTwist of particles. negative means no need to change."))
		,
		/* extra initializers */
		,
   		/* constructor */
		// initRun=true;
		,
		/* py */
		.def("_setWallEtaRoll", &DirectShearBox6D::setWallEtaRoll,(py::arg("er"),py::arg("bodyID")),
			"Set the rolling resistance coefficient (etaRoll) of the single body with bodyID to er, then refresh "
			"every contact so the new maxRollPl is propagated. etaRoll multiplies the contact radius to cap the "
			"rolling plastic moment.")
		.def("_setWallsEtaRoll", &DirectShearBox6D::setWallsEtaRoll,(py::arg("er"),py::arg("bodyIDs")),
			"Set the rolling resistance (etaRoll) of all bodies in the bodyIDs list and refresh their contacts.")
		.def("_setSpheresEtaRoll", &DirectShearBox6D::setSpheresEtaRoll,(py::arg("er")),
			"Set the rolling resistance (etaRoll) of all dynamic, non-clump-member spheres and refresh every contact.")
		.def("_setWallEtaTwist", &DirectShearBox6D::setWallEtaTwist,(py::arg("et"),py::arg("bodyID")),
			"Set the twisting resistance coefficient (etaTwist) of the single body with bodyID to et, then refresh "
			"every contact so the new maxTwistPl is propagated.")
		.def("_setWallsEtaTwist", &DirectShearBox6D::setWallsEtaTwist,(py::arg("et"),py::arg("bodyIDs")),
			"Set the twisting resistance (etaTwist) of all bodies in the bodyIDs list and refresh their contacts.")
		.def("_setSpheresEtaTwist", &DirectShearBox6D::setSpheresEtaTwist,(py::arg("et")),
			"Set the twisting resistance (etaTwist) of all dynamic, non-clump-member spheres and refresh every contact.")
		.def("setAllEta", &DirectShearBox6D::setAllEta,
			"Apply the whole rolling/twisting resistance plan at once: for every target*EtaRoll / target*EtaTwist "
			"attribute that is >= 0, push it into the corresponding wall/particle/pipe material, then call "
			"setInteractionsEta() once to recompute maxRollPl / maxTwistPl for all contacts.")
		.def("setExternal", &DirectShearBox6D::setExternal,
			"Same as DirectShearBox.setExternal(), plus it also applies the planned rolling/twisting resistances "
			"(setAllEta) so the 6D moment law is active before the external stress servo runs.")
		.def("showState", &DirectShearBox6D::showState,
			"Print a human-readable summary of the current simulation state (including the 6D etaRoll/etaTwist "
			"settings) to stdout.")
		)
	// clang-format on
	// DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(DirectShearBox6D);

} // namespace yade
