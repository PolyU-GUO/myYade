/*************************************************************************
*  Copyright (C) 2006 by luc Scholtes                                    *
*  luc.scholtes@hmg.inpg.fr                                              *
*  Original file: Law2_ScGeom_CapillaryPhys_Capillarity, from which this *
*  file is derived.                                                      *
*                                                                        *
*  Copyright (C) 2026 by Chang GUO                                       *
*  chang.guo@connect.polyu.hk                                            *
*  Modifications by Chang GUO (2026): the 6D / cohesive-frictional       *
*  variant of the law above (CapillaryCohPhys based); see the class doc. *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Law2_ScGeom6D_CapillaryCohPhys_Capillarity.hpp"
#include <pkg/common/ElastMat.hpp>
#include <pkg/dem/ScGeom.hpp>

#include <lib/base/Math.hpp>
#include <lib/high-precision/Constants.hpp>
#include <lib/high-precision/MathApprox.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/CapillaryPhys.hpp>
#include <pkg/dem/HertzMindlin.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

YADE_PLUGIN((Law2_ScGeom6D_CapillaryCohPhys_Capillarity));
CREATE_LOGGER(Law2_ScGeom6D_CapillaryCohPhys_Capillarity);

void Law2_ScGeom6D_CapillaryCohPhys_Capillarity::action()
{
	if (!scene) cerr << "scene not defined!";
	if (!capillary) postLoad(*this); //when the script does not define arguments, postLoad is never called
	shared_ptr<BodyContainer>& bodies      = scene->bodies;
	int                        sphereIndex = Sphere::getClassIndexStatic();

	// The 6D variant always works with CapillaryCohPhys (there is no CapillaryPhys /
	// MindlinCapillaryPhys switch as in the base class), so no contact model has to be detected
	// here: the loop below only invalidates the meniscii list, which must be re-initialized after
	// the first real contact has been created.
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (I->isReal()) {
			bodiesMenisciiList6D.initialized = false; //must be re-initialized after creation of first real contact in the model
			break;
		}
	}

	if (fusionDetection && !bodiesMenisciiList6D.initialized) bodiesMenisciiList6D.prepare(scene);//TODO

	FOREACH(const shared_ptr<Interaction>& interaction, *scene->interactions)
	{ // could be done in parallel as soon as OpenMPVector class (lib/base/openmp-accu.hpp) is extended See http://www.mail-archive.com/yade-dev@lists.launchpad.net/msg10842.html and msg11238.html
		/// interaction is real
		if (interaction->isReal()) {
			CapillaryCohPhys*        ContactPhysics = static_cast<CapillaryCohPhys*>(interaction->phys.get()); //the 6D variant always uses CapillaryCohPhys (cohesive-frictional model)

			unsigned int id1 = interaction->getId1();
			unsigned int id2 = interaction->getId2();

			/// definition of interacting objects (not necessarily in contact)
			ScGeom6D* currentContactGeometry = static_cast<ScGeom6D*>(interaction->geom.get());

			/// test of interacting bodies geometries since capillarity will be computed between spheres only
			int geometryIndex1 = (*bodies)[id1]->shape->getClassIndex();
			int geometryIndex2 = (*bodies)[id2]->shape->getClassIndex();
			if (!(geometryIndex1 == sphereIndex && geometryIndex2 == sphereIndex)) { // such interactions won't have a meniscus
				// thus we will ask for the interaction to be erased in this distant case, as we do w. sphere-sphere interactions below:
				if (currentContactGeometry->penetrationDepth < 0) scene->interactions->requestErase(interaction);
				continue;
			}
			/// Interacting Grains:
			// If you want to define a ratio between YADE sphere size and real sphere size
			Real alpha = 1;
			Real R1    = 0;
			R1         = alpha * std::min(currentContactGeometry->radius2, currentContactGeometry->radius1);
			Real R2    = 0;
			R2         = alpha * std::max(currentContactGeometry->radius2, currentContactGeometry->radius1);

			/// intergranular distance
			Real D = -alpha * currentContactGeometry->penetrationDepth;
			if (D <= 0 || createDistantMeniscii) { //||(scene->iter < 1) ) // a simplified way to define meniscii everywhere
				D = max(0.,
				        D); // defines fCap when spheres interpenetrate. D<0 leads to wrong interpolation as D<0 has no solution in the interpolation : this is not physically interpretable!! even if, interpenetration << grain radius.
				if (fusionDetection && !ContactPhysics->meniscus) bodiesMenisciiList6D.insert(interaction);
				ContactPhysics->meniscus = true;
			}
			Real Dinterpol = D / R2;

			/// Suction (Capillary pressure):
			Real Pinterpol = 0;
			Pinterpol = ContactPhysics->isBroken ? 0 : capillaryPressure * (R2 / surfaceTension);
			ContactPhysics->capillaryPressure = capillaryPressure;

			/// Capillary solution finder:
			if ((Pinterpol >= 0) && (ContactPhysics->meniscus)) {
				int* currentIndexes = ContactPhysics->currentIndexes;
				//If P=0, we use null solution
				MeniscusParameters solution(
				        Pinterpol ? capillary->interpolate(R1, R2, Dinterpol, Pinterpol, currentIndexes) : MeniscusParameters());
				/// If no bridge, delete the interaction if necessary and escape
				if (!solution.V) {
					if ((fusionDetection) || (ContactPhysics->isBroken))
						bodiesMenisciiList6D.remove(interaction);
					if (D > 0) {
						scene->interactions->requestErase(interaction);
						continue;
					} else if ((Pinterpol > 0) && (showError)) {
						LOG_ERROR(
						        "No meniscus found at a contact. capillaryPressure may be too large wrt. the loaded data files."); // V=0 at a contact reveals a problem if and only if uc* > 0
						showError = false;
					} //show error message once / avoid console spam
				}
				/// capillary adhesion force
				Real     Finterpol = solution.F;
				Vector3r fCap      = -Finterpol * (2 * Mathr::PI * (R2 / alpha) * surfaceTension) * currentContactGeometry->normal;
				ContactPhysics->fCap = fCap;
				/// meniscus volume
				Real Vinterpol = solution.V * pow(R2 / alpha, 3);
				ContactPhysics->vMeniscus = Vinterpol;
				if (Vinterpol != 0) ContactPhysics->meniscus = true;
				else
					ContactPhysics->meniscus = false;

				/// wetting angles
				ContactPhysics->Delta1 = min(
						max(solution.delta1, solution.delta2),
						90.0); // undesired values greater than 90 degrees (by few degrees) may be present in the capillary files for high r (~ 10) and very low suction and distance
				ContactPhysics->Delta2 = min(solution.delta1, solution.delta2);
				// nn11 and nn33
				ContactPhysics->nn11 = pow(R2 / alpha, 2) * solution.nn11;
				ContactPhysics->nn33 = pow(R2 / alpha, 2) * solution.nn33;
			}
			///interaction is not real	//If the interaction is not real, it should not be in the list
		} else if (fusionDetection)
			bodiesMenisciiList6D.remove(interaction);
	}
	if (fusionDetection) checkFusion();

#ifdef YADE_OPENMP
	const long size = scene->interactions->size();
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
	for (long i = 0; i < size; i++) {
		const shared_ptr<Interaction>& interaction = (*scene->interactions)[i];
#else
	FOREACH(const shared_ptr<Interaction>& interaction, *scene->interactions)
	{
#endif
		if (interaction->isReal()) {
			CapillaryCohPhys*        ContactPhysics = static_cast<CapillaryCohPhys*>(interaction->phys.get()); //the 6D variant always uses CapillaryCohPhys (cohesive-frictional model)

			if (ContactPhysics->meniscus) {
				if (fusionDetection) { //version with effect of fusion
					//BINARY VERSION : if fusionNumber!=0 then no capillary force
					short int& fusionNumber = ContactPhysics->fusionNumber;
					if (binaryFusion) {
						if (fusionNumber != 0) { //cerr << "fusion" << endl;
							ContactPhysics->fCap = Vector3r::Zero();
							continue;
						}
					}
					//LINEAR VERSION : capillary force is divided by (fusionNumber + 1) - NOTE : any decreasing function of fusionNumber can be considered in fact
					else if (fusionNumber != 0)
						ContactPhysics->fCap /= (fusionNumber + 1.);
				}
				scene->forces.addForce(interaction->getId1(), -(ContactPhysics->fCap));
				scene->forces.addForce(interaction->getId2(), ContactPhysics->fCap);
			}
		}
	}
}

void Law2_ScGeom6D_CapillaryCohPhys_Capillarity::checkFusion()
{
	//Reset fusion numbers
	FOREACH(const shared_ptr<Interaction>& interaction, *scene->interactions)
	{ // same remark for parallel loops, the most problematic part ?
		if (interaction->isReal()) static_cast<CapillaryCohPhys*>(interaction->phys.get())->fusionNumber = 0;
	}

	std::list<shared_ptr<Interaction>>::iterator firstMeniscus, lastMeniscus, currentMeniscus;
	Real                                         angle1 = -1.0;
	Real                                         angle2 = -1.0;

	for (int i = 0; i < bodiesMenisciiList6D.size(); ++i) { // i is the index (or id) of the body being tested
		CapillaryCohPhys*        cundallInteractionPhysics1 = NULL;
		CapillaryCohPhys*        cundallInteractionPhysics2 = NULL;
		if (!bodiesMenisciiList6D[i].empty()) {
			lastMeniscus = bodiesMenisciiList6D[i].end();
			for (firstMeniscus = bodiesMenisciiList6D[i].begin(); firstMeniscus != lastMeniscus;
			     ++firstMeniscus) { //FOR EACH MENISCUS ON THIS BODY...
				currentMeniscus = firstMeniscus;
				++currentMeniscus;
				cundallInteractionPhysics1 = YADE_CAST<CapillaryCohPhys*>((*firstMeniscus)->phys.get());
				if (i == (*firstMeniscus)->getId1()) angle1 = cundallInteractionPhysics1->Delta1; //get angle of meniscus1 on body i
				else
					angle1 = cundallInteractionPhysics1->Delta2;
				for (; currentMeniscus != lastMeniscus; ++currentMeniscus) { //... CHECK FUSION WITH ALL OTHER MENISCII ON THE BODY
					cundallInteractionPhysics2 = YADE_CAST<CapillaryCohPhys*>((*currentMeniscus)->phys.get());
					if (i == (*currentMeniscus)->getId1())
						angle2 = cundallInteractionPhysics2->Delta1; //get angle of meniscus2 on body i
					else
						angle2 = cundallInteractionPhysics2->Delta2;
					if (angle1 == 0 || angle2 == 0) cerr << "THIS SHOULD NOT HAPPEN!!" << endl;

					//cerr << "angle1 = " << angle1 << " | angle2 = " << angle2 << endl;

					Vector3r normalFirstMeniscus   = YADE_CAST<ScGeom6D*>((*firstMeniscus)->geom.get())->normal;
					Vector3r normalCurrentMeniscus = YADE_CAST<ScGeom6D*>((*currentMeniscus)->geom.get())->normal;

					Real normalDot = 0;
					if ((*firstMeniscus)->getId1() == (*currentMeniscus)->getId1()
					    || (*firstMeniscus)->getId2() == (*currentMeniscus)->getId2())
						normalDot = normalFirstMeniscus.dot(normalCurrentMeniscus);
					else
						normalDot = -(normalFirstMeniscus.dot(normalCurrentMeniscus));

					Real normalAngle = 0;
					if (normalDot >= 0) normalAngle = math::fastInvCos0(normalDot);
					else
						normalAngle = ((Mathr::PI)-math::fastInvCos0(-(normalDot)));

					if ((angle1 + angle2) * Mathr::DEG_TO_RAD > normalAngle) {
						++(cundallInteractionPhysics1->fusionNumber);
						++(cundallInteractionPhysics2->fusionNumber);
					}
				}
			}
		}
	}
}

//what is this function for? it is never called...
//TODO: remove?
BodiesMenisciiList6D::BodiesMenisciiList6D(Scene* scene)
{
	initialized = false;
	prepare(scene);
}

bool BodiesMenisciiList6D::prepare(Scene* scene)
{
	//cerr << "preparing bodiesInteractionsList" << endl;
	interactionsOnBody.clear();
	shared_ptr<BodyContainer>& bodies = scene->bodies;

	Body::id_t MaxId = -1;

	for (const auto& b : *bodies) {
		MaxId = max(MaxId, b->getId());
	}
	interactionsOnBody.resize(MaxId + 1);
	for (unsigned int i = 0; i < interactionsOnBody.size(); ++i) {
		interactionsOnBody[i].clear();
	}

	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{ // parallel version using Engine::ompThreads variable not accessible, this function is not one of the Engine..
		if (I->isReal()) {
			if (static_cast<CapillaryCohPhys*>(I->phys.get())->meniscus) insert(I);
		}
	}

	return initialized = true;
}

bool BodiesMenisciiList6D::insert(const shared_ptr<Interaction>& interaction)
{
	checkLengthBuffer(interaction);
	interactionsOnBody[interaction->getId1()].push_back(interaction);
	interactionsOnBody[interaction->getId2()].push_back(interaction);
	return true;
}


bool BodiesMenisciiList6D::remove(const shared_ptr<Interaction>& interaction)
{
	checkLengthBuffer(interaction);
	interactionsOnBody[interaction->getId1()].remove(interaction);
	interactionsOnBody[interaction->getId2()].remove(interaction);
	return true;
}

void BodiesMenisciiList6D::checkLengthBuffer(const shared_ptr<Interaction>& interaction)
{
	Body::id_t maxBodyId = std::max(interaction->getId1(), interaction->getId2());
	if (unsigned(maxBodyId) >= interactionsOnBody.size()) { interactionsOnBody.resize(maxBodyId + 1); }
}

std::list<shared_ptr<Interaction>>& BodiesMenisciiList6D::operator[](int index) { return interactionsOnBody[index]; }

int BodiesMenisciiList6D::size() { return interactionsOnBody.size(); }

void BodiesMenisciiList6D::display()
{
	std::list<shared_ptr<Interaction>>::iterator firstMeniscus;
	std::list<shared_ptr<Interaction>>::iterator lastMeniscus;
	for (unsigned int i = 0; i < interactionsOnBody.size(); ++i) {
		if (!interactionsOnBody[i].empty()) {
			lastMeniscus = interactionsOnBody[i].end();
			//cerr << "size = "<<interactionsOnBody[i].size() << " empty="<<interactionsOnBody[i].empty() <<endl;
			for (firstMeniscus = interactionsOnBody[i].begin(); firstMeniscus != lastMeniscus; ++firstMeniscus) {
				if (*firstMeniscus) {
					if (firstMeniscus->get()) cerr << "(" << (*firstMeniscus)->getId1() << ", " << (*firstMeniscus)->getId2() << ") ";
					else
						cerr << "(void)";
				}
			}
			cerr << endl;
		} else
			cerr << "empty" << endl;
	}
}

BodiesMenisciiList6D::BodiesMenisciiList6D() { initialized = false; }
} // namespace yade
