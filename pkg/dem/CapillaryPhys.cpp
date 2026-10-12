// Modifications by Chang GUO (2026): added Ip2_CohFrictMat_CohFrictMat_CapillaryCohPhys
// (CapillaryCohPhys, the capillary interaction physics for cohesive-frictional
// contacts) and the code moved here from Ip2_FrictMat_FrictMat_CapillaryPhys.cpp.
// Licensed under the GNU General Public License v2 or later; see file LICENSE.
#include <pkg/dem/CapillaryPhys.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.
using math::max;
using math::min;

YADE_PLUGIN((CapillaryPhys)(Ip2_FrictMat_FrictMat_CapillaryPhys)(CapillaryCohPhys)(Ip2_CohFrictMat_CohFrictMat_CapillaryCohPhys));

// The following code was moved from Ip2_FrictMat_FrictMat_CapillaryPhys.cpp

void Ip2_FrictMat_FrictMat_CapillaryPhys::go(
        const shared_ptr<Material>& b1 //FrictMat
        ,
        const shared_ptr<Material>& b2 // FrictMat
        ,
        const shared_ptr<Interaction>& interaction)
{
	ScGeom* geom = YADE_CAST<ScGeom*>(interaction->geom.get());
	if (geom) {
		if (!interaction->phys) {
			const shared_ptr<FrictMat>& sdec1 = YADE_PTR_CAST<FrictMat>(b1);
			const shared_ptr<FrictMat>& sdec2 = YADE_PTR_CAST<FrictMat>(b2);

			if (!interaction->phys) interaction->phys = shared_ptr<CapillaryPhys>(new CapillaryPhys());
			const shared_ptr<CapillaryPhys>& contactPhysics = YADE_PTR_CAST<CapillaryPhys>(interaction->phys);

			Real Ea = sdec1->young;
			Real Eb = sdec2->young;
			Real Va = sdec1->poisson;
			Real Vb = sdec2->poisson;
			Real Da = geom->radius1; // FIXME - multiply by factor of sphere interaction distance (so sphere interacts at bigger range that its geometrical size)
			Real Db = geom->radius2; // FIXME - as above
			Real fa = sdec1->frictionAngle;
			Real fb = sdec2->frictionAngle;
			Real Kn = 2 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); //harmonic average of two stiffnesses
			Real Ks = 2 * Ea * Da * Va * Eb * Db * Vb
			        / (Ea * Da * Va + Eb * Db * Va); //harmonic average of two stiffnesses with ks=V*kn for each sphere

			contactPhysics->tangensOfFrictionAngle = math::tan(math::min(fa, fb));
			contactPhysics->kn                     = Kn;
			contactPhysics->ks                     = Ks;
		}
	}
};
void Ip2_CohFrictMat_CohFrictMat_CapillaryCohPhys::go(
        const shared_ptr<Material>& b1 // CohFrictMat
        ,
        const shared_ptr<Material>& b2 // CohFrictMat
        ,
        const shared_ptr<Interaction>& interaction)
{
	CohFrictMat* sdec1 = static_cast<CohFrictMat*>(b1.get());
	CohFrictMat* sdec2 = static_cast<CohFrictMat*>(b2.get());
	ScGeom6D*    geom  = YADE_CAST<ScGeom6D*>(interaction->geom.get());

	//Create cohesive interactions only once
	if (setCohesionNow && cohesionDefinitionIteration == -1) cohesionDefinitionIteration = scene->iter;
	if (setCohesionNow && cohesionDefinitionIteration != -1 && cohesionDefinitionIteration != scene->iter) {
		cohesionDefinitionIteration = -1;
		setCohesionNow              = 0;
	}

	if (geom) {
		const auto normalAdhPreCalculated
		        = (normalCohesion) ? (*normalCohesion)(b1->id, b2->id) : math::min(sdec1->normalCohesion, sdec2->normalCohesion);
		const auto shearAdhPreCalculated = (shearCohesion) ? (*shearCohesion)(b1->id, b2->id) : math::min(sdec1->shearCohesion, sdec2->shearCohesion);

		if (!interaction->phys) {
			interaction->phys            = shared_ptr<CapillaryCohPhys>(new CapillaryCohPhys());
			CapillaryCohPhys* contactPhysics = YADE_CAST<CapillaryCohPhys*>(interaction->phys.get());
			Real          Ea             = sdec1->young;
			Real          Eb             = sdec2->young;
			Real          Va             = sdec1->poisson;
			Real          Vb             = sdec2->poisson;
			Real          Da             = geom->radius1;
			Real          Db             = geom->radius2;
			Real          fa             = sdec1->frictionAngle;
			Real          fb             = sdec2->frictionAngle;
			Real          Kn             = 2.0 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); //harmonic average of two stiffnesses
			Real          frictionAngle  = (!frictAngle) ? math::min(fa, fb) : (*frictAngle)(sdec1->id, sdec2->id, fa, fb);

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
			contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);

			if ((setCohesionOnNewContacts || setCohesionNow) && sdec1->isCohesive && sdec2->isCohesive) {
				contactPhysics->cohesionBroken = false;
				contactPhysics->normalAdhesion = normalAdhPreCalculated * pow(math::min(Db, Da), 2);
				contactPhysics->shearAdhesion  = shearAdhPreCalculated * pow(math::min(Db, Da), 2);
				geom->initRotations(*(Body::byId(interaction->getId1(), scene)->state), *(Body::byId(interaction->getId2(), scene)->state));
				contactPhysics->fragile = (sdec1->fragile || sdec2->fragile);
			}
			contactPhysics->kn = Kn;
			contactPhysics->ks = Ks;

			contactPhysics->maxRollPl         = min(sdec1->etaRoll * Da, sdec2->etaRoll * Db);
			contactPhysics->maxTwistPl        = min(sdec1->etaTwist * Da, sdec2->etaTwist * Db);
			contactPhysics->momentRotationLaw = (sdec1->momentRotationLaw && sdec2->momentRotationLaw);
		} else { // !isNew, but if setCohesionNow, all contacts are initialized like if they were newly created
			CapillaryCohPhys* contactPhysics = YADE_CAST<CapillaryCohPhys*>(interaction->phys.get());
			if ((setCohesionNow && sdec1->isCohesive && sdec2->isCohesive) || contactPhysics->initCohesion) {
				contactPhysics->cohesionBroken = false;
				contactPhysics->normalAdhesion = normalAdhPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);
				contactPhysics->shearAdhesion  = shearAdhPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);

				geom->initRotations(*(Body::byId(interaction->getId1(), scene)->state), *(Body::byId(interaction->getId2(), scene)->state));
				contactPhysics->fragile      = (sdec1->fragile || sdec2->fragile);
				contactPhysics->initCohesion = false;
			}
		}
	}
};
} // namespace yade
