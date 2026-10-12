//
// C++ Interface: Law2_ScGeom6D_CapillaryCohPhys_Capillarity
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

#pragma once

#include <core/GlobalEngine.hpp>
#include "Law2_ScGeom_CapillaryPhys_Capillarity.hpp"
namespace yade { // Cannot have #include directive inside.

/**
6D "capillary + cohesive" variant of Law2_ScGeom_CapillaryPhys_Capillarity.

It computes the same capillary (liquid bridge / meniscus) effects as the base class, but for
interactions whose physics is :yref:`CapillaryCohPhys<CapillaryCohPhys>` -- the cohesive-frictional
(CohFrictPhys based) capillary physics created by
:yref:`Ip2_CohFrictMat_CohFrictMat_CapillaryCohPhys` -- and whose geometry is a
:yref:`ScGeom6D<ScGeom6D>`. The 6D geometry is required because the cohesive part of the model
keeps the relative roll / twist orientations of the two bodies, hence a *6D* geometry functor
(e.g. :yref:`Ig2_Sphere_Sphere_ScGeom6D`) must be used for the interactions handled here.

Unlike the base class, this law performs no CapillaryPhys / MindlinCapillaryPhys switch: every
interaction is assumed to carry a CapillaryCohPhys, and meniscus overlaps (menisci fusion) are
tracked in the dedicated BodiesMenisciiList6D container declared below instead of
BodiesMenisciiList.

The control parameter is still the capillary pressure (or suction) Uc = Ugas - Uliquid, and the
liquid bridge properties (volume V, extent over the interacting grains delta1 and delta2) are read
from the same ascii files M(r=i) with i = R1/R2 as for the base class -- see the base class
documentation (Law2_ScGeom_CapillaryPhys_Capillarity.hpp) for how to obtain or regenerate them.
As in the base class, the formulation is valid only for pendular menisci involving two grains, and
menisci fusion uses the same assumptions: binary -> F_cap = 0 as soon as one overlap is detected,
linear -> F_cap = F_cap / (fusionNumber + 1).
*/

/// Interaction is used below in BodiesMenisciiList6D. MeniscusParameters, capillarylaw, Tableau
/// and TableauD need no declaration here: they are all fully defined in the base header
/// Law2_ScGeom_CapillaryPhys_Capillarity.hpp, which is included above.
class Interaction;

/// This container class is used to check if meniscii overlap, for the 6D / CapillaryCohPhys
/// variant of the law (the base class uses the equivalent BodiesMenisciiList). Wet interactions
/// are put in a series of lists, with one list per body.
class BodiesMenisciiList6D {
private:
	std::vector<std::list<shared_ptr<Interaction>>>
	        interactionsOnBody; // would require use of OpenMPVector (lib/base/openmp-accu.hpp) for parallel coding, see http://www.mail-archive.com/yade-dev@lists.launchpad.net/msg10842.html and msg11238.html

	//shared_ptr<Interaction> empty;

public:
	virtual ~BodiesMenisciiList6D() { };
	BodiesMenisciiList6D();
	BodiesMenisciiList6D(Scene*); //TODO: remove?
	bool                                prepare(Scene*);
	bool                                insert(const shared_ptr<Interaction>&);
	bool                                remove(const shared_ptr<Interaction>&);
	std::list<shared_ptr<Interaction>>& operator[](int);
	int                                 size();
	void                                display();
	void                                checkLengthBuffer(const shared_ptr<Interaction>&);

	bool initialized;
};

/// This is the constitutive law (the 6D "capillary + cohesive" variant; see the class comment above)
class Law2_ScGeom6D_CapillaryCohPhys_Capillarity : public Law2_ScGeom_CapillaryPhys_Capillarity {
public:
	virtual ~Law2_ScGeom6D_CapillaryCohPhys_Capillarity() { };
	void                     checkFusion() override;
	BodiesMenisciiList6D       bodiesMenisciiList6D;
	void action() override;
	// void postLoad(Law2_ScGeom6D_CapillaryCohPhys_Capillarity&);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(Law2_ScGeom6D_CapillaryCohPhys_Capillarity,Law2_ScGeom_CapillaryPhys_Capillarity,"This engine allows one to take into account capillary forces/effects between spheres coming from the presence of distinct interparticular liquid bridges (menisci) at existing interactions (between particle pairs, in the so-called pendular regime). In order to allow capillary forces between distant spheres, it is necessary to enlarge the bounding boxes using :yref:`Bo1_Sphere_Aabb::aabbEnlargeFactor` and make the Ig2 define define distant interactions via :yref:`interactionDetectionFactor<Ig2_Sphere_Sphere_ScGeom::interactionDetectionFactor>`. See also :yref:`createDistantMeniscii<Law2_ScGeom_CapillaryPhys_Capillarity.createDistantMeniscii>`. It is also necessary to disable interactions removal by the constitutive law (:yref:`Law2<Law2_ScGeom_FrictPhys_CundallStrack::neverErase>=True`). The interaction physics handled here is always :yref:`CapillaryCohPhys<CapillaryCohPhys>`, created by :yref:`Ip2_CohFrictMat_CohFrictMat_CapillaryCohPhys` (so this law goes with :yref:`CohFrictMat` and a cohesive-frictional contact law such as :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`), and the interaction geometry must be a :yref:`ScGeom6D<ScGeom6D>`\n\nThe control parameter for simulating liquid bridges is the :yref:`capillary pressure<Law2_ScGeom_CapillaryPhys_Capillarity::capillaryPressure>` (or suction) Uc = Ugas - Uliquid. Liquid bridges properties (volume V, extent over interacting grains delta1 and delta2) are computed as a result of the defined capillary pressure and of the interacting geometry (spheres radii and interparticular distance). For this purpose, ascii files M(r=i) with i=R1/R2 (and a possible :yref:`suffix<Law2_ScGeom_CapillaryPhys_Capillarity.suffCapFiles>`) are required, containing a set of results from the resolution of the Laplace-Young equation for different configurations of the interacting geometry, and can be downloaded from `yade-data/capillaryFiles <https://gitlab.com/yade-dev/yade-data/-/tree/main/capillaryFiles>`_ (direct download `here <https://gitlab.com/yade-dev/yade-data/-/archive/main/yade-data-main.zip?path=capillaryFiles>`_) for the case of a null wetting angle. They could also be regenerated, possibly for other conditions, by any user as per point II. in :ysrc:`examples/capillaryLaplaceYoung/README.md`.\n\nSee :ysrc:`examples/capillaryLaplaceYoung/` folder for example scripts.\n\nReferences: in english [Scholtes2009b]_ [Duriez2017b]_; more detailed, but in french [Scholtes2009d]_.",
	,,/*constructor*/
	showError = true;
	,
	 );
	// clang-format on
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(Law2_ScGeom6D_CapillaryCohPhys_Capillarity);

} // namespace yade
