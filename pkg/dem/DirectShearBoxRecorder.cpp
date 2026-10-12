/*************************************************************************
*  Copyright (C) 2026 by Chang GUO                                       *
*  chang.guo@connect.polyu.hk                                            *
*                                                                        *
*  DirectShearBoxRecorder: records the six-face stress / displacement    *
*  of the direct-shear box (added in V2.0.0).                            *
*                                                                        *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "DirectShearBoxRecorder.hpp"
#include <lib/high-precision/Constants.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/dem/DirectShearBox.hpp>
#include <preprocessing/dem/Shop.hpp>
#include <pkg/common/PeriodicEngines.hpp>

namespace yade { // Cannot have #include directive inside.

CREATE_LOGGER(DirectShearBoxRecorder);
DirectShearBoxRecorder::~DirectShearBoxRecorder() {};

// void DirectShearBoxRecorder::change2file2()
// {
// 	assert(!out.is_open());
// 	std::string fileTemp = file2;
// 	if (fileTemp.empty()) throw ios_base::failure(__FILE__ ": Empty filename.");
// 	out.open(fileTemp.c_str(), std::fstream::app);
// 	if (!out.good()) throw ios_base::failure(__FILE__ ": I/O error opening file `" + fileTemp + "'.");
// }


void DirectShearBoxRecorder::action()
{
	// Called every iterPeriod by the simulation loop.
	// 1. Find the (first) DirectShearBox / DirectShearBox6D engine in the scene and cache it.
	// 2. Make sure boundary stresses are fresh (call computeStress if it was not done this iter).
	// 3. Write one CSV line. Column order (also the header written on the first line):
	//    iteration, virtime, realtime,
	//    pvz            = pressor Z-velocity (wall_id_Z1),
	//    pz             = free height heightZ / scalingFactor,
	//    sheardisx/y/z  = shear-box displacement (Y divided by scalingFactor),
	//    refPz          = relative vertical displacement (heightZ - heightZ0)/scalingFactor,
	//    shearratey     = shear-box Y velocity,
	//    pressorgoal    = target pressor stress,
	//    normalstress   = measured pressor stress (stress[Z1][2]),
	//    shearstress    = -Fy(shearbox) / (depthY0*widthX0),
	//    porosity, d50 (currentd50), unb (unbalanced force), kineticE, meanstress.
	Scene* scene2;
	if (scene) scene2 = scene;
	else scene2 = Omega::instance().getScene().get();
	scene2->forces.sync();
	// at the beginning of the file; write column titles
	// out is std::ofstream class object, see in Recorder.hpp
	if (out.tellp() == 0) out << "iteration,virtime,realtime,pvz,pz,sheardisx,sheardisy,sheardisz,refPz,shearratey,pressorgoal,normalstress,shearstress,porosity,d50,unb,kineticE,meanstress" << endl;

	if (!dsb) {
		vector<shared_ptr<Engine>>::iterator itFirst = scene2->engines.begin();
		vector<shared_ptr<Engine>>::iterator itLast  = scene2->engines.end();
		for (; itFirst != itLast; ++itFirst) {
			const std::string& className = (*itFirst)->getClassName();
			if (className == "DirectShearBox" || className == "DirectShearBox6D") {
				LOG_DEBUG("DirectShearBox controller engine found");
				dsb = YADE_PTR_CAST<DirectShearBox>(*itFirst);
				break;
			}
		}
		if (!dsb) LOG_ERROR("DirectShearBoxe NOT found");
	}
	if (!(scene2->iter % dsb->computeStressInterval == 0)) dsb->computeStress();

	State* p_shear   = Body::byId(dsb->shearbox_id, scene2)->state.get();
	State* p_wall_z1 = Body::byId(dsb->wall_id_Z1, scene2)->state.get();

	out << boost::lexical_cast<string>(scene2->iter) << ","
		<< boost::lexical_cast<string>(scene2->time) << ","
		<< boost::lexical_cast<string>(getClock()-realtime0) << ","
		<< boost::lexical_cast<string>(p_wall_z1->vel[2]) << ","
		<< boost::lexical_cast<string>(dsb->heightZ/scalingFactor) << ","
		<< boost::lexical_cast<string>(p_shear->displ()[0]) << ","
		<< boost::lexical_cast<string>(p_shear->displ()[1]/scalingFactor) << ","
		<< boost::lexical_cast<string>(p_shear->displ()[2]) << ","
		<< boost::lexical_cast<string>((dsb->heightZ-dsb->heightZ0)/scalingFactor) << ","
		<< boost::lexical_cast<string>(p_shear->vel[1]) << ","
	    << boost::lexical_cast<string>(dsb->pressorGoal) << ","
	    << boost::lexical_cast<string>(dsb->stress[dsb->wallsZ1][2]) << ","
	    << boost::lexical_cast<string>(-(scene2->forces.getForce(dsb->shearbox_id)).y()/(dsb->depthY0*dsb->widthX0)) 
		<< ","
		<< boost::lexical_cast<string>(dsb->porosity) << ","
		<< boost::lexical_cast<string>(dsb->currentd50) << ","
	    << boost::lexical_cast<string>(Shop::unbalancedForce(false, scene2)) << "," 
	    << boost::lexical_cast<string>(Shop::kineticEnergy()) <<  ","
	    << boost::lexical_cast<string>(dsb->meanStress) 
		<<  endl;
}

YADE_PLUGIN((DirectShearBoxRecorder));

} // namespace yade
