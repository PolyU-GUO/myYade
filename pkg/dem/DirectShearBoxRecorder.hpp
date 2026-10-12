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


#pragma once

#include <pkg/common/Recorder.hpp>
#include <pkg/common/PeriodicEngines.hpp>
namespace yade { // Cannot have #include directive inside.

/*! \brief Record the state of a direct-shear simulation driven by DirectShearBox.

	This recorder reads the live state of the associated :yref:`DirectShearBox` engine
	(packing porosity, boundary stresses, shear-box displacement/velocity, applied
	frictions, total mass, unbalanced force, ...) and writes one line per record to the
	output text file. The header line lists every column so the format is self-describing
	(see DirectShearBoxRecorder::action() for the exact column order).
 */


class DirectShearBox;

class DirectShearBoxRecorder : public Recorder {
private:
	shared_ptr<DirectShearBox> dsb;
	Real realtime0=getClock();

public:
	virtual ~DirectShearBoxRecorder();
	void action() override;
	// void change2file2();

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(DirectShearBoxRecorder,Recorder,"Engine recording direct-shear variables of a :yref:`DirectShearBox` simulation (porosity, stresses, displacements, frictions, mass, unbalanced force). The first line of the output file lists all column names. Insert this recorder via DirectShearBox::setAxial().",
		((Real,scalingFactor,1,,"Scaling factor applied to the vertical stress (pz) and shear displacement, to account for particle-size scaling: simulated_dimension / real_dimension. Set it to (simulated d50 / real d50)."))
		// ((std::string,file2,"",,"File once axial action"))
		, 
		initRun=true;
		);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(DirectShearBoxRecorder);

} // namespace yade
