/*************************************************************************
*  Copyright (C) 2026 by Chang GUO                                       *
*  chang.guo@connect.polyu.hk                                            *
*                                                                        *
*  timeSaver: periodic snapshot-saving engine, added in V2.0.0.          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#pragma once

#include <pkg/common/PeriodicEngines.hpp>
#include <pkg/dem/NewtonIntegrator.hpp>

namespace yade { 
class timeSaver : public PeriodicEngine {
public:
	virtual ~timeSaver();
	void action() override;

	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(timeSaver,PeriodicEngine,
		"Engine that periodically saves the whole simulation to a .yade file, so a run interrupted "
		"by a crash/power loss can be resumed. By default it saves to savePath + \"timeSaver.yade.bz2\" "
		"every realPeriod seconds (1 hour if realPeriod is left at 0). It uses wall-clock time, not "
		"simulation time, so the period is independent of the time step."
		,
		/* attrs */
		((string,savePath,"log/",,"Directory where the timeSaver snapshot is written. You may add a descriptive sub-path; please verify it exists before running."))
		,
		/* extra initializers */
		,
   		/* constructor */
		if (!realPeriod) realPeriod=60*60;
		,
		/* py */
		)
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(timeSaver);

} // namespace yade
