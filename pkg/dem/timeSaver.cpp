/*************************************************************************
*  Copyright (C) 2026 by Chang GUO                                       *
*  chang.guo@connect.polyu.hk                                            *
*                                                                        *
*  timeSaver: periodic snapshot-saving engine, added in V2.0.0.          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "timeSaver.hpp"
#include <core/Omega.hpp>

namespace yade { // Cannot have #include directive inside.

CREATE_LOGGER(timeSaver);
YADE_PLUGIN((timeSaver));
timeSaver::~timeSaver() {};

void timeSaver::action()
{
	// Called every realPeriod seconds (wall-clock). Snapshot the whole scene so the
	// simulation can be resumed after an interruption. The previous snapshot is overwritten.
	Omega::instance().saveSimulation(savePath + "timeSaver.yade.bz2");
}

} // namespace yade
