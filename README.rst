===================================
YADE - Yet Another Dynamic Engine
===================================

Extensible open-source framework for discrete numerical models, focused on 
Discrete Element Method. The computation parts are written in c++ using flexible
object model, allowing independent implementation of new algorithms and interfaces. 
Python is used for rapid and concise scene construction, simulation control, 
postprocessing and debugging.

**Documentation:**

- https://www.yade-dem.org/doc/ (includes c++ and python reference)
- https://yade-dem.org/doc/Yade.pdf (PDF-file with documentation)
- https://yade-dem.org/doc/Yade.epub (EPUB-file with documentation)
 
**Downloads:**

- Packages for Ubuntu and Debian:
  
  - https://yade-dem.org/doc/installation.html#packages
  - External packages: https://launchpad.net/~yade-users/+archive/external

- Tarballs of different versions: 

  - https://gitlab.com/yade-dev/trunk/tags
  - Repository on gitlab.com: https://gitlab.com/yade-dev/trunk/-/archive/master/trunk-master.tar.gz
  - Repository on git.debian.org: http://git.debian.org/?p=debian-science/packages/yade.git

**Teams involved:**

- General users: https://launchpad.net/~yade-users
- Developers and advanced users: https://launchpad.net/~yade-dev
- Developers for distributed memory parallelization of Yade: https://launchpad.net/~yade-mpi
- Infrastructure maintainers: https://launchpad.net/~yade-pkg

**Miscellaneous:**

- `Wiki <https://yade-dem.org/wiki/>`_
- `Screenshots <https://www.yade-dem.org/wiki/Screenshots_and_videos>`_
- `Tutorial <https://yade-dem.org/doc/tutorial-examples.html>`_
- `Videos <https://yade-dem.org/doc/tutorial-more-examples-fast.html>`_
- `Yade page on Open Hub <http://www.openhub.net/p/yade/>`_
- `Yade on Debian PTS <https://tracker.debian.org/pkg/yade>`_
- `Yade on Debian salsa <https://salsa.debian.org/science-team/yade/commits/master>`_
- `Yade in Ubuntu <https://launchpad.net/ubuntu/+source/yade>`_
- `Launchpad page <http://www.launchpad.net/yade>`_


Copyright and license
=====================
This repository is a snapshot of YADE (Yet Another Dynamic Engine). The whole
upstream tree is included and its copyright notices and licence (GNU General
Public License, version 2 or later -- see ``LICENSE``) apply to it unchanged.

The components added on top of that tree -- the DirectShear suite
(``DirectShearBox`` / ``DirectShearBox6D`` / ``DirectShearBoxRecorder`` /
``timeSaver``) and the capillary model for the 6D (cohesive-frictional) contact
model (``CapillaryCohPhys`` with ``Law2_ScGeom6D_CapillaryCohPhys_Capillarity``,
plus the Shop / python helpers that go with them) -- are

::

  Copyright (C) 2026 by Chang GUO <chang.guo@connect.polyu.hk>

and are released under the same licence. Upstream files that were modified here
carry a notice of the change in their header; see ``COPYRIGHT`` for the list.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

Feature notes (for my own reference)
=====================================

DirectShearBox / DirectShearBox6D
----------------------------------
Engines to run a direct-shear (or interface-shear / soil-pipe) test.

* Two compaction modes, switched by ``internalCompaction``:
  - ``internalCompaction = true``: particles grow from ``initiald50`` toward
    ``targetd50`` (staged growth via r1..r3 / f0..f3), walls stay fixed.
  - ``internalCompaction = false``: the pressor (wall Z1) is servo-controlled to
    ``pressorGoal``; during ``sampling`` it steps through ``targetPressures``.
* After compaction, call ``setExternal()`` then ``setAxial()`` to start shearing.
  ``setAxial()`` inserts a ``DirectShearBoxRecorder`` automatically.
* ``autoSimulation = true`` lets the engine run the whole compaction -> shear
  workflow without external script stepping.
* DirectShearBox6D additionally controls rolling (etaRoll) / twisting (etaTwist)
  resistance on walls, particles and pipe.

DirectShearBoxRecorder
----------------------
Periodic recorder writing one CSV line per ``iterPeriod`` to ``Recorder::file``.
Columns (header line is written automatically):
iteration, virtime, realtime, pvz, pz, sheardisx, sheardisy, sheardisz, refPz,
shearratey, pressorgoal, normalstress, shearstress, porosity, d50, unb,
kineticE, meanstress.
Divide pz / sheardisy / refPz by ``scalingFactor`` (simulated d50 / real d50) to
recover real dimensions.

timeSaver
---------
Saves the whole scene to ``savePath + "timeSaver.yade.bz2"`` every ``realPeriod``
seconds of wall-clock time (default 1 hour), so a crashed run can be resumed.

Capillary suite
---------------
This release adds a capillary model which supports the 6D (cohesive-frictional)
contact model. The capillary law (``Law2_ScGeom_CapillaryPhys_Capillarity``) only
worked with ``CapillaryPhys`` / ``MindlinCapillaryPhys`` interactions -- no
cohesion, no 6D orientation data; this release adds ``CapillaryCohPhys`` and its
6D law ``Law2_ScGeom6D_CapillaryCohPhys_Capillarity``, so meniscus (liquid bridge)
effects can now be computed on ``ScGeom6D`` interactions as well.
(``CapillaryCohPhys`` is ``CohFrictPhys`` based, so the capillary model then also
carries cohesion and rolling / twisting resistance.)

The engines this release adds:

* ``CapillaryCohPhys`` -- capillary interaction physics (IPhys) for
  cohesive-frictional contacts. It is ``CohFrictPhys`` based, so besides the
  liquid-bridge fields (``meniscus``, ``isBroken``, ``capillaryPressure``,
  ``vMeniscus``, ``Delta1`` / ``Delta2``, ``fCap``, ``fusionNumber``,
  ``nn11`` / ``nn33``) it also carries normal / tangential cohesion and the
  rolling / twisting stiffness. It is created for ``CohFrictMat`` /
  ``CohFrictMat`` interactions by
  ``Ip2_CohFrictMat_CohFrictMat_CapillaryCohPhys``.
* ``Law2_ScGeom6D_CapillaryCohPhys_Capillarity`` -- the 6D counterpart of
  ``Law2_ScGeom_CapillaryPhys_Capillarity``. It expects ``CapillaryCohPhys``
  physics and ``ScGeom6D`` geometry, keeps its own ``BodiesMenisciiList6D`` for
  the menisci-overlap bookkeeping, and uses the same control parameter
  (``capillaryPressure``), the same ``surfaceTension`` and the same ``M(r=i)``
  files as the base law (the ``M(r=..)`` files must be reachable from the
  directory yade is run in). As for the base law, distant menisci require a
  matching ``Bo1_Sphere_Aabb::aabbEnlargeFactor`` /
  ``Ig2_Sphere_Sphere_ScGeom6D::interactionDetectionFactor`` pair and
  ``neverErase = True``.
* ``Shop.getCapillaryStressNew()`` -- Love-Weber capillary stress computed from
  ``CapillaryCohPhys`` (sphere-sphere interactions only). NOTE: this is a NEW
  function added beside the base ``Shop.getCapillaryStress()``: this release
  deliberately leaves ``getCapillaryStress()`` untouched and puts the
  CapillaryCohPhys aware version in ``getCapillaryStressNew()``.
* ``Shop.capillaryCalculation()`` -- the same aggregation in one pass: the six
  homogenised stress components, the capillary pressure and the water
  saturation (in percent).
* ``Shop.setFrictionEtaRoll()``, ``Shop.growParticlesWithInteractions()`` and
  ``Shop.growParticlesWith6DInteractions()`` -- companions used by the
  capillary / 6D workflows. All five helpers are exposed to python.

Version numbering scheme
------------------------
Releases are numbered with a leading ``V`` and three levels, e.g. ``V2.0.0``
(major.minor.patch): the git tag of a release carries that code, so this version
is tagged ``V2.0.0``.

A script's run command must match the value of ``v`` declared in the script's
"yade version" section, otherwise the script aborts with a Version error.

``V2.0.0`` is the release on this branch: it introduces the **DirectShear**
capability (``DirectShearBox6D`` / ``DirectShearBoxRecorder`` / ``timeSaver``)
and the **capillary model which supports the 6D contact model**
(``CapillaryCohPhys`` and its 6D law
``Law2_ScGeom6D_CapillaryCohPhys_Capillarity``), so that meniscus (liquid bridge)
effects can be computed on ``ScGeom6D`` (cohesive-frictional) interactions.
