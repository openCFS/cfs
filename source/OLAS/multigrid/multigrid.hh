// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

/**********************************************************
 * Modul header for module "multigrid"
 * Includes all the multigrid headers.
 * If INCLUDE_MULTIGRID_CC_FILES is defined, additionally
 * all source files (*.cc), that contain templated code, are
 * included, too.
 **********************************************************/

#include "OLAS/multigrid/ppflags.hh"

// include sources with implementation in *.cc
//#ifdef INCLUDE_MULTIGRID_CC_FILES

#include  "OLAS/multigrid/prematrix.cc"
#include  "OLAS/multigrid/depgraph.cc"
#include  "OLAS/multigrid/topology.cc"
#include  "OLAS/multigrid/agglomerate.cc"
#include  "OLAS/multigrid/transfer.cc"
#include  "OLAS/multigrid/hierarchylevel.cc"
#include  "OLAS/multigrid/smoother.cc"
#include  "OLAS/multigrid/gaussseidel.cc"
#include  "OLAS/multigrid/AFWsmoother.cc"
#include  "OLAS/multigrid/jacobi.cc"
//#include  "OLAS/multigrid/amg.cc"

#include  "OLAS/multigrid/prematrix.hh"
#include  "OLAS/multigrid/depgraph.hh"
#include  "OLAS/multigrid/topology.hh"
#include  "OLAS/multigrid/transfer.hh"
#include  "OLAS/multigrid/hierarchylevel.hh"
#include  "OLAS/multigrid/smoother.hh"
#include  "OLAS/multigrid/gaussseidel.hh"
#include  "OLAS/multigrid/AFWsmoother.hh"
#include  "OLAS/multigrid/jacobi.hh"
#include  "OLAS/multigrid/amg.hh"

//#endif // INCLUDE_MULTIGRID_CC_FILES
