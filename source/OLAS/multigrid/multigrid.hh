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

#include "multigrid/ppflags.hh"

// include sources with implementation in *.cc
//#ifdef INCLUDE_MULTIGRID_CC_FILES

//#include  "multigrid/prematrix.cc"
#include  "multigrid/depgraph.cc"
#include  "multigrid/topology.cc"
#include  "multigrid/transfer.cc"
//#include  "multigrid/hierarchylevel.cc"
#include  "multigrid/smoother.cc"
#include  "multigrid/gaussseidel.cc"
#include  "multigrid/jacobi.cc"
#include  "multigrid/amg.cc"

#include  "multigrid/preMatrix.hh"
#include  "multigrid/depgraph.hh"
#include  "multigrid/topology.hh"
#include  "multigrid/transfer.hh"
#include  "multigrid/hierarchylevel.hh"
#include  "multigrid/smoother.hh"
#include  "multigrid/gaussseidel.hh"
#include  "multigrid/jacobi.hh"
#include  "multigrid/amg.hh"

//#endif // INCLUDE_MULTIGRID_CC_FILES
