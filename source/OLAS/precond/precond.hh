// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_PRECOND_HEADERS_HH
#define OLAS_PRECOND_HEADERS_HH


#include "precond/baseprecond.hh"
#include "precond/bnprecond.hh"
#include "precond/generateprecond.hh"

#include "precond/idprecond.hh"
#include "precond/jacprecond.hh"
#include "precond/sbmdiagprecond.hh"
#include "precond/ilu0precond.hh"
#include "precond/ilutpprecond.hh"
#include "precond/mgprecond.hh"

#ifdef PARALLEL
#include "precond/parbnprecond.hh"
#include "precond/ddprecond.hh"
#include "precond/pargenerateprecond.hh"
#endif

#endif // OLAS_PRECOND_HEADERS_HH 
