// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_MATVEC_HEADERS_HH
#define OLAS_MATVEC_HEADERS_HH

#include "matvec/typedefs.hh"
#include "matvec/opdefs.hh"

#include "matvec/basematrix.hh"
#include "matvec/stdmatrix.hh"
#include "matvec/sbmmatrix.hh"
#include "matvec/basevector.hh"
#include "matvec/stdvector.hh"
#include "matvec/sbmvector.hh"
#include "matvec/vector.hh"
#include "matvec/generatematvec.hh"

// Specialized Matrices
#include "matvec/crs_matrix.hh"
#include "matvec/scrs_matrix.hh"
#include "matvec/diag_matrix.hh"

// simple sparse matrix container
#include "matvec/coordformat.hh"

// pattern management
#include "matvec/sparsitypatterns.hh"
#include "matvec/patternpool.hh"

#endif
