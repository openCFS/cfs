// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/*! \file generateidbchandler.hh

    This file provides functions that act as factories for the generation
    of instances of the templated IDBC_Handler class.
    \remark Note that the documentation of the functions in this file
    is not be found here, but in the description of the namespace
    OLAS to which they belong.
*/


#ifndef OLAS_GENERATE_IDBC_HANLDER_HH
#define OLAS_GENERATE_IDBC_HANLDER_HH

#include <set>

#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/basematrix.hh"

namespace CoupledField {

  class BaseGraphManager;
  // forward declaration
  class BaseIDBC_Handler;

  //! Function for generating an IDBC_Handler object for elimination approach
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( const std::set<FEMatrixType> usedFEMatrices,
                              BaseGraphManager *graphManager, UInt numPDEs,
                              UInt numIDBCs,
                              const BaseMatrix::EntryType eType,
                              bool sbmCase );

  //! Function for generating an IDBC_HandlerPenalty object for StdMatrices
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( UInt numIDBC, UInt blockSize,
                              const BaseMatrix::EntryType eType );

  //! Function for generating an IDBC_HandlerPenalty object for SBM_Matrices
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( UInt numIDBC, UInt numPDEs, UInt *bcOffsets,
                              const BaseMatrix::EntryType eType );
}

#endif

