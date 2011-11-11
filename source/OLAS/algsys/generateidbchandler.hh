#ifndef OLAS_GENERATE_IDBC_HANLDER_HH
#define OLAS_GENERATE_IDBC_HANLDER_HH
/*! \file generateidbchandler.hh

    This file provides functions that act as factories for the generation
    of instances of the templated IDBC_Handler class.
    \remark Note that the documentation of the functions in this file
    is not be found here, but in the description of the namespace
    OLAS to which they belong.
*/

#include <set>

namespace CoupledField {

  // forward declaration
  class BaseIDBC_Handler;
  class GraphManager;

  //! Function for generating an IDBC_Handler object for elimination approach
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( const std::set<FEMatrixType> usedFEMatrices,
                              GraphManager *graphManager, StdVector<UInt>& numIDBC,
                              const BaseMatrix::EntryType eType );

  //! Function for generating an IDBC_HandlerPenalty object for SBM_Matrices
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObjectPenalty( StdVector<UInt>& numIDBC,
                                     const BaseMatrix::EntryType eType );
}

#endif

