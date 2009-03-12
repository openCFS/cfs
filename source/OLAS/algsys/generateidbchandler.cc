// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "idbchandler.hh"
#include "idbchandlervoid.hh"
#include "idbchandlerpenalty.hh"
#include "generateidbchandler.hh"

// Required for using AssocType
#include "MatVec/typedefs.hh"


namespace CoupledField {


  // ******************************
  //   GenerateIDBC_HandlerObject
  // ******************************
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( const std::set<FEMatrixType> usedFEMatrices,
                              BaseGraphManager *graphManager, UInt numPDEs,
                              UInt numIDBCs, const BaseMatrix::EntryType eType,
                              bool sbmCase ) {

    // VERSION:
    //
    // - This version of the factory generates objects of type IDBC_Handler,
    //   i.e. objects that use the elimination approach to IDBC handling.
    //
    // - It works for both StandardSystem and SBM_System
    //
    // - If no IDBCs are present in the linear system an IDBC_HandlerVoid
    //   object will be generated instead.


    BaseIDBC_Handler *retVal = NULL;

    // If no inhomogeneous Dirichlet values are present
    // generate a pseudo idbc handler
    if ( numIDBCs == 0 ) {
      retVal = new IDBC_HandlerVoid();
      ASSERTMEM( retVal, sizeof(IDBC_HandlerVoid) );
      (*cla) << " GenerateIDBC_HandlerObject: Generated"
             << " IDBC_HandlerVoid" << std::endl;
    }

    // Generate a real idbc handler
    else {

      if ( eType == BaseMatrix::DOUBLE ) {
        retVal =
          new IDBC_Handler<Double>( usedFEMatrices, graphManager, numPDEs,
                                    sbmCase );
        ASSERTMEM( retVal, sizeof(IDBC_Handler<Double>) );
        (*cla) << " GenerateIDBC_HandlerObject: Generated"
               << " IDBC_Handler<Double>"
               << std::endl;
      }
      else if ( eType == BaseMatrix::COMPLEX ) {
        retVal =
          new IDBC_Handler<Complex>( usedFEMatrices, graphManager, numPDEs,
                                     sbmCase );
        ASSERTMEM( retVal, sizeof(IDBC_Handler<Complex>) );
        (*cla) << " GenerateIDBC_HandlerObject: Generated"
               << " IDBC_Handler<Complex>"
               << std::endl;
      }
      else {
        EXCEPTION( "GenerateIDBC_HandlerObject: Cannot generate "
                 << "IDBC_Handler<T> object with T = '"
                 << BaseMatrix::entryType.ToString( eType ) << "'");
      }

    }

    return retVal;
  }


  // ******************************
  //   GenerateIDBC_HandlerObject
  // ******************************
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( UInt numIDBC, UInt blockSize,
                              const BaseMatrix::EntryType eType ) {

    // VERSION:
    //
    // - This version of the factory generates objects of type
    //   IDBC_HandlerPenalty, i.e. objects that use the penalty
    //   approach to IDBC handling.
    //
    // - It works only for the StandardSystem class
    //
    // - If no IDBCs are present in the linear system an IDBC_HandlerVoid
    //   object will be generated instead.


    BaseIDBC_Handler *retVal = NULL;

    // If no inhomogeneous Dirichlet values are present
    // generate a pseudo idbc handler
    if ( numIDBC == 0 ) {
      retVal = new IDBC_HandlerVoid();
      ASSERTMEM( retVal, sizeof(IDBC_HandlerVoid) );
      (*cla) << " GenerateIDBC_HandlerObject: Generated"
             << " IDBC_HandlerVoid" << std::endl;
    }

    // Generate a real idbc handler
    else {

      if ( eType == BaseMatrix::DOUBLE ) {
        retVal = new IDBC_HandlerPenalty<Double>( numIDBC, blockSize );
        ASSERTMEM( retVal, sizeof(IDBC_HandlerPenalty<Double>) );
        (*cla) << " GenerateIDBC_HandlerObject: Generated"
               << " IDBC_HandlerPenalty<Double>"
               << std::endl;
      }
      else if ( eType == BaseMatrix::COMPLEX ) {
        retVal = new IDBC_HandlerPenalty<Complex>( numIDBC, blockSize );
        ASSERTMEM( retVal, sizeof(IDBC_HandlerPenalty<Complex>) );
        (*cla) << " GenerateIDBC_HandlerObject: Generated"
               << " IDBC_HandlerPenalty<Complex>"
               << std::endl;
      }
      else {
        EXCEPTION( "GenerateIDBC_HandlerObject: Cannot generate "
                 << "IDBC_HandlerPenalty<T> object with T = '"
                 << BaseMatrix::entryType.ToString( eType ) << "'");
      }

    }

    return retVal;
  }


  // ******************************
  //   GenerateIDBC_HandlerObject
  // ******************************
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( UInt numIDBC, UInt numPDEs,
                              UInt *bcOffsets, const BaseMatrix::EntryType eType ) {

    // VERSION:
    //
    // - This version of the factory generates objects of type
    //   IDBC_HandlerPenalty, i.e. objects that use the penalty
    //   approach to IDBC handling.
    //
    // - It works only for the SBM_System class


    BaseIDBC_Handler *retVal = NULL;

    if ( eType == BaseMatrix::DOUBLE ) {
      retVal = new IDBC_HandlerPenalty<Double>( numIDBC, numPDEs, bcOffsets );
      ASSERTMEM( retVal, sizeof(IDBC_HandlerPenalty<Double>) );
      (*cla) << " GenerateIDBC_HandlerObject: Generated"
             << " IDBC_HandlerPenalty<Double>"
             << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retVal = new IDBC_HandlerPenalty<Complex>( numIDBC, numPDEs, bcOffsets);
      ASSERTMEM( retVal, sizeof(IDBC_HandlerPenalty<Complex>) );
      (*cla) << " GenerateIDBC_HandlerObject: Generated"
             << " IDBC_HandlerPenalty<Complex>"
             << std::endl;
    }
    else {
      EXCEPTION( "GenerateIDBC_HandlerObject: Cannot generate "
               << "IDBC_HandlerPenalty<T> object with T = '"
               << BaseMatrix::entryType.ToString( eType ) << "'");
    }

    return retVal;
  }

}
