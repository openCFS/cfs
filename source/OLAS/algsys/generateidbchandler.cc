#include "IDBC_Handler.hh"
#include "IDBC_HandlerVoid.hh"
#include "IDBC_HandlerPenalty.hh"
#include "generateidbchandler.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

// define logging stream
DEFINE_LOG(genIdbc, "genIdbc")


  // ******************************
  //   GenerateIDBC_HandlerObject
  // ******************************
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( const std::set<FEMatrixType> usedFEMatrices,
                              GraphManager *graphManager, StdVector<UInt>& numIDBC,
                              const BaseMatrix::EntryType eType ) {

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
    bool hasIDBCs = false;
    
    for( UInt i = 0; i < numIDBC.GetSize(); ++i ) {
      if (numIDBC[i] > 0) {
        hasIDBCs = true;
        break;
      }
    }
    
    if ( !hasIDBCs ) {
      retVal = new IDBC_HandlerVoid();
      ASSERTMEM( retVal, sizeof(IDBC_HandlerVoid) );
      LOG_DBG(genIdbc) << " GenerateIDBC_HandlerObject: Generated"
                        << " IDBC_HandlerVoid";
    }

    // Generate a real idbc handler
    else {

      UInt numBlocks = numIDBC.GetSize();
      if ( eType == BaseMatrix::DOUBLE ) {
        retVal =
          new IDBC_Handler<Double>( usedFEMatrices, graphManager, numBlocks);
        ASSERTMEM( retVal, sizeof(IDBC_Handler<Double>) );
        LOG_DBG(genIdbc) << " GenerateIDBC_HandlerObject: Generated"
                         << " IDBC_Handler<Double>";
      }
      else if ( eType == BaseMatrix::COMPLEX ) {
        retVal =
          new IDBC_Handler<Complex>( usedFEMatrices, graphManager, numBlocks );
        ASSERTMEM( retVal, sizeof(IDBC_Handler<Complex>) );
        LOG_DBG(genIdbc) << " GenerateIDBC_HandlerObject: Generated"
                         << " IDBC_Handler<Complex>";
      }
      else {
        EXCEPTION( "GenerateIDBC_HandlerObject: Cannot generate "
                 << "IDBC_Handler<T> object with T = '"
                 << MatrixEntryTypeEnum.ToString( eType ) << "'");
      }

    }

    return retVal;
  }


  //********************************
  //    GenerateIDBC_HandlerObjectMH
  //********************************
BaseIDBC_Handler*
GenerateIDBC_HandlerObjectMH( const std::set<FEMatrixType> usedFEMatrices,
                            GraphManager *graphManager, StdVector<UInt>& numIDBC,
                            const BaseMatrix::EntryType eType, UInt M,UInt a ) {
  // VERSION:
  //
  // - This version of the factory generates objects of type IDBC_Handler,
  //   i.e. objects that use the elimination approach to IDBC handling
  //   for the multiharmonic case.
  //   Difference to standard:
  //    - Two more parameter (M,a)
  //    - Calls a different IDBC_handler(..., M ,a)
  //
  // - It works for both StandardSystem and SBM_System
  //
  // - If no IDBCs are present in the linear system an IDBC_HandlerVoid
  //   object will be generated instead.


  BaseIDBC_Handler *retVal = NULL;

  // If no inhomogeneous Dirichlet values are present
  // generate a pseudo idbc handler
  bool hasIDBCs = false;

  for( UInt i = 0; i < numIDBC.GetSize(); ++i ) {
    if (numIDBC[i] > 0) {
      hasIDBCs = true;
      break;
    }
  }

  if ( !hasIDBCs ) {
    retVal = new IDBC_HandlerVoid();
    ASSERTMEM( retVal, sizeof(IDBC_HandlerVoid) );
    LOG_DBG(genIdbc) << " GenerateIDBC_HandlerObject: Generated"
                      << " IDBC_HandlerVoid";
  }

  // Generate a real idbc handler
  else {
    UInt numBlocks = numIDBC.GetSize();
    if ( eType == BaseMatrix::COMPLEX ) {
      //here it calls a different IDBC_Handler
      retVal =
        new IDBC_Handler<Complex>( usedFEMatrices, graphManager, numBlocks , M, a );
      ASSERTMEM( retVal, sizeof(IDBC_Handler<Complex>) );
      LOG_DBG(genIdbc) << " GenerateIDBC_HandlerObject: Generated"
                       << " IDBC_Handler<Complex>";
    }
    else {
      EXCEPTION( "GenerateIDBC_HandlerObject: Can only generate "
               << "IDBC_HandlerMH<Complex>  not for T = '"
               << MatrixEntryTypeEnum.ToString( eType ) << "'");
    }

  }

  return retVal;
}
  // ************************************
  //   GenerateIDBC_HandlerObjectPenalty
  // ************************************
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObjectPenalty( StdVector<UInt>& numIDBC, 
                                     const BaseMatrix::EntryType eType ) {

    // VERSION:
    //
    // - This version of the factory generates objects of type
    //   IDBC_HandlerPenalty, i.e. objects that use the penalty
    //   approach to IDBC handling.
    //
    // - It works only for the SBM_System class


    BaseIDBC_Handler *retVal = NULL;

    if ( eType == BaseMatrix::DOUBLE ) {
      retVal = new IDBC_HandlerPenalty<Double>( numIDBC  );
      ASSERTMEM( retVal, sizeof(IDBC_HandlerPenalty<Double>) );
      LOG_DBG(genIdbc) << " GenerateIDBC_HandlerObject: Generated"
                       << " IDBC_HandlerPenalty<Double>";
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retVal = new IDBC_HandlerPenalty<Complex>( numIDBC );
      ASSERTMEM( retVal, sizeof(IDBC_HandlerPenalty<Complex>) );
      LOG_DBG(genIdbc)<< " GenerateIDBC_HandlerObject: Generated"
                      << " IDBC_HandlerPenalty<Complex>";
    }
    else {
      EXCEPTION( "GenerateIDBC_HandlerObject: Cannot generate "
               << "IDBC_HandlerPenalty<T> object with T = '"
               << MatrixEntryTypeEnum.ToString( eType ) << "'");
    }

    return retVal;
  }

}
