// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "algsys/idbchandler.cc"
#include "algsys/idbchandlervoid.hh"
#include "algsys/idbchandlerpenalty.cc"
#include "algsys/generateidbchandler.hh"

// Required for using assocType
#include "matvec/typedefs.cc"


namespace OLAS {


  // ******************************
  //   GenerateIDBC_HandlerObject
  // ******************************
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( const std::set<FEMatrixType> usedFEMatrices,
                              BaseGraphManager *graphManager, UInt numPDEs,
                              UInt numIDBCs, const MatrixEntryType eType,
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

    ENTER_FCN( "GenerateIDBC_HandlerObject" );

    BaseIDBC_Handler *retVal = NULL;

    // If no inhomogeneous Dirichlet values are present
    // generate a pseudo idbc handler
    if ( numIDBCs == 0 ) {
      retVal = new IDBC_HandlerVoid();
      AssertMem( retVal, sizeof(IDBC_HandlerVoid) );
      (*cla) << " GenerateIDBC_HandlerObject: Generated"
             << " IDBC_HandlerVoid" << std::endl;
    }

    // Generate a real idbc handler
    else {

      if ( eType == DOUBLE ) {
        retVal =
          New IDBC_Handler<Double>( usedFEMatrices, graphManager, numPDEs,
                                    sbmCase );
        AssertMem( retVal, sizeof(IDBC_Handler<Double>) );
        (*cla) << " GenerateIDBC_HandlerObject: Generated"
               << " IDBC_Handler<Double>"
               << std::endl;
      }
      else if ( eType == COMPLEX ) {
        retVal =
          New IDBC_Handler<Complex>( usedFEMatrices, graphManager, numPDEs,
                                     sbmCase );
        AssertMem( retVal, sizeof(IDBC_Handler<Complex>) );
        (*cla) << " GenerateIDBC_HandlerObject: Generated"
               << " IDBC_Handler<Complex>"
               << std::endl;
      }
      else {
        (*error) << "GenerateIDBC_HandlerObject: Cannot generate "
                 << "IDBC_Handler<T> object with T = '"
                 << Enum2String( eType ) << "'";
        Error( __FILE__, __LINE__ );
      }

      // Force instantiation of public member functions of IDBC_Handler class.
      // Note that the InstantiatePublicMethods() method itself should never
      // actually be called, but the compiler must not notice this. Thus, the
      // akward if condition.
      if ( usedFEMatrices.size() > MAX_NUM_FE_MATRICES ) {
        Error( "GenerateIDBC_HandlerObject: Internal error!", __FILE__,
               __LINE__ );
        if ( eType == DOUBLE ) {
          dynamic_cast< IDBC_Handler<Double>* >(retVal)
            ->InstantiatePublicMethods();
        }
        else {
          dynamic_cast< IDBC_Handler<Complex>* >(retVal)
            ->InstantiatePublicMethods();
        }
      }
    }

    return retVal;
  }


  // ******************************
  //   GenerateIDBC_HandlerObject
  // ******************************
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( UInt numIDBC, UInt blockSize,
                              const MatrixEntryType eType ) {

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

    ENTER_FCN( "GenerateIDBC_HandlerObject" );

    BaseIDBC_Handler *retVal = NULL;

    // If no inhomogeneous Dirichlet values are present
    // generate a pseudo idbc handler
    if ( numIDBC == 0 ) {
      retVal = new IDBC_HandlerVoid();
      AssertMem( retVal, sizeof(IDBC_HandlerVoid) );
      (*cla) << " GenerateIDBC_HandlerObject: Generated"
             << " IDBC_HandlerVoid" << std::endl;
    }

    // Generate a real idbc handler
    else {

      if ( eType == DOUBLE ) {
        retVal = New IDBC_HandlerPenalty<Double>( numIDBC, blockSize );
        AssertMem( retVal, sizeof(IDBC_HandlerPenalty<Double>) );
        (*cla) << " GenerateIDBC_HandlerObject: Generated"
               << " IDBC_HandlerPenalty<Double>"
               << std::endl;
      }
      else if ( eType == COMPLEX ) {
        retVal = New IDBC_HandlerPenalty<Complex>( numIDBC, blockSize );
        AssertMem( retVal, sizeof(IDBC_HandlerPenalty<Complex>) );
        (*cla) << " GenerateIDBC_HandlerObject: Generated"
               << " IDBC_HandlerPenalty<Complex>"
               << std::endl;
      }
      else {
        (*error) << "GenerateIDBC_HandlerObject: Cannot generate "
                 << "IDBC_HandlerPenalty<T> object with T = '"
                 << Enum2String( eType ) << "'";
        Error( __FILE__, __LINE__ );
      }

      // Force instantiation of public member functions of IDBC_HandlerPenalty
      // class. Note that the InstantiatePublicMethods() method itself should
      // never actually be called, but the compiler must not notice this.
      // Thus, the akward if condition.
      if ( blockSize == 0 ) {
        Error( "GenerateIDBC_HandlerObject: Internal error!", __FILE__,
               __LINE__ );
        if ( eType == DOUBLE ) {
          dynamic_cast< IDBC_HandlerPenalty<Double>* >(retVal)
            ->InstantiatePublicMethods();
        }
        else {
          dynamic_cast< IDBC_HandlerPenalty<Complex>* >(retVal)
            ->InstantiatePublicMethods();
        }
      }
    }

    return retVal;
  }


  // ******************************
  //   GenerateIDBC_HandlerObject
  // ******************************
  BaseIDBC_Handler*
  GenerateIDBC_HandlerObject( UInt numIDBC, UInt numPDEs,
                              UInt *bcOffsets, const MatrixEntryType eType ) {

    // VERSION:
    //
    // - This version of the factory generates objects of type
    //   IDBC_HandlerPenalty, i.e. objects that use the penalty
    //   approach to IDBC handling.
    //
    // - It works only for the SBM_System class

    ENTER_FCN( "GenerateIDBC_HandlerObject" );

    BaseIDBC_Handler *retVal = NULL;

    if ( eType == DOUBLE ) {
      retVal = New IDBC_HandlerPenalty<Double>( numIDBC, numPDEs, bcOffsets );
      AssertMem( retVal, sizeof(IDBC_HandlerPenalty<Double>) );
      (*cla) << " GenerateIDBC_HandlerObject: Generated"
             << " IDBC_HandlerPenalty<Double>"
             << std::endl;
    }
    else if ( eType == COMPLEX ) {
      retVal = New IDBC_HandlerPenalty<Complex>( numIDBC, numPDEs, bcOffsets);
      AssertMem( retVal, sizeof(IDBC_HandlerPenalty<Complex>) );
      (*cla) << " GenerateIDBC_HandlerObject: Generated"
             << " IDBC_HandlerPenalty<Complex>"
             << std::endl;
    }
    else {
      (*error) << "GenerateIDBC_HandlerObject: Cannot generate "
               << "IDBC_HandlerPenalty<T> object with T = '"
               << Enum2String( eType ) << "'";
      Error( __FILE__, __LINE__ );
    }

    // Force instantiation of public member functions of IDBC_HandlerPenalty
    // class. Note that the InstantiatePublicMethods() method itself should
    // never actually be called, but the compiler must not notice this.
    // Thus, the akward if condition.
    if ( numPDEs == 0 ) {
      Error( "GenerateIDBC_HandlerObject: Internal error!", __FILE__,
             __LINE__ );
      if ( eType == DOUBLE ) {
        dynamic_cast< IDBC_HandlerPenalty<Double>* >(retVal)
          ->InstantiatePublicMethods();
      }
      else {
        dynamic_cast< IDBC_HandlerPenalty<Complex>* >(retVal)
          ->InstantiatePublicMethods();
      }
    }

    return retVal;
  }

}
