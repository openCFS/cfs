#include "algsys/generateentrymanipulator.hh"
#include "algsys/entrymanipulatordouble.cc"
#include "algsys/entrymanipulatorcomplex.cc"

#ifdef __INTEL_COMPILER
#include "matvec/vector.cc"
#include "matvec/parvector.cc"
#include "matvec/communicator.cc"
#endif


namespace OLAS {


  // =======================================================
  //   Function for generation of EntryManipulator objects
  // =======================================================
  BaseEntryManipulator *GenerateEntryManipulatorObject( MatrixEntryType eType,
                                                        UInt bs ) {

    ENTER_FCN( "GenerateEntryManipulatorObject()" );

    BaseEntryManipulator *retVal = NULL;


    // **************
    //  Scalar cases
    // **************
    if ( bs == 1 ) {

      // ---------------
      //  Scalar Double
      // ---------------
      if ( eType == DOUBLE ) {
        retVal = New EntryManipulatorDouble();
        if ( retVal == NULL ) {
          (*error) << "GenerateEntryManipulatorObject: "
                   << "Failed to generate EntryManipulatorDouble!";
          Error( __FILE__, __LINE__ );
        }
        else {
          (*cla) << " GenerateEntryManipulatorObject:"
                 << " Generated EntryManipulatorDouble"
                 << std::endl;
        }
      }

      // ----------------
      //  Scalar Complex
      // ----------------
      else if ( eType == COMPLEX ) {
        retVal = New EntryManipulatorComplex();
        if ( retVal == NULL ) {
          (*error) << "GenerateEntryManipulatorObject: "
                   << "Failed to generate EntryManipulatorComplex!";
          Error( __FILE__, __LINE__ );
        }
        else {
          (*cla) << " GenerateEntryManipulatorObject:"
                 << " Generated EntryManipulatorComplex"
                 << std::endl;
        }
      }
    }


    // *************
    //  Block cases
    // *************

    else if ( bs > 1 ) {

      (*error) << " GenerateEntryManipulatorObject:"
               << "Block size = " << bs << "! Please re-compile with "
               << "USE_TVMET = yes to enable support for dof blocking!";
      Error( __FILE__, __LINE__ );
    }


    // ********************
    //  Strange parameters
    // ********************
    if ( retVal == 0 ) {
      (*error) << "GenerateEntryManipulatorObject(): "
               << "Generation of Entry Manipulator object failed for '"
               << "eType = " << Enum2String( eType )
               << "' and blockSize = " << bs;
      Error( __FILE__, __LINE__ );
    }

    // Return what we got
    return retVal;
  }

}
