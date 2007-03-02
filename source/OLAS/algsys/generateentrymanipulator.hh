#ifndef OLAS_GENERATEENTRYMANIPULATOR_HH
#define OLAS_GENERATEENTRYMANIPULATOR_HH


#include "utils/utils.hh"


/*! \file generateentrymanipulator.hh

    This file handles the creation of EntryManipulator objects,
    which is responsible for assembling the element matrices
    and hte Dirichlet values into the matrices and RHS.
    \remark Note that the documentation of the functions in this file
    cannot be found here, but in the description of the namespace
    OLAS to which they belong.
*/


namespace OLAS {


  class StdVector;
  class BaseEntryManipulator;

  //! Generates an object of a child of the BaseEntryManipulator class
  BaseEntryManipulator *GenerateEntryManipulatorObject( MatrixEntryType eType,
                                                        UInt bs );

}

#endif
