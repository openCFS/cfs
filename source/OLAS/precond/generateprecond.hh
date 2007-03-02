
/*! \file generateprecond.hh
 This module handles generation of preconditioner objects. It is also
 responsible for the instantiation of the templated preconditioners.
*/

#ifndef OLAS_GENERATEPRECOND_HH
#define OLAS_GENERATEPRECOND_HH

//!
namespace OLAS {

  //! Generate a standard preconditioner object

  //! This method will generate a standard preconditioner object that fits
  //! to the input matrix and return a pointer to that object.
  //! \param mat      %Matrix that is preconditioned
  //! \param ptype    Type of desired preconditoner
  //! \param myParams Pointer to parameter object for the preconditioner
  //! \param myReport Pointer to report object for the preconditioner
  BaseStdPrecond* GenerateStdPrecondObject( const StdMatrix &mat,
					    PrecondType ptype,
					    OLAS_Params *myParams,
					    OLAS_Report *myReport );

}

#endif
