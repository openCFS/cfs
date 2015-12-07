// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <set>

/*! \file generateprecond.hh
 This module handles generation of preconditioner objects. It is also
 responsible for the instantiation of the templated preconditioners.
*/

#ifndef OLAS_GENERATEPRECOND_HH
#define OLAS_GENERATEPRECOND_HH

#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/BaseMatrix.hh"
#include "OLAS/precond/BasePrecond.hh"


//!
namespace CoupledField {

  class BasePrecond;
  class BaseSBMPrecond;
  class StdMatrix;
  class SBM_Matrix;
  
  class SolStrategy;
  
  //! Generate a standard preconditioner object
  
  //! This method will generate a standard preconditioner object that fits
  //! to the input matrix and return a pointer to that object.
  //! \param mat       %Matrix that is preconditioned
  //! \param precondId Identifier token of preconditioner to ge generated
  //! \param xml       Pointer to ParamNode of <precondList> section
  //! \param olasInfo  Base below "OLAS" in info.xml
  BasePrecond* GenerateStdPrecondObject( const StdMatrix &mat,
                                         const std::string& precondId,
                                         PtrParamNode solverNode,
                                         PtrParamNode olasInfo );

  
  //! Generate a SBM preconditioner object
  
  //! This method will generate a SBM preconditioner object that fits
  //! to the input matrix and return a pointer to that object.
  //! \param mat       %Matrix that is preconditioned
  //! \param precondId Identifier token of preconditioner to ge generated
  //! \param xml       Pointer to ParamNode of <precondList> section
  //! \param olasInfo  Base below "OLAS" in info.xml
  BaseSBMPrecond* GenerateSBMPrecondObject( const SBM_Matrix &mat,
                                            const std::string& precondId,
                                            PtrParamNode solverNode,
                                            PtrParamNode olasInfo );
  
  //! Return list of compatible matrix types

  //! This method returns a list of all matrix storage types the solver can 
  //! handle. 
  //! \note In case the solver can handle any type of sparse matrix (i.e.
  //! Krylov based solvers, requiring just matrix-vector operations),
  //! the return set will be empty!
  std::set<BaseMatrix::StorageType> 
  GetPrecondCompatMatrixFormats(BasePrecond::PrecondType );
  
  //! Return if preconditioner is capable of solving a SBM system
  bool IsPreconndSBMCapable(BasePrecond::PrecondType );

}

#endif
