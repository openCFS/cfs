// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/stdmatrix.hh"
#include "MatVec/sbmmatrix.hh"

#include "MatVec/basevector.hh"
#include "MatVec/sbmvector.hh"

#include "baseprecond.hh"

namespace CoupledField {

  static EnumTuple precondTypeTuples[] = 
  {
    EnumTuple( BasePrecond::NOPRECOND, "noPrecond" ),
    EnumTuple( BasePrecond::ID, "Id" ),
    EnumTuple( BasePrecond::MG, "MG"),
    EnumTuple( BasePrecond::JACOBI, "Jacobi"),
    EnumTuple( BasePrecond::BLOCK_JACOBI, "BlockJacobi"),
    EnumTuple( BasePrecond::SSOR, "SSOR" ),
    EnumTuple( BasePrecond::ILU0, "ILU0" ),
    EnumTuple( BasePrecond::ILUTP, "ILUTP"),
    EnumTuple( BasePrecond::ILUK, "ILUK"),
    EnumTuple( BasePrecond::ILDL0, "ILDL0" ),
    EnumTuple( BasePrecond::ILDLK, "ILDLK" ),
    EnumTuple( BasePrecond::ILDLTP, "ILDLTP"),
    EnumTuple( BasePrecond::ILDLCN, "ILDLCN"),
    EnumTuple( BasePrecond::IC0, "IC0" ),
  };

  Enum<BasePrecond::PrecondType> BasePrecond::precondType = \
  Enum<BasePrecond::PrecondType>("Preconditioner Types",
      sizeof(precondTypeTuples) / sizeof(EnumTuple),
      precondTypeTuples); 


  void BaseStdPrecond::Apply( const BaseMatrix &sysmat, const BaseVector &r, 
                           BaseVector &z ) const {
    const StdMatrix& stdsysmat = dynamic_cast<const StdMatrix&>(sysmat);
    const SingleVector& stdr = dynamic_cast<const SingleVector&>(r);
    SingleVector& stdz = dynamic_cast<SingleVector&>(z);

    Apply(stdsysmat,stdr,stdz);
  }
  
  void BaseStdPrecond::Setup( BaseMatrix &sysMat ) {
    Setup( dynamic_cast<StdMatrix&>(sysMat) );
  }

  // ------------------------------------------------------------------------
  //  S B M   -  P R E C O N D I T I O N E R
  // ------------------------------------------------------------------------
  
  BaseSBMPrecond::BaseSBMPrecond( UInt numBlocks ) {
    readyToUse_ = false;
    numBlocks_ = numBlocks;
    stdPreconds_.Resize( numBlocks_ );
    stdPreconds_.Init(NULL);
  }
  
  BaseSBMPrecond::~BaseSBMPrecond() {
    
    for( UInt i = 0; i < numBlocks_; ++i) {
      if( stdPreconds_[i] != NULL ) {
        delete stdPreconds_[i];
      }
    }
    stdPreconds_.Clear();
  }
  
  void BaseSBMPrecond::SetPrecond(UInt blockNum, BaseStdPrecond* stdPrecond ) {
    
    // check block index
    if( blockNum > numBlocks_ ) {
      EXCEPTION("Can not assign a standard preconditioner for block #" << blockNum
               << " as the system has only size of " << numBlocks_ );
    }
    
    // check if preconditioner was already set for given block
    if( stdPreconds_[blockNum] != NULL ) {
      EXCEPTION("Can not set preconditioner for SBM block #" << blockNum 
                << " as there is already a preconditioner of type "
                << precondType.ToString(stdPreconds_[blockNum]->GetPrecondType())
                << " assigned." );
    }
    stdPreconds_[blockNum] = stdPrecond;
    
  }
  
  void BaseSBMPrecond::Apply( const BaseMatrix& sysmat, const BaseVector& r, 
                              BaseVector& z ) const {
      const SBM_Matrix& sbmsysmat = dynamic_cast<const SBM_Matrix&>(sysmat);
      const SBM_Vector& sbmr      = dynamic_cast<const SBM_Vector&>(r);
      SBM_Vector& sbmz            = dynamic_cast<SBM_Vector&>(z);
      Apply(sbmsysmat,sbmr,sbmz);
    }
  
  void BaseSBMPrecond::Apply( const SBM_Matrix &A, const SBM_Vector &r,
                              SBM_Vector &z ) const {
    
    // Loop over all rows
    for( UInt iRow = 0; iRow < numBlocks_; ++iRow ) {
    
      // If preconditioner for row is defined, apply it
      if( stdPreconds_[iRow] != NULL ) {
        const StdMatrix * stdMat = A.GetPointer(iRow,iRow);
        const SingleVector * rStd = r.GetPointer(iRow);
        SingleVector * zStd = z.GetPointer(iRow);
        stdPreconds_[iRow]->Apply(*stdMat, *rStd, *zStd );
      }
    }
  }

  void BaseSBMPrecond::Setup( BaseMatrix &A ) {
    SBM_Matrix& sbmA = dynamic_cast<SBM_Matrix&>(A);
    Setup(sbmA);
  }
  
  void BaseSBMPrecond::Setup( SBM_Matrix &A ) {
    for( UInt iRow = 0; iRow < numBlocks_; ++iRow ) {
      // If preconditioner for row is defined, apply it
      if( stdPreconds_[iRow] != NULL ) {
        stdPreconds_[iRow]->Setup(A(iRow,iRow));
      }
    }
  }

}
