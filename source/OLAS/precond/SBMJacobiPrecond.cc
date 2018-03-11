#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/precond/SBMJacobiPrecond.hh"
#include "MatVec/SBM_Matrix.hh"
#include "MatVec/SBM_Vector.hh"

namespace CoupledField {

  
  SBMJacobiPrecond::SBMJacobiPrecond( UInt numBlocks,
                                  PtrParamNode olasInfo,
                                  bool isMultiHarmonic)
  : BaseSBMPrecond(numBlocks, olasInfo) {
    if(isMultHarm_){
      stdPreconds_.Resize( 1 );
    }else{
      stdPreconds_.Resize( numBlocks_ );
    }
    stdPreconds_.Init(NULL);

    setMHPrecond_ = false;
    isMultHarm_ = isMultiHarmonic;
  }
  

  SBMJacobiPrecond::~SBMJacobiPrecond() {
    if(isMultHarm_){
      if(stdPreconds_[0] != NULL) delete stdPreconds_[0];
    }else{
      for( UInt i = 0; i < numBlocks_; ++i) {
          if( stdPreconds_[i] != NULL ) {
            delete stdPreconds_[i];
          }
        }
      }
    stdPreconds_.Clear();
  }
  

  void  SBMJacobiPrecond::SetPrecond(UInt blockNum, BasePrecond* stdPrecond ) {
    if(isMultHarm_){
      // In the multiharmonic case this method is only called for blockNum 0
      // and isn't allowed to be called one more time, therefore check this
      if(setMHPrecond_) EXCEPTION("SBMJacobiPrecond was already set for the multiharmonic case");
      stdPreconds_[blockNum] = stdPrecond;
      setMHPrecond_ = true;
    }else{
      // Non-multiharmonic case
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

  }

  void  SBMJacobiPrecond::Apply( const SBM_Matrix &A, const SBM_Vector &r,
                               SBM_Vector &z ) {
    
    assert( A.GetNumRows() <= numBlocks_ );
    // Loop over all rows
    for( UInt iRow = 0; iRow < A.GetNumRows(); ++iRow ) {

      if(isMultHarm_){
        // then we only have one stdPreconds_, since the
        // diagonal elements are all the same
        // If preconditioner for row is defined, apply it
        if( stdPreconds_[0] != NULL ) {
          const StdMatrix * stdMat = A.GetPointer(iRow,iRow);
          const SingleVector * rStd = r.GetPointer(iRow);
          SingleVector * zStd = z.GetPointer(iRow);
          stdPreconds_[0]->Apply(*stdMat, *rStd, *zStd );
        }
      }else{
        // If preconditioner for row is defined, apply it
        if( stdPreconds_[iRow] != NULL ) {
          const StdMatrix * stdMat = A.GetPointer(iRow,iRow);
          const SingleVector * rStd = r.GetPointer(iRow);
          SingleVector * zStd = z.GetPointer(iRow);
          stdPreconds_[iRow]->Apply(*stdMat, *rStd, *zStd );
        }
      }

      // If preconditioner for row is defined, apply it
      if( stdPreconds_[iRow] != NULL ) {
        const StdMatrix * stdMat = A.GetPointer(iRow,iRow);
        const SingleVector * rStd = r.GetPointer(iRow);
        SingleVector * zStd = z.GetPointer(iRow);
        stdPreconds_[iRow]->Apply(*stdMat, *rStd, *zStd );
      }
    }
  }

  void  SBMJacobiPrecond::Setup( SBM_Matrix &A ) {
    assert( A.GetNumRows() <= numBlocks_ );
    if(isMultHarm_){
      if( stdPreconds_[0] != NULL ) {
        stdPreconds_[0]->Setup(A(0,0));
      }
    }else{
      for( UInt iRow = 0; iRow < A.GetNumRows(); ++iRow ) {
        // If preconditioner for row is defined, apply it
        if( stdPreconds_[iRow] != NULL ) {
          stdPreconds_[iRow]->Setup(A(iRow,iRow));
        }
      }
    }

  }    

  void  SBMJacobiPrecond::GetPrecondSysMat( BaseMatrix& sysMat ) {
    if(isMultHarm_) EXCEPTION("SBMJacobiPrecond::GetPrecondSysMat Not adapted for multiharmonic yet!");
    SBM_Matrix& sbmA = dynamic_cast<SBM_Matrix&>(sysMat);
       // loop over all diagonal
       for( UInt iRow = 0; iRow < numBlocks_; ++iRow ) {
         for( UInt iCol = iRow; iCol < numBlocks_; ++iCol ) {

           if( iCol == iRow ) {
             if( stdPreconds_[iRow] != NULL ) {
               // get preconditioned sysmat from diagonal precond
               stdPreconds_[iRow]->GetPrecondSysMat( sbmA(iRow,iRow) );
             }
           } else { 
             sbmA(iRow,iCol).Init();
           }
         }
       }
  }


}
