#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/precond/SBMDiagPrecond.hh"
#include "MatVec/SBM_Matrix.hh"
#include "MatVec/SBM_Vector.hh"

namespace CoupledField {

  
  SBMDiagPrecond::SBMDiagPrecond( UInt numBlocks, 
                                  PtrParamNode olasInfo ) 
  : BaseSBMPrecond(numBlocks, olasInfo) {
    stdPreconds_.Resize( numBlocks_ );
      stdPreconds_.Init(NULL);
  }
  

  SBMDiagPrecond::~SBMDiagPrecond() {
    for( UInt i = 0; i < numBlocks_; ++i) {
      if( stdPreconds_[i] != NULL ) {
        delete stdPreconds_[i];
      }
    }
    stdPreconds_.Clear();
  }
  

  void  SBMDiagPrecond::SetPrecond(UInt blockNum, BasePrecond* stdPrecond ) {
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

  void  SBMDiagPrecond::Apply( const SBM_Matrix &A, const SBM_Vector &r,
                               SBM_Vector &z ) {
    
    assert( A.GetNumRows() <= numBlocks_ );
    // Loop over all rows
    for( UInt iRow = 0; iRow < A.GetNumRows(); ++iRow ) {

      // If preconditioner for row is defined, apply it
      if( stdPreconds_[iRow] != NULL ) {
        const StdMatrix * stdMat = A.GetPointer(iRow,iRow);
        const SingleVector * rStd = r.GetPointer(iRow);
        SingleVector * zStd = z.GetPointer(iRow);
        stdPreconds_[iRow]->Apply(*stdMat, *rStd, *zStd );
      }
    }
  }

  void  SBMDiagPrecond::Setup( SBM_Matrix &A ) {
    assert( A.GetNumRows() <= numBlocks_ );
    for( UInt iRow = 0; iRow < A.GetNumRows(); ++iRow ) {
      // If preconditioner for row is defined, apply it
      if( stdPreconds_[iRow] != NULL ) {
        stdPreconds_[iRow]->Setup(A(iRow,iRow));
      }
    }
  }    

  void  SBMDiagPrecond::GetPrecondSysMat( BaseMatrix& sysMat ) {
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
