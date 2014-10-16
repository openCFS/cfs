#ifndef OLAS_SBMDIAGPRECOND_HH
#define OLAS_SBMDIAGPRECOND_HH


#include "BasePrecond.hh"
#include "BNPrecond.hh"

namespace CoupledField {

  //! Preconditioner for SBM Matrices working on the diagonal sub-matrices

  //! This is a diagonal preconditioner for super-block matrices. It is
  //! diagonal in the sense that it will apply to each diagonal sub-matrix
  //! of the SBM matrix a preconditioner.
  class SBMDiagPrecond : public BaseSBMPrecond {

  public:

    //! Default constructor
    SBMDiagPrecond( UInt numBlocks, 
                    PtrParamNode olasInfo );
    
    //! Destructor
    virtual ~SBMDiagPrecond();

    //! Set preconditioner for standard matrices for block \blockNum

    //! This method allows to compose the SBM-preconditioner by single 
    //! standard-preconditioners for every row.
    virtual void SetPrecond(UInt blockNum, BasePrecond* stdPrecond );

    //! \copydoc BaseSBMPrecond(SBM_Matrix, SBM_Vector, SBM_Vector)
    virtual void Apply( const SBM_Matrix &A, const SBM_Vector &r,
                        SBM_Vector &z ) ;

    //! \copydoc BaseSBMPrecond::Setup(SBM_Matrix, PtrParamNode)
    virtual void Setup( SBM_Matrix &A );

    //! \copydoc BasePrecond::ExportPrecondSysMat
    virtual void GetPrecondSysMat( BaseMatrix& sysMat );

    //! \copydoc BasePrecond::GetPrecondType
    virtual PrecondType GetPrecondType() const {
      return SBM_DIAG;
    }

  private:
    //! Vector containing a preconditioner for every diagonal SBM block
       
    //! We store for every SBM-row a standard preconditioner
    StdVector<BasePrecond*> stdPreconds_;
  };

}

#endif
