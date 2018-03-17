#ifndef OLAS_SBMJACOBIPRECOND_HH
#define OLAS_SBMJACOBIPRECOND_HH


#include "BasePrecond.hh"
#include "BNPrecond.hh"

namespace CoupledField {

  //! Jacobi preconditioner for SBM Matrices

  //! This is a Jacobi preconditioner for super-block matrices.
  //! It will apply one specific preconditioner to every diagonal
  //! subblock (e.g. AMG, pardiso, LU, ...)
  class SBMJacobiPrecond : public BaseSBMPrecond {

  public:

    //! Default constructor
    SBMJacobiPrecond( UInt numBlocks,
                    PtrParamNode olasInfo);
    
    //! Destructor
    virtual ~SBMJacobiPrecond();

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

    bool setMHPrecond_;
  };

}

#endif
