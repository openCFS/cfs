#ifndef OLAS_SBMDIAGPRECOND_HH
#define OLAS_SBMDIAGPRECOND_HH

namespace OLAS {

  //! Preconditioner for SBM Matrices working on the diagonal sub-matrices

  //! This is a diagonal preconditioner for super-block matrices. It is
  //! diagonal in the sense that it will apply to each diagonal sub-matrix
  //! of the SBM matrix a preconditioner.
  class SBMDiagPrecond : public BaseSBMPrecond {

  public:

    //! Default constructor
    
    //! The default constructor does not do anything, but set the pointer of
    //! the preconditioner array to NULL.
    SBMDiagPrecond() : precond_(NULL) {};

    //! Applies the preconditioner by "solving" Az=r for z

    //! This method applies the preconditioner. Formally this means that for
    //! the given vectors r and z the linear system Az=r with the problem
    //! matrix A is solved for z. In the case of the SBMDiagPrecond
    //! preconditioner one preconditioning step consists of applying the
    //! individual preconditioners to the diagonal sub-matrices of A using the
    //! corresponding sub-vectors of r and z.
    //! \param A problem matrix
    //! \param r residual vector for current iteration step
    //! \param z output vector computed by the preconditioner
    virtual void Apply( const SBM_Matrix &A, const SBM_Vector &r,
			SBM_Vector &z ) const;

    //! A call of this method triggers the construction of the preconditioner.

    //! When this method is called the preconditioner will be constructed.
    //! In order to achieve this the method will trigger the setup phase of
    //! the different preconditioners for the diagonal sub-matrices of A.
    //! \param A problem matrix
    virtual void Setup( SBM_Matrix &A );

  private:

    //! Array containing the preconditioners for the diagonal sub-matrices
    BaseStdPrecond *precond_;

    //! Number of columns in the super-block matrix to which we are applied
    Integer Ncols_;

    //! Number of rows in the super-block matrix to which we are applied
    Integer Nrows_;

  };

}

#endif
