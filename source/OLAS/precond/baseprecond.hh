// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASEPRECOND_HH
#define OLAS_BASEPRECOND_HH

#include "General/environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  class BaseMatrix;
  class StdMatrix;
  class SBM_Matrix;
  class BaseVector;
  class SBM_Vector;
  class SingleVector;

  //! Generic preconditioner class
  class BasePrecond {

  public:
    //! Type of preconditioner

    //! This enumeration data type describes the type of Preconditioner which
    //! is applied before solving the algebraic system.
    //! The enumeration contains the following values:
    //! - NOPRECOND
    //! - ID
    //! - MG
    //! - JACOBI
    //! - SSOR
    //! - ILU0
    //! - ILUTP
    //! - ILUK
    //! - ILDLK
    //! - ILDLTP
    //! - ILDLCN
    typedef enum { NOPRECOND, ID, 

                  // === Classical Preconditioners ==
                  MG, JACOBI, BLOCK_JACOBI, SSOR, ILU0,ILUTP, 
                  ILUK, ILDL0, ILDLK, ILDLTP, ILDLCN, IC0,
                  
                  // === Solvers used as Preconditioners ==
                  DIRECT, RICHARDSON, CG, LANCZOS, QMR, GMRES,
                  MINRES, SYMMLQ, LAPACK_LU, LAPACK_LL, PARDISO,  
                  ILUPACK, LU_SOLVER, CHOLMOD, 
                  LDL_SOLVER, LDL_SOLVER2, DIAGSOLVER} PrecondType;
    static Enum<PrecondType> precondType;
    
  public:

    //! Default Constructor
    BasePrecond()  {
      readyToUse_ = false;
    };

    //! Default Destructor^
    virtual ~BasePrecond() {
    };

    //! A call of this method triggers the construction of the preconditioner.

    //! When this method is called the preconditioner will be constructed. This
    //! can involve a complex setup as e.g. in the case of an AMG
    //! preconditioner or be cheap as in the case of a Jacobi preconditioner.
    //! Note that only after this method has been called once, the precontioner
    //! can be applied.
    //! \param sysmat problem matrix
    //! \param analysis_id references to the info/analysis/progress/step(/substep) 
    //!                    element with the "analysis_id" attribute 
    virtual void Setup( BaseMatrix& sysmat,
                        PtrParamNode analysis_id ) = 0;

    //! Applies the preconditioner by "solving" Az=r for z

    //! This method applies the preconditioner. Formally this means that for
    //! the given vectors r and z the linear system Az=r with the problem
    //! matrix A is solved for z.
    //! \param sysmat problem matrix
    //! \param r residual vector for current iteration step
    //! \param z output vector computed by the preconditioner
    virtual void Apply(const BaseMatrix& sysmat, const BaseVector& r, 
                       BaseVector& z) = 0;

    //! Export precontitioned matrix \f[C^{-1}A\f]
    
    //! This methods tries to export the preconditioner matrix \f[C^{-1}\f]. 
    //! This is inly possible in some rare cases, where we can explicitly
    //! calculate the preconditioner in matrix form (e.g. Jacobi, Id).
    //! \note The argument #sysMat gets overwritten upon method return
    //! \note If the calculation is not implemented by the derived class,
    //! we issue a warning an silently return an empty matrix.
    virtual void GetPrecondSysMat( BaseMatrix& sysMat );

    
    //! Query type of preconditioner object
    virtual PrecondType GetPrecondType() const {
      return NOPRECOND;
    }

  protected:

    //! Before the preconditioner can be applied its setup phase must be
    //! performed. This is achieved by a call to Setup.
    bool readyToUse_;

    
    //! Pointer to parameter object

    //! This is a pointer to a parameter object containing the steering
    //! parameters for this preconditioner.
    PtrParamNode xml_;

    //! Pointer to report object

    //! This is a pointer to a report object which the preconditioner can use
    //! to store general information about its performance or setup phase.
    PtrParamNode infoNode_;

  };


//  //! Generic preconditioner class for standard matrices (not SBM matrices)
//  class BaseStdPrecond : public BasePrecond{
//
//  public:
//
//    //! Default constructor
//
//    //! The constructor has nothing to do but to set the attribute
//    //! readyToUse_ to false.
//    BaseStdPrecond() : readyToUse_(false) {
//    };
//
//    //! Default Destructor
//    virtual ~BaseStdPrecond() {
//    };
//
//    //! @copydoc BasePrecond::Apply(const BaseMatrix, const BaseVector, BaseVector)
//    virtual void Apply( const BaseMatrix &sysmat, const BaseVector &r, 
//                        BaseVector &z );
//    
//    //! Applies the preconditioner by "solving" Az=r for z
//
//    //! This version of the Apply method has an interface fitting to
//    //! StdMatrices and SingleVectors. It is purely virtual.
//    virtual void Apply( const StdMatrix &sysmat, const SingleVector &r,
//                        SingleVector &z) = 0;
//
//    //! A call of this method triggers the construction of the preconditioner.
//
//    //! This method implements the purely virtual Setup function inherited from
//    //! the BasePrecond class. It does this by dynamically down-casting the
//    //! input matrix and vectors from BaseMatrix/BaseVector type to the
//    //! StdMatrix/SingleVector type and calling the Setup method with the
//    //! corresponding interface. Thus, using this method with SBM matrices
//    //! or vectors will lead to a run-time error.
//    virtual void Setup( BaseMatrix &sysMat, PtrParamNode analysis_id );
//                        
//    //! A call of this method triggers the construction of the preconditioner.
//     
//    //! This version of the Setup method has an interface fitting to
//    //! StdMatrices and SingleVectors. It is purely virtual.
//    virtual void Setup( StdMatrix &sysMat, PtrParamNode anaylsis_id ) = 0;
//
//  protected:
//
//    //! This attribute keeps track on whether the preconditioner was set up
//
//    //! Before the preconditioner can be applied its setup phase must be
//    //! performed. This is achieved by a call to Setup.
//    bool readyToUse_;
//
//  };


  //! Generic Preconditioner for SBM Matrices
  
  //! A preconditioner for SBM-matrices / vectors can be considered as a
  //! composition of standard preconditioners, i.e. for every SBM-row,
  //! we can have a different standard preconditioner.
  class BaseSBMPrecond : public BasePrecond {

  public:

    //! Default constructor

    //! The constructor has nothing to do but to set the attribute
    //! readyToUse_ to false.
    //! \param numBlocks number of SBM-blocks
    BaseSBMPrecond(UInt numBlocks, 
                   PtrParamNode olasInfo);

    //! Default Destructor
    virtual ~BaseSBMPrecond();

    //! Return pointer to info object
    PtrParamNode GetInfoNode() { return infoNode_;}
    
    //! Set preconditioner for standard matrices for block \blockNum
    
    //! This method allows to compose the SBM-preconditioner by single 
    //! standard-preconditioners for every row.
    virtual void SetPrecond(UInt blockNum, BasePrecond* stdPrecond );

    //! Applies the preconditioner by "solving" Az=r for z

    //! This version of the Apply method has an interface fitting to
    //! SBM_Matrices and SBM_Vectors. It is purely virtual.
    virtual void Apply( const SBM_Matrix &A, const SBM_Vector &r,
                        SBM_Vector &z );

    //! Applies the preconditioner by "solving" Az=r for z

    //! This method implements the purely virtual Apply function inherited from
    //! the BasePrecond class. It does this by dynamically down-casting the
    //! input matrix and vectors from BaseMatrix/BaseVector type to the
    //! SBM_Matrix/SBM_Vector type and calling the Apply method with the
    //! corresponding interface. Thus, using this method with StdMatrices
    //! or SingleVectors will lead to a run-time error.
    virtual void Apply( const BaseMatrix& sysmat, const BaseVector& r, 
                        BaseVector& z );

    //! A call of this method triggers the construction of the preconditioner.

    //! This version of the Setup method has an interface fitting to
    //! SBM_Matrices and SBM_Vectors. It is purely virtual.
    virtual void Setup( SBM_Matrix &A, PtrParamNode analysis_id );

    //! A call of this method triggers the construction of the preconditioner.

    //! This method implements the purely virtual Setup function inherited from
    //! the BasePrecond class. It does this by dynamically down-casting the
    //! input matrix and vectors from BaseMatrix/BaseVector type to the
    //! SBM_Matrix/SBM_Vector type and calling the Setup method with the
    //! corresponding interface. Thus, using this method with SBM matrices
    //! or vectors will lead to a run-time error.
    virtual void Setup( BaseMatrix &A, PtrParamNode analysis_id );
    
    //! \copydoc BasePrecond::ExportPrecondSysMat
    virtual void GetPrecondSysMat( BaseMatrix& sysMat );
    
    virtual PrecondType GetPrecondType() const {return NOPRECOND;}

  private:
    
    //! Vector containing a preconditioner for every SBM-row
    
    //! We store for every SBM-row a standard preconditioner
    StdVector<BasePrecond*> stdPreconds_;

    //! Number of sbm rows/blocks
    UInt numBlocks_; 
    
    //! This attribute keeps track on whether the preconditioner was set up

    //! Before the preconditioner can be applied its setup phase must be
    //! performed. This is achieved by a call to Setup.
    bool readyToUse_;

  };

} // namespace

#endif
