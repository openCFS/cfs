#ifndef OLAS_BASEPRECOND_HH
#define OLAS_BASEPRECOND_HH

#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "OLAS/algsys/AlgebraicSys.hh"

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
    typedef enum { 
      NOPRECOND, ID, 

      // === Classical Preconditioners ==
      MG, JACOBI, BLOCK_JACOBI, SSOR, ILU0,ILUTP, 
      ILUK, ILDL0, ILDLK, ILDLTP, ILDLCN, IC0,
      SBM_DIAG, SBM_JACOBI,

      // === Solvers used as Preconditioners ==
      DIRECT, RICHARDSON, CG, LANCZOS, QMR, GMRES,
      MINRES, SYMMLQ, LAPACK_LU, LAPACK_LL, PARDISO_PRECOND, PARDISO_FACTREUSE,
      LU_SOLVER, CHOLMOD,
      LDL_SOLVER, LDL_SOLVER2, DIAGSOLVER} 
    PrecondType;
    static Enum<PrecondType> precondType;
    
  public:

    //! Default Constructor
    BasePrecond()  {
      readyToUse_ = false;
    };

    //! Default Destructor^
    virtual ~BasePrecond() {
    };

    //! Does constructor stuff only possible after child constructors are called 
    virtual void PostInit();
    
    //! A call of this method triggers the construction of the preconditioner.

    //! When this method is called the preconditioner will be constructed. This
    //! can involve a complex setup as e.g. in the case of an AMG
    //! preconditioner or be cheap as in the case of a Jacobi preconditioner.
    //! Note that only after this method has been called once, the precontioner
    //! can be applied.
    //! \param sysmat problem matrix
    virtual void Setup( BaseMatrix& sysmat) = 0;

    struct EdgeGeom{
      StdVector<Integer> eNodes; // edge nodes
      Double length;
    };

    //! A call of this method triggers the construction of the preconditioner, using
    //! algebraic multigrid (AMG).

    //! When this method is called the AMG-preconditioner will be constructed.
    //! This involves a complex setup (construction of hierarchy levels,
    //! transfer operators and solving the coarse system).
    //! Note that only after this method has been called once, the preconditioner
    //! can be applied.
    //! \param sysmat problem matrix
    //! \param auxmat auxiliary matrix
    //! \param amgType type of AMG-version (scalar, vectorial, edge)
    //! \param edgeIndNode connection of indices in the matrix and geometrical info
    //! \param nodeNumIndex connection of indices in the matrix and node-numbers
    virtual void SetupMG( BaseMatrix& sysmat,
                          BaseMatrix& auxmat,
                          const AMGType amgType,
                          const StdVector< StdVector< Integer> >& edgeIndNode,
                          const StdVector<Integer>& nodeNumIndex){}


    //! Applies the preconditioner by "solving" Az=r for z

    //! This method applies the preconditioner. Formally this means that for
    //! the given vectors r and z the linear system Az=r with the problem
    //! matrix A is solved for z.
    //! \param sysmat problem matrix
    //! \param r residual vector for current iteration step
    //! \param z output vector computed by the preconditioner
    virtual void Apply(const BaseMatrix& sysmat, const BaseVector& r, BaseVector& z) = 0;
    
    //! Mid-solve hook called by the host iterative solver, once per
    //! iteration, to ask whether it should abort and retry with a
    //! refreshed preconditioner. Returning true is a one-way commitment:
    //! the precond marks itself stale (next Setup() refactorises), and
    //! the host solver must abort, call Setup() on this precond, and
    //! restart the inner loop from the current iterate. A precond MUST
    //! return true at most once per host Solve() call.
    //! \param currentIter current iteration count of the host solver
    virtual bool ShouldAbortAndRefresh(UInt currentIter) {
      return false;
    }

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
    
    //! Return timer object for setup of preconditioner
    shared_ptr<Timer> GetSetupTimer() { return setupTimer_; }

    //! Return timer object for application of preconditioner
    shared_ptr<Timer> GetPrecondTimer() { return precondTimer_; }

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

    //! Pointer to timer object for setup of preconditioner
    shared_ptr<Timer> setupTimer_;

    //! Pointer to timer object for application of preconditioner
    shared_ptr<Timer> precondTimer_;

  };



  //! Generic Preconditioner for SBM Matrices
  
  //! General class for preconditioners acting on Super-Block matrices (SBM).
  //! In most cases the SBM-preconditioners are compound preconditioners, 
  //! which can consist of several instances of StdPreconditioners.
  class BaseSBMPrecond : public BasePrecond {

  public:

    //! Default constructor

    //! The constructor has nothing to do but to set the attribute
    //! readyToUse_ to false.
    //! \param numBlocks number of SBM-blocks
    //! \param olasInfo info paramNode for preconditioner 
    BaseSBMPrecond( UInt numBlocks, 
                    PtrParamNode olasInfo);

    //! Default Destructor
    virtual ~BaseSBMPrecond();

    //! Return pointer to info object
    PtrParamNode GetInfoNode() { return infoNode_;}
    

    //! Applies the preconditioner by "solving" Az=r for z

    //! This version of the Apply method has an interface fitting to
    //! SBM_Matrices and SBM_Vectors. It is purely virtual.
    virtual void Apply( const SBM_Matrix &A, const SBM_Vector &r,
                        SBM_Vector &z ) = 0;

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
    virtual void Setup( SBM_Matrix &A ) = 0;

    //! A call of this method triggers the construction of the preconditioner.

    //! This method implements the purely virtual Setup function inherited from
    //! the BasePrecond class. It does this by dynamically down-casting the
    //! input matrix and vectors from BaseMatrix/BaseVector type to the
    //! SBM_Matrix/SBM_Vector type and calling the Setup method with the
    //! corresponding interface. Thus, using this method with SBM matrices
    //! or vectors will lead to a run-time error.
    virtual void Setup( BaseMatrix &A );
    
    //! \copydoc BasePrecond::GetPrecondType
    virtual PrecondType GetPrecondType() const {return NOPRECOND;}

  protected:
   
    //! Number of sbm rows/blocks
    UInt numBlocks_; 
  };

} // namespace

#endif
