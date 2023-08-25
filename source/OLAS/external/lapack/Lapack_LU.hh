#ifndef LAPACK_LU_HH
#define LAPACK_LU_HH

#include "OLAS/solver/BaseSolver.hh"
#include "LapackGBMatrix.hh"

namespace CoupledField {

  //! This class implements the interface to LAPACK's LU decomposition.

  //! This class implements the interface to LAPACK's LU decomposition.
  //! It is designed to work with the LapackGBMatrix, which uses the banded
  //! storage format of LAPACK. The factorisation / solution process consists
  //! of the following steps:
  //! - Setup phase:
  //!   -# Scaling of the matrix using xGBEQU and xLAQyy. This step is
  //!      optional and will only be performed, if LAPACKLU_tryScaling is
  //!      set to 'true'.
  //!   -# An LU factorisation of the (scaled) matrix is computed using
  //!      xGBTRF.
  //! - Solution phase:
  //!   -# The solution of the linear system is computed by backward-forward
  //!      substitution using the factored matrix and xGBTRS.
  //!   -# If LAPACKLU_refineSol is set to 'true' the obtained solution is
  //!      improved by iterative refinement using xGBRFS.
  //!
  //! Here 'x' is either 'D' in the case of a solver for a problem with matrix
  //! entries of type F77REAL8 or 'Z' for F77COMPLEX16 entries.
  //! Consequently the LAPACKLU solver reads the following parameters from the
  //! myParams_ object.
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for LAPACKLU solver</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>LAPACKLU_tryScaling</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>If this is set to 'true' we try to scale the matrix entries
  //!         before performing the factorisation in order to equilibrate the
  //!         matrix and reduce its condition number.
  //!         </td>
  //!     </tr>
  //!     <tr>
  //!       <td>LAPACKLU_refineSol</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>If this is set to 'true' we preform some steps of iterative
  //!         refinement after the backward/forward substitution pair.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>LAPACKLU_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>If this is 'true' then the Solve method of the Lapack_LU class
  //!           will be verbose and write some information on the solution
  //!           process to the standard logfile (*cla)<b>(removed!)</b>.</td>
  //!     </tr>
  //!   </table>
  //! </center>
  //!
  //! \note
  //! - If scaling is performed (for LAPACKLU_tryScaling = 'true') the system
  //!   matrix is changed.
  //! - Note that this factorisation is <b>not done in place</b>.
  //!   The reason for this lies in the different storage requirements of the
  //!   LAPACK routines used for scaling and factoring. The LAPACK
  //!   factorisation routine requires additional super-diagonals (which the
  //!   scaling routines do not expect). Thus, we cannot perform the
  //!   factorisation in place, even if the original matrix offered the
  //!   additional space, since then scaling would not be possible.
  //! - The class stores its own copy of the factorised matrix. Thus, it is
  //!   safe to alter or destroy the system matrix after Setup has been called.
  //!   However, this also means that subsequent changes to the system matrix
  //!   will not affect the Solve unless Setup is called again. Note that the
  //!   above is only true, if no iterative refinement is performed. In the
  //!   latter case Solve requires the initial system matrix.
  //! - The amount of memory required to store the factorised version of the
  //!   \f$n\times n\f$ matrix with \f$w_u\f$ upper and \f$w_l\f$ lower
  //!   diagonals is \f$(2\,w_l + 1 + w_u)\cdot n\f$ entries.
  //!
  //! \htmlonly
  //! <center><img src="../AddDoc/lapack.jpg"></center>
  //! \endhtmlonly
  //! For further information see the
  //! <a href="http://www.netlib.org/lapack">LAPACK</a> pages at Netlib.
  class Lapack_LU : public BaseDirectSolver {

  public:

    //! Constructor

    //! The constructor does not do anything but initialise the vector of pivot
    //! indices and the copy of the system matrix to NULL and set the pointers
    //! to the internal communication objects. Note that, if the solver object
    //! is generated via this contructor, then the Setup method must be called
    //! before Solve can be used. Also note that the pointer to the report
    //! object is optional. If no object is assigned, then no report will be
    //! generated. This is the constructor which will be called by the
    //! GenerateSolverObject() factory.
    Lapack_LU( PtrParamNode solverNode, PtrParamNode olasInfo = PtrParamNode() );

    //! Alternate constructor

    //! The alternate constructor is designed to pass the given problem
    //! matrix directly on to the Setup routine, which will then perform
    //! the basic LU factorisation of the matrix. Note that the pointer to the
    //! report object is optional. If no object is assigned, then no report
    //! will be generated.
    Lapack_LU( BaseMatrix &mat, PtrParamNode solverNode, 
	       PtrParamNode olasInfo = PtrParamNode() );

    //! Default Destructor

    //! The default destructor needs to be deep in order to free the memory
    //! dynamically allocated for the vector of pivot indices and the copy of
    //! the system matrix.
    ~Lapack_LU();

    //! Computation of the LU factorisation

    //! The setup method takes care of the LU factorisation of the problem
    //! matrix. Since the factorisation is not done in place we can delegate
    //! the actual work to the private version of Setup() that expects a const
    //! reference.
    void Setup( BaseMatrix &sysmat);

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Direct solution of the linear system

    //! The solve method takes care of the solution of the linear system.
    //! The solution is computed via a backward/forward substitution pair
    //! using the LU factorisation of the problem matrix computed in the
    //! setup phase.
    //! Note that the sysmat input
    //! parameter will only be used, when an iterative refinement is performed.
    void Solve( const BaseMatrix &sysmat,	const BaseVector &rhs, BaseVector &sol);

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return LAPACK_LU
    SolverType GetSolverType() {
      return LAPACK_LU;
    }

  private:

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    Lapack_LU();

    //! Computation of the LU factorisation

    //! This setup method takes care of the LU factorisation of the problem
    //! matrix. Note that this factorisation is <b>not done in place</b>.
    //! The reason for this lies in the different storage requirements of the
    //! LAPACK routines used for scaling and factoring. The LAPACK
    //! factorisation routine requires additional super-diagonals (which the
    //! scaling routines do not expect). Thus, we cannot perform the
    //! factorisation in place, even if the original matrix offered the
    //! additional space, since then scaling would not be possible.
    void PrivateSetup( const BaseMatrix &sysmat );

    //! Factorisation for real matrices in double precision

    //! This method is called by the Setup routine for the case that we are
    //! dealing with a real matrix in double precision. It employ the
    //! DGB-versions of the LAPACK routines, most importantly DGBTRF.
    void FactoriseF77REAL8( const StdMatrix &stdmat );

    //! Factorisation for complex matrices in double precision

    //! This method is called by the Setup routine for the case that we are
    //! dealing with a complex matrix in double precision. It employ the
    //! ZGB-versions of the LAPACK routines, most importantly ZGBTRF.
    void FactoriseF77COMPLEX16( const StdMatrix &stdmat );

    //! Solution method for real matrices in double precision

    //! This method is called by the Solve routine for the case that we are
    //! dealing with a real matrix in double precision. It employs the
    //! DGB-versions of the LAPACK routines, most importantly DGBTRS.
    void SolveF77REAL8( const BaseVector &rhs, BaseVector &sol,
			const BaseMatrix *mat );

    //! Solution method for real matrices in double precision

    //! This method is called by the Solve routine for the case that we are
    //! dealing with a complex matrix in double precision. It employs the
    //! DGB-versions of the LAPACK routines, most importantly ZGBTRS.
    void SolveF77COMPLEX16( const BaseVector &rhs, BaseVector &sol,
			    const BaseMatrix *mat );

    //! Vector containing the pivot indices

    //! This vector contains the indices of the pivot elements as computed
    //! during the LU factorisation. Note that contrary to our usual policy
    //! this vector is 0-based. The reason is that we never acess it ourselves
    //! but only pass the starting address of the memory block to the LAPACK
    //! routines.
    int *pivots_;

    //! Factorised system matrix

    //! This vector contains the factorised system matrix. Since the LAPACK
    //! factorisation routine requires additional super-diagonals (which the
    //! scaling routines do not expect) we cannot perform the factorisation
    //! in place, even if the original matrix offered the additional space.
    LapackBaseMatrix *facmat_;

    //! Keep track on status of factorisation

    //! The factorisation routine will set this boolean to yes, once it has
    //! performed a factorisation. This gives us a confortable way of keeping
    //! track on the status of the solver and the internal memory allocation.
    bool amFactorised_;

    //@{
    //! Workspace for LAPACK routines

    //! This array contains dynamically allocated workspace required for the
    //! LAPACK routines that perform iterative refinement.
    double *workspaceF77REAL8_;
    std::complex<double> *workspaceF77COMPLEX16_;
    int *workspaceInt_;
    //@}

    //@{
    //! \name Attributes related to equilibration of system matrix

    //! This array contains the scaling factors for the rows of the matrix
    double *row_scalings_;

    //! This array contains the scaling factors for the columns of the matrix
    double *col_scalings_;

    //! Stores the type of scaling LAPACK performed
    char scalingType_;

    //@}

  };

}

#endif
