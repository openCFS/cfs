#ifndef LAPACK_LL_HH
#define LAPACK_LL_HH

#include "OLAS/solver/BaseSolver.hh"

namespace CoupledField {

  class StdMatrix;
  

  //! This class implements the interface to LAPACK's Cholesky decomposition.

  //! This class implements the interface to LAPACK's Cholesky decomposition.
  //! It is designed to work with the SCRS_Matrix class and makes use of the
  //! banded storage format of LAPACK. The factorisation / solution process
  //! consists of the following steps:
  //! - <em>Setup phase:</em>
  //!   A Cholesky factorisation of the system matrix is computed using xPBTRF.
  //! - <em>Solution phase:</em>
  //!   The solution of the linear system is computed by backward-forward
  //!   substitution using the factored matrix and xPBTRS.
  //!
  //! Here 'x' is either 'D' in the case of a solver for a problem with matrix
  //! entries of type DOUBLE or 'Z' for COMPLEX entries. This direct solver
  //! does not require any steering parameters for the factorisation process
  //! itself, but only for certain management tasks. Consequently the only two
  //! parameters the Lapack_LL solver reads from the myParams_ object are the
  //! following:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for LAPACKLL solver</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>LAPACKLL_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>If this is 'true' then the Setup() and Solve() methods of the
  //!           Lapack_LL class will be verbose and write some information to
  //!           the standard logfile (*cla) <b>(removed!)</b>.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>newMatrixPattern</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If Setup() is called, and an old factorisation exists, the
  //!           old factorisation is deleted. However, the method by default
  //!           assumes that the sparsity pattern of the matrix has not
  //!           changed, but only the values of the matrix entries have. If
  //!           this parameter is set to 'true', then also the bandwidth will
  //!           be re-computed and new memory for the band matrix will be
  //!           allocated, if the capacity of facMat_ should no longer be
  //!           sufficient.
  //!     </tr>
  //!   </table>
  //! </center>
  //!
  //! \note
  //! - Note that this factorisation is <b>not done in place</b>.
  //!   The reason for this lies in the fact that we must convert the matrix
  //!   from a sparse to a banded representation, before passing it on to
  //!   LAPACK.
  //! - The class stores its own copy of the factorised matrix. Thus, it is
  //!   safe to alter or destroy the system matrix after Setup() has been
  //!   called. However, this also means that subsequent changes to the system
  //!   matrix will not affect the Solve() unless Setup() is called again.
  //! - The amount of memory required to store the factorised version of the
  //!   \f$n\times n\f$ matrix with a bandwidth of \f$w\f$ is
  //!   \f$(w+1)n\f$ entries.
  //! - For the sake of efficiency Setup() will only compute the bandwidth
  //!   of a matrix, if there is no prior factorisation, or if the parameter
  //!   LAPACKLL_newMatrixPattern is true.
  //!
  //! \htmlonly
  //! <center><img src="../AddDoc/lapack.jpg"></center>
  //! \endhtmlonly
  //! For further information see the
  //! <a href="http://www.netlib.org/lapack">LAPACK</a> pages at Netlib.
  class Lapack_LL : public BaseDirectSolver {

  public:

    //! Constructor

    //! The constructor does not do anything but initialising the internal
    //! attributes to default values and setting the pointers to the internal
    //! communication objects. This is the constructor which will be called by
    //! the GenerateSolverObject() factory.
    Lapack_LL( PtrParamNode solverNode, PtrParamNode olasInfo = PtrParamNode() );

    //! Default Destructor

    //! The default destructor needs to be deep in order to free the memory
    //! dynamically allocated for the (factorised) system matrix and the
    //! right-hand side vector.
    ~Lapack_LL();

    //! Computation of the Cholesky factorisation

    //! The setup method takes care of the Cholesky factorisation of the
    //! problem matrix. Depending on the type of matrix entry the factorisation
    //! is delegated to the appropriate private method. Note that the
    //! factorisation over-writes the internal copy of the system matrix.
    void Setup( BaseMatrix &sysmat );

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Direct solution of the linear system

    //! The solve method takes care of the solution of the linear system.
    //! The solution is computed via a backward/forward substitution pair
    //! using the Cholesky factorisation of the problem matrix computed in the
    //! setup phase.
    //! Note that the method will neglect the sysMat input parameter.
    void Solve( const BaseMatrix &sysmat, 
		const BaseVector &rhs, BaseVector &sol );

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return LAPACK_LL
    SolverType GetSolverType() {
      return LAPACK_LL;
    }

  private:

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    Lapack_LL();

    //! Factorisation for real matrices in double precision

    //! This method is called by the Setup routine for the case that we are
    //! dealing with a real matrix in double precision. It employs the
    //! LAPACK routine DPBTRF for computing the Cholesky factorisation.
    void FactoriseReal( StdMatrix &stdMat );

    //! Factorisation for complex matrices in double precision

    //! This method is called by the Setup routine for the case that we are
    //! dealing with a complex matrix in double precision. It employs the
    //! LAPACK routine ZPBTRF for computing the Cholesky factorisation.
    void FactoriseComplex( StdMatrix &stdMat );

    //! Solution method for real matrices in double precision

    //! This method is called by the Solve routine for the case that we are
    //! dealing with a real matrix in double precision. It calls the LAPACK
    //! routine DPBTRS for the actual work.
    void SolveReal( const BaseVector &rhs, BaseVector &sol );

    //! Solution method for complex matrices in double precision

    //! This method is called by the Solve routine for the case that we are
    //! dealing with a complex matrix in double precision. It calls the LAPACK
    //! routine ZPBTRS for the actual work.
    void SolveComplex( const BaseVector &rhs, BaseVector &sol );

    //! Factorised system matrix

    //! This array contains the factorised system matrix. Since the LAPACK
    //! factorisation routine requires additional super-diagonals (which the
    //! scaling routines do not expect) we cannot perform the factorisation
    //! in place, even if the original matrix offered the additional space.
    Double *facmat_;

    //! This attribute keeps track of the length of the facMat_ array.

    //! This attribute keeps track of the length of the facMat_ array. If
    //! repeated factorisations are computed, the internal array is only
    //! re-allocated, if the capacity is no longer sufficient to hold the
    //! new band matrix.
    UInt facMatCapacity_;

    //! This attribute keeps track of the number of entries in the band matrix

    //! This attribute keeps track of the number of entries in the band version
    //! of the sparse matrix resp. in its factorisation.
    UInt facMatEntries_;

    //! This array is used as right-hand side vector for the LAPACK routine

    //! This right-hand side of the linear system is over-written by the
    //! solution routine in LAPACK. We therefore provide an array into which
    //! the right-hand side values are copied, before they are passed to LAPACK
    //! and from which the solution values are taken in turn.
    Double *lapackRHS_;

    //! This attribute keeps track of the length of the lapackRHS_ array.

    //! This attribute keeps track of the length of the lapackRHS_ array. If
    //! repeated factorisations/solutions are computed, the internal array is
    //! only re-allocated, if the capacity is no longer sufficient to hold the
    //! new right-hand side.
    UInt lapackRHSCapacity_;

    //! This attribute stores the bandwidth of the system matrix
    UInt sysMatBW_;

    //! This attribute keeps track of status of factorisation

    //! The factorisation routine will set this boolean to yes, once it has
    //! performed a factorisation. This gives us a comfortable way of keeping
    //! track of the status of the solver and the internal memory allocation.
    bool amFactorised_;

    //! Factorised matrix

    //! This array is used to store the Cholesky factorisation of the system
    //! matrix. We first generate a LAPACK band matrix from the system matrix
    //! in this array, which in turn is over-written by LAPACK's factorisation
    //! routine D/ZPBTRF. We use a double pointer for the interface. In the
    //! complex case it is of double length with each real part followed
    //! directly by the imaginary one.
    Double *facMat_;

  };

}

#endif
