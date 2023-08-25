#ifndef LUSOLVER_HH
#define LUSOLVER_HH



#include "OLAS/utils/math/CroutLU.hh"
#include "OLAS/utils/math/IterativeRefinement.hh"

#include "OLAS/solver/BaseSolver.hh"

namespace CoupledField {

  //! This class implements a sparse direct solver

  //! This class implements the direct solution of a sparse linear system of
  //! equations. The solution is obtained by first computing a LU decomposition
  //! of the system matrix, using a sparse Crout version of Gaussian
  //! Elimination, and then performing a pair of forward/backward
  //! substitutions.
  //! The LUSolver reads the following parameters from the myParams_ object:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for LUSolver object</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>CROUT_saveFacToFile</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If the value of this boolean parameter is true, then
  //!           factorisation matrix \f$F=L+U-I\f$ computed in the setup phase
  //!           of the solver, will be exported to a file in MatrixMarket
  //!           format. This IO operation will take place immediately after
  //!           the end of the setup phase.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>CROUT_facFileName</td>
  //!       <td align="center">string</td>
  //!       <td align="center">-</td>
  //!       <td>File name for exporting the factorisation matrix. This value
  //!           is only required, when CROUT_saveFacToFile is 'true'.
  //!           </td>
  //!     </tr>
  //!     <tr>
  //!       <td>LUSOLVER_itRefSteps</td>
  //!       <td align="center">&gt;= 0</td>
  //!       <td align="center">2</td>
  //!       <td>This parameter determines the number of iterative refinement
  //!           steps that are to be performed to improve the accuray of the
  //!           computed approximate solution.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>LUSOLVER_itRefVerbosity</td>
  //!       <td align="center">&gt;= 0</td>
  //!       <td align="center">2</td>
  //!       <td>This parameter determines the number of iterative refinement
  //!           steps that are to be performed to improve the accuray of the
  //!           computed approximate solution.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //!
  //! \note
  //! - The class stores its own copy of the factorised matrix. Thus, it is
  //!   safe to alter or destroy the system matrix after Setup has been called.
  //!   However, this also means that subsequent changes to the system matrix
  //!   will not affect the Solve unless Setup is called again.
  //! - Pivoting is currently not supported during the factorisation.
  //! - The class only works together with CRS_Matrices with scalar entries.
  template<typename T>
  class LUSolver : public BaseDirectSolver, public CroutLU<T> {

  public:

    //! Constructor
    LUSolver( PtrParamNode solverNode, PtrParamNode olasInfo = PtrParamNode() );

    //! Default Destructor

    //! The default destructor needs to be deep in order to free the memory
    //! dynamically allocated for the vector of pivot indices and the copy of
    //! the system matrix.
    ~LUSolver();

    //! Computation of the LU factorisation

    //! The setup method takes care of the LU factorisation of the problem
    //! matrix. Since the factorisation is not done in place we can delegate
    //! the actual work to the private version of Setup() that expects a const
    //! reference.
    void Setup( BaseMatrix& sysMat);

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Direct solution of the linear system

    //! The solve method takes care of the solution of the linear system.
    //! The solution is computed via a backward/forward substitution pair
    //! using the LU factorisation of the problem matrix computed in the
    //! setup phase.
    //! Note that the method will neglect the precond input parameter, since
    //! we perform a direct solution. Note also, that the sysmat input
    //! parameter will only be used, when an iterative refinement is performed.
    void Solve( const BaseMatrix &sysMat, const BaseVector &rhs, BaseVector &sol);

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return LU_SOLVER
    SolverType GetSolverType() {
      return LU_SOLVER;
    }

  private:

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    LUSolver();

    //! Dropping strategy for fill-in in an incomplete decomposition

    //! The CroutLU class requires that a derived class implements this method,
    //! which is used to drop fill-in and compute an incomplete LU
    //! decomposition. In the case of the LUSolver class we simply implement
    //! an empty method.
    //! \param k        index of current row / column
    //! \param vecZ     dense vector storing current row of factor U
    //! \param vecZFill vector containing indices of non-zero entries in vecZ
    //! \param vecW     dense vector storing current column of factor L
    //! \param vecWFill vector containing indices of non-zero entries in vecW
    void DropEntries( UInt k, std::vector<T> &vecZ,
		      std::vector<UInt> &vecZFill,
		      std::vector<T> &vecW,
		      std::vector<UInt> &vecWFill ) {};

    //! Compute maximal memory requirements of factorisation

    //! Blablubb Blablubb Blablubb Blablubb Blablubb Blablubb Blablubb Blablubb
    //! Blablubb Blablubb Blablubb Blablubb Blablubb Blablubb Blablubb Blablubb
    void ComputeRequiredMemory( UInt &maxColLengthL, UInt &maxRowLengthU,
				UInt &maxEntriesL, UInt &maxEntriesU,
				CRS_Matrix<T> &crsMat );

    //! Keep track on status of factorisation

    //! The Setup routine will set this boolean to true, once it has computed
    //! a factorisation. This gives us a confortable way of keeping track of
    //! the status of the solver and the internal memory allocation.
    bool amFactorised_;

    //! Object for performing iterative refinement
    IterativeRefinement iterativeRefiner_;
    
    UInt itRefSteps_;
    

  };

}

#endif
