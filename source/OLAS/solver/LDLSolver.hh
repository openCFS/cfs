#ifndef LDLSOLVER_HH
#define LDLSOLVER_HH



#include "OLAS/utils/math/LDLSystemSolve.hh"
#include "OLAS/utils/math/IterativeRefinement.hh"

#include "OLAS/solver/BaseSolver.hh"

namespace CoupledField {

  template<typename> class SCRS_Matrix;

  //! This class implements a sparse direct solver

  //! This class implements the direct solution of a sparse linear system of
  //! equations. The solution is obtained by first computing a \f$LDL^T\f$
  //! factorisation of the system matrix, using a sparse Crout version of
  //! Gaussian Elimination, and then performing a pair of forward/backward
  //! substitutions.
  //! The LDLSolver reads the following parameters from the myParams_ object:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for LDLSolver object</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>LDLSOLVER_saveFacToFile</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If the value of this boolean parameter is true, then the
  //!           factorisation matrix \f$F=D+L^T\f$ computed in the
  //!           factorisation phase of the solver, will be exported to a file
  //!           in MatrixMarket format. This IO operation will take place
  //!           immediately after the end of the factorisation phase.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>LDLSOLVER_facFileName</td>
  //!       <td align="center">string</td>
  //!       <td align="center">-</td>
  //!       <td>File name for exporting the factorisation matrix. This value
  //!           is only required, when LDLSOLVER_saveFacToFile is 'true'.
  //!           </td>
  //!     </tr>
  //!     <tr>
  //!       <td>LDLSOLVER_facPatternOnly</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If the value of this boolean parameter is true, then only the
  //!           sparsity pattern of the  factorisation matrix \f$F=D+L^T\f$
  //!           will be exported at the end of the analyse phase.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>LDLSOLVER_itRefSteps</td>
  //!       <td align="center">&gt;= 0</td>
  //!       <td align="center">2</td>
  //!       <td>This parameter determines the number of iterative refinement
  //!           steps that are to be performed to improve the accuray of the
  //!           computed approximate solution.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>LDLSOLVER_itRefVerbosity</td>
  //!       <td align="center">&gt;= 0</td>
  //!       <td align="center">2</td>
  //!       <td>This parameter determines the number of iterative refinement
  //!           steps that are to be performed to improve the accuray of the
  //!           computed approximate solution.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>LDLSOLVER_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>This parameter is used to control the verbosity of the solver.
  //!           If it is set to 'true', then the object will log some
  //!           information on the solution process to the standard %OLAS
  //!           report stream (*cla) ---removed logging---.
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
  class LDLSolver : public BaseDirectSolver, public LDLSystemSolve<T> {

  public:

    //! Constructor
    LDLSolver( PtrParamNode solverNode, PtrParamNode olasInfo );

    //! Default Destructor

    //! The default destructor needs to be deep in order to free the memory
    //! dynamically allocated for the vector of pivot indices and the copy of
    //! the system matrix.
    ~LDLSolver();

    //! Computation of the LDL factorisation

    //! The setup method takes care of the LDL factorisation of the problem
    //! matrix.
    void Setup( BaseMatrix &sysMat );

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Direct solution of the linear system

    //! The solve method takes care of the solution of the linear system.
    //! The solution is computed via a backward/forward substitution pair
    //! using the LDL factorisation of the problem matrix computed in the
    //! setup phase.
    //! Note that the sysmat input
    //! parameter will only be used, when an iterative refinement is performed.
    void Solve( const BaseMatrix &sysMat,
                const BaseVector &rhs, BaseVector &sol );

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return LDL_SOLVER
    SolverType GetSolverType() {
      return LDL_SOLVER;
    }

  private:

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    LDLSolver();

    //! Analyse

    void Analyse( const SCRS_Matrix<T> &sysMat );

    //! Factorise

    void Factorise( const SCRS_Matrix<T> &sysMat );

    //! Keep track on status of factorisation

    //! The Setup routine will set this boolean to true, once it has computed
    //! a factorisation. This gives us a confortable way of keeping track of
    //! the status of the solver and the internal memory allocation.
    bool amFactorised_;

    //! Dense auxilliary vector for factorisation
    T *denseVec_;

    //! Length of dense auxilliary vectors
    UInt auxVecSize_;

    //! Dimension of problem

    //! This attribute stores the dimension of the problem, i.e. the order
    //! of the problem matrix
    UInt sysMatDim_;

    //! Auxilliary index vector for factorisation

    //! In this index vector we store for every row of the matrix the index
    //! (in the cidxU_/dataU_ array) of the first non-zero entry with column
    //! index larger or equal to the current value of loop counter of the
    //! factorisation.
    UInt *firstU_;

    //! List with indices of rows to be scanned

    //! This list contains the indices of those rows that must be scanned
    //! when updating the auxilliary index vectors after each elimination step.
    //! The entries in the list are in ascending row index order.
    Integer *scanList_;

    //! List with indices of rows contributing in current elimination step

    //! List with indices of rows contributing in current elimination step.
    //! The entries in the list are in ascending row index order and
    //! activeList_ is a subset of scanList_.
    Integer *activeList_;

    //! Inverses of values of factor \f$D\f$

    //! In this vector the entries of the diagonal factor \f$D\f$ are stored.
    //! To be more precise, it is their inverses, since for performance reasons
    //! we already compute the inverses in the factorisation step. Divisions
    //! can be much more costly than mutliplications on some architectures.
    T* dataD_;

    //! Values of factor \f$U=L^T\f$
    T* dataU_;

    //! Values of factor \f$U=L^T\f$
    UInt* cidxU_;

    //! Values of factor \f$U=L^T\f$
    UInt* rptrU_;

    //! Object for performing iterative refinement
    IterativeRefinement iterativeRefiner_;

    //! Export LDL factorisation matrix to a file in MatrixMarket format

    //! The method will export the factorisation matrix to an ascii file
    //! according to the MatrixMarket specifications. By factorisation matrix
    //! we understand the matrix \f$F=D+L^T\f$.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname       name of output file
    //! \param patternOnly if true, only the sparsity pattern is exported
    void ExportFactorisation( const char *fname, bool patternOnly = false );

    UInt itRefSteps_;
  };

}

#endif
