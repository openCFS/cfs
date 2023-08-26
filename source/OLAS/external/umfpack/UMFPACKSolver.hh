#ifndef UMFPACK_SOLVER_HH
#define UMFPACK_SOLVER_HH

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/StdVector.hh"
#include "OLAS/solver/BaseSolver.hh"

namespace CoupledField {

  //! This class implements the interface to PARDISO's LU decomposition
  //! and solving of a sparse system of linear equations.
  //! There is no need to initialise it, the Setup-Method does in fact
  //! nothing at all.
  //! After Solve is called the BaseMatrix (expected to be either a CRS
  //! or SCRS-matrix) will be reordered using the Nested Dissection or
  //! the Minimum Degree Algorithm, then it will be LU-factorized and
  //! finally solved by backward-forward substitution.
  //!
  //! The following parameters are read from the myParams-object:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for PARDISO solver</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>PARDISO_definite</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If this is set to 'true' the matrix is expected to be
  //!         positive definit.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>PARDISO_hermitian</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If this is set to 'true' the matrix is expected to be
  //!         complex and hermitian.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>PARDISO_symstructure</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If this is set to 'true' the matrix is expected to be
  //!         symmetric in structure.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>PARDISO_ordering</td>
  //!       <td align="center">NESTED_DISSECTION, MINUMUM_DEGREE or
  //!         NOREORDERING</td>
  //!       <td align="center">NESTED_DISSECTION</td>
  //!       <td>This parameter determines the strategy for computing the
  //!           re-ordering of unknowns in the symbolic factorisation.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>PARDISO_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If 'true' the %UMFPACKSolver object will log some minimal
  //!           status reports to the standard %OLAS logfile.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>PARDISO_stats</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>Setting this parameter to 'true' will triggers the pardiso
  //!           library to print statistical information on the solution
  //!           process to the screen (msgLvl_ = 1).
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //!
  //! \note
  //! - The statistical information that Pardiso prints to the screen for
  //!   PARDISO_stats = true, is also contained in the output parts of iparm_.
  //!   Thus, we could also write our own method that re-directs the statistics
  //!   to the standard logfile.
  template<typename T>
  class UMFPACKSolver : public BaseDirectSolver {

  public:

    //! Constructor

    //! The constructor does not do anything but set the pointers to the
    //! internal communciation objects.
    //! Note that the pointer to the report object is optional.
    UMFPACKSolver (PtrParamNode solverNode, PtrParamNode olasInfo);

    //! Default Destructor

    //! The destructor is not used at all.
    ~UMFPACKSolver();

    //! Reordering and factorization of the linear system

    //! After Setup is called the BaseMatrix (expected to be either a CRS
    //! or SCRS-matrix) will be reordered using the Nested Dissection or
    //! the Minimum Degree Algorithm and then it will be LU-factorised.
    void Setup( BaseMatrix &sysmat);

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Direct solution of the linear system

    //! After Solve is called the matrix (which has already to be factorised
    //! by a call of Setup) is finally solved by backward-forward substitution.
    void Solve( const BaseMatrix &sysmat,
                const BaseVector &rhs, BaseVector &sol );

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return PARDISO
    SolverType GetSolverType() {
      return UMFPACK;
    }

  private:

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    UMFPACKSolver();

    //! Handle for symbolic factorization data of UMFPACK
    void *Symbolic;

    //! Handle for numeric factorization data of UMFPACK
    void *Numeric ;

    //@{
    //! Vectors containing the information about the CRS- or SCRS-matrix

    //! This pointer is used to hold the address of a part of the internal
    //! (S)CRS matrix structures. The related memory segment must not
    //! be altered of deleted. Therefore the pointer is const!
    Integer *Ap, *Ai;

    //! This array contains the entries of the problem matrix of the linear
    //! system that is to be solved. We implicitely assume that a Complex can
    //! directly be cast into an array of two Doubles.
    Double* Ax;

    //! Index arrays for system matrix in compressed column storage (CSC)
    //! format. Output from umfpack_*i_transpose of Ap, Ai.  
    std::vector<Integer> Rp, Ri;

    //! Data array for system matrix in compressed column storage (CSC)
    //! format. Output from umfpack_*i_transpose of Ax.  
    std::vector<Double> Rx;

    //@}

    //! Stored information about the storage type and entry type of the matrix
    BaseMatrix::StorageType stype;
    BaseMatrix::EntryType etype;

    //! Printing level of UMFPACK
    Integer printingLevel_;

    //! Print statistics while solving.
    bool printStats_;

    Double tol_;    
    Double symtol_;
    Double irstep_;

    //! The type of the matrix in a special encoding used by Pardiso
    int mType_;

    //! Dimension of the linear system
    int probDim_;

    //! The number of right hand sides Pardiso should solve the system for
    //! at one pass
    int nrhs_;

    //! Array containing entries of problem matrix

    //! A flag specifying if Setup is being called for the first time
    bool firstCall_;

    //! Do we solve a system with complex or double entries?
    bool isComplex_;

    //! number of non zero entries
    UInt nnz_;

    Integer status;

    std::vector<double> Info;
    std::vector<double> Control;
    
    Double strategy_;
    Double ordering_;
    
  };

}

#endif
