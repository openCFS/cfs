// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef PARDISO_SOLVER_HH
#define PARDISO_SOLVER_HH

#include <string>

#include "matvec/matvec.hh"
#include "solver/solver.hh"

namespace OLAS {

  //! This class implements the interface to PARDISO's LU decomposition
  //! and solving of a sparse system of linear equations.

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
  //!       <td>If 'true' the %PardisoSolver object will log some minimal
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
  class PardisoSolver : public BaseDirectSolver {

  public:

    //! Constructor

    //! The constructor does not do anything but set the pointers to the
    //! internal communciation objects.
    //! Note that the pointer to the report object is optional.
    PardisoSolver (OLAS_Params *myParams, OLAS_Report *myReport = NULL);

    //! Default Destructor

    //! The destructor is not used at all.
    ~PardisoSolver();

    //! Reordering and factorization of the linear system

    //! After Setup is called the BaseMatrix (expected to be either a CRS
    //! or SCRS-matrix) will be reordered using the Nested Dissection or
    //! the Minimum Degree Algorithm and then it will be LU-factorised.
    void Setup( BaseMatrix &sysmat);

    //! Direct solution of the linear system

    //! After Solve is called the matrix (which has already to be factorised
    //! by a call of Setup) is finally solved by backward-forward substitution.
    //! Note that the method will neglect the precond input parameter, since
    //! we perform a direct solution.
    void Solve( const BaseMatrix &sysmat, const BasePrecond &precond,
                const BaseVector &rhs, BaseVector &sol );

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return PARDISO
    SolverType GetSolverType() {
      return PARDISO;
    }

  private:

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    PardisoSolver();

    //! enum of PARIDO's error codes
    typedef enum {
      NO_ERROR = 0,
      INPUT_INCONSISTENT = -1,
      NOT_ENOUGH_MEMORY = -2,
      REORDERING_PROBLEM = -3,
      ZERO_PIVOT = -4,
      INTERNAL_ERROR = -5,
      PREORDERING_FAILED = -6,
      DIAGONAL_MATRIX = -7,
      INT_OVERFLOW = -8
    } PardisoError;

    //! Returns a string describing the given error code
    std::string GetErrorString(int err_code);

    //@{
    //! Vectors containing the information about the CRS- or SCRS-matrix

    //! This pointer is used to hold the address of a part of the internal
    //! (S)CRS matrix structures. The related memory segment must not
    //! be altered of deleted. Therefore the pointer is const!
    const Integer *rowPtr_;
    const Integer *colPtr_;
    const T *datPtr_;
    //@}

    //! Stored information about the storage type and entry type of the matrix
    MatrixStorageType stype;
    MatrixEntryType etype;

    //! A working array for the Pardiso-Routines

    //! This is the internal memory address array used by Pardiso. It must
    //! contain 64 entries, which are 4-byte integers on a 32-bit architecture
    //! and 8-byte integers on a 64-bit architecture.
    int pt_[64];

    //! Array for passing information to and from Pardiso

    //! This integer array is used to communicate parameters to and from
    //! Pardiso. In order to not confuse the indices, we use a static array
    //! going from iparm_[0] to iparm_[64], where iparm_[0] is unused and
    //! pass the address (iparm_ + 1) to Pardiso. Thus, we have the same
    //! indices for the parameters in PardisoSolver/C++ as in Pardiso/F77.
    int iparm_[65];

    //! The type of the matrix in a special encoding used by Pardiso
    int mType_;

    //! Maximal number of factors with identical nonzero-structure Pardiso
    //! should keep in memory at the same time.
    int maxfct;

    //! The number of the matrix (out of the maxfct ones) that should be used
    //! for the solution process
    int mnum;

    //! Dimension of the linear system
    int probDim_;

    //! The number of right hand sides Pardiso should solve the system for
    //! at one pass
    int nrhs;

    //! Specifies verbosity of Pardiso itself

    //! This attribute stores a number specifying the verbosity of the Pardiso
    //! solver itself. For a value of 0 no output is generated, will setting
    //! the message level to 1 results in statistical information being printed
    //! to the standard output. The attribute's value is set in the Setup()
    //! method corresponding to the PARDISO_stats parameter.
    int msgLvl_;

    //! Integer value of zero for passing to Pardiso
    int zeroINT_;

    //! Floating point value of zero for passing to Pardiso
    Double zeroDBL_;

    //! Array containing entries of problem matrix

    //! This array contains the entries of the problem matrix of the linear
    //! system that is to be solved. We imnplicitely assume that a Complex can
    //! directly be cast into an array of two Doubles.
    Double *theMatrix_;

    //! A flag specifying if Setup is being called for the first time
    bool firstCall_;

    //! Array with identity reordering

    //! This pointer is either NULL or points to a one-based array containing
    //! an identity re-ordering. The latter is used to force Pardiso to use
    //! the original matrix pattern, in the case that NOREORDERING is specified
    //! by the user via the PARDISO_ordering parameter.
    int *idPerm_;

    //! Size of identity permutation array
    int idPermSize_;

  };

}

#endif
