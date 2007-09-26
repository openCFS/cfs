// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ILDL_PRECOND_HH
#define ILDL_PRECOND_HH

#include "matvec/matvec.hh"
#include "solver/solver.hh"
#include "utils/math/ldlsystemsolve.hh"
#include "precond/ILDLPrecond/baseildlfactoriser.hh"
#include "precond/ILDLPrecond/ildlkfactoriser.hh"


namespace OLAS {


  //! This is a base class for all incomplete LDL preconditioners

  //! This is a base class for all incomplete LDL preconditioners
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
  //!         <b>Parameters for ILDLPrecond object</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPrecond_subType</td>
  //!       <td align="center">ILDLK, ILDLT, ILDLG </td>
  //!       <td align="center">-</td>
  //!       <td>This parameter describes the variant of ILDL preconditioning
  //!           that is to be used in the factorisation. A description of the
  //!           parameters for these subType can be found in the documentation
  //!           for the corresponding factoriser classes: ILDLKFactoriser,
  //!           ILDLTPFactoriser and ILDLCNFactoriser</td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPRECOND_saveFacToFile</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If the value of this boolean parameter is true, then the
  //!           factorisation matrix \f$F=D+L^T\f$ computed in the
  //!           factorisation phase of the preconditioner, will be exported to
  //!           a file in MatrixMarket format. This IO operation will take
  //!           place immediately after the end of the factorisation phase.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPRECOND_facFileName</td>
  //!       <td align="center">string</td>
  //!       <td align="center">-</td>
  //!       <td>File name for exporting the factorisation matrix. This value
  //!           is only required, when ILDLPRECOND_saveFacToFile is 'true'.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPRECOND_facPatternOnly</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If the value of this boolean parameter is true, then only the
  //!           sparsity pattern of the  factorisation matrix \f$F=D+L^T\f$
  //!           will be exported at the end of the analyse phase. This value
  //!           is only required, when ILDLPRECOND_saveFacToFile is 'true'.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPRECOND_logging</td>
  //!       <td align="center">0, 1, 2</td>
  //!       <td align="center">1</td>
  //!       <td>This parameter is used to control the verbosity of the solver.
  //!           For 0 no logging is done (besides the progress report of the
  //!           factorisation). For 1 information on the factorisation process
  //!           is logged and for 2 also some information on the application of
  //!           the preconditioner is written to the standard %OLAS
  //!           report stream (*cla).
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //!
  //! \note
  //! - The class stores its own copy of the factorised matrix. Thus, it is
  //!   safe to alter or destroy the system matrix after Setup has been called.
  //!   However, this also means that subsequent changes to the system matrix
  //!   will not affect the Apply method unless Setup is called again.
  //! - Pivoting is currently not supported during the factorisation.
  //! - The class only works together with SCRS_Matrices with scalar entries.
  template<class T>
  class ILDLPrecond : public BNPrecond<ILDLPrecond<T>, SCRS_Matrix<T>, T >,
                      public LDLSystemSolve<T> {

  public:

    // ILDLKFactoriser is a close friend of ours
    // template <class A> 
    // friend class ILDLKFactoriser;

    //! Constructor
    ILDLPrecond( const StdMatrix &stdMat, OLAS_Params *myParams,
                 OLAS_Report *myReport = NULL );

    //! Default Destructor

    //! The default destructor needs to be deep in order to free the memory
    //! dynamically allocated for the vector of pivot indices and the copy of
    //! the system matrix.
    ~ILDLPrecond();

    //! Computation of the incomplete LDL factorisation

    //! The setup method takes care of the incomplete LDL factorisation of the
    //! problem matrix. The actual factorisation is delegated to the internal
    //! factoriser_ object.
    void Setup( StdMatrix &sysMat );

    //! Direct solution of the linear system

    //! The solve method takes care of the solution of the linear system.
    //! The solution is computed via a backward/forward substitution pair
    //! using the LDL factorisation of the problem matrix computed in the
    //! setup phase.
    //! Note that the method will neglect the precond input parameter, since
    //! we perform a direct solution. Note also, that the sysmat input
    //! parameter will only be used, when an iterative refinement is performed.
    void Apply( const StdMatrix &stdMat, const SparseVector &,
                SparseVector &z ) const;

    //! Query type of this preconditioner.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type PrecondType.
    //! \return The return value depends on the speficic ILDL variant this
    //! object uses. The latter is stored in the myVariant_ attribute.
    PrecondType GetPrecondType() const {
      return myVariant_;
    }

  private:

    //! Object implementing the specific variant of ILDL factorisation

    //! This attribute stores a pointer to a BaseILDLFactoriser object. The
    //! object itself is generated by the GenerateFactoriser() method depending
    //! on the value of the ILDLPRECOND_subType parameter.
    BaseILDLFactoriser<T> *factoriser_;

    //! Variant of ILDL factorisation
    PrecondType myVariant_;

    //! Keep track on status of factorisation

    //! The Setup routine will set this boolean to true, once it has computed
    //! a factorisation. This gives us a confortable way of keeping track of
    //! the status of the preconditioner and the internal memory allocation.
    bool amFactorised_;

    //! Dimension of problem

    //! This attribute stores the dimension of the problem, i.e. the order
    //! of the problem matrix
    UInt sysMatDim_;

    //@{ \name Administration of matrix factors
    //! The matrix factor \f$U=L^T\f$ is stored in CRS storage format without
    //! the diagonal entries, which we know are all one. Within each row, the
    //! non-zero entries are stored in lexicographic ordering with respect to
    //! their column indices starting from the diagonal.

    //! Values of non-zero entries of factor \f$U=L^T\f$
    std::vector<T> dataU_;

    //! Column indices of non-zero entries factor \f$U=L^T\f$
    std::vector<UInt> cidxU_;

    //! Pointers to row starts in contiguous storage of factor \f$U=L^T\f$
    std::vector<UInt> rptrU_;

    //! Inverses of values of factor \f$D\f$

    //! In this vector the entries of the diagonal factor \f$D\f$ are stored.
    //! To be more precise, it is their inverses, since for performance reasons
    //! we already compute the inverses in the factorisation step. Divisions
    //! can be much more costly than mutliplications on some architectures.
    std::vector<T> dataD_;

    //@}

    //! Default constructor

    //! The default constructor is private, since we want the solver object to
    //! be initialised with a pointer to a parameter object right at
    //! instantiation.
    ILDLPrecond();

    //! Generate factoriser (obviously :(
    void GenerateFactoriser();

    //! Export ILDL factorisation matrix to a file in MatrixMarket format

    //! The method will export the factorisation matrix to an ascii file
    //! according to the MatrixMarket specifications. By factorisation matrix
    //! we understand the matrix \f$F=D+L^T\f$.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname       name of output file
    //! \param patternOnly if true, only the sparsity pattern is exported
    void ExportFactorisation( const char *fname, bool patternOnly = false );

  };

}

#endif
