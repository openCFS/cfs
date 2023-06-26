#ifndef OLAS_CGSOLVER_HH
#define OLAS_CGSOLVER_HH

#include <iostream>
#include <fstream>



#include "BaseSolver.hh"

namespace CoupledField {

  class BasePrecond;

  //! Class for a Preconditioned Conjugate Gradient Solver

  //! This class represents an iterative solver of Conjugate Gradient type. It
  //! implements the classical Hennes-Stiefel CG approach for symmetric (resp.
  //! hermitian), positive definite matrices and in combination with a
  //! preconditioning of the linear system. With respect to preconditioning
  //! the class implements the untransformed variant of the CG method.
  //! The template parameter T relates
  //! to the scalar entry type of the system matrix in the linear system the
  //! solver is applied to.\n
  //! \n
  //! The behaviour of the CGSolver can be tuned via the following
  //! parameters which are read from the parameter object myParams_.
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for preconditioned CG solver</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>CG_maxIter</td>
  //!       <td align="center">&gt; 0</td>
  //!       <td align="center">50</td>
  //!       <td>Maximal number of iteration steps allowed for the CG method.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>CG_epsilon</td>
  //!       <td align="center">&gt; 0</td>
  //!       <td align="center">1e-6</td>
  //!       <td>This is the threshold for the stopping test. The CG iterations
  //!           terminated once the following holds
  //!           \f[
  //!           \frac{\|z^{(k)}\|_2}{\|z^{(0)}\|_2} < \varepsilon
  //!           \f]
  //!           where \f$z^{(k)}=M^{-1}(b-Ax^{(k)})\f$ is the preconditioned
  //!           residual of the k-th approximation.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>CG_resDirectly</td>
  //!       <td align="center">&gt; 0</td>
  //!       <td align="center">50</td>
  //!       <td>This is the interval for performing a full residual
  //!           computation. See also #resDirectly_ below.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>CG_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>This parameter is used to control the verbosity of the method.
  //!           If it is set to 'yes', then information on the CG convergence,
  //!           like e.g. the norm of the residual per iteration step will
  //!           be logged to the standard %OLAS report stream (*cla) ---removed logging--.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //! \n\n The above parameters are read from myParams_ every time a solve
  //! of the linear system is triggered. Once the system has been solved, the
  //! CGSolver writes three values into the myReport_ object:\n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="2" align="center">
  //!         <b>Values the preconditioned CG solver writes into myReport_</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Report</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>numIter</td>
  //!       <td>This is the number of iterations steps the CG method has
  //!         performed.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>finalNorm</td>
  //!       <td>This is final norm of the relative (preconditioned) residual:
  //!         \f$\displaystyle\frac{\|z^{(k)}\|_2}{\|z^{(0)}\|_2}\enspace.\f$
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>solutionIsOkay</td>
  //!       <td>This boolean is set to true, if the iteration was terminated
  //!           because the stopping criterion was fulfilled and to false
  //!           otherwise.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  template<typename T>
  class CGSolver : public BaseIterativeSolver {

  public:

    /** The CG constructor initialized the variables but does not
     * start any calculation or allocate huge memory.
     * The dafaults are handlended in the impementation itself.
     * @param xml_ the ParamNode which might contain the cg element - can be NULL
     * @param infoNode_ report object for storing general information on solution process */
    CGSolver(PtrParamNode xml, PtrParamNode olasInfo )
        : r_(NULL),
          s_(NULL),
          d_(NULL),
          q_(NULL) {
      xml_       = xml;
      infoNode_ = olasInfo->Get("cg");
      resDirectly_ = 0;
    };


    //! Default destructor

    //! The default destructor needs to be deep to free dynamically allocated
    //! memory.
    ~CGSolver();

    //! Solve the linear system sysmat*sol=rhs for sol

    //! \param sysmat System matrix of the linear system to be solved
    //! \param rhs Right-hand side vector of the linear system
    //! \param sol Solution vector of linear system
    void Solve(const BaseMatrix& sysmat,  
	       const BaseVector& rhs, BaseVector& sol );

    //! Dummy setup method

    //! This method implements the pure virtual setup function defined in the
    //! BaseSolver class. In the case of the CG solver, there is nothing to
    //! do actually.
    //! \note Depending on how the BaseSolver interface develops this method
    //! might be removed again.
    void Setup( BaseMatrix &sysmat );

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return CG
    SolverType GetSolverType() {
      return CG;
    }

  private:

    //@{
    //! Auxilliary vector

    //! This is an auxilliary vector. It is dynamically allocated by the
    //! class and destroyed by the destructor.
    BaseVector *r_, *s_, *d_, *q_;
    //@}

    //! Constant for determining when to compute residual directly

    //! In the CG method the residual vector is not computed directly, but
    //! via a recurrence relation. This is done for the sake of efficiency.
    //! However, due to rounding errors, the residual vector computed via the
    //! recurrence relation can start to deviate from the real residual vector.
    //! Thus, the residual is computed directly via the formula \f$r = b-Ax \f$
    //! every resDirectly_ steps.
    UInt resDirectly_;

    //! Default Constructor

    //! The default constructor is private in order to disallow its use. The
    //! reason is that we need pointers to the parameter and report objects to
    //! perform correct setup of an object of this class.
    CGSolver(){};

  };

} //namespace


#endif // OLAS_CGSOLVER_HH
