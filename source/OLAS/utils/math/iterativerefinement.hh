// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ITERATIVE_REFINEMENT_HH
#define ITERATIVE_REFINEMENT_HH

#include "General/defs.hh"

namespace CoupledField {

  class BaseDirectSolver;
class BaseMatrix;
  class BasePrecond;
class BaseVector;


  //! This processing class implements iterative refinement

  //! This class implements iterative refinement. The latter is a process
  //! to improve an approximate solution \f$x^\ast\f$ of a linear system
  //! \f$Ax=b\f$ typically computed by some kind of direct method
  //! \f$\mathcal{M}\f$. A single sweep of iterative refinement consists of
  //! the three steps:
  //!
  //! -# Compute the residual of the approximate solution \f$r=b-Ax^\ast\f$
  //! -# Solve the system \f$A\delta x=r\f$ for the update \f$\delta x\f$
  //! -# Update approximate solution \f$x^\ast \leftarrow x^\ast+\delta x\f$
  //!
  //! The solution in step 2 is typically computed using the same method
  //! \f$\mathcal{M}\f$, which in the case of a direct method is comparatively
  //! cheap, since the existing factorisation of \f$A\f$ can be used for this.
  //! The above steps can be repeated until the solution is good enough.
  //! \note
  //! - This class implements <em>fixed precision</em> iterative refinement,
  //!   i.e. the computation of the residual is performed using the same
  //!   floating point data type as that used for computation of \f$x^\ast\f$.
  //! - For a long time there was a general believe that iterative refinement
  //! would only be sensible, if the residual computation was performed in
  //! enhanced precision arithmetic, e.g. using double when \f$x^\ast\f$
  //! had been computed using floats. However, by means of rounding error
  //! analysis it can be shown that also <em>fixed precision</em> iterative
  //! refinement can be advantageous. For a detailed consideration of the
  //! numerical aspects of iterative refinement see e.g.
  //! <em>Accuracy and Stability of Numerical Algorithms</em>, Nicholas J.
  //! Higham, SIAM, 2nd ed.
  //! - The class is a pure processing class, i.e. it does not store any
  //!   data it operates on (the exception being two auxilliary vectors).
  class IterativeRefinement {

  public:

    // Default constructor
    IterativeRefinement();

    // Destructor
    ~IterativeRefinement();

    //! Perform iterative refinement

    //! \param mySolver solver used to compute solution \f$x^\ast\f$
    //! \param sysMat   matrix \f$A\f$ of the linear system
    //! \param sol      approximate solution \f$x^\ast\f$ of the linear system
    //! \param rhs      right hand side \f$b\f$ of the linear system
    //! \param numSteps number of refinement steps
    //! \param logLevel determines verbosity of the refinement process
    //! we currently distinguish three levels\n
    //! <center>
    //!   <table border="1" width="80%" cellpadding="10">
    //!    <tr>
    //!     <td align="center"><b>logLevel</b></td>
    //!     <td align="center"><b>meaning</b></td>
    //!    </tr>
    //!    <tr>
    //!     <td align="center">0</td>
    //!     <td>no reports to %OLAS log-file</td>
    //!    </tr>
    //!    <tr>
    //!     <td align="center">1</td>
    //!     <td>reports general startup information</td>
    //!    </tr>
    //!    <tr>
    //!     <td align="center">&gt; 1</td>
    //!     <td>prints residual norm for each step</td>
    //!    </tr>
    //!  </table>
    //! </center>
    void Refine( BaseDirectSolver &mySolver, BaseMatrix const &sysMat,
                 BaseVector &sol, BaseVector const &rhs, UInt& numSteps,
                 UInt logLevel );

  private:

    //! Generates auxilliary vectors
    void GenerateAuxilliaryVectors( const BaseMatrix &sysMat );

    //! Auxilliary vector for computing/storing the residual
    BaseVector *residual_;

    //! Auxilliary vector for computing/storing the update \f$\delta x\f$
    BaseVector *update_;

    //! Dummy preconditioner
    BasePrecond *dummyPrecond_;

  };

}

#endif
