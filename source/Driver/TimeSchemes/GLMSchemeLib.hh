// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GLMSchemeLib.hh
 *       \brief    Class declaration for all GLM schemes available
 *
 *       \date     03/01/2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef GLMSCHEMELIB_HH_
#define GLMSCHEMELIB_HH_

#include "MatVec/Matrix.hh"
#include "Domain/Domain.hh"
namespace CoupledField {

/*! \class GLMScheme
 *    \brief The base class defining all variable necessary to define a GLM
 *    @author A.Hueppe
 *    @date 01/2012
 *
 *  Every attribute is public and can be directly altered. For adding a scheme one has to
 *  fill the parameters according to the schemes definition and has to implement the
 *  ComputeCoefficients method in which the tableau for the scheme is given for each solution
 *  derivative order. If some orders cannot be given e.g. solderivorder 0 with an explicit RK4 it is the
 *  programmers duty to give an error.
 *  A detailed information about the ideas behind a GLM can be found in the developer manual.
 */
class GLMScheme{
  public:

  Double dtCurrent_      = -1.0;
  Double dtPrev1_        = -1.0;
  Double dtPrev2_        = -1.0;
  Double prev_dtCurrent_ = -1.0;
  Double prev_dtPrev1_   = -1.0;
  Double prev_dtPrev2_   = -1.0;

  bool   adaptiveEnabled_ = false;  // adaptive timestepping active (BDF2 or Newmark)
  Double local_error_  = 0.0;
  bool   initialized_  = false;
  bool   coefChanged_  = false;  // true when a0/dt changed; triggers system matrix rebuild

  /// Enumeration for each GLM scheme available
  typedef enum{
    TRAPEZOIDAL = 1,
    NEWMARK = 2,
     BDF2 = 3,
     RK4 = 4
  } SchemeType;

  
    GLMScheme();

    virtual ~GLMScheme();

    /*!
     * This method creates the tableau for the current timescheme based on the
     * requested solution order, i.e. formerly known as effective mass or stiffness and
     * a given timestep
     * @param solDerivOrder Which order of time derivative has the solution
     * @param deltaT Current timestepsize
     */
    virtual void ComputeCoefficients(UInt solDerivOrder,Double deltaT)=0;

    //! Adaptive: advance dt history one step (called once per attempt from BeginStep).
    //! Default no-op; BDF2 overrides it (its tableau depends on the step ratio dtCurrent_/dtPrev1_).
    virtual void AdvanceAdaptiveStep(Double /*newDt*/) {}

    /*!
     * Transforms a given BC value according to the current scheme formulation.
     * E.g. if the user specifies a dirichlet condition on the unknown but we solve for
     * the first time derivative we need to trans form the BC value. For Trapezoidal,
     * Newmark and HHT method we give here a general algorithm which goes like
     *
     * \ol Transform the coefficient matrix to valDerivOrder
     * \ol Multiply GLM Vector with the transformed tableau coefficients on the correct row
     * \ol update the value
     *
     * For other schemes we need to overwrite this method
     *
     * @param glm Reference to the current glm vector
     * @param value The initial BC value
     * @param valDerivOrder The requested time derivative of the BC value
     * @param eqnNumber The equation number of the BC value for accessing the GLM
     * @return The transformed BC value
     */
    virtual Double TransformBC(const StdVector< SingleVector* > & glm, Double value,
                                  UInt valDerivOrder, Integer eqnNumber);

    ///Sets some additional things depending on the scheme
    virtual void PrepareStage(UInt i,Double aTime, Domain* domain){

    };

    //! Get type of scheme
    virtual SchemeType GetType() const = 0;
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++
    //Define Scheme Formulation
    //++++++++++++++++++++++++++++++++++++++++++++++++++
    ///Define maximum time derivative order for the scheme
    UInt maxDerivOrder_;

    ///Define the order of derivative of the solution of the equation system
    ///e.g. for Tapezoidal: EffMass->solDerivOrder_=1 or EffStiff->solderivOrder_=0
    UInt solDerivOrder_;

    ///Define number of stages
    UInt numStages_;
    

    /*!
     Store the matrix of coefficients defining the scheme for a given solution order
     The number of rows is (maxDerivOrder_+1)*(solDerivOrder_) + sizeGLMVec_
     The number of cols is numStages + sizeGLMVec_
     */
    Matrix<Double> schemeCoefs_;

    //++++++++++++++++++++++++++++++++++++++++++++++++++
    // DEFINE DIMENSIONS OF GLM VECTOR
    //++++++++++++++++++++++++++++++++++++++++++++++++++
    ///Define number of old solutions to be stored
    UInt numOldSols_;

    ///Define number of first order time derivatives of solutions to be stored in GLM vector
    UInt numSol1stDerivs_;

    ///Define number of second order time derivatives of solution to be stored in GLM vector
    UInt numSol2ndDerivs_;

    ///Just for convenience we store the total size of the GLM vector
    UInt sizeGLMVec_;

    ///Store the current time step size
    Double curTStepSize_;

    //++++++++++++++++++++++++++++++++++++++++++++++++++
    // Optimization flags
    //++++++++++++++++++++++++++++++++++++++++++++++++++
    ///If this is true we do not need an update step for the solution
    ///e.g. trapezoidal eff. stiff.
    bool lastStageIsSolution_;

    ///For some schemes we can use the computed right hand side values for the update step
    bool usePredictors_;

  private:


};


/*! \class Trapezoidal
 *    \brief Defines coefficients for Trapezoidal timestepping
 *    @author A.Hueppe
 *    @date 01/2012
 *
 *  The tableau for trapezoidal time stepping in effective mass formulation reads as
 *  <table border=0  cellpadding="5"  cellspacing="0" >
 *  <tr>
 *     <td  style="border-right: 1px solid black; padding: 5px;">\f$ \gamma_\mathrm{P} \Delta t\f$ </td>
 *     <td> -1 </td>
 *     <td> \f$ (\gamma_\mathrm{P}-1) \Delta t\f$</td>
 *  </tr><tr>
 *     <td style="border-right: 1px solid black; padding: 5px;">1</td>
 *     <td> 0 </td>
 *     <td> 0</td>
 *  </tr>
 *  <tr><td colspan="3"><hr></td></tr>
 *  <tr>
 *     <td style="border-right: 1px solid black; padding: 5px;">\f$ \gamma_\mathrm{P} \Delta t\f$ </td>
 *     <td> 1 </td>
 *     <td> \f$ (1 - \gamma_\mathrm{P}) \Delta t\f$</td>
 *  </tr><tr>
 *     <td style="border-right: 1px solid black; padding: 5px;"> 1 </td>
 *     <td> 0 </td>
 *     <td> 0 </td>
 *  </tr>
 *  </table>
 *
 *  for gamma_ = 0.5 the scheme is second order accurate
 *
 */

class Trapezoidal : public GLMScheme{
  public:

    Trapezoidal(Double gamma);

    //! \copydoc GLMScheme::ComputeCoefficients(UInt,Double)
    virtual void ComputeCoefficients(UInt solDerivOrder,Double deltaT);

    //! \copydoc GLMScheme::GetType
    virtual SchemeType GetType() const {
      return TRAPEZOIDAL;
    }
    
  private:
    /*!
     * parameter for switching between implicit and explicit scheme
     *  for gamma = 0.5 the scheme is second order accurate
     */
    Double gamma_;

};

/*! \class Newmark
 *    \brief Defines coefficients for Newmark timestepping
 *    @author A.Hueppe
 *    @date 01/2012
 *
 *  The tableau for Newmark time stepping in effective mass form reads as
 *  <table border=0  cellpadding="5"  cellspacing="0" >
 *  <tr>
 *     <td  style="border-right: 1px solid black; padding: 5px;">\f$ \beta \Delta t^2\f$ </td>
 *     <td> -1 </td>
 *     <td> \f$ - \Delta t\f$</td>
 *     <td> \f$ (\beta - 0.5) \Delta t^2\f$</td>
 *  </tr><tr>
 *     <td  style="border-right: 1px solid black; padding: 5px;">\f$ \gamma \Delta t\f$ </td>
 *     <td> 0.0 </td>
 *     <td> \f$ -1.0\f$</td>
 *     <td> \f$ (\gamma - 1.0) \Delta t\f$</td>
 *  </tr>
 *  <tr>
 *     <td  style="border-right: 1px solid black; padding: 5px;">1 </td>
 *     <td> 0.0 </td>
 *     <td> 0.0 </td>
 *     <td> 0.0 </td>
 *  </tr>
 *  <tr><td colspan="4"><hr></td></tr>
 *  <tr>
 *     <td  style="border-right: 1px solid black; padding: 5px;">\f$ \beta \Delta t^2\f$ </td>
 *     <td> 1 </td>
 *     <td> \f$ \Delta t\f$</td>
 *     <td> \f$ (0.5 - \beta) \Delta t^2\f$</td>
 *  </tr><tr>
 *     <td  style="border-right: 1px solid black; padding: 5px;">\f$ \gamma \Delta t\f$ </td>
 *     <td> 0.0 </td>
 *     <td> \f$ 1.0\f$</td>
 *     <td> \f$ (1.0 - \gamma) \Delta t\f$</td>
 *  </tr>
 *  <tr>
 *     <td  style="border-right: 1px solid black; padding: 5px;">1 </td>
 *     <td> 0.0 </td>
 *     <td> 0.0 </td>
 *     <td> 0.0 </td>
 *  </tr>
 *  </table>
 *
 *  For gamma = 0.5 and beta = 0.25, the scheme is second order accurate if we have no damping i.e. no dependence on first order
 *  time derivative.
 */
class Newmark : public GLMScheme{
  public:

    Newmark(Double gamma,Double beta,Double alpha=0.0);

    //! \copydoc GLMScheme::GetType
    virtual SchemeType GetType() const {
      return NEWMARK;
    }
    
    //! \copydoc GLMScheme::ComputeCoefficients(UInt,Double)
    virtual void ComputeCoefficients(UInt solDerivOrder,Double deltaT);


    virtual void PrepareStage(UInt i,Double aTime, Domain* domain);

    Double GetGamma() const { return gamma_; }
    Double GetBeta()  const { return beta_; }
    Double GetAlpha() const { return alpha_; }

  private:
    /*!parameter for switching between implicit and explicit scheme
     * for gamma = 0.5 the scheme is second order accurate */
    Double gamma_;

    /*!parameter for switching between implicit and explicit scheme
     * for beta_ = 0.25 the scheme is second order accurate in the
     * absence of damping matrix*/
    Double beta_;

    /*!parameter for alpha method retains second order accuracy
     * in presence of a damping matrix
     * alpha_ = 0 corresponds to the standard newmark scheme*/
    Double alpha_;

};

/*! \class BDF2
 *    \brief Defines coefficients for BDF2 timestepping
 *    @author A.Hueppe
 *    @date 02/2014
 */
class Bdf2 : public GLMScheme{
  public:

  Bdf2();

    //! \copydoc GLMScheme::GetType
    virtual SchemeType GetType() const {
      return BDF2;
    }

    //! \copydoc GLMScheme::ComputeCoefficients(UInt,Double)
    virtual void ComputeCoefficients(UInt solDerivOrder,Double deltaT);

    //! \copydoc GLMScheme::AdvanceAdaptiveStep(Double)
    virtual void AdvanceAdaptiveStep(Double newDt) override;

    virtual void PrepareStage(UInt i,Double aTime, Domain* domain){
     /// domain->GetMathParser()->SetValue( MathParser_GLOB_HANDLER,
     ///                                    "t", aTime+(alpha_*curTStepSize_) );
    }
  private:

  Double w_;
};


class RungeKutta4 : public GLMScheme{
  public:

    RungeKutta4();

    //! \copydoc GLMScheme::GetType
    virtual SchemeType GetType() const {
      return RK4;
    }

    //! \copydoc GLMScheme::ComputeCoefficients(UInt,Double)
    virtual void ComputeCoefficients(UInt solDerivOrder,Double deltaT);

    Double TransformBC(const StdVector< SingleVector* > & glm,
                       Double value,
                       UInt valDerivOrder,
                       Integer eqnNumber);

    virtual void PrepareStage(UInt i,Double aTime);

};



}

#endif /* GLMSCHEMELIB_HH_ */
