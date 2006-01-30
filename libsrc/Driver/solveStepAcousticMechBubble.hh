#ifndef FILE_SOLVESTEPACOUSTICMECHBUBBLE
#define FILE_SOLVESTEPACOUSTICMECHBUBBLE

#include "stdSolveStep.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODEDescr/Gilmoredimlos.hh"
#include "ODESolve/ODESolver_RKF45.hh"
#include "ODESolve/ODESolver_Rosenbrock.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepAcousticMechBubble : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepAcousticMechBubble(DirectCoupledPDE& apde, BubbleDynType bubbleDynType);

    //! Destructor
    virtual ~SolveStepAcousticMechBubble();


    //----------------------- TRANSIENT---------------------------------------
    //! base method for solving one transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void SolveStepTrans( const Boolean reset ) ;


    //! solves for one nonlinear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransBubble( const Boolean reset );

    //!
    void ComputeBubbleRHS();

    //!  return pointer to vector with first derivative of solution
    const StdVector<Double>& GetResultData(std::string resultType);

  private:

    //! Pointer to acoustic PDE
    SinglePDE * acouPDE_;

    //! Pointer to mechanic PDE
    SinglePDE * mechPDE_;

    //-------------for bubble dynamics-------------------------------------
    //!method for ...
    void SetupBubbleDynamics();

    //! Pointer to class of bubbledynamical method
    StdVector<BubbleODE*> ptBubble_;

    //! Pointer to the solver class for ode's
    //! Later there could be a choice which solver should be used, then
    //! pointer needs to point to base class
    BaseODESolver *ptODESolver_;

    //! Vector contains on input the initial data of radius und velocity
    //! for the bubbledynamic, on the output the solution 
    StdVector<Double> bubbleValues_;

    //! Vector contains radius of bubbles for each element
    StdVector<Double> radiusWork_;
    StdVector<Double> radius_;

    //! Vector contains velocity of bubble walls for each element
    StdVector<Double> velocityWork_;
    StdVector<Double> velocity_;

    // In actual cavitation model: bubbledensity is constant in time, 
    // but may vary in space; in complexer models it can vary in time too.
    // Attention: for now assumption that density is constant in space. 
    //! Density of number of bubbles per unit volume
    Double bubbleDensity_;

    //! Attribute describing model for bubble dynamics
    BubbleDynType bubbleDynType_;


    //! Parameters for the bubble dynamics
    //! \param initRadius_ Initial radius of bubble
    //! \param initVel_ Initial velocity of bubble
    //! \param density_ Density of fluid
    //! \param sonicVel_ Sonic velocity of fluid
    //! \param pStatic_ Static pressure of fluid
    //! \param pVapour_ Vapour pressure of fluid
    //! \param surfaceTension_ Surface tension of fluid
    //! \param polytrop_ Polytropic index of fluid
    //! \param viscosity_ Viscosity of fluid
    Double initRadius_;
    Double initVel_;
    Double density_;
    Double sonicVel_;
    Double pStatic_;
    Double pVapour_;
    Double surfaceTension_;
    Double polytrop_;
    Double viscosity_;



    //! Parameter for the dimensionless case
    StdVector<Double> yNoDim_;
    Double tNoDim_;
    Double pressureNoDim;
    Double presDerivNoDim;  

    //! Suggested step size for ODESolver
    Double hTry_;
   


    DirectCoupledPDE* directPDE_;

  };


#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SolveStepAcousticBubble
  //! This class is used, if simulation acoustics
  //! computed in bubbly-liquid, it couples the 
  //! acoustic with the bubbledynamics  
  //! 
  //! \purpose
  //! 
  //! \collab
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! So far one has to comment in or out the way he
  //! wishes to compute, dimensionless or not, choice of
  //! solver, model for bubbledynamics in case of Gilmore
  //! dimensionless


#endif

} // end of namespace

#endif

