#ifndef FILE_SOLVESTEPACOUSTICBUBBLE
#define FILE_SOLVESTEPACOUSTICBUBBLE

#include "stdSolveStep.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODESolve/ODESolver_RKF45.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepAcousticBubble : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepAcousticBubble(StdPDE& apde, BubbleDynType bubbleDynType);

    //! Destructor
    virtual ~SolveStepAcousticBubble();


    //----------------------- TRANSIENT---------------------------------------
    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    void SolveStepTrans(const UInt kstep, const Double asteptime,
                        const Boolean reset) ;


    //! solves for one nonlinear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */   
    void StepTransBubble(const UInt kstep, const Double asteptime,
                         const Boolean reset);

    //!
    void ComputeBubbleRHS();

    //!  return pointer to vector with first derivative of solution
    const StdVector<Double>& GetResultData(std::string resultType);

  private:

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


  };

} // end of namespace

#endif

