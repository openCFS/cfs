#ifndef FILE_BUBBLEDRIVER_2004
#define FILE_BUBBLEDRIVER_2004

#include "singleDriver.hh"
#include "ODEDescr/BubbleODE.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODEDescr/Gilmoredimlos.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODESolve/ODESolver_RKF45.hh"
#include "ODESolve/ODESolver_ExplEuler.hh"
#include "ODESolve/ODESolver_Rosenbrock.hh"

#include <cstdio>      // needed for test phase output

namespace CoupledField {
  
  // forward class declarations
  class TimeFunc;


  //! Driver for bubbledynamical problems.
  class BubbleDriver : public SingleDriver {

  public:
    //! Constructor
    //! \param adomain pointer to class Domain
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time
    //! \param driverTag tag for current driver section
    //! \param true, if driver is part of  multiSequence
    BubbleDriver(Domain      *adomain,
                 UInt        stepOffset = 0,
                 Double      timeOffset = 0.0,
                 std::string driverTag = "anyTag",
                 Boolean     isPartOfSequence = FALSE);

    //! Destructor 
    ~BubbleDriver();

    //! This method constitutes the actual driving method which controls the
    //! solution process for the stand-alone bubbledynamics.
    void SolveProblem();

  private:
    //! \param numstep_ Number of steps
    //! \param isavebegin_ First step to save
    //! \param isaveincr_  Increment for saving steps
    //! \param isaveend_   Last step to save
    UInt numstep_,isavebegin_,isaveincr_,isaveend_;

    //! \param firstdt_ Increment stepsize for calling ODE-Solver
    Double firstdt_;

    //! Pointer to class of bubbledynamical method
    BubbleODE *ptBubble_;

    //! Vector contains on input the initial data
    //! of radius und velocity for the bubbledynamic,
    //! on the output the solution 
    StdVector<Double> y_;

    //! Pointer to the solver class for ode's
    //! Later there could be a choice of solver
    BaseODESolver *ptODESolver_;

    //! Pressure values needed for bubbledynamics,
    //! if not computed in bubbledynamics or given 
    Double pressure_;

    //! Values of pressure derivative needed for bubbledynamics
    //! if not computed in bubbledynamics or given 
    Double dpresdt_;

    //! \param pressureAmpl_ Amplitude of pressure, needed only
    //! if pressure is directly computed
    //! \param frequency_ Frequency of pressure, needed only
    //! if pressure is directly computed
    Double pressureAmpl_;
    Double frequency_;

    //! Needed for printing of bubblevalue- results
    FILE *fp;

    //! Attribute describing model for bubble dynamics
    BubbleDynType bubbleDynType_;

    //! pointer to time functions
    TimeFunc * ptTimeFunc_;  
    Char *fncFileName_; 

    // For computation of dimensionless varibles, 
    // these must be available in solve
    //! \param initRadius_ Initial radius of bubble
    //! \param initVel_ Initial velocity of bubble
    //! \param density_ Density of fluid
    //! \param sonicVel_ Sonic velocity of fluid
    //! \param pStatic_ Static pressure of fluid
    //! \param pVapour_ Vapour pressure of fluid
    //! \param surfaceTension_ Surface tension of fluid
    //! \param polytrop_ Polytropic index of fluid
    //! \param viscosity_ Viscosity of fluid
    Double initRadius_, initVel_, density_, sonicVel_, 
      pStatic_, pVapour_, surfaceTension_, 
      polytrop_, viscosity_;
  };


#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class BubbleDriver
  //! 
  //! \purpose 
  //! This class serves as driver class for the stand-alone
  //! computation of the bubbledynamics.
  //! 
  //! \collab 
  //! The Bubbledriver calls
  //! either ODESolver_RKF45::Solve or 
  //! ODESolver_Rosenbrock::Solve;
  //! which then call the computation-method of one of 
  //! the bubbledynamical models, Gilmoredimlos::CompDeriv
  //! or Gilmore::CompDeriv or KellerMiksis::CompDeriv
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! So far one has to choose the options of the following points
  //! by directly commenting in and out: 
  //! a) the choice of the solver class,
  //! b) wether or not the dimensionless bubbledynamical equations 
  //!    has to be used, 
  //! c) the way the pressure is determined,
  //! d) wether or not the solver routine should be called at 
  //!    predetermined steps or rather called only once 
  //! At least the first two choices should be included in the XML-file

#endif


}

#endif // FILE_BUBBLEDRIVER
