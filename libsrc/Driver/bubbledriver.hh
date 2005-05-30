#ifndef FILE_BUBBLEDRIVER_2004
#define FILE_BUBBLEDRIVER_2004

#include "singleDriver.hh"
#include "ODEDescr/BubbleODE.hh"

// IMPORTANT: here all bubble methods need to be included
#include "ODESolve/ODESolver_RKF45.hh"
#include "ODESolve/ODESolver_ExplEuler.hh"

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
                 UInt     stepOffset = 0,
                 Double      timeOffset = 0.0,
                 std::string driverTag = "anyTag",
                 Boolean     isPartOfSequence = FALSE);

    //! Destructor 
    ~BubbleDriver();

    //! This method constitutes the actual driving method which controls the
    //! solution process for the bubbledynamics.
    //! test phase: this is only the bubbledynamics
    //! Later it will be the coupled solution process of acoustical and 
    //! bubbeldynamical behaviour
    void SolveProblem();

  private:
    //! 
    UInt numstep_,isavebegin_,isaveincr_,isaveend_;

    //! 
    Double firstdt_;

    //! Pointer to class of bubbledynamical method
    BubbleODE *ptBubble_;

    //! Vector contains on input the initial data of radius und velocity
    //! for the bubbledynamic, on the output the solution 
    //! Needs to be changed if bubbledynamics is coupled to acoustics
    StdVector<Double> y_;

    //! Pointer to the solver class for ode's
    //! Later there could be a choice which solver should be used, then
    //! pointer needs to point to base class
    BaseODESolver *ptODESolver_;

    //! Pressure values needed for bubbledynamics
    //! Will be changed when coupled to acoustics
    Double pressure_;
    Double pressureAmpl_;
    //! Values of pressure derivative needed for bubbledynamics
    //! Will be changed when coupled to acoustics
    Double dpresdt_;

    Double frequency_;

    //!  Needed for print of bubblevalues, will bechanged later
    FILE *fp;

    //! Attribute describing model for bubble dynamics
    BubbleDynType bubbleDynType_;

    TimeFunc * ptTimeFunc_;   //!< pointer to time functions
    Char *fncFileName_; 
  };

}

#endif // FILE_BUBBLEDRIVER
