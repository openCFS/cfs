#ifndef FILE_SOLVESTEPELEC
#define FILE_SOLVESTEPELEC

#include "StdSolveStep.hh"
#include "Materials/Models/Hysteresis.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Electrostatic-PDE

  class SolveStepElec : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepElec(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepElec();


    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepStatic( );

    //! base method for solving one static step 
    void SolveStepStatic();

    //! solves for one nonlinear static step 
    void StepStaticNonLinEpsDiff();

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initializations before execution the SolveStep-method
    void PreStepTrans( );

    //! base method for solving one transient step 
    void SolveStepTrans()
    {SolveStepStatic();};

    //! solves for one linear transient step 
    void StepTransLin()
    {StepStaticLin();};

    //! routine for actions after the SolveStep-method
     void PostStepTrans()
    {PostStepStatic();};

    //! computes differential permittivity
    void SetPreviousVals4Hyst();

  private:

    bool doInit_;

#if 0 // Clang does not need this
    Hysteresis * hyst_;

    Double Esat_;
    Double Psat_;
    Double Ec_;
    UInt Dir_;
    bool isVirgin_;
#endif

    //for differential permittivity
    Vector<Double> Eprevious_;
    Vector<Double> Dprevious_;
    Matrix<Double> epsDiff_;
  };

} // end of namespace

#endif

