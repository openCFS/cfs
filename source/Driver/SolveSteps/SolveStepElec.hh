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
    void SolveStepStatic(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);

    //! solves for one nonlinear static step 
    void StepStaticNonLinEpsDiff(PtrParamNode analysis_id);

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initializations before execution the SolveStep-method
    void PreStepTrans( );

    //! base method for solving one transient step 
    void SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL)
    {SolveStepStatic(analysis_id, adjointParams);};

    //! solves for one linear transient step 
    void StepTransLin(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL)
    {StepStaticLin(analysis_id, adjointParams);};

    //! routine for actions after the SolveStep-method
     void PostStepTrans()
    {PostStepStatic();};

    //! computes differential permittivity
    void SetPreviousVals4Hyst();

  private:

    bool doInit_;

    Hysteresis * hyst_;

    Double Esat_;
    Double Psat_;
    Double Ec_;
    UInt Dir_;
    bool isVirgin_;

    //for differential permittivity
    Vector<Double> Eprevious_;
    Vector<Double> Dprevious_;
    Matrix<Double> epsDiff_;
  };

} // end of namespace

#endif

