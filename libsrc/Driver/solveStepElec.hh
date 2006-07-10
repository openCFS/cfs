#ifndef FILE_SOLVESTEPELEC
#define FILE_SOLVESTEPELEC

#include "stdSolveStep.hh"
#include "Utils/hysteresis.hh"

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
    void SolveStepStatic( );

    //! solves for one nonlinear static step 
    void StepStaticNonLinEpsDiff( );

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans( )
    {PreStepStatic();};

    //! base method for solving one transient step 
    void SolveStepTrans( )
    {SolveStepStatic();};

    //! solves for one linear transient step 
    void StepTransLin( )
    {StepStaticLin();};

    //! routine for actions after the SolveStep-method
     void PostStepTrans()
    {PostStepStatic();};

    //! compute polarization and add the term to RHS
    void AddPolarizationToRHS();

    //! update the hysteresis values
    void DoUpdateHyst();

    //! compute constant part of RHS for differential permittivity formulation 
    void ComputeConstPartRHS();

    //! computes differential permittivity
    void ComputeDiffEpsilon();

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

