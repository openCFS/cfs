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
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepStatic( const Boolean reset );

    //! base method for solving one static step 
    //!\param reset TRUE: perfrom new assembly, etc
    void SolveStepStatic( const Boolean reset );

    //! solves for one nonlinear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepStaticNonLin( const Boolean reset );

    //! solves for one nonlinear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepStaticNonLinEpsDiff( const Boolean reset );

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepTrans( const Boolean reset )
    {PreStepStatic(reset);};

    //! base method for solving one transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void SolveStepTrans( const Boolean reset )
    {SolveStepStatic(reset);};

    //! solves for one linear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransLin( const Boolean reset )
    {StepStaticLin(reset);};

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

    Boolean doInit_;

    Hysteresis * hyst_;

    Double Esat_;
    Double Psat_;
    Double Ec_;
    UInt Dir_;
    Boolean isVirgin_;

    //for differential permittivity
    Vector<Double> Eprevious_;
    Vector<Double> Dprevious_;
    Matrix<Double> epsDiff_;
  };

} // end of namespace

#endif

