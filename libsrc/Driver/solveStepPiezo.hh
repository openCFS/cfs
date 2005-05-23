#ifndef FILE_SOLVESTEPPIEZO
#define FILE_SOLVESTEPPIEZO

#include "stdSolveStep.hh"


namespace CoupledField
{

  //! Base class for solution of a single step: Piezotrostatic-PDE

  class SolveStepPiezo : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepPiezo(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepPiezo();


    //----------------------- TRANSIENT---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans(const Integer kstep, const Double asteptime,
                              const Boolean reset);

    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
                                const Boolean reset);

    //! solves for one nonlinear transient step (with hysteresis) 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepTransNonLinEpsDiff(const Integer kstep, const Double asteptime,
                                        const Boolean reset);

    //! update the hysteresis values
    void DoUpdateHyst();

    //! computes differential permittivity
    void ComputeDiffEpsilon();

  private:

    Boolean doInit_;

    //for differential permittivity
    Vector<Double> Eprevious_;
    Vector<Double> Dprevious_;
    Matrix<Double> epsDiff_;
  };

} // end of namespace

#endif

