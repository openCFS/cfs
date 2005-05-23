#ifndef FILE_SOLVESTEPACOUSTIC
#define FILE_SOLVESTEPACOUSTIC

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepAcoustic : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepAcoustic(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepAcoustic();


    //----------------------- TRANSIENT---------------------------------------
    //! solves for one nonlinear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */   
    void StepTransNonLin(const Integer kstep, const Double asteptime,
                         const Boolean reset);

    //! compute nonlinear part of RHS
    void AddNonLinRHS();

  private:


  };

} // end of namespace

#endif

