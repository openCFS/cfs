#ifndef FILE_SOLVESTEPACOUSTIC
#define FILE_SOLVESTEPACOUSTIC

#include "baseSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepAcoustic : public BaseSolveStep
  {

  public:

    //! Constructor
    SolveStepAcoustic(BasePDE& apde);

    //! Destructor
    virtual ~SolveStepAcoustic();


    //----------------------- TRANSIENT---------------------------------------
    //! solves for one nonlinear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */   
    void StepTransNonLin(const Integer kstep, const Double asteptime,
			 const Integer level, const Boolean reset);

    //! compute nonlinear part of RHS
    void AddNonlinearRHS();

  private:


  };

} // end of namespace

#endif

