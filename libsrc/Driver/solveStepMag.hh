#ifndef FILE_SOLVESTEPMAG
#define FILE_SOLVESTEPMAG

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepMag :  public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMag(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMag();

    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset);

    //! solves for one linear static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepStaticNonLin(const Integer kstep, const Double aTime,
				  const Integer level, const Boolean reset);
    
    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepStatic(const Integer level);


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


  private:


  };

} // end of namespace

#endif

