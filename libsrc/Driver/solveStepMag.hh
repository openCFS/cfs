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
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic(const UInt kstep, const Double asteptime,
                               const Boolean reset);

    //! solves for one linear static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepStaticNonLin(const UInt kstep, const Double aTime,
                                  const Boolean reset);
    
    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
    */  
    virtual void PostStepStatic();


    //----------------------- TRANSIENT---------------------------------------
    //! solves for one nonlinear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    void StepTransNonLin(const UInt kstep, const Double asteptime,
                         const Boolean reset);


  private:


  };

} // end of namespace

#endif

