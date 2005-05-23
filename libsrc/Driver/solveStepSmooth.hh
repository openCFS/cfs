#ifndef FILE_SOLVESTEPSMOOTH
#define FILE_SOLVESTEPSMOOTH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step

  class SolveStepSmooth : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepSmooth(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepSmooth();

    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic(const Integer kstep, const Double asteptime,
                               const Boolean reset);

    //! solves for one nonlinear static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepStaticNonLin(const Integer kstep, const Double asteptime,
                                  const Boolean reset);

    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
    */  
    virtual void PostStepStatic(const Integer kstep, const Double asteptime);


    //----------------------- TRANSIENT---------------------------------------  
    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Integer kstep, const Double steptime, 
                                const Boolean updatesysmat)
    {SolveStepStatic(kstep,steptime,updatesysmat);};

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans(const Integer kstep, const Double asteptime,
                              const Boolean reset)
    {PreStepStatic(kstep,asteptime,reset);};
  
    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
    */  
    virtual void PostStepTrans(const Integer kstep, const Double asteptime)
    {PostStepStatic(kstep,asteptime);};
  
    //! solves for one nonlinear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepTransNonLin(const Integer kstep, const Double asteptime,
                                 const Boolean reset)
    {StepStaticNonLin(kstep,asteptime,reset);};


  private:


  };

} // end of namespace

#endif

