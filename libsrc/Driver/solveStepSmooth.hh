#ifndef FILE_SOLVESTEPSMOOTH
#define FILE_SOLVESTEPSMOOTH

#include "baseSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step

  class SolveStepSmooth : public BaseSolveStep
  {

  public:

    //! Constructor
    SolveStepSmooth(BasePDE& apde);

    //! Destructor
    virtual ~SolveStepSmooth();

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

    //! solves for one nonlinear static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepStaticNonLin(const Integer kstep, const Double asteptime,
				  const Integer level, const Boolean reset);

     //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepStatic(const Integer kstep, const Double asteptime,
				const Integer level);


   //----------------------- TRANSIENT---------------------------------------  
    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Integer kstep, const Double steptime, 
				const Integer level, 
				const Boolean updatesysmat)
    {SolveStepStatic(kstep,steptime,level,updatesysmat);};

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean reset)
    {PreStepStatic(kstep,asteptime,level,reset);};
  
    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepTrans(const Integer kstep, const Double asteptime,
			       const Integer level)
    {PostStepStatic(kstep,asteptime,level);};
  
    //! solves for one nonlinear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepTransNonLin(const Integer kstep, const Double asteptime,
				 const Integer level, const Boolean reset)
    {StepStaticNonLin(kstep,asteptime,level,reset);};


  private:


  };

} // end of namespace

#endif

