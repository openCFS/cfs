#ifndef FILE_SOLVESTEPELEC
#define FILE_SOLVESTEPELEC

#include "baseSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Electrostatic-PDE

  class SolveStepElec : public BaseSolveStep
  {

  public:

    //! Constructor
    SolveStepElec(BasePDE& apde);

    //! Destructor
    virtual ~SolveStepElec();


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

    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepStatic(const Integer kstep, const Double asteptime,
				const Integer level);

    //----------------------- TRANSIENT---------------------------------------
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

    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
				const Integer level, const Boolean reset)
    {SolveStepStatic(kstep,asteptime,level,reset);};

    //! solves for one linear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepTransLin(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean reset)
    {StepStaticLin(kstep,asteptime,level,reset);};

    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepTrans(const Integer kstep, const Double asteptime,
			       const Integer level)
    {PostStepStatic(kstep,asteptime,level);};



  private:


  };

} // end of namespace

#endif

