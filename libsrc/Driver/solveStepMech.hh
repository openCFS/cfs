#ifndef FILE_SOLVESTEPMECH
#define FILE_SOLVESTEPMECH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step

  class SolveStepMech : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMech(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMech();

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
//     virtual void StepStaticLin(const Integer kstep, const Double asteptime,
// 			       const Integer level, const Boolean reset);

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
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
//     virtual void PreStepTrans(const Integer kstep, const Double asteptime,
// 			      const Integer level, const Boolean reset);

    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
//     virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
// 				const Integer level, const Boolean updatesysmat);

    //! solves for one linear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
//     virtual void StepTransLin(const Integer kstep, const Double asteptime,
// 			      const Integer level, const Boolean updatesysmat);


    //! solves for one linear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    void StepTransNonLin(const Integer kstep, const Double asteptime,
			 const Integer level, const Boolean reset);

    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
//     virtual void PostStepTrans(const Integer kstep, const Double asteptime,
// 			       const Integer level);


  private:


  };

} // end of namespace

#endif

