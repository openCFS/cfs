#ifndef FILE_BASESOLVESTEP
#define FILE_BASESOLVESTEP

#include "Utils/StdVector.hh"
#include "General/environment.hh"

namespace CoupledField
{

  //! Base class for solution of a single step

  class BaseSolveStep
  {

  public:

    

    //! Destructor
    virtual ~BaseSolveStep();


     //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset) = 0;
 
    //! base method for solving one static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepStatic(const Integer kstep, const Double asteptime,
				 const Integer level, const Boolean reset) = 0;

    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepStatic(const Integer kstep, const Double asteptime,
				const Integer level)  = 0;



    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean reset) = 0;

    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
				const Integer level, const Boolean reset) = 0;

    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepTrans(const Integer kstep, const Double asteptime,
			       const Integer level) = 0;

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */   
    virtual void PreStepHarmonic(const Integer freqStep, const Double frequency, 
				 Integer level, const Boolean reset) = 0;

    //!  base method for solving one harmonic step 
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepHarmonic(const Integer freqStep, const Double frequency, 
				   Integer level, const Boolean reset) = 0;

    //!  routine for actions after the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void PostStepHarmonic(const Integer freqStep, const Double frequency, 
				  Integer level, const Boolean reset) = 0;


    virtual const StdVector<Double>& GetResultData(std::string resultType)
    {Error("BaseSolveStep: GetResultData not implemented");}

  protected:

    //! Constructor
    BaseSolveStep();

  };


} // end of namespace

#endif

