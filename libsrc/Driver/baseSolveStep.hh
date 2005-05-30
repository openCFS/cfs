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
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic(const UInt kstep, const Double asteptime,
                               const Boolean reset) = 0;
    
    //! base method for solving one static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepStatic(const UInt kstep, const Double asteptime,
                                 const Boolean reset) = 0;

    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
    */  
    virtual void PostStepStatic(const UInt kstep, const Double asteptime)  = 0;



    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans(const UInt kstep, const Double asteptime,
                              const Boolean reset) = 0;

    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const UInt kstep, const Double asteptime,
                                const Boolean reset) = 0;

    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
    */  
    virtual void PostStepTrans(const UInt kstep, 
                               const Double asteptime) = 0;
    
    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param reset TRUE: perfrom new assembly, etc
    */   
    virtual void PreStepHarmonic(const UInt freqStep, const Double frequency, 
                                 const Boolean reset) = 0;

    //!  base method for solving one harmonic step 
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepHarmonic(const UInt freqStep, const Double frequency, 
                                   const Boolean reset) = 0;

    //!  routine for actions after the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void PostStepHarmonic(const UInt freqStep, const Double frequency, 
                                  const Boolean reset) = 0;

    //! Set the current time step value
    virtual void SetTimeStep( Double dt ) = 0;
    

  protected:

    //! Constructor
    BaseSolveStep();

  };


} // end of namespace

#endif

