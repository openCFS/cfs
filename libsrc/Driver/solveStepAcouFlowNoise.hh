#ifndef FILE_SOLVESTEPACOUFLOWNOISE
#define FILE_SOLVESTEPACOUFLOWNOISE

#include "baseSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepAcouFlowNoise : public BaseSolveStep
  {

  public:

    //! Constructor
    SolveStepAcouFlowNoise(BasePDE& apde);

    //! Destructor
    virtual ~SolveStepAcouFlowNoise();


    //----------------------- TRANSIENT---------------------------------------
    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param steptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    void SolveStepTrans(const Integer kstep, const Double steptime, const Integer level, 
			const Boolean reset);

  private:


  };

} // end of namespace

#endif

