#ifndef FILE_SOLVESTEPACOUFLOWNOISE
#define FILE_SOLVESTEPACOUFLOWNOISE

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepAcouFlowNoise : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepAcouFlowNoise(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepAcouFlowNoise();


    //----------------------- TRANSIENT---------------------------------------
    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param steptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    void SolveStepTrans(const UInt kstep, const Double steptime, 
                        const Boolean reset);

  private:


  };

} // end of namespace

#endif

