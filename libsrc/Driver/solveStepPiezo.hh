#ifndef FILE_SOLVESTEPPIEZO
#define FILE_SOLVESTEPPIEZO

#include "baseSolveStep.hh"
#include "Utils/hysteresis.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Piezoelectrics

  class SolveStepPiezo : public BaseSolveStep
  {

  public:

    //! Constructor
    SolveStepPiezo(BasePDE& apde);

    //! Destructor
    virtual ~SolveStepPiezo();


    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean reset);

    //! solves for one nonlinear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */   
    void StepTransNonLin(const Integer kstep, const Double asteptime,
			 const Integer level, const Boolean reset);

    //! compute polarization and add the term to RHS
    void AddPolarizationToRHS();

    //! update the hysteresis values
    void DoUpdateHyst();

  private:

    Boolean doInit_;

    Hysteresis * hyst_;

    Double Esat_;
    Double Psat_;
    Double Ec_;
    Integer Dir_;
    Boolean isVirgin_;

  };

} // end of namespace

#endif

