#ifndef FILE_ITERSOLVESTEP
#define FILE_ITERSOLVESTEP

#include "baseSolveStep.hh"
//#include "CoupledPDE/itercoupledpde.hh"
//#include "CoupledPDE/pdecoupling.hh"

namespace CoupledField
{
  // forward class declarations
  class IterCoupledPDE;
  class PDECoupling;

  //! Derived class for step-wise solving of iterative coupled StdPDEs
  class IterSolveStep : public BaseSolveStep
  {

  public:

    //! Constructor
    IterSolveStep(IterCoupledPDE& apde);

    //! Destructor
    virtual ~IterSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset)  {;};
 
    //! base method for solving one static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepStatic(const Integer kstep, const Double asteptime,
				 const Integer level, const Boolean reset);

    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepStatic(const Integer kstep, const Double asteptime,
				const Integer level) {;};



    //----------------------- TRANSIENT---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean reset) {;};


    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
				const Integer level, const Boolean updatesysmat);
    
    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepTrans(const Integer kstep, const Double asteptime,
			       const Integer level) {;};

    //----------------------- HARMONIC---------------------------------------
    
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */   
    virtual void PreStepHarmonic(const Integer freqStep, const Double frequency, 
				 Integer level, const Boolean reset) {;};


    //!  base method for solving one harmonic step 
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepHarmonic(const Integer freqStep, const Double frequency, 
				   Integer level, const Boolean reset);


    //!  routine for actions after the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void PostStepHarmonic(const Integer freqStep, const Double frequency, 
				  Integer level, const Boolean reset) {;};



  protected:

    //! calculates the norm of a vector
    Double CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval); 
    

    //! reference to PDE
    IterCoupledPDE &rPDE_;

    //! reference to coupling
    StdVector<PDECoupling*> & rCouplings_;

    //! reference to current level of solution
    Integer &actlevel_;
    
  };

} // end of namespace

#endif

