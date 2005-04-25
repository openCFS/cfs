#ifndef FILE_SOLVESTEPELEC
#define FILE_SOLVESTEPELEC

#include "stdSolveStep.hh"
#include "Utils/hysteresis.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Electrostatic-PDE

  class SolveStepElec : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepElec(StdPDE& apde);

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

    //! base method for solving one static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepStatic(const Integer kstep, const Double asteptime,
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

    //! solves for one nonlinear static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepStaticNonLinEpsDiff(const Integer kstep, const Double asteptime,
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

    //! compute polarization and add the term to RHS
    void AddPolarizationToRHS();

    //! update the hysteresis values
    void DoUpdateHyst();

    //! compute constant part of RHS for differential permittivity formulation 
    void ComputeConstPartRHS();

    //! computes differential permittivity
    void ComputeDiffEpsilon();

  private:

    Boolean doInit_;

    Hysteresis * hyst_;

    Double Esat_;
    Double Psat_;
    Double Ec_;
    Integer Dir_;
    Boolean isVirgin_;

    //for differential permittivity
    Vector<Double> Eprevious_;
    Vector<Double> Dprevious_;
    Matrix<Double> epsDiff_;
  };

} // end of namespace

#endif

