// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_TIMESTEPPING_2001
#define FILE_TIMESTEPPING_2001

#include <map>
#include <ostream>
#include <utility>

#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/vector.hh"

namespace CoupledField {


  // forward class declaration
  class BaseSystem;

  // typedef of timesteps: the higher the number, the older the time step
  // these are also used for the coefficient vectors. 
  // @warning: sol_timeStepVec_ may have different number of entries as
  // coefficients
  typedef enum {NO_TIMESTEPTYPE = -1, TIMESTEP_0 = 0, TIMESTEP_1 = 1, \
    TIMESTEP_2 = 2, \
    PREDICTOR_1, PREDICTOR_2, \
    CORRECTOR_1, CORRECTOR_2, \
    COEFFRHS}
  TIMEStepType;
  typedef enum {NO_DERIVTYPE = -1, FIRST_DERIV = 1, SECOND_DERIV = 2}
  DERIVType;


  //! Base class for time stepping algorithms

  class TimeStepping
  {
  public:
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 

    TimeStepping( BaseSystem * algebraicsystem );

    //! deconstructor
    virtual ~TimeStepping();
  
    //! initialization
    //! \param rhsSIze total number of entries in the rhs vector
    virtual void Init( Double dt, UInt rhsSize ) = 0;
    
    //! Reinitialization (setting all remembered vectors to zero)    
    virtual void ReInit(){
      EXCEPTION("Not implemented here!");
    }

    virtual void setSubSteps( UInt subSteps )
    {EXCEPTION("Error not implemented!");};

    virtual void resetDeltaT( )
    {EXCEPTION("Error not implemented!");};

    //! get matrix factors for effective systems matrix
    virtual const std::map<FEMatrixType,Double>  &
    GetEffSysMatFactors( ) const;

    //! perform predictor step
    virtual void Predictor(Vector<Double>& solold)=0;
  
    //! perform corrector step
    virtual void Corrector(Vector<Double>& solnew)=0;
  
    //! perform an update to RHS
    virtual void UpdateRHS()=0;

    //! Substract Stiffness from RHS
    virtual void SubstractStiffnessFromRHS(Vector<Double>& actSol)
    {EXCEPTION("Error not implemented!");};

    //! perform an update to RHS with actual solution (for nonlin calculation)
    virtual void UpdateRHS(Vector<Double>& actSol)
    {EXCEPTION("Error not implemented!");};

    //! perform calculations at end of timestep different from above methods 
    //! e.g. when fractional damping is used, store current in memory
    virtual void AdvanceTimestep(Vector<Double>& solnew){;};

    /**
     * Sets the time step timeStepType to the value timeStep_vec
     * @param timeStep_vec The time step which should be set
     * @param timeStepType Which time step is to be set
     */
    virtual void SetOld(const Vector<Double> & timeStep_vec, TIMEStepType timeStepType)
    {
      sol_timeStepVec_[timeStepType] = timeStep_vec;
    }

    /**
     * Returns the time step timeStepType as a vector
     * @param timeStepType The time step which is demanded
     * @retrun The vector of a particular time step
     * @throws exception if time step has not been set
     */
    virtual const Vector<Double>& GetOld(TIMEStepType timeStepType)
    {
      std::map<TIMEStepType, Vector<Double> >::iterator iter = \
        sol_timeStepVec_.find(timeStepType);
      if (iter == sol_timeStepVec_.end())
      {
        EXCEPTION("Tried to access timeStepType: " << timeStepType << " which is not set!");
      }
      return iter->second;
    }

    /**
     * Sets the time step timeStepType to the value timeStep_vec
     * @param timeStep_vec The time step which should be set
     * @param timeStepType Which time step is to be set
     */
    virtual void SetDeriv(const Vector<Double> & deriv_vec, DERIVType derivType)
    {
      solDeriv_vec_[derivType] = deriv_vec;
    }

    /**
     * Returns the dervative as a vector
     * @param derivType The time step which is demanded
     * @retrun The vector of a particular derivative
     * @throws exception if derivative has not been set
     */
    virtual const Vector<Double>& GetDeriv(DERIVType derivType)
    {
      std::map<DERIVType, Vector<Double> >::iterator iter = \
        solDeriv_vec_.find(derivType);
      if (iter == solDeriv_vec_.end())
      {
        EXCEPTION("Tried to access derivType: " << derivType << " which is not set!");
      }
      return iter->second;
    }

    /**
     * Checks if time step exists
     * @param timeStepType the time step you want
     * @return true or false
     */
    virtual bool is_SolTimeStep_set(TIMEStepType timeStepType)
    {
      std::map<TIMEStepType, Vector<Double> >::iterator iter = \
        sol_timeStepVec_.find(timeStepType);
      if (iter == sol_timeStepVec_.end())
      {
        return false;
      }
      return true;
    }

    /**
     * Checks if derivative exists
     * @param derivType the time step you want
     * @return true or false
     */
    virtual bool is_Deriv_set(DERIVType derivType)
    {
      std::map<DERIVType, Vector<Double> >::iterator iter = \
        solDeriv_vec_.find(derivType);
      if (iter == solDeriv_vec_.end())
      {
        return false;
      }
      return true;
    }

    /**
     * Gives back the coefficients needed for the time stepping algorithm.
     * The coefficient vector has length zero if it has not been set.
     * @return returns a map of the time stepping coefficients 
     * @warning The coefficients are not set in every time stepping algorithms
     */
    virtual const std::map<TIMEStepType, Double>& GetCoefficients() const
    {
      return sol_timeStepCoeff_;
    }

    /**
     * Gives back the derivatives
     * @return returns a map of the derivatives
     */
    virtual std::map<TIMEStepType, Vector<Double> >& GetTimeStepMap()
    {
      return sol_timeStepVec_;
    }
    /**
     * Gives back the derivatives
     * @return returns a map of the derivatives
     */
    virtual std::map<DERIVType, Vector<Double> >& GetDeriveMap()
    {
      return solDeriv_vec_;
    }

    //! set the time step
    void SetTimeStep(Double dt) 
    { dt_ = dt;};

    //! get the time step size
    inline const Double& GetTimeStep() const
    { return dt_;};

    //! get beta coefficient from Newmark time stepping scheme
    virtual Double GetNewmarkBeta() { 
      EXCEPTION( "Not implemented here" );
      return -1.0;
    }
    
    //! get gamma coefficient from Newmark time stepping scheme
    virtual Double GetNewmarkGamma() {
      EXCEPTION("Not implemented here");
      return(-1.0);
    }

    //! set gamma coefficient from Trapezoidal time stepping scheme
    virtual void SetTrapezoidalGamma(Double val) { 
      EXCEPTION( "Not implemented here" );
    }

    //! Dirichlet boundary condition has to be adapted
    virtual Double DirichletBC4EffMassMatrix(Double val, Integer eq) {
      EXCEPTION( "DirichletBC4EffMassMatrix not implemented" );
      return -1.0;
    }
    /**
     * returns the solution type of the wanted derivative
     * @param solType The solution type from which the name of the first/second
     * drivative is wanted
     * @param derivType The derivative of the solution type from which we want
     * the name (solution type)
     * @return Returns the name of the derivative of a solution type
     */
    SolutionType mapDerivToSolutionType(const SolutionType solType, \
        const DERIVType derivType) const;
    /**
     * Checks if a solution type is actually a derivative of a nother solution.
     * @param solType The type to check if is a derivative
     * @return True or false if the solution type is actually a derivative of a
     * calculated solution type
     * @warning The function may return false negatives, meaning it maybe that a
     * derivative may not be recognised as such. So check it if it supports
     * your derivative.
     */
    bool isDeriv(const SolutionType solType) const;

  protected:

    //! Checks if a given FE-Matrixtype is defined
    bool FeMatrixPresent( FEMatrixType type);

    //! pointer to algebraic system
    BaseSystem * algsys_;

    //! total number of entries in the rhs vector
    UInt rhsSize_;
    
    //! time step size
    Double dt_;     
    
    //! Flag for omitting time step initialization (i.e. first predictor step)
    bool omitFirstPredictor_;

    //! matrix factors
    std::map<FEMatrixType,Double> matrix_factors_;

    //! Derivative and time steps stored in a map
    std::map<TIMEStepType, Vector<Double> > sol_timeStepVec_;
    std::map<DERIVType, Vector<Double> > solDeriv_vec_;
    //! coefficient vector for time stepping
    //! @param sol_timeStepCoeff_, [TIMESTEP_0] is the last calculated time step, [TIMESTEP_1] is the
    //! time step before that and so forth
    std::map<TIMEStepType, Double> sol_timeStepCoeff_;
    
  private:
   

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class TimeStepping
  //! 
  //! \purpose 
  //! This is a base class for all classes, which performa a timestepping
  //! i.e. solve a simple ODE. These are for exmaples trapezoidal rule
  //! and newmark schemes.
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif // FILE_TIMESTEPPING
