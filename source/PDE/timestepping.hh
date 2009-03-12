// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_TIMESTEPPING_2001
#define FILE_TIMESTEPPING_2001

#include <map>

#include "General/environment.hh"
#include "Utils/nodestoresol.hh"

namespace CoupledField {


  // forward class declaration
  class BaseSystem;

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

    //! set vector with first derivative
    virtual void SetDeriv1(const Vector<Double> & deriv1) {
      solderiv1_ = deriv1;
      isDeriv1Set_ = true;}
  
    //! set vector with second derivative
    virtual void SetDeriv2(const Vector<Double> & deriv2) {
      solderiv2_ = deriv2;
      isDeriv2Set_ = true; }

    //! set vector with second derivative
    virtual void SetOld1(const Vector<Double> & tn_1) {
      sol_tn_1_ = tn_1;
      isSolTN1Set_ = true; }

    //!  return pointer to vector with first derivative of solution
    virtual const Vector<Double>& GetDeriv1() const { return solderiv1_;}
  
    //! return pointer to vector with second derivative of solution
    virtual const Vector<Double>& GetDeriv2() const { return solderiv2_;}

    virtual const Vector<Double>& GetOld1() const { return sol_tn_1_;}

    //! set the time step
    void SetTimeStep(Double dt) 
    { dt_ = dt;};

    //! get the time step size
    Double GetTimeStep() 
    { return dt_;};

    //! get beta coefficient from Newmark time stepping scheme
    virtual Double GetNewmarkBeta() { 
      EXCEPTION( "Not implemented here" );
      return -1.0;
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

  protected:

    //! Checks if a given FE-Matrixtype is defined
    bool FeMatrixPresent( FEMatrixType type);

    //! pointer to algebraic system
    BaseSystem * algsys_;

    //! total number of entries in the rhs vector
    UInt rhsSize_;
    
    //! time step size
    Double dt_;     

    //! matrix factors
    std::map<FEMatrixType,Double> matrix_factors_;

    //! first and second time derivative of solution
    Vector<Double> solderiv1_, solderiv2_, sol_tn_1_, sol_tn_2_;

    //! Flag indicating if 1st derivative was already set from outside
    bool isDeriv1Set_;

    //! Flag indicating if 2nd derivative was already set from outside
    bool isDeriv2Set_;

    //! Flag indicating if soltution of last time step is already set due to a restart file
    bool isSolTN1Set_;

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
