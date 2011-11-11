#ifndef FILE_NEWMARK_2001
#define FILE_NEWMARK_2001

#include "TimeStepping.hh"

namespace CoupledField
{
  class ParamNode;

  //! class for time stepping of hyperbolic PDE: method is Newmark
  class Newmark: public TimeStepping
  {
  public:
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    Newmark( AlgebraicSys * algebraicsystem, PtrParamNode systemNode );

    //! destructor
    virtual ~Newmark();
  
    //! initilization
    //! \param rhsSize total number of entries in the rhs vector
    void Init( Double dt, UInt rhsSize );
    
    //! Reinitialize (set all vectors to zero)
    void ReInit();

    void setSubSteps( UInt subSteps );
    void resetDeltaT( );

    //! perform predictor step
    void Predictor(SBM_Vector& solold);

    //! perform corrector step
    void Corrector(SBM_Vector& solnew);

    //! perform an update to RHS
    void UpdateRHS();
    
    //! perform an update to RHS with actual solution (for nonlin calculation)
    virtual void UpdateRHS(SBM_Vector& actSol);  
    
    //! substracts -Ku from RHS
    virtual void SubstractStiffnessFromRHS(SBM_Vector& actSol);
    
    //! perform calculations at end of timestep (set the new timestep)
    void AdvanceTimestep(SBM_Vector& solnew);

    Double GetNewmarkBeta(){
      return beta_;
    }
    
    Double GetNewmarkGamma(){
      return gamma_;
    }


  private:

    //! compute parameters for multiplication
    void CalcParameters(Double dt);
    Double globalDeltaT_;
    
    //@{
    //! integration parameters
    Double alpha_, gamma_, beta_, nu_;  
    //@}

    //! predictor for nodal solution
    SBM_Vector solpred_, solpredSAVE_;

    //! predictor for derivative of solution
    SBM_Vector solderiv1pred_, solderiv1predSAVE_;

  };




  class NewmarkEffMass: public TimeStepping
  {
  public:
    //! constructor
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    NewmarkEffMass( AlgebraicSys * algebraicsystem, PtrParamNode systemNode,
                    bool intExplicit = false );

    //! destructor
    virtual ~NewmarkEffMass();
  
    //! initilization
    //! \param rhsSIze total number of entries in the rhs vector
    void Init( Double dt, UInt rhsSize );
    
    //! perform predictor step
    void Predictor(SBM_Vector& solold);

    //! perform corrector step
    void Corrector(SBM_Vector& solnew);

    //! perform an update to RHS
    void UpdateRHS();

    //! Dirichlet boundary condition has to be adapted
    Double DirichletBC4EffMassMatrix(Double val, Integer eq);

    //! set the time step
    void SetTimeStep(Double dt) 
    { dt_ = dt;};
  

  private:

    //! compute parameters for multiplication
    void CalcParameters(Double dt);
   
    //! time step
    Double dt_;
    
    //@{
    //! integration parameters
    Double alpha_, gamma_, beta_, nu_;
    //@}

    //! nodal solution
    SBM_Vector sol_;
    
    //! predictor for nodal solution
    SBM_Vector solpred_;
    
    //! predictor for derivative of solution
    SBM_Vector solderiv1pred_;
  
    // if true, we perform explicit time stepping
    bool intExplicit_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class Newmark
  //! 
  //! \purpose 
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

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class NewmarkEffMass
  //! 
  //! \purpose 
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




}

#endif // FILE_NEWMARK
