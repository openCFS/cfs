// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NEWMARK_2001
#define FILE_NEWMARK_2001

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "MatVec/vector.hh"
#include "timestepping.hh"

namespace CoupledField {
class BaseSystem;
}  // namespace CoupledField

namespace CoupledField
{

  //! class for time stepping of hyperbolic PDE: method is Newmark
  class Newmark: public TimeStepping
  {
  public:
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    Newmark( BaseSystem * algebraicsystem, PtrParamNode systemNode );

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
    void Predictor(Vector<Double>& solold);

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    void UpdateRHS();
    
    //! perform an update to RHS with actual solution (for nonlin calculation)
    virtual void UpdateRHS(Vector<Double>& actSol);  
    
    //! substracts -Ku from RHS
    virtual void SubstractStiffnessFromRHS(Vector<Double>& actSol);
    
    //! perform calculations at end of timestep (set the new timestep)
    void AdvanceTimestep(Vector<Double>& solnew);

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
    Vector<Double> solpred_, solpredSAVE_;

    //! predictor for derivative of solution
    Vector<Double> solderiv1pred_, solderiv1predSAVE_;

    //! solution and first derivative of previous time step;
    Vector<Double> solPrevious_, solderiv1Previous_;
  };




  class NewmarkEffMass: public TimeStepping
  {
  public:
    //! constructor
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    NewmarkEffMass( BaseSystem * algebraicsystem, PtrParamNode systemNode,
                    bool intExplicit = false );

    //! destructor
    virtual ~NewmarkEffMass();
  
    //! initilization
    //! \param rhsSIze total number of entries in the rhs vector
    void Init( Double dt, UInt rhsSize );
    
    //! perform predictor step
    void Predictor(Vector<Double>& solold);

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

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
    Vector<Double> sol_;
    
    //! predictor for nodal solution
    Vector<Double> solpred_;
    
    //! predictor for derivative of solution
    Vector<Double> solderiv1pred_;
  
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
