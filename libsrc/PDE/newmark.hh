#ifndef FILE_NEWMARK_2001
#define FILE_NEWMARK_2001

#include "timestepping.hh"


namespace CoupledField
{

  //! class for time stepping of hyperbolic PDE: method is Newmark

  class Newmark: public TimeStepping
  {
  public:
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    //! \param rhsSIze total number of entries in the rhs vector
    Newmark( BaseSystem * algebraicsystem, UInt rhsSize );

    //! destructor
    virtual ~Newmark();
  
    //! initilization
    void Init( std::map<FEMatrixType,Double> & matrix_factors,
               Double dt );

    //! perform predictor step
    void Predictor(Vector<Double>& solold);

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    void UpdateRHS();

    //! perform an update to RHS with actual solution (for nonlin calculation)
    virtual void UpdateRHS(Vector<Double>& actSol);  



  private:

    //! compute parameters for multiplication
    void CalcParameters(Double dt);

    //@{
    //! integration parameters
    Double alpha_, gamma_, beta_;  
    //@}

    //@{
    //! coefficients from Newmark method
    Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_;
    //@}

    //! flag indicating if damping matrix is needed
    Boolean damping_;

    //! predictor for nodal solution
    Vector<Double> solpred_;

    //! predictor for derivative of solution
    Vector<Double> solderiv1pred_;
  };




  class NewmarkEffMass: public TimeStepping
  {
  public:
    //! constructor
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    //! \param rhsSIze total number of entries in the rhs vector
    NewmarkEffMass( BaseSystem * algebraicsystem, UInt rhsSize );

    //! destructor
    virtual ~NewmarkEffMass();
  
    //! initilization
    void Init( std::map<FEMatrixType,Double> & matrix_factors, 
               Double dt );
    
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
    Double alpha_, gamma_, beta_;
    //@}
    
    //@{
    //! coefficients from Newmark method
    Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_; 
    //@}

    //! flag indicating if damping matrix is needed
    Boolean damping_;
    
    //! nodal solution
    Vector<Double> sol_;
    
    //! predictor for nodal solution
    Vector<Double> solpred_;
    
    //! predictor for derivative of solution
    Vector<Double> solderiv1pred_;
  

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
