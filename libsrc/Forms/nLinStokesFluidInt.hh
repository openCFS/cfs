#ifndef FILE_NLINSTOKESFLUIDINT
#define FILE_NLINSTOKESFLUIDINT

#include <Elements/basefe.hh>
#include <DataInOut/MaterialData.hh>

#include <Forms/stokesFluidInt.hh>
#include <General/environment.hh>

namespace CoupledField
{

  // =============================================================================
  // base class for nonlinear stokes fluid
  // =============================================================================

  /// base class for calculation of nonlinear stokes fluid flows
class nLinStokesFluidInt : public StokesFluidInt
{
public:
  /// Constructor
  //nLinStokesFluidInt(BaseFE * aptelem, Double density, Double dynamicViscosity);

  /// Constructor
  nLinStokesFluidInt(Double density, Double dynamicViscosity);
  
  /// Destructor
  virtual ~nLinStokesFluidInt();  

  /// in nonlinear calculations, the actual flow velocities of the element is needed
  /*!
    \param velocity (input) Matrix with velocities v of all nodes of actual element
    \f[ \left( \begin{array}{ccc} 
    v_{x1} &  v_{x2} &  v_{x3} \\
    v_{y1} &  v_{y2} &  v_{y3} \\
    \end{array}\right) \f]         
  */
  void setActElemVelocity(Matrix<Double>& velocity) {elemVelocity_ = velocity;};  

  /// in nonlinear calculations, the actual velocities of the element is needed
  /*!
    \param velocity (input) Matrix with velocities v of all nodes of actual element
    \f[ \left( \begin{array}{ccc} 
    v_{x1} &  v_{x2} &  v_{x3} \\
    v_{y1} &  v_{y2} &  v_{y3} \\
    \end{array}\right) \f]         
  */
  virtual void SetActElemSol(Matrix<Double>& velocity) {elemVelocity_ = velocity;};


protected:    
  /// Calculation of stiffmess matrix
  virtual void CalcElementMatrix(Matrix<Double> & ptCoord, 
                                 Matrix<Double> & elemMat) {
    Error( "nLinStokesFluidInt::CalcElementMatrix() not correctly overwritten!",
             __FILE__, __LINE__);
  };

  /// element matrix containing the convective term (v*dv/dx)
//  virtual void calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord);

  /// displacement of all nodes of actual element
  Matrix<Double> elemVelocity_;


  char * className;
};
  
  // =============================================================================
  // 3D nonlinear stokes fluid
  // =============================================================================

 /// class for calculation of 3d convective nonlinear stokes fluid
//  nonlinear convective term
class nLinStokesFluid3DInt_Convective : public nLinStokesFluidInt
{
public:

  /// Constructor
//  nLinStokesFluid3DInt_Convective(BaseFE * aptelem, 
//                                  Double density, Double dynamicViscosity);

  /// Constructor
  nLinStokesFluid3DInt_Convective(Double density, Double dynamicViscosity);

  
  /// Destructor
  virtual ~nLinStokesFluid3DInt_Convective();  
  
protected:  
/// Calculation of stiffmess matrix
  virtual void CalcElementMatrix(Matrix<Double> & ptCoord,
                                 Matrix<Double> & elemMat); 
                                 
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 8;};
};


// =============================================================================
// nonlinear plane stokes fluid (2D)
// =============================================================================

/// class for calculation of plane strain nonlinear elasticity
// first part: nonlinear B-matrix
class nLinStokesFluidPlaneInt_Convective : public nLinStokesFluidInt
{
public:

  /// Constructor
  nLinStokesFluidPlaneInt_Convective(Double density, Double dynamicViscosity);

  
  /// Destructor
  virtual ~nLinStokesFluidPlaneInt_Convective();  
  
protected:  
/// Calculation of stiffmess matrix
  virtual void CalcElementMatrix(Matrix<Double> & ptCoord,
                                 Matrix<Double> & elemMat); 
                                 
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 4;};  
};



// =============================================================================
// nonlinear axisymmetrical stokes fluid
// =============================================================================

/// class for calculation of axisymmetric nonlinear elasticity
// first part: nonlinear B-matrix
class nLinStokesFluidAxiInt_Convective : public nLinStokesFluidInt
{
public:

  /// Constructor
  nLinStokesFluidAxiInt_Convective(Double density, Double dynamicViscosity);

  
  /// Destructor
  virtual ~nLinStokesFluidAxiInt_Convective();  
  
protected:  
/// Calculation of stiffmess matrix
  virtual void CalcElementMatrix(Matrix<Double> & ptCoord,
                                 Matrix<Double> & elemMat); 
                                 
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 4;};  
};

/// class for calculation of right-hand-side of nonlinear stokes fluid
class nLinStokesFluid_linFormInt : public nLinStokesFluidInt
{
public:
  /// constructor
  nLinStokesFluid_linFormInt(Double density,
                             Double dynamicViscosity,
                             Boolean axi=FALSE);

  /// destructor
  virtual ~nLinStokesFluid_linFormInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result);

  /// in nonlinear calculations, the actual velocity of the element is needed
  /*!
    \param velocity (input) Matrix with velocities v of all nodes of actual element
    \f[ \left( \begin{array}{ccc} 
    v_{x1} &  v_{x2} &  v_{x3} \\
    v_{y1} &  v_{y2} &  v_{y3} \\
    \end{array}\right) \f]         
  */
  void setActElemVelocity(Matrix<Double>& velocity) {elemVelocity_ = velocity;};  


  /// in nonlinear calculations, the actual velocity of the element is needed
  /*!
    \param velocity (input) Matrix with velocities v of all nodes of actual element
    \f[ \left( \begin{array}{ccc} 
    v_{x1} &  v_{x2} &  v_{x3} \\
    v_{y1} &  v_{y2} &  v_{y3} \\
    \end{array}\right) \f]         
  */
  virtual void SetActElemSol(Matrix<Double>& velocity) {elemVelocity_ = velocity;};
 
  
protected:
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 3;};

  /// material data
  Double density_;
  Double dynamicViscosity_;

  /// velocities of all nodes of actual element
  Matrix<Double> elemVelocity_;
};

} //end namespace

#endif // FILE_XXX

