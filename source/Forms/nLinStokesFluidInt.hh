// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINSTOKESFLUIDINT
#define FILE_NLINSTOKESFLUIDINT

#include <Elements/basefe.hh>

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

  /// element matrix containing the convective term (v*dv/dx)
//  virtual void calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord);

  /// displacement of all nodes of actual element
  Matrix<Double> elemVelocity_;

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
  nLinStokesFluid3DInt_Convective(Double density, Double dynamicViscosity);

  
  /// Destructor
  virtual ~nLinStokesFluid3DInt_Convective();  
  
protected:  
  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
                                 
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
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
                                 
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
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
                                 
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
                             bool axi=false);

  /// destructor
  virtual ~nLinStokesFluid_linFormInt();

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Vector<Double> & elemVec, EntityIterator& it);

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

