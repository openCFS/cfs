// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CURLCURLNODEINT
#define FILE_CURLCURLNODEINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
class CurlCurlNode2DInt : public BaseForm
{
public:

  /// Constructor
  CurlCurlNode2DInt( BaseMaterial* matData, bool axi=false, 
                     bool coordUpdate = false );

  /// 
  virtual ~ CurlCurlNode2DInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );

  // returns the B operator matrix
  void calcBMat( Matrix<Double> &rotMat,
                 UInt ip, Matrix<Double> &ptCoord );
  
protected: 
private:
  /// multiplicative value for laplace integration 
  Double matVal_;
};



  /// Class for calculation  element stiffness matrix of curl-curl operator
class CurlCurlNode3DInt : public BaseForm
{
public:

  //! Constructor 
  CurlCurlNode3DInt( BaseMaterial* matData, bool coordUpdate = false );

  //! Destructor 
  virtual ~CurlCurlNode3DInt();

  //! Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
  
  //! returns curl and div - matrix
  void calcBMat( Matrix<Double> &bMatCurl, Matrix<Double> &bMatDiv,
		 UInt ip, Matrix<Double> &ptCoord );

protected: 
  UInt nrDofs_;

private:

  //! multiplicative value for curl-curl operator 
  Double matVal_;

  //! true, if we have an orthotropic material
  bool isOrthotropic_;

  //! contains the orthotropic reluctivities
  Vector<Double> reluctivityVec_;

};


  /// Class for calculation  of coupling matrix between vector and scalar potential
class MagCoupVectorScalarPotentialInt : public BaseForm
{
public:

  /// Constructor
 MagCoupVectorScalarPotentialInt(Double elecCond, bool coordUpdate = false );

  /// 
  virtual ~ MagCoupVectorScalarPotentialInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
  
protected: 
private:

  /// multiplicative value for curlc-curl operator 
  Double matVal_;

  UInt nrDofsVec_;
};


}

#endif // FILE_CURLCURLNODEINT
