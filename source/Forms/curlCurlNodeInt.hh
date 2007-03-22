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
  CurlCurlNode2DInt(Double laplVal, bool axi=false, 
                    bool coordUpdate = false );

  /// 
  virtual ~ CurlCurlNode2DInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
  
protected: 
private:
  /// multiplicative value for laplace integration 
  Double matVal_;
};


}

#endif // FILE_CURLCURLNODEINT
