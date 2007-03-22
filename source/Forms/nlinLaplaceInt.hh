// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINLAPLACEINT
#define FILE_NLINLAPLACEINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a nonlinear laplace operator
  /// the multiplicative factor is a function of time

class nLinLaplaceInt : public BaseForm
{
public:

  /// Constructor
  nLinLaplaceInt(Double laplVal, bool axi=false, bool isSpeedVariable = false, 
                 bool coordUpdate = false );

  /// 
  virtual ~nLinLaplaceInt();
  
  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );



private:
  /// multiplicative value for laplace integration 
  Double laplVal_;
};


}

#endif // FILE_LAPLACEINT
