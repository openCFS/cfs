// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LAPLACEXYZINT
#define FILE_LAPLACEXYZINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a splitted laplace operator
class LaplaceXYZInt : public BaseForm
{
public:


  /// Constructor
  LaplaceXYZInt(Double laplVal, const UInt nrDofsPerNode, bool axi=false);

  /// 
  virtual ~LaplaceXYZInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );


private:

  //!
  UInt nrDofsPerNode_;

  /// multiplicative value for laplace integration 
  Double laplVal_;
};


}

#endif // FILE_LAPLACEXYZINT
