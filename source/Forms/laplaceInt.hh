// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LAPLACEINT
#define FILE_LAPLACEINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
class LaplaceInt : public BaseForm
{
public:

  LaplaceInt(Double laplVal, bool axi=false, bool coordUpdate = false);

  /** alternative constructor with MaterialDescriptor, laplVal_ is then not used!
   * @see MassInt() */
  LaplaceInt(BaseMaterial* mat, const MaterialDescriptor& md, bool axi=false, bool coordUpdate = false);

  /// 
  virtual ~LaplaceInt();
  
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
