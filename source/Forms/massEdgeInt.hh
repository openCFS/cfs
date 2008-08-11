// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MASSEDGEINT_1
#define FILE_MASSEDGEINT_1

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix
class MassEdgeInt : public BaseForm
{
public:
  /// Constructor
  MassEdgeInt( Double acond, 
               bool scaleByEdgeSize,
               bool coordUpdate = false );

  /// Destructor
  virtual ~MassEdgeInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );


private:
  // multiplicative value for mass integrator
  Double conductivity_;
  
  // flag indicating use of scaling by edge size
  bool scaleByEdgeSize_;
};

}

#endif // FILE_MASSEDGEINT
