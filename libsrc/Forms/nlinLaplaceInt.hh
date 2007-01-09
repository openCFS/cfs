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
  nLinLaplaceInt(Double laplVal, bool axi=false, 
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
