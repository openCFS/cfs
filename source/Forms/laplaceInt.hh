#ifndef FILE_LAPLACEINT
#define FILE_LAPLACEINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
class LaplaceInt : public BaseForm
{
public:

  /// Constructor
  LaplaceInt(Double laplVal, bool axi=false, 
             bool coordUpdate = false );

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
