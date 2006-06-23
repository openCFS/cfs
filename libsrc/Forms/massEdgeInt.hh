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
  MassEdgeInt(Double acond);

  /// Destructor
  virtual ~MassEdgeInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );


private:
  // multiplicative value for mass integrator
  Double conductivity_;
  
};

}

#endif // FILE_MASSEDGEINT
