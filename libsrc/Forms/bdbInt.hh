#ifndef FILE_BDBINT
#define FILE_BDBINT

#include "baseForm.hh" 

namespace CoupledField
{

  /// abstract class for calculation of bdb forms
class BDBInt : public BaseForm
{

public:
  /// Constructor with pointer to BaseElem
  BDBInt(BaseFE * aptelem);


  /// Destructor
  virtual ~BDBInt();


  /// Function for calculation bdb matrix 
  virtual void CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & );

protected:
  /// returns B - matrix for BDB
  virtual void getBMat(Matrix<Double> & bMat);

  /// returns D - matrix for BDB
  virtual void getDMat(Matrix<Double> & dMat);
  

};

}

#endif // FILE_MASSINT
