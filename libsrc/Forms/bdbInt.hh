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
  BDBInt(BaseFE * aptelem, MaterialData & matData);


  /// Destructor
  virtual ~BDBInt();


  /// Function for calculation bdb matrix 
  virtual void CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> &);

protected:
  /// returns B - matrix for BDB
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)=0;

  /// returns D - matrix for BDB
  virtual void calcDMat(Matrix<Double> & dMat)=0;

  /// returns dimension of D matrix
  virtual Integer getDimD()=0;

  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs()=0;

};

}

#endif // FILE_MASSINT
