#ifndef FILE_BASEFORMINT_1
#define FILE_BASEFORMINT_1

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
