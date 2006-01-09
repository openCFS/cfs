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
  LaplaceXYZInt(BaseFE * aptelem, Double laplVal, Boolean axi=FALSE);

  /// Constructor
  LaplaceXYZInt(Double laplVal, const UInt nrDofsPerNode, Boolean axi=FALSE);

  /// 
  virtual ~LaplaceXYZInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);


  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

  
  virtual void SetActElemSol(Matrix<Double>& disp){};


private:

  //!
  UInt nrDofsPerNode_;

  /// multiplicative value for laplace integration 
  Double laplVal_;
};


}

#endif // FILE_LAPLACEXYZINT
