#ifndef FILE_CURLCURLNODEINT
#define FILE_CURLCURLNODEINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
class CurlCurlNode2DInt : public BaseForm
{
public:
  /// Constructor
 CurlCurlNode2DInt(BaseFE * aptelem, Double laplVal, Boolean axi=FALSE);

  /// Constructor
 CurlCurlNode2DInt(Double laplVal, Boolean axi=FALSE);

  /// 
  virtual ~ CurlCurlNode2DInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected: 
private:
  /// multiplicative value for laplace integration 
  Double matVal_;
};


}

#endif // FILE_CURLCURLNODEINT
