#ifndef FILE_NLINCURLCURLNODEINT
#define FILE_NLINCURLCURLNODEINT

#include "Utils/ApproxData.hh"
#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
class nLinCurlCurlNode2DInt : public BaseForm
{
public:
  /// Constructor
  nLinCurlCurlNode2DInt(BaseFE * aptelem, Double laplVal, Boolean axi=FALSE);

  /// Constructor
  nLinCurlCurlNode2DInt(ApproxData *nlinFnc, Double startVal, Boolean axi=FALSE);

  /// 
  virtual ~nLinCurlCurlNode2DInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //
  virtual void SetNonLinMethod(std::string atype);
  
  //
  virtual void SetActElemSol(Matrix<Double>& magPot) 
  { magPot.ConvertToVec_AppendRows(magPot_);};

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected: 
private:
  /// start value for reluctivity
  Double startmatVal_;
  ApproxData *nlinFnc_;
  Vector<Double> magPot_;
  NonLinMethod nonLinType_;

};


}

#endif // FILE_NLINCURLCURLNODEINT
