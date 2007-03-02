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
  nLinCurlCurlNode2DInt(ApproxData *nlinFnc, Double startVal, 
                        bool axi=false, bool coordUpdate = false );

  /// 
  virtual ~nLinCurlCurlNode2DInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix( Matrix<Double>& stiffMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
    
  //
  virtual void SetNonLinMethod(std::string atype);

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
