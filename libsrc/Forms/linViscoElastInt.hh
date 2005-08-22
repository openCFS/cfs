#ifndef FILE_LINVISCOELASTINT
#define FILE_LINVISCOELASTINT

#include "linElastInt.hh"

namespace CoupledField
{

  /// class for calculation of the modified stiffness matrix and the material dependent matrix on the right hand side for the viscoelastic material behavior (relaxation and creep) 
class LinViscoElastInt : public linElastInt
{  
public:
  /// Constructor
  LinViscoElastInt(BaseFE * aptelem, MaterialData & matDat, std::string geomType, std::string matrixType,Double timeStep);

  /// Constructor
  LinViscoElastInt(MaterialData & matDat, std::string geomType, std::string matrixType, Double timeStep);

  /// Destructor
  virtual ~LinViscoElastInt();

  
protected:
  
  /// returns D - matrix for linViscoElastInt
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns dimension of D matrix
  virtual UInt getDimD();
  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs();  

private:
  void  GetModifiedMaterialMat(Matrix<Double>& cMat);

  void CalcAInvMat(Matrix<Double>& aInvMat);


  BaseFE * aptelem_;
  std::string geomType_;
  std::string matrixType_;
  Double timeStep_;


private:
  Double dampAlpha_, dampBeta_, fracDeriv_, timeStepPowerFracDeriv_,elastModule_;

};

} //end namespace

#endif //FILE_LINVISCOELASTINT
