#ifndef FILE_PIEZOSTRESSSTRAIN_04
#define FILE_PIEZOSTRESSSTRAIN_04

#include <Elements/basefe.hh>
#include <DataInOut/MaterialData.hh>

#include <Forms/linPiezoInt.hh>
#include <General/environment.hh>

namespace CoupledField
{
  
//! class for calculation of mechanical stresses of piezos
class PiezoStressStrain : public linPiezoInt
{
public:

  /// Constructor
  PiezoStressStrain(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  PiezoStressStrain(MaterialData & matData);
  
  /// Destructor
  virtual ~PiezoStressStrain();  
  
  /// calculates Piola-Kirchoff-stresses (vector notation)
  void CalcStressVec(Vector<Double>& StressVec, Integer ip, Matrix<Double> & ptCoord);  

  /// in stress calculations, the actual displacement of the element is needed
  /*!
  \param disp (input) Matrix with displacement d of all nodes of actual element
  \f[ \left( \begin{array}{ccc} 
             d_{x1} &  d_{x2} &  d_{x3} \\
             d_{y1} &  d_{y2} &  d_{y3} \\
	     \end{array}\right) \f]	    
  */
  virtual void SetActElemSol(Matrix<Double>& disp) {
    ENTER_FCN( "PiezoStressStrain::SetActElemSol" );
    elemDisp_ = disp;};

protected:  
  /// returns D - matrix (material matrix)
  //  virtual void calcDMat(Matrix<Double> & dMat);
  
  /// returns B 
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)=0;

  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs()=0;  

  /// returns nr. of stress components
  virtual Integer getStressComp()=0;  

  /// displacement of all nodes of actual element
  Matrix<Double> elemDisp_;

};
  

//! 3D Cauchy stresses, linear strains
class PiezoStressStrain3D : public PiezoStressStrain
{
public:

  /// Constructor
  PiezoStressStrain3D(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  PiezoStressStrain3D(MaterialData & matData);
  
  /// Destructor
  virtual ~PiezoStressStrain3D();  
  
protected:  
  /// returns D - matrix (material matrix)
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns B 
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 4;};

  /// returns the size of the material d-matrix
  virtual Integer getDimD(){return 9;};

  /// returns nr. of stress components
  virtual Integer getStressComp() { return 6;};  

};
  



//linear plane strain stress
class PiezoStressStrainPlaneStrain : public PiezoStressStrain
{
public:
  /// Constructor
  PiezoStressStrainPlaneStrain(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  PiezoStressStrainPlaneStrain(MaterialData & matData);
  
  /// Destructor
  virtual ~PiezoStressStrainPlaneStrain();  
  
protected:  
  /// returns D - matrix (material matrix)
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns B 
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};  

  /// returns the size of the material d-matrix
  virtual Integer getDimD(){return 5;};

 /// returns nr. of stress components
 virtual Integer getStressComp() { return 3;};  
};
  

//linear stress for axisymmetric case
class PiezoStressStrainAxi : public PiezoStressStrain
{
public:
  /// Constructor
  PiezoStressStrainAxi(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  PiezoStressStrainAxi(MaterialData & matData);
  
  /// Destructor
  virtual ~PiezoStressStrainAxi();  
  
protected:  
  /// returns D - matrix (material matrix)
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns B 
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};  

  /// returns the size of the material d-matrix
  virtual Integer getDimD(){return 6;};

 /// returns nr. of stress components
 virtual Integer getStressComp() { return 4;};  
};
  

} //end namespace

#endif // FILE_XXX


