#ifndef FILE_NLINELASTINT_03
#define FILE_NLINELASTINT_03

#include <Elements/basefe.hh>
#include <DataInOut/MaterialData.hh>

#include <Forms/linElastInt.hh>

namespace CoupledField
{
  

  /// base class for calculation of nonlinear linear elasticity
class nLinElastInt : public linElastInt
{
public:
  /// Constructor
  nLinElastInt(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinElastInt(MaterialData & matData);
  
  /// Destructor
  virtual ~nLinElastInt();  

  /// in nonlinear calculations, the actual displacement of the element is needed
  /*!
  \param disp (input) Matrix with displacement d of all nodes of actual element
  \f[ \left( \begin{array}{ccc} 
             d_{x1} &  d_{x2} &  d_{x3} \\
             d_{y1} &  d_{y2} &  d_{y3} \\
	     \end{array}\right) \f]	    
  */
  void setActElemDispl(Matrix<Double>& disp) {elemDisp_ = disp;};  

  /// in nonlinear calculations, the actual displacement of the element is needed
  /*!
  \param disp (input) Matrix with displacement d of all nodes of actual element
  \f[ \left( \begin{array}{ccc} 
             d_{x1} &  d_{x2} &  d_{x3} \\
             d_{y1} &  d_{y2} &  d_{y3} \\
	     \end{array}\right) \f]	    
  */
  virtual void SetActElemSol(Matrix<Double>& disp) {elemDisp_ = disp;};


protected:    

  /// returns B - matrix for BDB
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);


  /// displacement of all nodes of actual element
  Matrix<Double> elemDisp_;
};
  
  
 


/// class for calculation of 3d nonlinear linear elasticity
// derived from nLinElastInt with 3d material matrix 
// first part: nonlinear B-matrix
class nLinMech3dInt_BNonLin : public nLinElastInt
{
public:

  /// Constructor
  nLinMech3dInt_BNonLin(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinMech3dInt_BNonLin(MaterialData & matData);

  
  /// Destructor
  virtual ~nLinMech3dInt_BNonLin();  
  
protected:  
  /// returns D - matrix for BDB
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns dimension of D matrix
  virtual Integer getDimD(){return 6;};
  
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};  
};
  



/// class for calculation of 3d nonlinear linear elasticity
// derived from nLinMech3dInt_BNonLin (3d material matrix is needed)
// second part: regarding internal stresses
class nLinMech3dInt_PiolaStress : public nLinMech3dInt_BNonLin
{
public:
  friend class nLinMech_linFormInt;
  

  /// Constructor
  nLinMech3dInt_PiolaStress(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinMech3dInt_PiolaStress(MaterialData & matData);
  
  /// Destructor
  virtual ~nLinMech3dInt_PiolaStress();  
  
protected:  
  /// returns D - matrix for BDB (size: 9x9, contains the 2. Piola-Kirchhoff-Stress tensor!!)
  virtual void calcDMat(Matrix<Double> & dMat, Integer ip, Matrix<Double> & ptCoord);
  
  /// returns B - matrix for BDB
  virtual void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

  /// returns dimension of D matrix
  virtual Integer getDimD(){return piolaDimD_;};
  
  /// returns nr. of degrees of freedom
  virtual Integer getNrDofs(){return 3;};  


protected:
  /// sets the size of d-matrix (needed for lin and nonlin B-matrix)
  virtual void setPiolaDimD(Integer actDim){piolaDimD_ = actDim;};

  /// calculates Piola-Kirchoff-stresses (vector notation)
  void CalcStressVec(std::vector<Double>& piolaStressVec, Integer ip, Matrix<Double> & ptCoord);  

  /// returns linear B - matrix
  virtual void calcLinBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

  /// returns nonlinear B - matrix
  virtual void calcNonLinBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord);

  /// returns material D-matrix for 3d mechanics
  virtual void calcMaterialDMat(Matrix<Double> & dMat);


private:
  /// conversion of stress vector to stress tensor
  void convertStressVecToTensor(Matrix<Double>& stressTensor, std::vector<Double>& piolaStress);
  
  /// dimension of d-matrix (has to be changed for some dirty implementation features ... )
  Integer piolaDimD_;
  
};
  



/// class for regarding 3d prestress
class PreStressInt : public nLinMech3dInt_PiolaStress
{
public:
  // preStressLinFormInt uses calcDMat from this class
  friend class PreStressLinFormInt;
  
  /// Constructor
  PreStressInt(BaseFE * aptelem, MaterialData & matData, Double aPreStressVal, Directions stressDir);
  
  
  /// Destructor
  virtual ~PreStressInt();  
  
protected:
  /// calculates pre-stresses (vector notation)
  void CalcStressVec(std::vector<Double>& piolaStressVec, Integer ip, Matrix<Double> & ptCoord);


private: 
  /// 
  Double preStressVal_;

  /// direction of stresses (enumeration type)
  Directions preStressDir_;
};




  
} //end namespace

#endif // FILE_XXX


