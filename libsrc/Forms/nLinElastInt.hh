#ifndef FILE_NLINELASTINT_03
#define FILE_NLINELASTINT_03

#include <Elements/basefe.hh>
#include <DataInOut/MaterialData.hh>

#include <Forms/linElastInt.hh>
#include <General/environment.hh>

namespace CoupledField
{
  


  // =============================================================================
  // base class for nonlinear mechanics
  // =============================================================================

  /// base class for calculation of nonlinear elasticity
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
  virtual void calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord);

  /// calcs the material matrix for the 2d case
  virtual void Calc2DMaterialMatrix(Matrix<Double> & dMat, enum orientation2D actOrientation);


  /// displacement of all nodes of actual element
  Matrix<Double> elemDisp_;


  char * className;
};
  
  
 

  // =============================================================================
  // 3D nonlinear mechanics
  // =============================================================================


  /// class for calculation of 3d nonlinear elasticity
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
  virtual UInt getDimD(){return 6;};
  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 3;};
};




/// class for calculation of 3d nonlinear elasticity
// second part: regarding internal stresses
class nLinMechInt_PiolaStress : public nLinElastInt
{
public:
  friend class nLinMech_linFormInt;
  

  /// Constructor
  nLinMechInt_PiolaStress(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinMechInt_PiolaStress(MaterialData & matData);
  
  /// Destructor
  virtual ~nLinMechInt_PiolaStress();  
  
protected:  
  /// returns D - matrix for BDB (size: 9x9, contains the 2. Piola-Kirchhoff-Stress tensor!!)
  virtual void calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord);
  
  /// returns B - matrix for BDB
  virtual void calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord);

  /// returns dimension of D matrix
  virtual UInt getDimD(){return piolaDimD_;};
  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs()=0;  


protected:
  /// sets the size of d-matrix (needed for lin and nonlin B-matrix)
  virtual void setPiolaDimD(UInt actDim){piolaDimD_ = actDim;};

  /// calculates Piola-Kirchoff-stresses (vector notation)
  virtual void CalcStressVec(Vector<Double>& piolaStressVec, UInt ip, Matrix<Double> & ptCoord);  

  /// returns linear B - matrix
  virtual void calcLinBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord);

  /// returns nonlinear B - matrix
  virtual void calcNonLinBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord);

  /// returns material D-matrix for 3d mechanics
  virtual void calcMaterialDMat(Matrix<Double> & dMat);

  /// returns the size of the material d-matrix
  virtual UInt getMaterialDMatSize()=0;

  /// returns the size of the full piola d-matrix
  virtual UInt getFullPiolaDMatSize()=0;

  /// conversion of stress vector to stress tensor
  virtual void convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress);
  
private:
  
  /// dimension of d-matrix (has to be changed for some dirty implementation features ... )
  UInt piolaDimD_;
};
  


/// class for calculation of 3d nonlinear elasticity
// second part: regarding internal stresses
class nLinMech3dInt_PiolaStress : public nLinMechInt_PiolaStress
{
public:

  /// Constructor
  nLinMech3dInt_PiolaStress(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinMech3dInt_PiolaStress(MaterialData & matData);
  
  /// Destructor
  virtual ~nLinMech3dInt_PiolaStress();  
  
protected:  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 3;};


  /// returns the size of the material d-matrix
  virtual UInt getMaterialDMatSize(){return 6;};


  /// returns the size of the full piola d-matrix
  virtual UInt getFullPiolaDMatSize(){return 9;};
};
  





// =============================================================================
// nonlinear plane strain mechanics (2D)
// =============================================================================

/// class for calculation of plane strain nonlinear elasticity
// first part: nonlinear B-matrix
class nLinMechPlaneStrainInt_BNonLin : public nLinElastInt
{
public:

  /// Constructor
  nLinMechPlaneStrainInt_BNonLin(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinMechPlaneStrainInt_BNonLin(MaterialData & matData);

  
  /// Destructor
  virtual ~nLinMechPlaneStrainInt_BNonLin();  
  
protected:  
  /// returns D - matrix for BDB
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns dimension of D matrix
  virtual UInt getDimD(){return 3;};
  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 2;};  
};




/// class for calculation of 2d nonlinear elasticity
// second part: regarding internal stresses
class nLinMechPlaneStrainInt_PiolaStress : public nLinMechInt_PiolaStress
{
public:
  /// Constructor
  nLinMechPlaneStrainInt_PiolaStress(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinMechPlaneStrainInt_PiolaStress(MaterialData & matData);
  
  /// Destructor
  virtual ~nLinMechPlaneStrainInt_PiolaStress();  
  
protected:  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 2;};  

  /// conversion of stress vector to stress tensor
  virtual void convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress);

  /// returns material D-matrix for 3d mechanics
  virtual void calcMaterialDMat(Matrix<Double> & dMat);

  /// returns the size of the material d-matrix
  virtual UInt getMaterialDMatSize(){return 3;};

  /// returns the size of the full piola d-matrix
  virtual UInt getFullPiolaDMatSize(){return 4;};
};
  







// =============================================================================
// nonlinear axisymmetrical mechanics
// =============================================================================

/// class for calculation of axisymmetric nonlinear elasticity
// first part: nonlinear B-matrix
class nLinMechAxiInt_BNonLin : public nLinElastInt
{
public:

  /// Constructor
  nLinMechAxiInt_BNonLin(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinMechAxiInt_BNonLin(MaterialData & matData);

  
  /// Destructor
  virtual ~nLinMechAxiInt_BNonLin();  
  
protected:  
  /// returns D - matrix for BDB
  virtual void calcDMat(Matrix<Double> & dMat);

  /// returns dimension of D matrix
  virtual UInt getDimD(){return 4;};
  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 2;};  
};




/// class for calculation of 2d axisymmetric nonlinear elasticity
// second part: regarding internal stresses
class nLinMechAxiInt_PiolaStress : public nLinMechInt_PiolaStress
{
public:
  /// Constructor
  nLinMechAxiInt_PiolaStress(BaseFE * aptelem, MaterialData & matData);

  /// Constructor
  nLinMechAxiInt_PiolaStress(MaterialData & matData);
  
  /// Destructor
  virtual ~nLinMechAxiInt_PiolaStress();  
  
protected:  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 2;};  

  /// conversion of stress vector to stress tensor
  virtual void convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress);

  /// returns material D-matrix for 3d mechanics
  virtual void calcMaterialDMat(Matrix<Double> & dMat);

protected:
  /// returns the size of the material d-matrix
  virtual UInt getMaterialDMatSize(){return 4;};

  /// returns the size of the full piola d-matrix
  virtual UInt getFullPiolaDMatSize(){return 5;};  
};
  






// =============================================================================
// prestress
// =============================================================================



/// class for regarding 3d prestress
class PreStressInt : public nLinMechInt_PiolaStress
{
public:
  // preStressLinFormInt uses calcDMat from this class
  friend class PreStressLinFormInt;
  
  /// Constructor
  PreStressInt(BaseFE * aptelem, MaterialData & matData, Vector<Double> aPreStressVal);
  
  //! Constructor
  PreStressInt(MaterialData & matData, Vector<Double> aPreStressVal);
  
  /// Destructor
  virtual ~PreStressInt();  

  //! computation of D-matrix (consists of stress values)
  void calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord);
  
  //!
  void calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord);

  //!
  void convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {Error("convertStressVecToTensor not implemented");}

protected:
  /// calculates pre-stresses (vector notation)
  void CalcStressVec(Vector<Double>& piolaStressVec, UInt ip, Matrix<Double> & ptCoord);

  /// returns the size of the full piola d-matrix
  virtual UInt getFullPiolaDMatSize(){
    Error( "Not implemented here", __FILE__, __LINE__ );
    return 0;
  }

  /// returns dimension of D matrix
  virtual UInt getDimD() {
    Error( "Not implemented here", __FILE__, __LINE__ );
    return 0;
  }
  
  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){
    Error( "Not implemented here", __FILE__, __LINE__ );
    return 0;
  }

  /// returns the size of the material d-matrix
  virtual UInt getMaterialDMatSize() {
    Error( "Not implemented here", __FILE__, __LINE__ );
    return 0;
  }


private: 
  /// 
  Vector<Double> preStressVal_;

};


/// class for regarding 3d prestress
class PreStressInt3D : public PreStressInt
{
public:
  /// Constructor
  PreStressInt3D(BaseFE * aptelem, MaterialData & matData, Vector<Double> aPreStressVal);
  
  //! Constructor
  PreStressInt3D(MaterialData & matData, Vector<Double> aPreStressVal);
  
  /// Destructor
  virtual ~PreStressInt3D();  

  //!
  void convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress);

protected:
  /// returns the size of the full piola d-matrix
  virtual UInt getFullPiolaDMatSize(){return 9;};

  /// returns dimension of D matrix
  virtual UInt getDimD(){return 6;};

  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 3;};

  /// returns the size of the material d-matrix
  virtual UInt getMaterialDMatSize(){return 6;};

};


/// class for regarding 2d prestress in plane strain case
class PreStressIntPlaneStrain : public PreStressInt
{
public:
  /// Constructor
  PreStressIntPlaneStrain(BaseFE * aptelem, MaterialData & matData, Vector<Double> aPreStressVal);
  
  //! Constructor
  PreStressIntPlaneStrain(MaterialData & matData, Vector<Double> aPreStressVal);
  
  /// Destructor
  virtual ~PreStressIntPlaneStrain();  

  //!
  void convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress);

protected:
  /// returns the size of the full piola d-matrix
  virtual UInt getFullPiolaDMatSize(){return 4;};

  /// returns dimension of D matrix
  virtual UInt getDimD(){return 3;};

  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 2;};

  /// returns the size of the material d-matrix
  virtual UInt getMaterialDMatSize(){return 3;};

};


/// class for regarding 2d axi prestress
class PreStressIntAxi : public PreStressInt
{
public:
  /// Constructor
  PreStressIntAxi(BaseFE * aptelem, MaterialData & matData, Vector<Double> aPreStressVal);
  
  //! Constructor
  PreStressIntAxi(MaterialData & matData, Vector<Double> aPreStressVal);
  
  /// Destructor
  virtual ~PreStressIntAxi();  

  //!
  void convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress);

protected:
  /// returns the size of the full piola d-matrix
  virtual UInt getFullPiolaDMatSize(){return 5;};

  /// returns dimension of D matrix
  virtual UInt getDimD(){return 4;};

  /// returns nr. of degrees of freedom
  virtual UInt getNrDofs(){return 2;};

  /// returns the size of the material d-matrix
  virtual UInt getMaterialDMatSize(){return 4;};

};

  
} //end namespace

#endif // FILE_XXX


