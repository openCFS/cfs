#ifndef MATERIAL_DATA
#define MATERIAL_DATA

/**************************************************************************/
/* File:   MaterialData.hh                                            */
/* Author: Fred Hofer & Michael Schinnerl                                 */
/* Date:   26.06.1998                                                     */ 
/*         Rev. 16.12.1999                                                */
/*                                                                        */
/* Stores characteristics of linear materials. The                        */
/* information is read from a file                                        */
/*                                                                        */
/**************************************************************************/

#include <General/environment.hh>
#include <Matrix/matrix.hh>

namespace CoupledField
{

  //! Class for Material Data
  /*! 
    Interface class for getting material data
  */

class MaterialData
{
private:
  Double density;
  Double compressibility;
  Double damp_alfa;
  Double damp_beta;
  Double BoverA;
  Double eModule;
  Double nu;
  Double LameLambda;
  Double LameMu;
//  Double permeability;
//  Double conductivity;
  Double permMx, permMy, permMz;   // permanent magnetization
  Integer scaledMatDat;
  char * name; 
  //char name[stringLength]; 
  std::string bhCurveFile_; //!< name of BH-Curve datafile

  //const ARRAY <NonlinSpline*> & magneticSpline;
  //  ARRAY <NonlinSpline*> * magneticSpline;

  /// contains the stiffnes matrix, the piezelectric coefficients and the permitivity matrix
  Matrix<Double> * piezoMatrix;
  Matrix<Double> * piezoMatrixC;

  Matrix<Double> * permeaMatrix;
  Matrix<Double> * conducMatrix;

  Integer matNr;
  Integer nonlin;
  
public:
  //  MaterialData(const Integer& MatNr, const Integer& Nonlin, const ARRAY<NonlinSpline*> & MagneticSpline);
  
  MaterialData();

  MaterialData(const MaterialData& mat);

  ~MaterialData();
  
  /// set the material number 
  void SetMatNr(const Integer& MatNr){matNr = MatNr; };

  // defines material either as linear or nonlinear and sets the permeability variable if necesarry
  //  void DefLin(const Integer& isLin);

  // in the linear case MagneticSpline has 1 NonlinSpline with 1 factor = const mue !!
  // the splines must express the magntic field strength as function of the induction B!!!
  // $H = a_0 + a_1B + a_2B^2 + ... !
  //  void SetNonLinSpline(ARRAY<NonlinSpline*> * MagneticSpline){magneticSpline = MagneticSpline;};

  /// set the e-module
  void SetEModule(const Double& EModule){eModule = EModule;};

  /// set nu
  void SetNu(const Double& Nu){nu = Nu;};


  /// set density of the material
  void SetDensity(const Double& Density){density = Density;};

  /// set Lame-Paramter lambda
  void SetLambda(const Double& lambda){LameLambda = lambda;};

  /// set Lame-Paramter mu
  void SetMu(const Double& mu){LameMu = mu;};

  /// set permeability of the material
//  void SetPermeability(const Double& aPerm){permeability = aPerm;}

  /// set nonlinearity datafile for BH-curve for magnetic material
  void SetBHCurveFileName( const Char *filename) {
    bhCurveFile_.assign( filename );
  }

  /// set compressibility of the material
  void SetCompressibility(const Double& compr){compressibility = compr;}

  /// set damping coefficients
  void SetDampingCoeffs(const Double& Damp_alfa, const Double& Damp_beta)
    {damp_alfa=Damp_alfa; damp_beta=Damp_beta;};


  /// set BoverA (nonlinearity parameter in nonlinear acoustics)
  void SetBoverA(const Double& BA){BoverA = BA;};


  /// set values of permanent magnetization
  void SetPermMag(const Double& mX, const Double& mY, const Double& mZ)
    {permMx=mX, permMy=mY, permMz=mZ;};

  /// get values of permanent magnetization
  void GetPermMag( Double& mX,  Double& mY, Double& mZ)
  {mX=permMx, mY=permMy, mZ=permMz;};

  /// get nonlinearity datafile name for BH-curve for magnetic material
  std::string& GetBHCurveFileName() { return bhCurveFile_; }

  /// set one value of the data-matrix on position (i,j)
  void SetPiezoMatrixData(const Integer& i, const Integer& j, const Double& value)
    {(*piezoMatrix)(i,j) = value;};

  /// set one value of the data-matrix on position (i,j)
  void SetPiezoMatrixDataC(const Integer& i, const Integer& j, const Double& value)
    {(*piezoMatrixC)(i,j) = value;};

  /// get the value of the data-matrix on position (i,j)
  void GetPiezoMatrixData(const Integer& i, const Integer& j, Double& value)
    {value = (*piezoMatrix)(i,j);};

  /// get the value of the data-matrix on position (i,j)
  void GetPiezoMatrixDataC(const Integer& i, const Integer& j, Double& value)
    {value = (*piezoMatrixC)(i,j);};

  /// return a pointer to the data-matrix
  Matrix<Double> * GetMatrixC(){return piezoMatrixC;};

  /// return a pointer to the data-matrix
  Matrix<Double> * GetMatrix(){return piezoMatrix;};

  /// Rotates piezo Material Matrix. Input are the three solid angels (radian measure)

  //! Rotates the piezoelectric material matrix. Input are the three solid angels (radian measure).
  //! Special case: The choice of 1 for one of the angels, rotates piezo matrix in that way
  //! that it is polarized in the given direction. The choice of 0 does not performs any rotation
  //! Poling in z - direction (a3=1) is given by default.
  void RotateMaterialMatrix(const Double& a1, const Double& a2, const Double& a3);

  //! set one value of the permeability-matrix on position (i,j)
  void SetPermeability(const Integer& i, const Integer& j, const Double& value);

  //! get the value of the permeability-matrix on position (i,j)
  void GetPermeability(const Integer& i, const Integer& j, Double &value);

  //! return a pointer to the permeability matrix
  Matrix<Double> * GetPermeaMatrix(){return permeaMatrix;};


  //! set one value of the permeability-matrix on position (i,j)
  void SetConductivity(const Integer& i, const Integer& j, const Double& value);

  //! get the value of the permeability-matrix on position (i,j)
  void GetConductivity(const Integer& i, const Integer& j, Double &value);

  //! return a pointer to the permeability matrix
  Matrix<Double> * GetConducMatrix(){return conducMatrix;};

  /// set size of data-matrix in x and y direction to nrElems3d x nrElems3d
  /// this matrix includes the stiffness, piezoelectric coupling and permitivity matrix
  void DefFull3dMatrix(){
    piezoMatrix = new Matrix<Double>; 
    piezoMatrix->Resize(GetNrElems3d(), GetNrElems3d() );
    piezoMatrixC = new Matrix<Double>; 
    piezoMatrixC->Resize(GetNrElems3d(), GetNrElems3d() );};

  /// set conductivity of the material
  void SetConductivity(const Double& Conductivity);

  /// set name of the material
  void SetName(const char* Name);

  /// 
  void SetScaledFlag(){scaledMatDat = 1;};

  Integer IsMatDatScaled(){return scaledMatDat;};

  // is material nonlinear?
  //  Integer GetNonlin() const { return nonlin; };

  /// get nr of elems of the matrix in the 3d case
  Integer GetNrElems3d() const { return 9; };

  /// get nr of elems of the matrix in the 2d case
  Integer GetNrElems2d() const { return 5; };

  /// get number of the matNr
  Integer GetMatNr() const { return matNr; };

  /// get E-Module
  Double GetEModule() const {return eModule; };

  /// get density
  Double GetDensity() const {return density;};

  /// get LameLambda
  Double GetLameLambda() const {return LameLambda;}
  
  /// get LameMu
  Double GetLameMu() const {return LameMu;}

  /// get permittivity
  Double GetPermittivity(Integer i, Integer j) const {return (*piezoMatrix)[i+6][j+6];};

  /// get compressibility
  Double GetCompressibility() const {return compressibility;};

  /// get permeability of the material
//  Double GetPermeability() const {return permeability;};

  /// get alfa damping coefficient
  Double GetDampingAlfa() const {return damp_alfa;};

  /// get beta damping coefficient
  Double GetDampingBeta() const {return damp_beta;};  

  /// get BoverA for nonlinear acoustics
  Double GetBoverA() const {return BoverA;}; 

  /// get nu
  Double GetNu() const { return nu; };

  /// get conductivity
//  Double GetConductivity() const {return conductivity; };

  /// get name of the material
  char * GetMaterialName(); 

  //! get density and compressity for acoustic equation from material data-file
   /*!
	\param density out: value of density from the material file
	\param compress out: value of compressity from the material file
	\param matnum material number
	\param keyword name of material in the  material file
  */
  // for compatibility with elena
//   void ReadDensityAndCompressity(Double & density, Double & compress, const Integer matnum, const std::string keyword)
//   { Error("Not implemented",__FILE__,__LINE__);}


  //! get dielectric term for an electrostatic equation from the material file
  /*!
	\param dielectr out: value of dielectric term from the material file
	\param matnum in: material number
  */
  // for compatibility with elena
//   void ReadDielectricTerms(Double & dielectr,const Integer matnum)
//   { Error("Not implemented",__FILE__,__LINE__);}

  // get value for error functional (see Institusbericht Michael Schinnerl: Behandlung gekoppelter
  // Systeme mit Mehrgitterverfahren)
  // Double GetErrorFunctional(const Double& MagB) const;  
};

} // end of namespace
#endif
