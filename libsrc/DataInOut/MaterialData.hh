#ifndef MATERIAL_DATA
#define MATERIAL_DATA

/**************************************************************************/
/* File:   LoadMaterialData.hh                                            */
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

class MaterialData
{
private:
  Double density;
  Double damp_alfa;
  Double damp_beta;
  Double eModule;
  Double nu;
  Double constPerm;
  Double conductivity;
  Integer scaledMatDat;
  char * name; 
  //char name[stringLength]; 

  //const ARRAY <NonlinSpline*> & magneticSpline;
  //  ARRAY <NonlinSpline*> * magneticSpline;

  /// contains the stiffnes matrix, the piezelectric coefficients and the permitivity matrix
  Matrix<Double> * piezoMatrix;

  Integer matNr;
  Integer nonlin;
  
public:
  //  MaterialData(const Integer& MatNr, const Integer& Nonlin, const ARRAY<NonlinSpline*> & MagneticSpline);
  
  MaterialData();
  
  /// set the material number 
  void SetMatNr(const Integer& MatNr){matNr = MatNr; };

  // defines material either as linear or nonlinear and sets the constPerm variable if necesarry
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

  /// set damping coefficients
  void SetDampingCoeffs(const Double& Damp_alfa, const Double& Damp_beta)
    {damp_alfa=Damp_alfa; damp_beta=Damp_beta;};

  /// set one value of the data-matrix on position (i,j)
  void SetMatrixData(const Integer& i, const Integer& j, const Double& value)
    {(*piezoMatrix)(i,j) = value;};

  /// get the value of the data-matrix on position (i,j)
  void GetMatrixData(const Integer& i, const Integer& j, Double& value)
    {value = (*piezoMatrix)(i,j);};

  /// return a pointer to the data-matrix
  Matrix<Double> * GetMatrix(){return piezoMatrix;};

  /// set size of data-matrix in x and y direction to nrElems3d x nrElems3d
  /// this matrix includes the stiffness, piezoelectric coupling and permitivity matrix
  void DefFull3dMatrix(){piezoMatrix = new Matrix<Double>; piezoMatrix->Resize(GetNrElems3d(), GetNrElems3d() );};

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

  /// get alfa damping coefficient
  Double GetDampingAlfa() const {return damp_alfa;};

  /// get beta damping coefficient
  Double GetDampingBeta() const {return damp_beta;};  

  /// get nu
  Double GetNu() const { return nu; };

  /// get conductivity
  Double GetConductivity() const {return conductivity; };

  // get nonlinear permiability
  //Double GetPermiability(const Double& MagB) const;

  /// get const permiability - if material is linear!
  Double GetConstPermiability() const {return constPerm;};

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
  void ReadDensityAndCompressity(Double & density, Double & compress, const Integer matnum, const std::string keyword)
  { Error("Not implemented",__FILE__,__LINE__);}


  //! get dielectric term for an electrostatic equation from the material file
  /*!
	\param dielectr out: value of dielectric term from the material file
	\param matnum in: material number
  */
  // for compatibility with elena
  void ReadDielectricTerms(Double & dielectr,const Integer matnum)
  { Error("Not implemented",__FILE__,__LINE__);}

  // get value for error functional (see Institusbericht Michael Schinnerl: Behandlung gekoppelter
  // Systeme mit Mehrgitterverfahren)
  // Double GetErrorFunctional(const Double& MagB) const;  
};

} // end of namespace
#endif
