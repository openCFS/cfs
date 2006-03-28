#ifndef MATERIAL_DATA
#define MATERIAL_DATA

/**************************************************************************/
/* File:   MaterialData.hh                                                */
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

namespace CoupledField {

  //! Class for Material Data
  /*! 
    Interface class for getting material data
  */

  class MaterialData {

  private:

    // material type: fluid
    Double density;
    Double compressibility;
    Double damp_alfa;
    Double damp_beta;
    Double BoverA;
    // material type: thermic
    Double heatCapacity;
    Double thermalConductivity;

    Double eModule;
    Double nu;
    Double LameLambda;
    Double LameMu;
    Double permMx, permMy, permMz;   // permanent magnetization

    Double Esat;  //electric field intensity value for satuaration
    Double Psat;  //electric polarization value for satuaration
    Double aJiles_;
    Double alphaJiles_;
    Double kJiles_;
    Double cJiles_;
    UInt dirPol; //direction of polarization
    std::string hystType_;

    Boolean scaledMatDat;
    char * name; 
    //char name[stringLength]; 
    std::string bhCurveFile_; //!< name of BH-Curve datafile

    //const ARRAY <NonlinSpline*> & magneticSpline;
    //  ARRAY <NonlinSpline*> * magneticSpline;

    /// contains the stiffnes matrix, the piezelectric coefficients and the permitivity matrix
    Matrix<Double> * piezoMatrix_;
    Matrix<Double> * piezoMatrixC_;

    Matrix<Double> * permeaMatrix;
    Matrix<Double> * conducMatrix;

    Matrix<Complex> * complexPiezoMatrix_;
    Matrix<Double>  stiffnessMatrix_;
    Matrix<Complex> complexStiffnessMatrix_;

    UInt matNr;
    Integer nonlin;
  
  public:

    //! Default constructor
    MaterialData();

    //! Copy constructor
    MaterialData( const MaterialData &mat );

    ~MaterialData();

    void GetStiffnessMatrix(Matrix<Double> &stiffMat){
      stiffMat =  stiffnessMatrix_;};
    void GetStiffnessMatrix(Matrix<Complex> & stiffMat){
      stiffMat = complexStiffnessMatrix_;};

  
    /// set the material number 
    void SetMatNr(const UInt& MatNr){matNr = MatNr; };

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
    void SetBoverA( const Double &BA )
    {BoverA = BA;}

    /// set parameters for heat conduction
    void SetThermic(const Double& HeatCapacity, const Double ThermalConductivity)
    {heatCapacity = HeatCapacity; thermalConductivity = ThermalConductivity;}

    /// set electric field value for saturation
    void SetEsat(Double& val) {Esat = val;};

    /// set electric polarization value for saturation
    void SetPsat(Double& val) {Psat = val;};

    /// set direction of electric polarization
    void SetDirPol(UInt& val) {dirPol = val;};

    /// set type of hysteresis
    void SetHysteresisType(std::string& atype) {hystType_ = atype;};

    /// 
    void SetJiles_a(Double& val) {aJiles_ = val;};

    /// 
    void SetJiles_alpha(Double& val) {alphaJiles_ = val;};

    /// 
    void SetJiles_k(Double& val) {kJiles_ = val;};

    /// 
    void SetJiles_c(Double& val) {cJiles_ = val;};

    /// set values of permanent magnetization
    void SetPermMag(const Double& mX, const Double& mY, const Double& mZ)
    {permMx=mX, permMy=mY, permMz=mZ;};

    /// get values of permanent magnetization
    void GetPermMag( Double& mX,  Double& mY, Double& mZ)
    {mX=permMx, mY=permMy, mZ=permMz;};

    /// get nonlinearity datafile name for BH-curve for magnetic material
    std::string& GetBHCurveFileName() { return bhCurveFile_; }

    /// set one value of the data-matrix on position (i,j)
    void SetPiezoMatrixData(const UInt& i, const UInt& j, const Double& value)
    {(*piezoMatrix_)(i,j) = value;
      (*complexPiezoMatrix_)[i][j]=Complex(value,(*complexPiezoMatrix_)[i][j].imag());
    };

    /// set one value of the data-matrix on position (i,j)
    void SetPiezoMatrixDataC(const UInt& i, const UInt& j, const Double& value)
    {(*piezoMatrixC_)(i,j) = value;
      (*complexPiezoMatrix_)[i][j]=Complex((*complexPiezoMatrix_)[i][j].real(),value);
    };

    /// get the value of the data-matrix on position (i,j)
    void GetPiezoMatrixData(const UInt& i, const UInt& j, Double& value)
    {value = (*piezoMatrix_)(i,j);};

    /// get the value of the data-matrix on position (i,j)
    void GetPiezoMatrixDataC(const UInt& i, const UInt& j, Double& value)
    {value = (*piezoMatrixC_)(i,j);};

    /// returns a pointer to the imaginary part of data-matrix
    Matrix<Double> * GetMatrixC(){return piezoMatrixC_;};

    /// returns a pointer to the real part of data-matrix
    Matrix<Double> * GetMatrix(){return piezoMatrix_;};

    /// returns pointer to complex material matrix
    Matrix<Complex> * GetComplexMaterialMatrix(){
      return complexPiezoMatrix_;
    };

    /// Rotates piezo Material Matrix. 
    /// Input are the three solid angels (radian measure)

    //! Rotates the piezoelectric material matrix. 
    //! Input are the three solid angels (radian measure).
    //! Special case: The choice of 1 for one of the angels,
    //!  rotates piezo matrix in that way
    //! that it is polarized in the given direction. 
    //! The choice of 0 does not performs any rotation
    //! Poling in z - direction (a3=1) is given by default.
    void RotateMaterialMatrix(const Double& a1, const Double& a2, const Double& a3);

    //! set one value of the permeability-matrix on position (i,j)
    void SetPermeability(const UInt& i, const UInt& j, const Double& value);

    //! get the value of the permeability-matrix on position (i,j)
    void GetPermeability(const UInt& i, const UInt& j, Double &value);

    //! return a pointer to the permeability matrix
    Matrix<Double> * GetPermeaMatrix(){return permeaMatrix;};

    //! set one value of the permeability-matrix on position (i,j)
    void SetConductivity(const UInt& i, const UInt& j, const Double& value);

    //! get the value of the permeability-matrix on position (i,j)
    void GetConductivity(const UInt& i, const UInt& j, Double &value);

    //! return a pointer to the permeability matrix
    Matrix<Double> * GetConducMatrix(){return conducMatrix;};

    /// set size of data-matrix in x and y direction to nrElems3d x nrElems3d
    /// this matrix includes the stiffness, piezoelectric coupling and permitivity matrix
    void DefFull3dMatrix(){
      piezoMatrix_ = new Matrix<Double>; 
      piezoMatrix_->Resize(GetNrElems3d(), GetNrElems3d() );
      piezoMatrixC_ = new Matrix<Double>; 
      piezoMatrixC_->Resize(GetNrElems3d(), GetNrElems3d() );
      complexPiezoMatrix_ = new Matrix<Complex>;
      complexPiezoMatrix_->Resize(GetNrElems3d(), GetNrElems3d() );};

    /// set conductivity of the material
    void SetConductivity(const Double& Conductivity);

    /// set name of the material
    void SetName(const char* Name);

    /// 
    void SetScaledFlag(){scaledMatDat = 1;};

    Boolean IsMatDatScaled(){return scaledMatDat;};

    // is material nonlinear?
    //  UInt GetNonlin() const { return nonlin; };

    /// get nr of elems of the matrix in the 3d case
    UInt GetNrElems3d() const { return 9; };

    /// get nr of elems of the matrix in the 2d case
    UInt GetNrElems2d() const { return 5; };

    /// get number of the matNr
    UInt GetMatNr() const { return matNr; };

    /// get E-Module
    Double GetEModule() const {return eModule; };

    /// get density
    Double GetDensity() const {return density;};

    /// get LameLambda
    Double GetLameLambda() const {return LameLambda;}
  
    /// get LameMu
    Double GetLameMu() const {return LameMu;}

    /// get permittivity
    Double GetPermittivity(UInt i, UInt j) const {return (*piezoMatrix_)[i+6][j+6];};

    /// get compressibility
    Double GetCompressibility() const {return compressibility;};

    /// get permeability of the material
    //  Double GetPermeability() const {return permeability;};

    /// get alfa damping coefficient
    Double GetDampingAlfa() const {return damp_alfa;};

    /// get beta damping coefficient
    Double GetDampingBeta() const {return damp_beta;};  

    /// get BoverA for nonlinear acoustics
    Double GetBoverA() const {return BoverA;}

    /// get heatConduction
    Double GetHeatCapacity() const {return heatCapacity;}

    /// get thermalConductivity
    Double GetThermalConductivity() const {return thermalConductivity;}

    /// get nu
    Double GetNu() const { return nu; };

    /// get saturated electric field value
    Double GetEsat() const {return Esat;};

    /// get saturated electric polarization value
    Double GetPsat() const {return Psat;};

    /// get direction of polarization
    UInt GetDirPol() const {return dirPol;};

    /// get hysteresis type
    std::string GetHysteresisType() const {return hystType_;};

    /// 
    Double GetJiles_a() const {return aJiles_;};

    /// 
    Double GetJiles_alpha() const {return alphaJiles_;};

    /// 
    Double GetJiles_k() const {return kJiles_;};

    /// 
    Double GetJiles_c() const {return cJiles_;};

    /// get conductivity
    //  Double GetConductivity() const {return conductivity; };

    /// get name of the material
    char * GetMaterialName(); 

    //! Check whether material has imaginary part
    bool IsComplex() {
      std::string matName( name );
      UInt tail = matName.size() - 5;
      bool retVal = true;
      if ( matName.size() <= 5 || matName.substr(tail) != "-imag" ) {
        retVal = false;
      }
      return retVal;
    }

    //! get density and compressity for acoustic equation from material data-file
    /*!
      \param density out: value of density from the material file
      \param compress out: value of compressity from the material file
      \param matnum material number
      \param keyword name of material in the  material file
    */
    // for compatibility with elena
    //   void ReadDensityAndCompressity(Double & density, Double & compress, const UInt matnum, const std::string keyword)
    //   { Error("Not implemented",__FILE__,__LINE__);}


    //! get dielectric term for an electrostatic equation from the material file
    /*!
      \param dielectr out: value of dielectric term from the material file
      \param matnum in: material number
    */
    // for compatibility with elena
    //   void ReadDielectricTerms(Double & dielectr,const UInt matnum)
    //   { Error("Not implemented",__FILE__,__LINE__);}

    // get value for error functional (see Institusbericht Michael Schinnerl: Behandlung gekoppelter
    // Systeme mit Mehrgitterverfahren)
    // Double GetErrorFunctional(const Double& MagB) const;  
  };

} // end of namespace
#endif
