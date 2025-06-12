#ifndef BASEMATERIAL_DATA
#define BASEMATERIAL_DATA

#include <map>
#include <set>

#include "General/Environment.hh"

#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"

#include "Utils/ApproxData.hh"

#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField {

  // forward class declarations
  class MathParser;
  class CoordSystem;
  class Hysteresis;
  class ElemList;
  //class PiezoMicroModelHF;
  class PiezoMicroModelBK;

  //! Class for Material Data
  /*! 
    Base class for handling material data
  */

  class BaseMaterial {

  public:
    
    //! Helper struct, defining a single nonlinear dependency
    struct MatDescriptorNl {
      
      //! Constructor
      MatDescriptorNl();
      
      //! Destructor
      ~MatDescriptorNl();
      
      //! Type of approximation (smooth spline etc.)
      ApproxCurveType approxType;
      
      //! File containing the nonlinearity
      std::string fileName;
      
      //! Accuracy of measured values (needed for approximation)
      Double measAccuracy;
      
      //! Maximum value of approximation (limit value)
      Double maxVal;
      
      //! Angle, for which this curve is defined (0 for no angular dependency)
      Double angle;

      //! Temperature, for which this curve is defined
      Double temperature;
      
      //! Scaling factor of anisotropic behavior in z-direction 
      //! Explanation: we assume the same material behavior in z-direction as in xy-plane 
      //! therefore the same curves with an optional scaling factor are used for calculation
      Double zScaling;

      //! Pointer to interpolation class
      ApproxData* approxData;

      //! Scale factor material data read from file
      Double factor;

      //! string for analytic expression
      std::string analyticExpr;
      
      //! string for derivative
      std::string analyticExprDeriv;

      //! string for derivative w.r.t parameter 1
      std::string analyticExprDerivP1;

      //! string for derivative w.r.t parameter 2
      std::string analyticExprDerivP2;

      //! string for derivative w.r.t parameter 3
      std::string analyticExprDerivP3;            

      //! string for derivative w.r.t parameter 4
      std::string analyticExprDerivP4;      
    };
    
    //@{ \name public typedefs
    typedef std::map<MaterialType, std::string> StringMap;
    typedef std::map<MaterialType, Integer> IntegerMap;
    typedef std::map<MaterialType, PtrCoefFct> CoefMap;
    typedef std::map<MaterialType, MatDescriptorNl> NonLinIsoMap;
    typedef std::map<MaterialType, StdVector<MatDescriptorNl> > NonLinAnisoMap;
    typedef std::map<MaterialType, StdVector<MatDescriptorNl> > NonLinIsoMapTempDependBHcurves;

    //! Denote the symmetry type 
    typedef enum {
      NOSYMMETRY, GENERAL, ISOTROPIC, ORTHOTROPIC, TRANS_ISOTROPIC
    } SymmetryType;
    static Enum<SymmetryType> SymmetryTypeEnum;

    //! Denote the material model
    typedef enum {
      MAT_MODEL_LINEAR, MAT_MODEL_LINEARVISCOELASTIC
    } MaterialModel;
    static Enum<MaterialModel> MaterialModelEnum;

    //@}
    

    //! Default constructor
    
    //! Default constructor
    //! \param mp Pointer to MathParser
    //! \param defaultCoosy Pointer to default coordinate system
    BaseMaterial( MaterialClass matClass, MathParser* mp, CoordSystem * defaultCoosy);

    //! Destructor
    virtual ~BaseMaterial();


    //! Trigger finalization of material
    virtual void Finalize() {};

    /** Print the material data which is actually read and stored in isSet */
    void ToInfo(PtrParamNode in, SubTensorType stt = NO_TENSOR,
                const Vector<double>* rot = NULL);

    //! set the name of the material set
    void SetName(const std::string &name);

    //! returns the name of the material
    std::string GetName() const;

    //! returns the class of the material
    MaterialClass GetClass() const;

    // ======================================================================
    //  Coefficient Function Related Methods
    // ======================================================================
    //@{ \name Coefficient Function Related Method

    //! Return a tensor valued coefficient function (linear)
    virtual PtrCoefFct GetTensorCoefFnc(MaterialType matType,
                                        SubTensorType type,
                                        Global::ComplexPart matDataType,
                                        bool transposed = false ) const;

    //! Return vector-valued coefficient function (linear)
    virtual PtrCoefFct GetVectorCoefFnc(MaterialType matType,
                                        Global::ComplexPart matDataType) const;

    //! Return scalar-valued coefficient function (linear)
    virtual PtrCoefFct GetScalCoefFnc(MaterialType matType,
                                      Global::ComplexPart matDataType) const;

    //! Return a specific sub-tensor as coefficient function (linear)
    virtual PtrCoefFct GetSubTensorCoefFnc( MaterialType matType, 
                                            SubTensorType tensorType,
                                            Global::ComplexPart matDataType,
                                            bool transposed = false) const;

    //! Return a sub-vector in Voigt notation
    virtual PtrCoefFct GetSubVectorCoefFnc( MaterialType matType,
                                            SubTensorType tensorType,
                                            Global::ComplexPart matDataType
                                                = Global::COMPLEX) const;

    //! Return tensor-valued coefficient function for nonlinear function
    virtual PtrCoefFct GetTensorCoefFncNonLin( MaterialType matType,
                                               SubTensorType type,
                                               Global::ComplexPart matDataType,
                                               PtrCoefFct dependency );

    //! Return scalar-valued coefficient function for derivative w.r.t. parameter 
    virtual PtrCoefFct GetScalCoefFncNonLinDerivParam(MaterialType matType,
                                                      Global::ComplexPart matDataType,
                                                      PtrCoefFct fluxCoef) {
      EXCEPTION("GetScalCoefFncNonLinDerivParam not implemented");
    }                                                                                                                                                        
    
    //! Return scalar-valued coefficient function for a matrial model
    //virtual PtrCoefFct GetScalCoefFncModel(shared_ptr<CoefFunction> matModel);

    //! Return tensor-valued coefficient function for a matrial model (e.g., vector hysteresis)
    //virtual PtrCoefFct GetTensorCoefFncModel(shared_ptr<CoefFunction> matModel);

    //! Return vector-valued coefficient function for a matrial model (e.g., vector hysteresis)
    //virtual PtrCoefFct GetVectorCoefFncModel(shared_ptr<CoefFunction> matModel);

    //! Return scalar-valued coefficient function for nonlinear function
    virtual PtrCoefFct GetScalCoefFncNonLin(MaterialType matType,
                                            Global::ComplexPart matDataType,
                                            PtrCoefFct dependency,
                                            PtrCoefFct temperature_dependency = NULL );
                                            
    //! Return scalar-valued coefficient function for nonlinear function (only for magstrict nu)
    virtual PtrCoefFct GetScalCoefFncNonLin_MagStrict(MaterialType matType,
                                            Global::ComplexPart matDataType,
                                            PtrCoefFct dependency );
					    
    //! Return scalar-valued coefficient function for nonlinear function
    //! where the value is calculated depending on the value of
    //! \param temperatureCoef and \param elecPotCoef on \param regs.
    virtual PtrCoefFct GetScalCoefFncMultivariateNonLin(
        MaterialType matType,
        NonLinType nlType,
        Global::ComplexPart matDataType,
        StdVector<PtrCoefFct> dependencies,
        StdVector<RegionIdType> & regs);

    //! Set a coefficient function
    virtual void SetCoefFct( MaterialType matType, PtrCoefFct coef );

    //@}
    
    //! Query if a given parameter is set
    virtual bool IsSet( MaterialType matType ) const;

    //! set the symmetry type
    void SetSymmetryType(MaterialType matType, SymmetryType symType);

    //! get the symmetry type
    SymmetryType GetSymmetryType(MaterialType matType) const;

    //! Returns the material model selected for this material
    MaterialModel GetModel() const;

    //! Set the material model
    void SetModel(MaterialModel model);

    //! set a scalar string material parameter
    virtual void SetString(const std::string& param, MaterialType matType);

    //! set a scalar integer material parameter
    void SetScalar(Integer param, MaterialType matType);

    //! Set a scalar real constant material parameter
    virtual void SetScalar(Double param, MaterialType matType,
                           Global::ComplexPart dataType );

    //! set a scalar complex constant material parameter
    virtual void SetScalar(Complex param, MaterialType matType,
                           Global::ComplexPart dataType );
    
    //! set a real constant material vector
    virtual void SetVector(const Vector<Double>& param, MaterialType matType,
                           Global::ComplexPart dataType);

    //! set a complex constant material vector
    virtual void SetVector(const Vector<Complex>& param, MaterialType matType,
                           Global::ComplexPart dataType);

    //! set a real constant material tensor
    virtual void SetTensor(const Matrix<Double>& param, MaterialType matType,
                           Global::ComplexPart dataType);

    //! set a complex constant material tensor
    virtual void SetTensor(const Matrix<Complex>& param, MaterialType matType,
                           Global::ComplexPart dataType );
    
    //! Set a nonlinear isotropic approximation
    virtual void SetNonLinMatIso( MaterialType matType, MatDescriptorNl& data );
    
    //! Set a nonlinear anisotropic approximation
    virtual void SetNonLinMatAniso( MaterialType matType, StdVector<MatDescriptorNl>& data );
    
    //! Set a nonlinear isotropic approximation for temperature dependent BH-curves
    virtual void SetNonLinMatIsoTempDependBH( MaterialType matType, StdVector<MatDescriptorNl>& data );
    
    //! get a string material parameter
    virtual void GetString( std::string& param, MaterialType matType) const;
    
    //! Return a constant scalar-valued coefficient function.

    //! This method returns a constant, scalar-valued material data. If
    //! the data is not constant (e.g. time/freq/spatial dependency),
    //! an exception is thrown.
    virtual void GetScalar( Double& param, MaterialType matType,
                            Global::ComplexPart matDataType = Global::REAL) const;
    virtual void GetScalar( Complex& param, MaterialType matType,
                            Global::ComplexPart matDataType = Global::COMPLEX) const;

    //! get a integer material parameter
    virtual void GetScalar( Integer& param, MaterialType matType) const;

    //! Get a constant vector
    virtual void GetVector( Vector<Double>& param, MaterialType matType,
                            Global::ComplexPart matDataType = Global::REAL ) const;
    virtual void GetVector( Vector<Complex>& param, MaterialType matType,
                            Global::ComplexPart matDataType = Global::COMPLEX ) const;

    //! get a real material tensor
    virtual void GetTensor( Matrix<Double>& param, MaterialType matType, 
                            Global::ComplexPart dataType = Global::REAL,
                            SubTensorType subTensor = FULL ) const;
     
    //! get a complex material tensor
    virtual void GetTensor( Matrix<Complex>& param, MaterialType matType,
                            Global::ComplexPart dataType = Global::COMPLEX,
                            SubTensorType subTensor = FULL ) const;
    
    //! Rotate a material tensor by coordinates
    virtual void RotateTensorByRotationAngles( const Vector<Double> &coord, 
                                               MaterialType matType,
                                               bool persistent = false );

    //! Rotate all tensor material parameters by given rotation angle
    virtual void RotateAllTensorsByRotationAngles(  const Vector<Double>& rotAngle, 
                                                    bool persistent = false );

    //! Rotate a material tensor by 3x3 rotation matrix
    virtual void RotateTensorByRotationMatrix( const Matrix<Double>& rotMat, 
                                               MaterialType matType,
                                               bool persistent = false );

    //! Rotate all tensor material parameters by rotation matrix
    virtual void RotateAllTensorsByRotationMat(  const Matrix<Double>& rotMatrix,
                                                 bool persistent = false );

    //! Rotate all tensor material parameters by given coordinate system
    //! and axis mapping.
    virtual void RotateAllTensorsByCoordSys(CoordSystem *coordSys,
                                            const StdVector<std::string> &axisMap,
                                            const StdVector<Double> &axisFactors,
                                            bool persistent = false);

    //! Rotates the tensor in a way that is represents the attached
    //! coordinate system behaviour (cartesian, cylindrical, spherical)
    //! in this point
    virtual void RotateTensorByPointCoord( const Vector<Double> &coord, 
                                           MaterialType matType );

    //! Pass coordinate system to material
    void SetCoordSys( CoordSystem* system );

    //! Get coordinate system from material
    CoordSystem* GetCoordSys();

    // ======================================================================
    //  Hysteresis Related Information
    // ======================================================================
    //@{ \name Hysteresis Related Information

    //Initialize hysteresis
    //virtual void InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
    //                       bool isInverse = false, bool computeInverse = false );

    //Initialize hysteresis
    // calls either Preisach or VectorPreisach depending on the dimensions
    virtual void InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                           bool isInverse = false, bool computeInverse = false, UInt dim = 1);

    //Initialize hysteresis
    virtual void InitVecHyst( UInt numElemSD, shared_ptr<ElemList> actSDList, 
                              UInt dim );

    //! get hysteresis operator
    Hysteresis* getHysteresis() {
      return hyst_;
    };

    //! get hysteresis operator
    Hysteresis** getVecHysteresis() {
      return vecHyst_;
    };

    bool IsSetHysteresis () {
      return isHysteresis_;
    };

    //set values for differential material approach: scalar hyst
    virtual void SetPreviousHystVal( UInt nrElem, Double Xval ) {
      EXCEPTION( "SetPreviousHystVal not implemented" );
    };

    //set values for differential material approach: vector-hyst
    virtual void SetPreviousHystVal( UInt nrElem, Vector<Double>& Xval ) {     
      EXCEPTION( "SetPreviousHystVal not implemented" );
    };

    //! computes the scalar hysteresis value
    virtual Double ComputeScalarHystVal( UInt nrElem, Double Xval ) {
      EXCEPTION( "ComputeScalarHystVal not implemented" );
      return 1.0;
    };

    //! compute scalar differential parameter
    virtual Double ComputeScalarDiffVal( UInt nrElem, Double Xval ) {
      EXCEPTION( "ComputeScalarDiffValue not implemented" );
      return 1.0;
    };

    //! compute the vector hysteresis values
    virtual void ComputeVectorHystVal( UInt nrElem, Vector<Double>& Xin, 
                                       Vector<Double>& Yout ) {
      EXCEPTION( "ComputeVectorHystVal not implemented" );
    };

    //! get scalar hysteresis value
    virtual Double GetScalarHystVal( UInt nrElem );

    //! get previous scalar hysteresis value
    virtual Double GetScalarHystPrevVal( UInt nrElem );

    //! get vector of hysteresis value
    virtual void GetVectorHystVal( UInt nrElem, Vector<Double>& Val ) {
      EXCEPTION( "ComputeVectorHystVal not implemented" );
    };

    //! Set an anhysteretic material model
    virtual void SetAnhystMagModel( const std::string name );
    std::string GetAnhystMagModel(){return anhystereticModel_;};

    //@}

    // ======================================================================
    //  Micro-Piezoelectric Model
    // ======================================================================
    //@{ \name Micro-Piezoelectric Model

    //Initialize piezoelectric-micro-model
    virtual void InitPiezoMicro( UInt numElemSD, shared_ptr<ElemList> actSDList, 
                                 BaseMaterial* mechMat, BaseMaterial* elecMat,
                                 SubTensorType tensorType, Double dt);

    //! returns the material tensors
    void GetEffectiveTensors( Matrix<Double>& matMechC,
                              Matrix<Double>& matMechS,
                              Matrix<Double>& matElec,
                              Matrix<Double>& matPiezo,
                              Vector<Double>& stress, 
                              Vector<Double>& elecField,
                              UInt elemIdx, 
                              bool recompute,
                              bool previous );

    //! returns ireversible strain and polarization
    void GetEffectiveIrreversibleValues( Vector<Double>& Pirr,
                                         Vector<Double>& Sirr,
                                         UInt elemIdx,
                                         bool recompute,
                                         bool previous );
    
    //!
    void ComputeEffectiveCouplingTensor(Matrix<Double>& dmatEff, 
                                        Vector<Double>& elecFieldAct,
                                        Vector<Double>& elecFieldPrev,
                                        UInt elemIdx);

    //! get micro-piezo-object
    PiezoMicroModelBK* GetMicroPiezoModel() {
      return piezoMicroModel_;
    };

    //@}
    //========================== micro-piezoelectric-model:end ===================

    //! Set the mathParser strings for the Rayleigh damping coefficients alpha and beta
    //! Takes care of the different cases of computing the Rayleigh damping coefficients depending on
    //! the material parameters set in the material xml
    //! \param alpha String for the Rayleigh damping coefficient alpha
    //! \param beta String for the Rayleigh damping coefficient beta
    void GetRayleighCoeffStrings(std::string &alpha, std::string &beta);

    //! Create and return the mathParser strings for the Rayleigh damping coefficients alpha and beta,
    //! when the damping is specified in terms of a frequency-adapted or frequency-constant lossTangensDelta.
    //! \param alpha String for the Rayleigh damping coefficient alpha
    //! \param beta String for the Rayleigh damping coefficient beta
    void GetFreqAdaptedRayleighCoeffStrings(std::string &alpha, std::string &beta);

    /** converts MaterialClass to the corresponding MaterialType tensor. Extend for your needs */
    static MaterialType ConvertMaterialClass(MaterialClass mc);

    //! Rotate 2nd order tensor in Voigt notation
    template<typename T>
    static void PerformRotationVoigt( const Matrix<Double>& rotMatrix,
                                      Vector<T>& rotMatTensor,
                                      const Vector<T>& origMatTensor );

  protected:

    //! Error for material type not defined
    void matTypeNotAllowed(MaterialType matType, const std::string& dim ) const;

    //! data type not allowed in set/get-function
    void dataTypeNotAllowed4SetGet(Global::ComplexPart datType, const std::string& msg ) const;

    //! Error for data type not allowed
    void dataTypeNotAllowed(Global::ComplexPart datType, MaterialType matType ) const;

    //! Error for material type not in file
    void matTypeNotInDataBase(MaterialType matType, const std::string& dim ) const;

    //! Error for data type not allowed in set-function
    void setMakesNoSense(Global::ComplexPart datType, const std::string& msg ) const;

    //! Error for not available subtype of tensor
    static void subTensorNotAvailable(MaterialType matType, SubTensorType subTensor);

    //! Error for material data is not constant, but was queried via Get{Scalar,Vector,Tensor}
    void matDataNotConstant(MaterialType matType) const;

    //! Rotate tensorial coefficient function
    virtual void PerformRotation( const Matrix<Double>& rotMatrix, PtrCoefFct& rotatedCoef,
                                  PtrCoefFct origCoef);
    //! Rotate tensorial coefficient function (only for constant CoefFct)
    virtual void PerformRotationConst( const Matrix<Double>& rotMatrix, PtrCoefFct& rotatedCoef,
                                       PtrCoefFct origCoef);

    //! Compute a full 3x3 tensor from different representations
    virtual void CalcFull3x3Tensor( MaterialType isoPProp,
                                    MaterialType* orthoProp,
                                    MaterialType tensorProp );

    /** helper for ToInfo(). If the  imaginary part is zero, only the real part is printed  */
    void StoreTensor(PtrParamNode in, PtrCoefFct tensorFunc);


    //! name of material
    std::string name_;

    //! class of material database
    const MaterialClass class_;

    //! set, which knows about the allowed material parameters for a material class
    std::set<MaterialType> isAllowed_;

    //! set, which knows, which material parameters have been set
    std::set<MaterialType> isSet_;

    //! map, which knows about the material data defined by strings
    StringMap stringParams_;

    //! map, which knows about the material data defined by integers
    IntegerMap integerParams_;

    //! map storing the isotropic nonlinear material parameters
    NonLinIsoMap nonlinIsoParams_;
    
    //! map storing the anisotropic nonlinear material parameters
    NonLinAnisoMap nonlinAnisoParams_;

    //! map storing the temperature-dependent nonlinear material parameters
    NonLinIsoMapTempDependBHcurves nonlinIsoTempDependBHParams_;

    //! name of anhysteretic model version
    std::string anhystereticModel_;

    //! name of anhysteretic model version
    std::string hessiantype_;

    // ========================================================
    //  New coefficient based material representation
    // ========================================================
    //! Tensor coefficients
    CoefMap tensorCoef_;

    //! Original / unrotated tensor coefficients
    CoefMap tensorOrigCoef_;
    
    //! Vector coefficients
    CoefMap vectorCoef_;
    
    //! Scalar coefficients
    CoefMap scalarCoef_;
    
    //! Pointer to math parser instance
    MathParser*  mp_;
    
    //! Pointer to attaches coordinate system
    CoordSystem* coosy_;

    //! Store symmetryType per material data type
    std::map<MaterialType, SymmetryType> symmetryType_;

    //! Material model used for this material
    MaterialModel model_;

    //! hysteresis object
    Hysteresis* hyst_;

    //! hysteresis object
    Hysteresis* hystY_;

    //! hysteresis object
    Hysteresis* hystZ_;

    //! hysteresis object
    Hysteresis** vecHyst_;

    //!
    UInt dimVecHyst_;

    //!
    bool isHystInverse_;

    //!
    bool computeHystInverse_;

    //! hysteresis set 
    bool isHysteresis_;

    Vector<Double> Xprevious_; //! previous Xval in hysteresis
    Vector<Double> Yprevious_; //! previous Yval in hysteresis

    //! for vector version
    //! note that we do not use the matrices below due to different sorting at the moment
    Vector<Double>* XpreviousVEC_;
    Vector<Double>* YpreviousVEC_;

    // dimension of hystersis: 1 = Preisach, 2,3 = VectorPreisach
    // independent from SimplePreisach which uses dimVecHyst_
    UInt dim_;

    Matrix<Double> vecXprevious_; //! previous Xval in hysteresis
    Matrix<Double> vecYprevious_; //! previous Yval in hysteresis

    Matrix<Double> actDiffVal4VecHyst_;
    Matrix<Double> previousDiffVal4VecHyst_;

    std::map<UInt, UInt> globalElem2Local_;

    //! yes, piezoelectric micro-model is switched on
    bool isPiezoMicroModel_;	

    //! object for piezo-micro-modeling
    PiezoMicroModelBK* piezoMicroModel_;

  public:
    typedef enum {ALPHA_BETA, ADAPTED_LOSS_TANGENS} RayleighDampType;
    void SetRayleighType(RayleighDampType raylDampType) {raylDampType_ = raylDampType;};

  private:
    RayleighDampType raylDampType_;
  };

} // end of namespace

#endif
