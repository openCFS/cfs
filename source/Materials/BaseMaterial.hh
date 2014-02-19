#ifndef BASEMATERIAL_DATA
#define BASEMATERIAL_DATA

#include <map>
#include <set>

#include "General/Environment.hh"

#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"

#include "Utils/ApproxData.hh"
#include "Utils/mathParser/mathParser.hh"

#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField {

  // forward class declarations
  class CoordSystem;
  class Hysteresis;
  class ElemList;
  class PiezoMicroModelHF;
  class PiezoMicroModelBK;
  class CoefFunction;
  class BaseBOperator;
  class BaseFeFunction;

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

      //! Pointer to interpolation class
      ApproxData* approxData;
      
    };
    
    //@{ \name public typedefs
    typedef std::map<MaterialType, Matrix<Complex> > tensorMap;
    typedef std::map<MaterialType, Vector<Complex> > vectorMap;
    typedef std::map<MaterialType, Complex > scalarMap;
    typedef std::map<MaterialType, std::string > stringMap;
    typedef std::map<MaterialType, Integer > integerMap;
    typedef std::map<MaterialType, MathParser::HandleType > handleMap;
    typedef std::map<MaterialType, PtrCoefFct> CoefMap;
    typedef std::map<MaterialType, MatDescriptorNl> nonLinIsoMap;
    typedef std::map<MaterialType, StdVector<MatDescriptorNl> > nonLinAnisoMap;

    //! Denote the symmetry type 
    typedef enum {GENERAL, ISOTROPIC, ORTHOTROPIC, TRANS_OTHOTROP } SymmetryType;
    
    //@}
    

    //! Default constructor
    
    //! Default constructor
    //! \param mp Pointer to MathParser
    //! \param defaultCoosy Pointer to default coordinate system
    BaseMaterial( MathParser* mp,
                  CoordSystem * defaultCoosy);

    //! Destructor
    virtual ~BaseMaterial();


    //! Trigger finalization of material (calculation of rotated matrices)
    virtual void Finalize() {};

    /** Print the material data which is actually read and stored in isSet */
    void ToInfo(PtrParamNode in);

    //! set the name of the material set
    void SetName(const char* name) {
      matFileName_.assign( name ) ;
    }

    // returns the name of the material file
    std::string& GetName() {
      return matFileName_;
    }

    //! returns the name of the material database
    std::string& GetMaterialDatabaseName() {
      return materialDatabaseName_;
    }

    // ======================================================================
    //  Coefficient Function Related Methods
    // ======================================================================
    //@{ \name Coefficient Function Related Method

    //! Return a tensor valued coefficient function (linear)
    virtual PtrCoefFct GetTensorCoefFnc(MaterialType matType,
                                        SubTensorType type,
                                        Global::ComplexPart matDataType,
                                        bool transposed = false );

    //! Return vector-valued coefficient function (linear)
    virtual PtrCoefFct GetVectorCoefFnc(MaterialType matType,
                                        Global::ComplexPart matDataType);

    //! Return scalar-valued coefficient function (linear)
    virtual PtrCoefFct GetScalCoefFnc(MaterialType matType,
                                      Global::ComplexPart matDataType);

    //! Return a specific sub-tensor as coefficient function (linear)
    virtual PtrCoefFct GetSubTensorCoefFnc( MaterialType matType, 
                                            SubTensorType tensorType,
                                            bool transposed  ) ;

    //! Return tensor-valued coefficient function for nonlinear function
    virtual PtrCoefFct GetTensorCoefFncNonLin( MaterialType matType,
                                               SubTensorType type,
                                               Global::ComplexPart matDataType,
                                               PtrCoefFct dependency );
    
    //! Return scalar-valued coefficient function for nonlinear function
    virtual PtrCoefFct GetScalCoefFncNonLin(MaterialType matType,
                                            Global::ComplexPart matDataType,
                                            PtrCoefFct dependency );
					    
    //! Return scalar-valued coefficient function for nonlinear function
    //! where the value is calculated depending on the value of \param temperatureCoef and \param elecPotCoef on \param regs. 
    virtual PtrCoefFct GetScalCoefFncMultivariateNonLin(MaterialType matType,
                                                   NonLinType nlType,
                                                   Global::ComplexPart matDataType,
                                                   StdVector<PtrCoefFct> dependencies,
 					           StdVector<RegionIdType> & regs) 
						   {EXCEPTION("not implemented in base class"); PtrCoefFct dummy; return dummy;}

    //@}
    
    //! get info, which material parameter is set
    std::set<MaterialType> GetIsSetInfo() const
    { return isSet_;};

    //! get info, which material parameter is complex
    std::set<MaterialType> GetIsComplexInfo() const
    { return isComplex_;};

    //! get the tensors
    tensorMap GetTensorParams() const
    { return tensorParams_;};

    //! get the scalar material parameters
    scalarMap GetScalarParams() const
    { return scalarParams_;};

    //! get the material parameters defined by a string
    stringMap GetStringParams() const
    { return stringParams_;};

    //! get the integer material parameters 
    integerMap GetIntegerParams() const
    { return integerParams_;};

    //! Query if a given parameter is set
    virtual bool IsSet( MaterialType matType ) const;

    //! set the symmetry type
    void SetSymmetryType(MaterialType matType, SymmetryType symType) {
      symmetryType_[matType]=symType; 
    };

    //! get the symmetry type
    SymmetryType GetSymmetryType(MaterialType matType) const {
      SymmetryType ret = GENERAL;
      std::map<MaterialType, SymmetryType>::const_iterator it;
      it = symmetryType_.find(matType);
      if( it == symmetryType_.end() ) {
        EXCEPTION("Could not find symmetry type for material '"
            << MaterialTypeEnum.ToString(matType) << "'");
      } else {
        ret = it->second;
      }
         
      return ret; 
    };

   //! set a scalar string material parameter
    virtual void SetScalar(const std::string& param, MaterialType matType);

    //! set a scalar integer material parameter
    virtual void SetScalar(int param, MaterialType matType);

    //! set a scalar real material parameter
    virtual void SetScalar(double param, MaterialType matType, Global::ComplexPart dataType ) 
    {
     EXCEPTION("not implemented for " << materialDatabaseName_);
    };


    //! set a scalar complex material parameter
    virtual void SetScalar(Complex param, MaterialType matType, Global::ComplexPart dataType )
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }
    
    //! set a scalar parameter in string representation (gets later parsed by MathParser
    virtual void SetScalar(const std::string& param, MaterialType matType, Global::ComplexPart dataType )
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }


    //! set a real vector
    virtual void SetVector(const Vector<Double>& param, MaterialType matType, Global::ComplexPart dataType)
    {
      EXCEPTION("not implemented");
    }

    //! set a real material tensor
    virtual void SetTensor(const Matrix<Double>& param, MaterialType matType, Global::ComplexPart dataType)
    {
      EXCEPTION("not implemented");
    }
    

    //! set a complex material tensor
    virtual void SetTensor(const Matrix<Complex>& param, MaterialType matType, Global::ComplexPart dataType ) 
    {
      EXCEPTION("not implemented");
    }
    
    //! Set a nonlinear isotropic approximation
    virtual void SetNonLinMatIso( MaterialType matType, MatDescriptorNl& data );
    
    //! Set a nonlinear anisotropic approximation
    virtual void SetNonLinMatAniso( MaterialType matType, StdVector<MatDescriptorNl>& data );
    
    //! Set a coefficient function
    virtual void SetCoefFct( MaterialType matType, PtrCoefFct coef );


    //! get a string material parameter
    virtual void GetScalar( std::string& param, MaterialType matType) const;
    
    //! get a integer material parameter
    virtual void GetScalar( Integer& param, MaterialType matType) const
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }
    
    void GetScalar( Integer& param, MaterialType matType, Global::ComplexPart dataType) const;

    //! get a scalar real material parameter
    virtual void GetScalar( Double& param, MaterialType matType, Global::ComplexPart dataType ) const 
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }

    //! get a scalar complex material parameter
    virtual void GetScalar( Complex& param, MaterialType matType, Global::ComplexPart dataType ) const
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }

    //! get a real vector
    virtual void GetVector( Vector<Double>& param, MaterialType matType, Global::ComplexPart dataType ) const
    {
      EXCEPTION("not implemented");      
    }

    //! get a complex vector
     virtual void GetVector( Vector<Complex>& param, MaterialType matType, Global::ComplexPart dataType ) const
     {
       EXCEPTION("not implemented");
     }

    //! get a real material tensor
    virtual void GetTensor( Matrix<Double>& param, MaterialType matType, 
                            Global::ComplexPart dataType, SubTensorType = FULL ) const
    {
      EXCEPTION("not implemented");      
    }
     
    //! get a complex material tensor
    virtual void GetTensor( Matrix<Complex>& param, MaterialType matType, Global::ComplexPart dataType, SubTensorType = FULL ) const
    {
      EXCEPTION("not implemented");      
    }
    
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

    //! Rotates the tensor in a way that is represents the attached
    //! coordinate system behaviour (cartesian, cylindri, spherical)
    //! in this point
    virtual void RotateTensorByPointCoord( const Vector<Double> &coord, 
                                           MaterialType matType );

    //! Pass coordinate system to material
    void SetCoordSys( CoordSystem* system ) {coosy_ = system;}

    //! Get coordinate system from material
    CoordSystem* GetCoordSys() { return coosy_; }

    // ======================================================================
    //  Hysteresis Related Information
    // ======================================================================
    //@{ \name Hysteresis Related Information

    //Initialize hysteresis
    virtual void InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList, 
                           bool isInverse = false, bool computeInverse = false );

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


    //! Compute Rayleigh parameters
    void ComputeRayleighDamping( std::string& alpha, std::string& beta,
                                 Double dampFreq, Double RatioDeltaF,
                                 bool adjustDamping, bool isHarmonic );


  protected:

    //! Error for material type not defined
    void matTypeNotAllowed(MaterialType matType, const std::string& dim ) const;

    //! data type not allowed in set/get-function
    void dataTypeNotAllowed4SetGet(Global::ComplexPart datType, const std::string& msg ) const;

    //! Error for data type not allowed
    void dataTypeNotAllowed(Global::ComplexPart datType, MaterialType matType ) const;

    //! Error for material type not in file
    void matTypeNotInDataBase(MaterialType matType, const std::string& dim ) const;

    //! data type not allowed in set-function
    void setMakesNoSense(Global::ComplexPart datType, const std::string& msg ) const;

    //! Error for not available subtype of tensor
    static void subTensorNotAvailable(MaterialType matType, SubTensorType subTensor);

    //! Rotate tensorial coefficient function
    virtual void PerformRotation( const Matrix<Double>& rotMatrix, PtrCoefFct& rotatedCoef,
                                  PtrCoefFct origCoef);

    //! name of material database
    std::string materialDatabaseName_;

    //! name of material file
    std::string matFileName_;

    //! set, which knows about the allowed material parameters for a material class
    std::set<MaterialType> isAllowed_;

    //! set, which knows, which material parameters have been set
    std::set<MaterialType> isSet_;

    //! set, which knows about if the material data is complex or just real
    std::set<MaterialType> isComplex_;

    //! map, which knows about the material data defined by strings
    stringMap stringParams_;

    //! map, which knows about the material data defined by strings
    integerMap integerParams_;

    //! map, which knows about the scalar material parameters being set during read in 
    scalarMap scalarParams_;
    
    //! map for real part of scalar material data (string representation)
    stringMap scalarStringParamsReal_;
    
    //! map for mathParser Handles of real part of scalar material parameters
    handleMap scalarStringHandlesReal_;
    
    //! map for imaginary part of scalar material data (string representation)
    stringMap scalarStringParamsImag_;

    //! map for mathParser Handles of real part of scalar material parameters
    handleMap scalarStringHandlesImag_;
        
    //! map, which knows about the actual tensorial material parameters 
    vectorMap vectorParams_;

    //! map, which knows about the actual tensorial material parameters 
    tensorMap tensorParams_;

    //! map, which knows about the original tensorial material parameters before being rotated
    tensorMap tensorParamsOrig_;

    //! map storing the isotropic nonlinear material parameters
    nonLinIsoMap nonlinIsoParams_;
    
    //! map storing the anisotropic nonlinear material parameters
    nonLinAnisoMap nonlinAnisoParams_;
    
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
    MathParser *  mp_;
    
    //! Pointer to attaches coordinate system
    CoordSystem * coosy_;

    //! Store symmetryType per material data type
    std::map<MaterialType, SymmetryType> symmetryType_;

    //! hysteresis object
    Hysteresis * hyst_;

    //! hysteresis object
    Hysteresis * hystY_;

    //! hysteresis object
    Hysteresis * hystZ_;

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

    Matrix<Double> vecXprevious_; //! previous Xval in hysteresis
    Matrix<Double> vecYprevious_; //! previous Yval in hysteresis

    Matrix<Double> actDiffVal4VecHyst_;
    Matrix<Double> previousDiffVal4VecHyst_;

    std::map<UInt, UInt> globalElem2Local_;

    //! yes, piezoelectric micro-model is switched on
    bool isPiezoMicroModel_;	

    //! object for piezo-micro-modeling
    PiezoMicroModelBK* piezoMicroModel_;

  };

} // end of namespace

#endif
