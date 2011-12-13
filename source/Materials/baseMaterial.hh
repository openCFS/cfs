// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef BASEMATERIAL_DATA
#define BASEMATERIAL_DATA

#include <stddef.h>
#include <map>
#include <set>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
//#include "Utils/piezoMicroModel.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {

class ApproxData;
  // forward class declarations
  class CoordSystem;
  class ElemList;
  class Hysteresis;
  class PiezoMicroModelBK;


  struct nlMatDescriptor {
    ApproxCurveType approxType;
    std::string fileName;
    Double measAccuracy;
    Double maxVal;
  };


  //! Class for Material Data
  /*! 
    Base class for handling material data
  */

  class BaseMaterial {

  public:
    typedef std::map<MaterialType, Matrix<Complex> > tensorMap;
    typedef std::map<MaterialType, Vector<Complex> > vectorMap;
    typedef std::map<MaterialType, Complex > scalarMap;
    typedef std::map<MaterialType, std::string > stringMap;
    typedef std::map<MaterialType, Integer > integerMap;
    typedef std::map<MaterialType, MathParser::HandleType > handleMap;
    typedef std::map<MaterialType, nlMatDescriptor> nonLinMatInfoMap;

    typedef enum {GENERAL, ISOTROPIC, ORTHOTROPIC, TRANS_OTHOTROP } SymmetryType;

    //! Default constructor
    BaseMaterial();

    //! Destructor
    virtual ~BaseMaterial();


    //! Trigger finalization of mataterial (calculation of rotated matrices)
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

    //! get info, which material parameter is set
    std::set<MaterialType> GetIsSetInfo() const
    { return isSet_;};

    //! get infor, which material parameter is complex
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

    //! set, which material parameter has to be approx. /interp. 
    //! and is needed by bi- and linear-forms
    void NeedApproxMatCurve(  MaterialType type );

    virtual ApproxData* GetNonlinFnc( MaterialType matType ) {
      EXCEPTION("BaseMaterial: GetNlinFnc() not implemented");
      return NULL;
    };

    //! Query if a given parameter is set
    virtual bool IsSet( MaterialType matType ) const;

    //! set the symmetry type
    void SetSymmetryType(SymmetryType symType) {
      symmetryType_=symType; 
    };

    //! get the symmetry type
    SymmetryType GetSymmetryType() const {
      return symmetryType_; 
    };

   //! set a scalar string material parameter
    virtual void SetScalar(const std::string& param, MaterialType matType);

    //! set a scalar integer material parameter
    virtual void SetScalar(int param, MaterialType matType);

    //! set a scalar real material parameter
    virtual void SetScalar(double param, MaterialType matType, Global::ComplexPart dataType ) = 0;


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

    //! set info for the nonliear approx. / interp. of a material
    virtual void SetNonLinMat(MaterialType matType, nlMatDescriptor& data ) {
      nonLinMatInfo_[matType] = data;
    }

    //! get a string material parameter
    virtual void GetScalar( std::string& param, MaterialType matType) const;
    
    //! get a integer material parameter
    virtual void GetScalar( Integer& param, MaterialType matType) const
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }
    
    void GetScalar( Integer& param, MaterialType matType, Global::ComplexPart dataType) const;

    //! get a scalar real material parameter
    virtual void GetScalar( Double& param, MaterialType matType, Global::ComplexPart dataType ) const = 0;

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

    //! get a real material tensor
    virtual void GetTensor( Matrix<Double>& param, MaterialType matType, Global::ComplexPart dataType, SubTensorType = FULL ) const
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

    //Initialize approximations of nonlinear curves
    virtual void InitApproxCurves() {;};

    //======================================= Hysteresis methods ============================
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


//     //! compute scalar differential parameter: scalar hyst
//     virtual Double ComputeScalarDiffVal( UInt nrElem, Double Xval ) {
//       EXCEPTION( "ComputeScalarDiffValue not implemented" );
//       return 1.0;
//     };

//     //! compute scalar differential parameter: vector hyst
//     virtual Double ComputeScalarDiffVal( UInt nrElem, Vector<Double>& Xval );

    //! computes the scalar hystereis value
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

    //! get scalar hystereis value
    virtual Double GetScalarHystVal( UInt nrElem );

    //! get previous scalar hystereis value
    virtual Double GetScalarHystPrevVal( UInt nrElem );

    //! get vector of hysteresis value
    virtual void GetVectorHystVal( UInt nrElem, Vector<Double>& Val ) {
      EXCEPTION( "ComputeVectorHystVal not implemented" );
    };


    //========================== micro-piezoelectric-model: start==================

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

    //! rotate a tensor
    virtual void PerformRotation( Matrix<Complex>& rotMatrix,  Matrix<Complex>& matMatrix,
				  const Matrix<Complex>& origMatMatrix);

    //! name of material database
    std::string materialDatabaseName_;

    //! name of material file
    std::string matFileName_;

    //! contains all info to a nonlinear material parameter
    //! e.g., approximation Type, accuracy, etc.
    nonLinMatInfoMap nonLinMatInfo_;

    //! set, which knows the material parameter to be approximated
    std::set<MaterialType> needApproxMatCurves_;

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
    
    //! map for imagomary part of scalar material data (string representation)
    stringMap scalarStringParamsImag_;

    //! map for mathParser Handles of real part of scalar material parameters
    handleMap scalarStringHandlesImag_;
        
    //! map, which knows about the actual tensorial material parameters 
    vectorMap vectorParams_;

    //! map, which knows about the actual tensorial material parameters 
    tensorMap tensorParams_;

    //! map, which knows about the original tensorial material parameters before being rotated
    tensorMap tensorParamsOrig_;

    //! Pointer to math parser instance
    MathParser *  mp_;
    
    //! Pointer to attaches coordinate system
    CoordSystem * coosy_;

    SymmetryType symmetryType_;

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

    //! object for piezo-micor-modeling
    PiezoMicroModelBK* piezoMicroModel_;

  };

} // end of namespace

#endif
