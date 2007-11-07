// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef BASEMATERIAL_DATA
#define BASEMATERIAL_DATA

#include <map>
#include <set>

#include <General/environment.hh>
#include <Matrix/matrix.hh>
#include <Utils/vector.hh>
#include <Utils/ApproxData.hh>

namespace CoupledField {

  // forward class declarations
  class CoordSystem;
  class Hysteresis;
  class ElemList;

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

    typedef enum {GENERAL, ISOTROPIC, ORTHOTROPIC, TRANS_OTHOTROP } SymmetryType;

    //! Default constructor
    BaseMaterial();

    //! Destructor
    virtual ~BaseMaterial();


    //! Trigger finalization of mataterial (calculation of rotated matrices)
    virtual void Finalize() {};

    //! set the name of the material set
    void SetName(const Char* name) {
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

    //! set file name containing the nonlinear data
    void SetNonlinFileName( const Char *filename) {
      nonlinFileName_.assign( filename );
    }

    //! get nonlinear file name
    std::string& GetNonlinFileName() {
      return nonlinFileName_;
    }

    //! set, which nonlinear curves are needed by forms
    void NeedApproxMatCurve( ApproxMaterialCurves type );

    virtual ApproxData* GetNonlinFncBH() {
      EXCEPTION("BaseMaterial: GetNlinFncBH() not implemented");
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
    virtual void SetScalar(const std::string& param, MaterialType matType)
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }

    //! set a scalar integer material parameter
    virtual void SetScalar(int param, MaterialType matType);

    //! set a scalar real material parameter
    virtual void SetScalar(double param, MaterialType matType, DataType dataType ) = 0;


    //! set a scalar complex material parameter
    virtual void SetScalar(Complex param, MaterialType matType, DataType dataType )
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }


    //! set a real vector
    virtual void SetVector(const Vector<Double>& param, MaterialType matType, DataType dataType)
    {
      EXCEPTION("not implemented");
    }

    //! set a real material tensor
    virtual void SetTensor(const Matrix<Double>& param, MaterialType matType, DataType dataType)
    {
      EXCEPTION("not implemented");
    }


    //! set a complex material tensor
    virtual void SetTensor(const Matrix<Complex>& param, MaterialType matType, DataType dataType ) 
    {
      EXCEPTION("not implemented");
    }


    //! get a string material parameter

    virtual void GetScalar( std::string& param, MaterialType matType) const
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }
    
    
    //! get a integer material parameter
    virtual void GetScalar( Integer& param, MaterialType matType) const
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }
    
    void GetScalar( Integer& param, MaterialType matType, DataType dataType) const;

    //! get a scalar real material parameter
    virtual void GetScalar( Double& param, MaterialType matType, DataType dataType ) const = 0;

    //! get a scalar complex material parameter
    virtual void GetScalar( Complex& param, MaterialType matType, DataType dataType ) const
    {
      EXCEPTION("not implemented for " << materialDatabaseName_);
    }

    //! get a real vector
    virtual void GetVector( Vector<Double>& param, MaterialType matType, DataType dataType ) const
    {
      EXCEPTION("not implemented");      
    }

    //! get a real material tensor
    virtual void GetTensor( Matrix<Double>& param, MaterialType matType, DataType dataType, SubTensorType = FULL ) const
    {
      EXCEPTION("not implemented");      
    }
     
    //! get a complex material tensor
    virtual void GetTensor( Matrix<Complex>& param, MaterialType matType, DataType dataType, SubTensorType = FULL ) const
    {
      EXCEPTION("not implemented");      
    }
    

    //! rotate a material tensor by rotation angles given in degree
    virtual void RotateTensorByRotationAngles( const Vector<Double>& rotAngle, MaterialType matType, bool persistent = false );

    //! Rotate all tensor material parameters by given rotation angle
    virtual void RotateAllTensorsByRotationAngles(  const Vector<Double>& rotAngle, bool persistent = false );

    //! Rotates the tensor in a way that is represents the attached
    //! coordinate system behaviour (cartesian, cylindri, spherical)
    //! in this point
    virtual void RotateTensorByPointCoord( const Vector<Double> &coord, MaterialType matType );

    //! Pass coordinate system to material
    void SetCoordSys( CoordSystem* system ) {coosy_ = system;}

    //! Get coordinate system from material
    CoordSystem* GetCoordSys() { return coosy_; }

    //Initialize approximations of nonlinear curves
    virtual void InitApproxCurves() {;};

    //Initialize hysteresis
    virtual void InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList, 
                           bool isInverse = false, bool computeInverse = false );

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

    virtual Double ConvertVec2ScalarWithSign( Vector<Double>& vecVal );

    //set values for differential material approach
    virtual void SetPreviousHystVal( UInt nrElem, Double Xval ) {
      Error( "SetPreviousHystVal not implemented", __FILE__, __LINE__ );
    };

    //set values for differential material approach
  virtual void SetPreviousHystVal( UInt nrElem, Vector<Double>& Xval );

    //! compute scalar differential parameter
    virtual Double ComputeScalarDiffVal( UInt nrElem, Double Xval ) {
      Error( "ComputeScalarDiffValue not implemented", __FILE__, __LINE__ );
      return 1.0;
    };

    //! compute scalar differential parameter
    virtual Double ComputeScalarDiffVal( UInt nrElem, Vector<Double>& Xval );

    //! compute scalar differential parameters
    virtual void ComputeScalarDiffValues( UInt nrElem, Vector<Double>& in,
                                            Vector<Double>& scalarValues ) {
      Error( "ComputeScalarDiffValues not implemented", __FILE__, __LINE__ );
    };

    //! computes the scalar hystereis value
    virtual Double ComputeScalarHystVal( UInt nrElem, Double Xval ) {
      Error( "ComputeScalarHystVal not implemented", __FILE__, __LINE__ );
      return 1.0;
    };

    //! computes the scalar hystereis value
    virtual Double ComputeScalarHystVal( UInt nrElem, Vector<Double>& Xval );

    //! compute the vector hysteresis values
    virtual void ComputeVectorHystVal( UInt nrElem, Vector<Double>& Xin, 
                                       Vector<Double>& Yout ) {
      Error( "ComputeVectorHystVal not implemented", __FILE__, __LINE__ );
    };

    //! computes the scalar hystereis value
    virtual Double GetScalarHystVal( UInt nrElem );

    //! computes the scalar hystereis value
    virtual Double GetScalarHystPrevVal( UInt nrElem );

    //! get vector of hysteresis model
    virtual void GetVectorHystVal( UInt nrElem, Vector<Double>& Val ) {
      Error( "GetVectorHystVal not implemented", __FILE__, __LINE__ );
    };

    //! Compute and set damping parameters alpha and beta 
    void ComputeRayleighDamping(Double dampFreq, Double RatioDeltaF);

  protected:

    //! Error for material type not defined
    void matTypeNotAllowed(MaterialType matType, const std::string& dim ) const;

    //! data type not allowed in set/get-function
    void dataTypeNotAllowed4SetGet(DataType datType, const std::string& msg ) const;

    //! Error for data type not allowed
    void dataTypeNotAllowed(DataType datType, MaterialType matType ) const;

    //! Error for material type not in file
    void matTypeNotInDataBase(MaterialType matType, const std::string& dim ) const;

    //! data type not allowed in set-function
    void setMakesNoSense(DataType datType, const std::string& msg ) const;

    //! Error for not available subtype of tensor
    void subTensorNotAvailable(MaterialType matType, SubTensorType subTensor) const;

    //! rotate a tensor
    virtual void PerformRotation( Matrix<Complex>& rotMatrix,  Matrix<Complex>& matMatrix,
				  const Matrix<Complex>& origMatMatrix);

    //! name of material database
    std::string materialDatabaseName_;

    //! name of material file
    std::string matFileName_;

    //! name of file containing the nonlinear data
    std::string nonlinFileName_;

    //! set, which knows, which material parameters have been set
    std::set<ApproxMaterialCurves> needApproxMatCurves_;

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

    //! map, which knows about the actual tensorial material parameters 
    vectorMap vectorParams_;

    //! map, which knows about the actual tensorial material parameters 
    tensorMap tensorParams_;

    //! map, which knows about the original tensorial material parameters before being rotated
    tensorMap tensorParamsOrig_;

    //! Pointer to attaches coordinate system
    CoordSystem * coosy_;

    SymmetryType symmetryType_;

    //! hysteresis object
    Hysteresis * hyst_;

    //! hysteresis object
    Hysteresis** vecHyst_;

    //!
    bool isHystInverse_;

    //!
    bool computeHystInverse_;

    //! hysteresis set 
    bool isHysteresis_;

    Vector<Double> Xprevious_; //! previous Xval in hysteresis
    Vector<Double> Yprevious_; //! previous Yval in hysteresis

    std::map<UInt, UInt> globalElem2Local_;
  };

  std::ostream& operator << ( std::ostream & , const  BaseMaterial &);

} // end of namespace

#endif
