#ifndef BASEMATERIAL_DATA
#define BASEMATERIAL_DATA

#include <map>
#include <set>

#include <General/environment.hh>
#include <Matrix/matrix.hh>
#include <Utils/vector.hh>

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

    //! Query if a given parameter is set
    virtual bool IsSet( MaterialType matType ) const;

   //! set a scalar string material parameter
    virtual void SetScalar( std::string& param, const MaterialType& matType) {
      Error("SetScalar not implemented",__FILE__,__LINE__); };

    //! set a scalar integer material parameter
    virtual void SetScalar( Integer& param, const MaterialType& matType);

    //! set a scalar real material parameter
    virtual void SetScalar( Double& param, const MaterialType& matType,
			    const DataType& dataType ) {
      Error("SetScalar not implemented",__FILE__,__LINE__); };


    //! set a scalar complex material parameter
    virtual void SetScalar( Complex& param, const MaterialType& matType, 
			    const DataType& dataType ) {
      Error("SetScalar not implemented",__FILE__,__LINE__); };


    //! set a real material tensor
    virtual void SetTensor( Matrix<Double>& param, const MaterialType& matType,
			    const DataType& dataType ){
      Error("SetTensor not implemented",__FILE__,__LINE__); };


    //! set the symmetry type
    void SetSymmetryType(SymmetryType symType) {
      symmetryType_=symType; 
    };

    //! get the symmetry type
    SymmetryType GetSymmetryType() const {
      return symmetryType_; 
    };

    //! set a complex material tensor
    virtual void SetTensor( Matrix<Complex>& param, const MaterialType& matType,
			    const DataType& dataType ) {
      Error("SetTensor not implemented",__FILE__,__LINE__); };

    //! get a string material parameter
    void GetScalar( std::string& param, const MaterialType& matType) const {
     Error("GetScalar not implemented",__FILE__,__LINE__); };

    void GetScalar( Integer& param, const MaterialType& matType, 
		    const DataType& dataType) const;

    //! get a scalar real material parameter
    virtual void GetScalar( Double& param, const MaterialType& matType, 
			    const DataType& dataType ) const {
      Error("GetScalar not implemented",__FILE__,__LINE__); };


    //! get a scalar complex material parameter
    virtual void GetScalar( Complex& param, const MaterialType& matType, 
			    const DataType& dataType ) const {
      Error("GetScalar not implemented",__FILE__,__LINE__); };


    //! get a real material tensor
    virtual void GetTensor( Matrix<Double>& param, const MaterialType& matType,
			    const DataType& dataType, 
			    const SubTensorType = FULL ) const {
      Error("GetTensor not implemented",__FILE__,__LINE__); };

    //! get a complex material tensor
    virtual void GetTensor( Matrix<Complex>& param, const MaterialType& matType,
			    const DataType& dataType,
			    const SubTensorType = FULL ) const {
      Error("GetTensor not implemented",__FILE__,__LINE__); };

    //! rotate a material tensor by rotation angles given in degree
    virtual void RotateTensorByRotationAngles( Vector<Double>& rotAngle, 
					       const MaterialType& matType);

    //! Rotates the tensor in a way that is represents the attached
    //! coordinate system behaviour (cartesian, cylindri, spherical)
    //! in this point
    virtual void RotateTensorByPointCoord( Vector<Double> coord,
                                           const MaterialType& matType );

    //! Pass coordinate system to material
    void SetCoordSys( CoordSystem* system ) {coosy_ = system;}

    //! Get coordinate system from material
    CoordSystem* GetCoordSys() { return coosy_; }

    //Initialize hysteresis
    void InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList);

  protected:

    //! Error for material type not defined
    void matTypeNotAllowed( const MaterialType& matType, 
			    std::string dim ) const;

    //! data type not allowed in set/get-function
    void dataTypeNotAllowed4SetGet( const DataType& datType, 
				    const std::string msg ) const;

    //! Error for data type not allowed
    void dataTypeNotAllowed( const DataType& datType, 
			     const MaterialType& matType ) const;

    //! Error for material type not in file
    void matTypeNotInDataBase( const MaterialType& matType, std::string dim ) const;

    //! data type not allowed in set-function

    void setMakesNoSense( const DataType& datType, std::string msg ) const;

    //! Error for not available subtype of tensor
    void subTensorNotAvailable( const MaterialType& matType, 
			      const SubTensorType subTensor ) const;

    //! rotate a tensor
    virtual void PerformRotation( Matrix<Complex>& rotMatrix,  Matrix<Complex>& matMatrix,
				  const Matrix<Complex>& origMatMatrix);

    //! name of material database
    std::string materialDatabaseName_;

    //! name of material file
    std::string matFileName_;

    //! name of file containing the nonlinear data
    std::string nonlinFileName_;

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
    tensorMap tensorParams_;

    //! map, which knows about the original tensorial material parameters before being rotated
    tensorMap tensorParamsOrig_;

    //! Pointer to attaches coordinate system
    CoordSystem * coosy_;

    SymmetryType symmetryType_;

    //! hysteresis object
    Hysteresis * hyst_;

    std::map<UInt, UInt> globalElem2Local_;
  };

  std::ostream& operator << ( std::ostream & , const  BaseMaterial &);

} // end of namespace

#endif
