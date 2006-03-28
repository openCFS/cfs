#ifndef FILE_BASEFORM_
#define FILE_BASEFORM_

#include "Utils/StdVector.hh"
#include "Utils/vector.hh"
#include "Elements/basefe.hh"
#include "Domain/surfElem.hh"
#include "Materials/baseMaterial.hh"

namespace CoupledField
{

  //! Base class for calculation of bdb element matrices
  class BaseForm 
  {
  public:

    //! Constructor
    BaseForm(BaseFE * aptelem, BaseMaterial* matData, 
	     SubTensorType subTensor = FULL );

    //! Constructor
    BaseForm(BaseMaterial* matData, SubTensorType subTensor = FULL);

    //! Constructor
    BaseForm(BaseFE * aptelem);

    //! Constructor
    BaseForm();

    //! Deconstructor
    virtual ~BaseForm();

    //! Virtual function, implemented in derived classes
    virtual void CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & StiffMat);

    //! Virtual function, implemented in derived classes
    virtual void CalcComplexElementMatrix(Matrix<Double>& ptCoord,
                                          Matrix<Complex> & StiffMat,
                                          Double & beta, Double & omega);
  
    /// Calculation of vector of right hand side 
    virtual void CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & result)
    {Error("CalcElemVector not implemented!",__FILE__,__LINE__);};

    virtual void CalcElemVector4Dip(Matrix<Double>& ptCoord, 
                                    const StdVector<UInt> & connecth, 
                                    Vector<Double> & Result, 
                                    const Vector<Double> gradN_x_P)
    { Error(" CalcElemVector4Dip is not implemented for this class",__FILE__,__LINE__);};

    /// Calculation of vector of right hand side given from quadrupole contribution
    virtual void CalcElemVector4Quad(Matrix<Double>& ptCoord,
                                     const StdVector<UInt> & connecth,
                                     const Matrix<Double> & FlowData, 
                                     Vector<Double> & Result)
    { Error(" CalcElemVector4Quad is not implemented for this class",__FILE__,__LINE__);};

    /// Extraction of element velocity values from total flowdata matrix to a matrix (connecth, dim)
    virtual void GetQttiesOfElement(Matrix<Double>& elVec,
                                    const Matrix<Double>& FlowData,
                                    const StdVector<UInt>& connecth, 
                                    UInt matrixRow)
    { Error(" GetQttiesOfElement is not implemented for this class",__FILE__,__LINE__);};

    //! Prints the bilinear form
    virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

    //
    virtual void SetNonLinMethod(std::string atype) {;};

    //! needed for fractional damping model
    virtual void SetFracDamping() {;};

    //! needed for fractional damping model
    virtual void UnsetFracDamping() {;};

    //! needed for fractional damping model
    virtual Boolean IsFracDamping()
    {return isFracDamping_;};

    //! needed for fractional damping model
    virtual void SetRaylDamping() 
    { isRaylDamping_ = TRUE; };


    //! needed for fractional damping model
    virtual Boolean IsRaylDamping()
    {return isRaylDamping_;};

    //! set additional multiplicative factor for matrix
    virtual void SetFactor(Double factor) {;};

    //! set second multiplicative factor for matrix
    virtual void SetSecondFactor(Double factor) {;};

    //! sets pointer to actual element
    void SetElemPtr(BaseFE * elemPtr){ptelem = elemPtr;};

    //! sets pointer to actual material
    void SetMaterial( BaseMaterial* matPtr );

    //! query pointer to actual material
    BaseMaterial* GetMaterial() {
      ENTER_FCN( "BaseForm::GetMaterial" );
      return ptMaterial;
    }

    //! get base type of biliearform: STIFFNESS, DAMPING, MASS
    virtual FEMatrixType GetBaseType() 
    { return baseType_; };

    //! sets actual element solution
    virtual void SetActElemSol(Matrix<Double>& disp)
    {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

    //! sets actual element solution
    virtual void SetActElemSol(CFSMatrix & disp)
    {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

    //! sets actual element solution
    virtual void SetActElemSol(Vector<Double>& disp)
    {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

    //! sets actual first time derivative of element solution
    virtual void SetActElemSolDeriv1(Matrix<Double>& disp)
    {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

    //! sets actual first time derivative of element solution
    virtual void SetActElemSolDeriv1(Vector<Double>& disp)
    {Error("SetActElemSolDeriv1 not implemented!",__FILE__,__LINE__);};

    //! sets actual 2nd time derivative of element solution
    virtual void SetActElemSolDeriv2(Vector<Double>& disp)
    {Error("SetActElemSolDeriv2 not implemented!",__FILE__,__LINE__);};

    //! reads the values y(x) out of the file with name fncName  
    void ReadNlinFunc( std::string fncName, Vector<double> &xval,
                       Vector<Double> &yval ) {
      Error("ReadNlinFunc not implemented!", __FILE__, __LINE__ );
    };

    //! set the integration point
    void SetIntPoint(Vector<Double> point)
    { intPoint_ = point; isSetIntPoint_ = TRUE;};

    //!
    void SetDofZero(UInt posdof)
    {dofzero_ = posdof; };

    //!
    void SetPiezoMaterialType(piezoMaterialType & pMatType)
    { piezoMatType_=pMatType; };

    //!
    void SetMaterialArray(Matrix<Double>* mat)
    { materialArray_ = mat; };

    //! set softening type for forms
    void SetSofteningModel(std::string type) {
      softeningModel_ = type;
    }

    //!
    void SetSubdomain(UInt sd)
    {actSD_ = sd; };

    //!
    void SetElemNr(UInt nr)
    {actElemNr_ = nr; };

    //! set min/max of x,y,z coordinates form where PML starts
    virtual void SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer) {;};

    void SetFrequency(Double freq)
    {frequency_ = freq;};

  protected:

    //! pointer to reference element
    BaseFE  * ptelem;   
    
    //! pointer to material data
    BaseMaterial* ptMaterial ;

    //! true for axisymmetric setup
    Boolean isaxi_;

    //! type of tensor (plabeStrain, ..)
    SubTensorType subTensorType_;

    //
    Vector<Double> intPoint_;
    //
    Boolean isSetIntPoint_;

    piezoMaterialType piezoMatType_;     //! default = realMaterialParamter, piezoMatType_ = imagMaterialParamter if we consider complex-valued material Paramter;

    Boolean isFracDamping_;   //!< if true Assemble::AssembleMatrices will retrieve an additional multiplicative factor
    Boolean isRaylDamping_;   //!< if true Assemble::AssembleMatrices will retrieve an additional multiplicative factor
    UInt dofzero_;   //!< for multidof-handling, where one dof is zero (e.g. piezoelectric PDE)

    FEMatrixType baseType_;  // base type: STIFFNESS, DAMPING, MASS

    Matrix<Double>* materialArray_;

    UInt actSD_;
    UInt actElemNr_;

    // specifies model of softening
    std::string softeningModel_;

    //! softening part (bending or shear) of softening Model
    std::string softeningPart_;

    //maximal length of an edge within an element
    Double maxEdgeLength_;

    //minimal length of an edge within an element
    Double minEdgeLength_;

    //! current frequency
    Double frequency_;

  private:

    //! Should we delete the material data object?
    bool delMatDataAtEnd_;

  };


  //! Base class for surface integrators

  //! This class defines an abstract interface for all kinds of surface
  //! inetgrators. Since surface elements have no own material, they need
  //! information about ther one or two volume neighbours and their materials.
  //! Additionally, often the normal of the surface element is needed in order
  //! to computa a surface integral.
  class SurfForm : public BaseForm {

  public: 
    
    //! standard constructor
    SurfForm();

    //! standard destructor
    virtual ~SurfForm();
    
    //! Set pointer to surface element
    void SetSurfElem( SurfElem * ptSurfElem);

    //! Set normal pointing out of first volume element
    void SetFirstVoluNormal( Vector<Double> & n );

    //! Set information of first interface side
    void SetFirstVoluInfo( const std::string & name,
                           const StdVector<RegionIdType> & regionIds,
                           const StdVector<BaseMaterial*>& materials );

    //! Set information for second interface side
    void SetSecondVoluInfo( const std::string & name,
                            const StdVector<RegionIdType> & regionIds,
                            const StdVector<BaseMaterial*>& materials );

    //! set additional multiplicative factor for matrix
    void SetFactor(Double factor); 

    //!
    void SetFormulation(SolutionType aformulation) 
    { formulation_ = aformulation;};

  protected:

    //! Current surface element
    SurfElem * actElem_;

    //! Normal pointing out of first volume element
    Vector<Double> normal_;

    //! Name of PDE on first interface side
    std::string firstPDEName_;

    //! Name of PDE on second interface side
    std::string secondPDEName_;

    //! Region Ids on first interface side
    StdVector<RegionIdType> firstRegionIds_;

    //! Region ids on second interface side
    StdVector<RegionIdType> secondRegionIds_;

    //! Materials on first interface side
    StdVector<BaseMaterial*> firstMaterials_;

    //! Materials on second interface side
    StdVector<BaseMaterial*> secondMaterials_;

    //! Multiplicative factor for matrix
    Double factor_;

    //! formulation type
    SolutionType formulation_;


  };

} //end namespace

#endif // FILE_BASEFORM
