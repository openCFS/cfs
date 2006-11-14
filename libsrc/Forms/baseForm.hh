#ifndef FILE_BASEFORM_
#define FILE_BASEFORM_

#include "Utils/StdVector.hh"
#include "Utils/vector.hh"
#include "Elements/basefe.hh"
#include "Domain/surfElem.hh"
#include "Materials/baseMaterial.hh"
#include "Domain/entityList.hh"

#ifndef INTEGLIB
#include "Utils/mathParser/mathParser.hh"
#endif

namespace CoupledField
{
  
  // forward class declaration
  class AnsatzFct;

  //! Base class for calculation of bdb element matrices
  class BaseForm 
  {
  public:

    //! Constructor
    BaseForm(BaseMaterial* matData, SubTensorType subTensor = FULL,
             bool coordUpdate = false );

    //! Deconstructor
    virtual ~BaseForm();



    //! Return the name of this (bi)linearform
    const std::string& GetName() const { return name_;}

    //! Return true if element vector/matrix depends on solution 
    //! or its derivatives
    bool IsSolDependent() { return isSolDependent_; }

#ifndef INTEGLIB
    //! Virtual function
    virtual void CalcElementMatrix( Matrix<Double>& stiffMat,
                                    EntityIterator& ent1, 
                                    EntityIterator& ent2 ){
      
      Error( "Not implemented here", __FILE__, __LINE__ );}
    
    //! Virtual function, implemented in derived classes
    virtual void CalcComplexElementMatrix( Matrix<Complex> & StiffMat,
                                           EntityIterator& ent1, 
                                           EntityIterator& ent2,
                                           Double & beta, Double & omega) {
      Error( "Not implemented here", __FILE__, __LINE__ );}

#endif

    //
    virtual void SetNonLinMethod(std::string atype) {;};

    //! needed for fractional damping model
    virtual void SetFracDamping() {;};

    //! needed for fractional damping model
    virtual void UnsetFracDamping() {;};

    //! needed for fractional damping model
    virtual bool IsFracDamping()
    {return isFracDamping_;};

    //! set additional multiplicative factor for matrix
    virtual void SetFactor(Double factor) {;};

    //! set second multiplicative factor for matrix
    virtual void SetSecondFactor(Double factor) {;};

    //! sets pointer to actual material
    void SetMaterial( BaseMaterial* matPtr );

    //!
    void SetAnsatzFct( shared_ptr<AnsatzFct> actFct1,
                       shared_ptr<AnsatzFct> actFct2 = shared_ptr<AnsatzFct>() );

    bool IsCoordUpdate() { return coordUpdate_; }
    
    //! query pointer to actual material
    BaseMaterial* GetMaterial() {
      ENTER_FCN( "BaseForm::GetMaterial" );
      return ptMaterial;
    }

    //! Set solution class for non-linear integrators
    void SetSolution( NodeStoreSol<Double>& sol ) {
      sol_ = &sol; }

    //! Set first time derivative of solution for non-linear integrators
    void SetSolDeriv1( NodeStoreSol<Double>& solDeriv1 ) {
      solDeriv1_ = &solDeriv1; }

    //! Set second time derivative of solution for non-linear integrators
    void SetSolDeriv2( NodeStoreSol<Double>& solDeriv2 ) {
      solDeriv2_ = &solDeriv2; }

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
    { intPoint_ = point; isSetIntPoint_ = true;};

    //!
    void SetMatDataType(DataType & pMatType)
    { matDataType_=pMatType; };

    //! set softening type for forms
    void SetSofteningModel(std::string type) {
      softeningModel_ = type;
    }

    //! set min/max of x,y,z coordinates form where PML starts
    virtual void SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer) {;};
      
#ifndef INTEGLIB
    //! Get reference element and coordinates from element iterator
    virtual void ExtractElemInfo( EntityIterator& it);
    
#endif

  protected:

    //! pointer to reference element
    BaseFE  * ptelem;   

#ifndef INTEGLIB
    //! current entity iterators
    EntityIterator it1_;
    EntityIterator it2_;

    //! pointer to ansatz fct
    shared_ptr<AnsatzFct> ansatzFct1_;
    shared_ptr<AnsatzFct> ansatzFct2_;
#endif

    //! flag indicating updated lagrangian formulation
    bool coordUpdate_;

    //! flag indicating of value of element matrix/vector depends on
    //! its solution or its derivatives
    bool isSolDependent_;

    //! matrix with coordinates of current element
    Matrix<Double> ptCoord_;
    
    //! pointer to material data
    BaseMaterial* ptMaterial;

    //! true for axisymmetric setup
    bool isaxi_;

    //! type of tensor (plabeStrain, ..)
    SubTensorType subTensorType_;

    //! name of (bi)linearform
    std::string name_;

    //
    Vector<Double> intPoint_;
    //
    bool isSetIntPoint_;

    DataType matDataType_;     //! default = realMaterialParamter, piezoMatType_ = imagMaterialParamter if we consider complex-valued material Paramter;

    bool isFracDamping_;   //!< if true Assemble::AssembleMatrices will retrieve an additional multiplicative factor

    FEMatrixType baseType_;  // base type: STIFFNESS, DAMPING, MASS

    // specifies model of softening
    std::string softeningModel_;

    //! softening part (bending or shear) of softening Model
    std::string softeningPart_;

    //maximal length of an edge within an element
    Double maxEdgeLength_;

    //minimal length of an edge within an element
    Double minEdgeLength_;

#ifndef INTEGLIB
    //! Handle for MathParser object
    MathParser::HandleType mHandle_;

    //! Pointer to MathParser object
    MathParser * mParser_;
#endif

    //! solution vector
    NodeStoreSol<Double>* sol_;

    //! first derivative of solution
    NodeStoreSol<Double>* solDeriv1_;

    //! second derivative of solution
    NodeStoreSol<Double>* solDeriv2_;

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
                           const std::map<RegionIdType, BaseMaterial*>
                           & materials );

    //! Set information for second interface side
    void SetSecondVoluInfo( const std::string & name,
                            const std::map<RegionIdType, BaseMaterial*>
                            & materials );

    //! set additional multiplicative factor for matrix
    void SetFactor(Double factor); 

    //!
    void SetFormulation(SolutionType aformulation) 
    { formulation_ = aformulation;};

 
#ifndef INTEGLIB
    //! Get reference element and coordinates from element iterator
    void ExtractElemInfo( EntityIterator& it);
#endif

  protected:

    //! Current surface element
    const SurfElem * actElem_;

    //! Normal pointing out of first volume element
    Vector<Double> normal_;

    //! Name of PDE on first interface side
    std::string firstPDEName_;

    //! Name of PDE on second interface side
    std::string secondPDEName_;

    //! Materials on first interface side
    std::map<RegionIdType, BaseMaterial *> firstMaterials_;

    //! Materials on second interface side
    std::map<RegionIdType, BaseMaterial *> secondMaterials_;
    
    //! Multiplicative factor for matrix
    Double factor_;

    //! formulation type
    SolutionType formulation_;

  };

} //end namespace

#endif // FILE_BASEFORM
