// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEFORM_
#define FILE_BASEFORM_

#include "Utils/StdVector.hh"
#include "MatVec/vector.hh"
#include "Elements/basefe.hh"
#include "Domain/surfElem.hh"
#include "Materials/baseMaterial.hh"
#include "Domain/entityList.hh"
#include "Optimization/Design/DesignElement.hh"
#include "PDE/timestepping.hh"

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

    /** For bimaterial optimization it is good to be sure in what 'mode' the forms are initialized.
     * This is important for MassInt and other scalar forms where the material data is given
     * as scalar in the constructor.
     * This class allows an optional form constructor which describes the material such that we can be
     * sure that we are in a mode we expect. To use the class is optional!
     * The scalar material properties are then obtained each time. Eventually BaseMaterial needs a boost.
     * If the material descriptor is used, the scalard material value from the other constructor is not
     * used -> see MassInt.
     * Does not work for tensors but could be extended.
     * It might be generally a good idea to use the MaterialDescriptor */
     class MaterialDescriptor
     {
     public:
       /** NOT_SET: the material descriptor is not used
        *  SCALAR: material->GetScalar(mat_1)
        *  MINUS_SCALAR: -1 * SCALAR in acoustic-mech coupling case when not pressure formulation
        *  MAT_1_MAT_1_BY_MAT_2: material->GetScalar(mat_1)^2 / material->GetScalar(mat_2) standard acoustic mass
        *  MINUS_MAT_1_MAT_1_BY_MAT_2 additional -1 in mech coupling case when not pressure formulation */
       enum Type { NOT_SET = -1, SCALAR = 0, MINUS_SCALAR, MAT_1_MAT_1_BY_MAT_2, MINUS_MAT_1_MAT_1_BY_MAT_2};

       /** The default constructor sets to NOT_SET */
       MaterialDescriptor();
       MaterialDescriptor(Type type, MaterialClass mat_class, MaterialType mat_1, Global::ComplexPart dataType);
       MaterialDescriptor(Type type, MaterialClass mat_class, MaterialType mat_1, MaterialType mat_2, Global::ComplexPart dataType);

       /** is the material descriptor set? */
       bool Enabled() const { return type != NOT_SET; }

       /** extract the value from material by type if not NOT_SET, convenience function */
       double GetScalar(BaseMaterial* bm);

       /** Applies the ersatz material factor if it applies, eventually with bimaterial, otherwise the default.
        * This is NOT the factor but the applied factor!
        * @param default_mat_value if the material descriptor is not initialized. */
       double GetErsatzMaterial(BaseForm* form, const Elem* elem, double default_mat_value);

       Type type;
       MaterialType mat_1;
       MaterialType mat_2;

       MaterialClass mat_class;

       Global::ComplexPart data_type;
     private:
       void Init(Type type, MaterialClass mat_class, MaterialType mat_1, MaterialType mat_2, Global::ComplexPart data_type);
     };


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

    /** allow maniputlation fo the flag, used for optimization loops */
    void SetSolDependent(bool value) { isSolDependent_ = value;  } 

    //! Return true if element vector/matrix is complex
    bool IsComplex() { return isComplex_; }

    /** This is a SIMP optimization helper. It works for ADB and BDB forms and returns the material tensor
     * @param factor this scales the material tensor, it might be 1.0
     * @param derivative only interesting in the bimat case
     * @param bimat optional (otherwise NULL) bimaterial interpolation
     * @param out result */
    void GetScaledMaterial(Double factor, bool derivative, BaseMaterial* bimat, Matrix<Double>& out);

#ifndef INTEGLIB
    //! Virtual function
    virtual void CalcElementMatrix( Matrix<Double>& stiffMat,
                                    EntityIterator& ent1, 
                                    EntityIterator& ent2 ){
      
      EXCEPTION( "Not implemented here" );}

    //! Virtual function
    /** this version is for convenience, it gets called from ErsatzMaterial, 
     *  it is overwritten in linElastInt which does care about parametric material optimization
     *  all other desccendents don't have it implemented, so they are called without direction */
    virtual void CalcElementMatrix( Matrix<Double>& stiffMat,
                                    EntityIterator& ent1, 
                                    EntityIterator& ent2,
                                    const DesignElement::Type ){
      CalcElementMatrix( stiffMat, ent1, ent2 );
    }

    //! Virtual function
    virtual void CalcElementMatrix( Matrix<Complex>& stiffMat,
                                    EntityIterator& ent1, 
                                    EntityIterator& ent2 ){
      
      EXCEPTION( "Not implemented here" );}
    
#endif

    virtual void CalcMinMaxStrain( EntityIterator& ent1, 
                                   EntityIterator& ent2 ){;};

    virtual void ResetMinMaxStrain(){;};

    
    //! define diagonal mass matrix
    virtual void SetDiagMass() {
      EXCEPTION( "Not implemented here" );
    };

    //
    virtual void SetNonLinMethod(NonLinMethodType atype) {;};

    //! set second multiplicative factor for matrix
    virtual void SetSecondFactor( const std::string& factor) {;};

    //! sets pointer to actual material
    void SetMaterial( BaseMaterial* matPtr );

    //!
    void SetAnsatzFct( shared_ptr<AnsatzFct> actFct1,
                       shared_ptr<AnsatzFct> actFct2 = shared_ptr<AnsatzFct>() );

    //! Query if bilinearform uses updated lagrangian formulation
    bool IsCoordUpdate() { return coordUpdate_; }

    //! Set whether to use Updated Coordinates (ShapeOptimization does this)
    void SetUseCoordUpdate(bool coordUpdate){
      coordUpdate_ = coordUpdate;
    }

    //! Query if resulting element matrix is symmetric
    bool IsSymmetric() { return isSymmetric_; }
    
    //! query pointer to actual material
    BaseMaterial* GetMaterial() {
      return ptMaterial;
    }

    //! Set solution class for non-linear integrators
    void SetSolution( NodeStoreSol<Double>& sol ) {
      sol_ = &sol; }

    //! Set solution class for non-linear integrators
    void SetGridSolution( NodeStoreSol<Double>& gridSol ) {
      gridSol_ = &gridSol; }


    // ! In case of coupled PDEs set sol1
   void SetSolution1(NodeStoreSol<Double> & sol){
     sol1_= &sol;
   }

    // ! In case of coupled PDEs set sol2
    void SetSolution2(NodeStoreSol<Double> & sol){
      sol2_= &sol;
    }
    
    //! Set first time derivative of solution for non-linear integrators
    void SetSolDeriv1( NodeStoreSol<Double>& solDeriv1 ) {
      solDeriv1_ = &solDeriv1; }

    //! Set second time derivative of solution for non-linear integrators
    void SetSolDeriv2( NodeStoreSol<Double>& solDeriv2 ) {
      solDeriv2_ = &solDeriv2; }

    //! sets actual element solution
    virtual void SetActElemSol(Matrix<Double>& disp)
    {EXCEPTION("SetActElemSol not implemented!");};

    //! sets actual element solution
    virtual void SetActElemSol(DenseMatrix & disp)
    {EXCEPTION("SetActElemSol not implemented!");};

    //! sets actual element solution
    virtual void SetActElemSol(Vector<Double>& disp)
    {EXCEPTION("SetActElemSol not implemented!");};

    //! sets actual first time derivative of element solution
    virtual void SetActElemSolDeriv1(Matrix<Double>& disp)
    {EXCEPTION("SetActElemSol not implemented!");};

    //! sets actual first time derivative of element solution
    virtual void SetActElemSolDeriv1(Vector<Double>& disp)
    {EXCEPTION("SetActElemSolDeriv1 not implemented!");};

    //! sets actual 2nd time derivative of element solution
    virtual void SetActElemSolDeriv2(Vector<Double>& disp)
    {EXCEPTION("SetActElemSolDeriv2 not implemented!");};

    //! reads the values y(x) out of the file with name fncName  
    void ReadNlinFunc( std::string fncName, Vector<double> &xval,
                       Vector<Double> &yval ) {
      EXCEPTION("ReadNlinFunc not implemented!");
    };

    //! set the integration point
    void SetIntPoint(Vector<Double> point)
    { intPoint_ = point; isSetIntPoint_ = true;}

    //! unset integration point
    void UnsetIntPoint() { isSetIntPoint_ = false;}

    //!
    void SetMatDataType(Global::ComplexPart & pMatType)
    { matDataType_=pMatType; };

    //! set softening type for forms
    void SetSofteningModel(std::string type) {
      softeningModel_ = type;
    }

    //! set min/max of x,y,z coordinates form where PML starts
    virtual void SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer,
                           const std::string& coordSysId ) {;};
      
#ifndef INTEGLIB
    //! Get reference element and coordinates from element iterator
    virtual void ExtractElemInfo( EntityIterator& it);
    
#endif

    /**
     * Set Timestepping for non linear solvers
     * @param ts_alg the time stepping algorithm
     */
    void setTimeStepping(TimeStepping* const ts_alg) { TS_alg_ = ts_alg;};

    /** MassInt and LaplaceInt have optional constructors where the material descriptors can be set.*/
    const MaterialDescriptor& GetMaterialDescriptor() const { return md_; };

    /** Convenience function if one requires the B mat for the given element.
     * @param ip 1-based :(
     * @param geo_nonlin shall the element coordinates obtained linear or nonlinear? */
    void CalcBMatOnly(Matrix<Double> &bMat, UInt ip, Elem* elem, bool geo_nonlin = false);

    /** @see other CalcBMatOnly() */
    void CalcBMatOnly(Matrix<Double> &bMat, Vector<Double>& intPoint, Elem* elem, bool geo_nonlin = false);

    /** @see other CalcBMatOnly(). This is just a variant in the parameters */
    void CalcBMatOnly(Matrix<Double> &bMat, UInt ip, BaseFE* elem, Matrix<Double> &ptCoord);

  protected:

    /** Gets the factor for dMat to perform the ersatz material ansatz.
     * This is a pseudo density only, see SIMP optimization. 
     * It includes any potential exponent, hence it doesn't need to be the density 
     * This is done e.g. by SIMP or level set. Implemented for mechanic and piezo.
     * The region why the childs have to reimplement this method is, that
     * Domain::GetErsatzMaterial() identifies the proper penaltization via RTTI (typid)
     * @param elem the element
     * @return 1.0 if nothing is to be done or a factor */ 
    virtual Double GetErsatzMaterialFactor(const Elem* elem); 

    /** Some derived classes have a natural tensor e.g. PIEZO_TENSOR, MECH_STIFFNESS_TENSOR. Extend if you need it */
    virtual MaterialType getDMaterialType() { EXCEPTION("not implemented"); }


    /** Computes the discretized differential operator at the given integration point.
     * Possibly all forms overwrite this method. */
    virtual void CalcBMat(Matrix<Double>& bMat, UInt ip, const Matrix<Double>& ptCoord) {
      EXCEPTION("not implemented");
    }

    /**
     * Get Timestepping for non linear solvers
     * @return
     */
    const TimeStepping* getTimeStepping() const
    {
      if (TS_alg_ == NULL)
      {
        EXCEPTION("Time stepping is not set!");
      }
      return TS_alg_;
    };

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

    //! flag indicating if element matrix is symmetric
    bool isSymmetric_;

    //! flag indicating updated lagrangian formulation
    bool coordUpdate_;

    //! flag indicating of value of element matrix/vector depends on
    //! its solution or its derivatives
    bool isSolDependent_;

    //! true, if element matrix is complex
    bool isComplex_;

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

    Global::ComplexPart matDataType_;     //! default = realMaterialParamter, piezoMatType_ = imagMaterialParamter if we consider complex-valued material Paramter;

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

    //! solution vector, sol1 and sol2 are solutions in case id direct couplings
    NodeStoreSol<Double>* sol_, *sol1_, *sol2_;

    //! grid solution vector
    NodeStoreSol<Double>* gridSol_;

    //! first derivative of solution
    NodeStoreSol<Double>* solDeriv1_;

    //! second derivative of solution
    NodeStoreSol<Double>* solDeriv2_;

    /** to be used for bimaterial optimization problems or for a more detailed description of the used material.
     * Makes sense for scalar integrator. For most cases not used
     * @see MassInt() */
    MaterialDescriptor md_;

  private:

    /**
     * Needed for some integrators
     */
    TimeStepping* TS_alg_;
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
    
  };

} //end namespace

#endif // FILE_BASEFORM
