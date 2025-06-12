#ifndef CFS_FORMS_CONTEXT_HH
#define CFS_FORMS_CONTEXT_HH


#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

namespace CoupledField
{

  //! Forward class declarations
  class EntityList;
  class SinglePDE;
  class EqnMap;
  class EntityIterator;
  struct ResultInfo;
  class FeSpace;
  class BaseFeFunction;
  class LinearForm;
  class MathParser;




  //! Base class for wrapping a (bi)linearform

  //! This class is used as a wrapper around a bilinear-form (assembled onto
  //! matrices) or a linear-form (assembled onto the right hand side).
  //! As the class BaseForm itself serves only as a calculation scheme,
  //! this class (and its derivatives) store the following information:
  //! - destination matrix (bilinearform)
  //! - information about non-linearity (re-assembling flag, dependent values)
  //! - associated equation numbers for each element matrix


  //! Class for wrapping a bilinearform
  class BiLinFormContext  {

  public:



    //! Constructor
    //! \param biLinForm pointer to the bilinearform to be wrapped
    //! \param destMat destination Matrix (STIFFNESS, MASS, ...) of the
    //!                bilinearform
    BiLinFormContext( BiLinearForm* biLinForm, FEMatrixType destMat );

    //! Destructor
    virtual ~BiLinFormContext();

    // ======================================================
    //  MATRIX ASSEMBLING INFORMATION
    // ======================================================

    //! Returns true if a non-linear dependency (geometry,
    //! solution) is present for the wrapped bilinearform
    bool IsNonLin();

    //! bilinearform is part of Newton tangential matrix
    bool IsNewtonBiLinearForm() {
      return integrator_->IsNewtonBiLinearForm();
    }

    static void SetEnums();

    //! Get destination matrix
    FEMatrixType GetDestMat() const { return destMat_; }

    //! Set destination matrix
    void SetDestMat(FEMatrixType destMat) { destMat_ = destMat; }

    //! Defines a secondary destination for the element matrix
    void SetSecDestMat( FEMatrixType aSecMat,
                        std::string aSecMatFac ); 

    //! Returns matrix type of the secondary matrix
    FEMatrixType GetSecDestMat() const { return secDestMat_; }

    //! Returns the factor the secondary matrix gets multiplied with (string representation)
    std::string GetSecMatFac() const;

    //! Returns the current value of the secondary matrix factor (evaluated, number representation)
    Double EvalSecMatFac() const;
    
    //! Returns the integrator
    BiLinearForm * GetIntegrator() {return integrator_; };

    //! Return entry type of matrix (real/imag part)
    Global::ComplexPart GetEntryType() {return entryType_;};

    //! Set entrytype for matrix (real/imag part)
    void SetEntryType( Global::ComplexPart &pEntryType ){
      entryType_ = pEntryType;};

    //! Set eqn evaluation to volume for A operator
    void SetUseVolEqnA( bool useVolEqn ){
      useVolEqnA_ = useVolEqn;
    }

    //! Set eqn evaluation to volume for B operator
    void SetUseVolEqnB( bool useVolEqn ){
      useVolEqnB_ = useVolEqn;
    }

    // ======================================================
    //  MAPPING METHODS
    // ======================================================

    //! Map equations for bilinear form for combination of two given entities
    virtual void MapEqns( EntityIterator& it1,
                          EntityIterator& it2,
                          StdVector<Integer>& eqnVec1,
                          StdVector<Integer>& eqnVec2,
                          FeFctIdType& id1, FeFctIdType& id2 );

    // ======================================================
    // ENTITIES / RESULTS
    // ======================================================

    //! Set the result types and entities the bilinearform is working on
    void SetEntities( shared_ptr<EntityList> list1,
                      shared_ptr<EntityList> list2 );

    //! Note: In the case of coupling two different PDEs via Nitsche formulation, fct1 should be the unknown
    //! from Master side and fct2 should be the unknown from slave side
    void SetFeFunctions( shared_ptr<BaseFeFunction> fct1, 
                         shared_ptr<BaseFeFunction> fct2); 

    //! Return first set of current entities
    shared_ptr<EntityList> GetFirstEntities() { return ent1_; }

    //! Return second set of current entities
    shared_ptr<EntityList> GetSecondEntities() { return ent2_; }

    //! Returns pointer to first pde
    SinglePDE* GetFirstPde ()  { return ptPde1_; }

    //! Returns pointer to second pde
    SinglePDE* GetSecondPde () { return ptPde2_; }

    //! Returns information about first result info
    shared_ptr<ResultInfo> GetFirstResultInfo() { return result1_; }

    //! Returns information about second result info
    shared_ptr<ResultInfo> GetSecondResultInfo() { return result2_; }

    //! Returns information about first feFunction
    weak_ptr<BaseFeFunction> GetFirstFeFunction() { return feFct1_; }

    //! Returns information about second feFunction
    weak_ptr<BaseFeFunction> GetSecondFeFunction() { return feFct2_; }

    //! Set function for SetCounterPart
    void SetCounterPart( bool setCounterPart ) {
      setCounterPart_ = setCounterPart;
    }

    //! Set function for setNegate
    void SetNegate(bool setNegate) {
      negateEntries_ = setNegate;
    }

    //! Check, if element matrix has to be assembled to
    //! upper and lower part of global matrix
    bool IsSetCounterPart() const {
    	return setCounterPart_;
    }

    //! to check if we need to negate the integrator
    bool IsSetNegate() const {
      return negateEntries_;
    }

    //! Check if this is a diagonal bilinear form

    //! This method returns, if the bilinear form gets assembled on a
    //! diagonal block, i.e. it has the same primary and secondary
    //! fundtion Id, as well as entity lists.
    virtual bool isDiagonal() {
      return ((result1_ == result2_) && (ent1_ == ent2_));
    }


    //! create human readable debug output */
    std::string ToString();

  protected:

    //! Pointer to bilinearform
    BiLinearForm * integrator_;

    //! Destination matrix type
    FEMatrixType destMat_;

    //! Secondary destination matrix
    FEMatrixType secDestMat_;

    //! Handle for secondary matrix factor
    unsigned int secMatFacHandle_;
    
    //! Pointer to math parser instance
    MathParser* mathParser_;

    //! Entry type of matrix (real/imag part)
    Global::ComplexPart entryType_;

    // Flag indicating assembling of the integrator
    // in the counterpart of the pde location
    bool setCounterPart_;

    // Flag indicating negating the entries of the element matrix
    bool negateEntries_;

    // Flag to indicate if the number of functions shall be aquired from the volume or surface element
    // This is needed for e.g. a gradient evaluated at a surface, since we perform the evaluation at the
    // integration points of the surface, but we need all DOFs of the element for the calculation of the gradient
    // We can set this for A- and B operator seperately
    bool useVolEqnA_;

    bool useVolEqnB_;


    // ======================================================
    //  MAPPING DATA
    // ======================================================

    //! Pointer to first pde
    SinglePDE* ptPde1_;

    //! Pointer to second pde
    SinglePDE* ptPde2_;

    //! Pointer to first result type
    shared_ptr<ResultInfo> result1_;

    //! Pointer to second result type
    shared_ptr<ResultInfo> result2_;

    //! Pointer to first entity list
    shared_ptr<EntityList> ent1_;

    //! Pointer to second entity list
    shared_ptr<EntityList> ent2_;

    //! Pointer to first FeFunction
    weak_ptr<BaseFeFunction> feFct1_;

    //! Pointer to second FeFunction
    weak_ptr<BaseFeFunction> feFct2_;

  };



  // -------------------------------------------------------------------------

  //! Class for wrapping a linearform (assembled to the right hand side)
  class LinearFormContext {

  public:

    //! Constructor
    LinearFormContext( LinearForm* linearForm );

    //! Destructor
    virtual ~LinearFormContext();

    //! Return integrator
    LinearForm* GetIntegrator() { return integrator_; }

    //! Returns true if a non-linear dependency (geometry,
    //! solution) is present for the wrapped linearform
    bool IsNonLin();

    // ======================================================
    //  MAPPING METHODS
    // ======================================================

    /** Map equations of linearform for a given entitylist (e.g. NODE_LIST)
      * @param eqnVec returned equations. Might be empty, e.g. if nodes  are not part of simulation */
    void MapEqns( EntityIterator& it,
                  StdVector<Integer>& eqnVec,
                  FeFctIdType& id );

    // ======================================================
    // ENTITIES / RESULTS
    // ======================================================

    //! Set pointer to pde where the form is defined from
    void SetPtPde(SinglePDE* ptPde );

    //! Set pointer to pde where the form is defined from
    void SetFeFunction(shared_ptr<BaseFeFunction>  fct );

    //! Set the result types and entities the linearform is defined on
    void SetEntities( shared_ptr<EntityList> list );

    //! Return first set of current entities
    shared_ptr<EntityList> GetEntities() { return ent_; }

    //! returns pointer to first pde
    SinglePDE* GetPde () { return ptPde_; }
    
    //! Returns information about result info
    shared_ptr<ResultInfo> GetResultInfo() { return result_; }

    std::string ToString() const;

  protected:

    //! Pointer to linearform
    LinearForm* integrator_;

    // ======================================================
    //  MAPPING DATA
    // ======================================================

    //! Pointer to pde
    SinglePDE* ptPde_;

    //! Pointer to result type
    shared_ptr<ResultInfo> result_;

    //! Pointer to entity list
    shared_ptr<EntityList> ent_;

    //! Pointer to equation map
    //shared_ptr<EqnMap> map_;

    //! Pointer to FeFunction
    weak_ptr<BaseFeFunction> feFct_;

  }; // class LinearFormContext
  
  //! Specialized context for non-conforming interfaces (NcSurfABInt)
  class NcBiLinFormContext : public BiLinFormContext {
    
  public:
    
    //! Constructor
    //! \param biLinForm pointer to the bilinearform to be wrapped
    //! \param destMat destination Matrix (STIFFNESS, MASS, ...) of the
    //!                bilinearform
    NcBiLinFormContext(BiLinearForm *biLinForm, FEMatrixType destMat)
        : BiLinFormContext( biLinForm, destMat ) {

      isMoving_ = false;
    };

    //! Destructor
    virtual ~NcBiLinFormContext(){};

    // ======================================================
    //  MAPPING METHODS
    // ======================================================

    //! Map equations for bilinear form for combination of two given entities
    virtual void MapEqns( EntityIterator& it1,
                          EntityIterator& it2,
                          StdVector<Integer>& eqnVec1,
                          StdVector<Integer>& eqnVec2,
                          FeFctIdType& id1, FeFctIdType& id2 );

    //! Returns all equations of this context
    virtual void GetEqns( StdVector<Integer>& eqnVec1,
                             StdVector<Integer>& eqnVec2,
                             FeFctIdType& id1, FeFctIdType& id2 ) const;
    
    //! set the moving flag
    virtual void SetMotion(bool moving){
      isMoving_ = moving;
    }

    //! obtain moving flag
    virtual bool GetMotion(){
      return isMoving_;
    }

    //! Does the Context needs a fully populated matrix
    virtual bool NeedsFullMatrix(){
      return isMoving_;
    }



  protected:

    //flag indicating if we are moving
    bool isMoving_;

  }; // class NcBiLinFormContext

  class SurfaceBiLinFormContext : public NcBiLinFormContext {
   public:


     SurfaceBiLinFormContext( BiLinearForm* biLinForm, FEMatrixType destMat, BiLinearForm::CouplingDirection currentDirection);

     //! Destructor
     virtual ~SurfaceBiLinFormContext();

     virtual void MapEqns( EntityIterator& it1,
                               EntityIterator& it2,
                               StdVector<Integer>& eqnVec1,
                               StdVector<Integer>& eqnVec2,
                               FeFctIdType& id1, FeFctIdType& id2 );

     //! Returns all equations of this context
     virtual void GetEqns( StdVector<Integer>& eqnVec1,
                              StdVector<Integer>& eqnVec2,
                              FeFctIdType& id1, FeFctIdType& id2 ) const;


     //! Check if this is a diagonal bilinear form

     //! This method returns, if the bilinear form gets assembled on a
     //! diagonal block, i.e. it has the same primary and secondary
     //! fundtion Id, as well as entity lists.
     virtual bool isDiagonal() {
       return ((result1_ == result2_) && (ent1_ == ent2_) &&
               (currentDirection_== BiLinearForm::PRIM_PRIM ||
                currentDirection_== BiLinearForm::SEC_SEC ) );
     }

     //! Does the Context needs a fully populated matrix
     virtual bool NeedsFullMatrix(){
       bool direction = (currentDirection_ == BiLinearForm::PRIM_SEC ||
                         currentDirection_ == BiLinearForm::SEC_PRIM);
       return (direction && isMoving_) ;
     }
   protected:
     BiLinearForm::CouplingDirection currentDirection_;
   };
} // end of namespace

#endif
