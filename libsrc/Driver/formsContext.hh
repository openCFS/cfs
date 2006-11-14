#ifndef CFS_FORMS_CONTEXT_HH
#define CFS_FORMS_CONTEXT_HH


#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Utils/dampLayer.hh"

namespace CoupledField
{
  
  //! Forward class declarations
  class BaseForm;
  class LinearForm;
  class EntityList;
  class SinglePDE;
  class EqnMap;
  class EntityIterator;
  class ResultDof;

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
    BiLinFormContext( BaseForm* biLinForm, FEMatrixType destMat );

    //! Destructor
    virtual ~BiLinFormContext();

    // ======================================================
    //  MATRIX ASSEMBLING INFORMATION
    // ======================================================    
    
    //! Returns true if a non-linear dependency (geometry,
    //! solution) is present for the wrapped bilinearform
    bool IsNonLin();

    //! Get destination matrix
    FEMatrixType GetDestMat() const { return destMat_; }
    
    //! Set destination matrix
    void SetDestMat(FEMatrixType destMat) {
      destMat_ = destMat; }
    
    //! Defines a secondary destination for the element matrix
    void SetSecDestMat( FEMatrixType aSecMat, 
                          Double aSecMatFac ) {
      secDestMat_ = aSecMat;
      secMatFac_ = aSecMatFac; }

    //! initialize object for damping layer
    void SetDampLayer(std::string& dampingTypeFnc, 
		      Vector<Double>& mPoint, 
		      Double& dampFactor, 
		      Double& dampFactorMax, 
		      Double& startRadius, 
		      Double& endRadius);
    
    //! Returns matrix type of the secondary matrix 
    FEMatrixType GetSecDestMat() const { return secDestMat_; } 

    //! Returns the factor the secondary matrix gets multiplied with
    Double GetSecMatFac() const {return secMatFac_;} 

    //! Returns the integrator
    BaseForm * GetIntegrator() {return integrator_; };

    //! Return type of used materialtype (real/complex)
    DataType GetMatDataType(){return matDataType_;};

    //! Set type of used materialtype (real/complex)
    void SetMatDataType(DataType &pMatType){
      matDataType_ = pMatType;};

    // ======================================================
    //  MAPPING METHODS
    // ======================================================

    //! Map equations for bilinear form for combination of two given entities
    void MapEqns( EntityIterator& it1, 
                  EntityIterator& it2,
                  StdVector<Integer>& eqnVec1, 
                  StdVector<Integer>& eqnVec2,
                  PdeIdType& id1, PdeIdType& id2 );
    
    // ======================================================
    // ENTITIES / RESULTS
    // ======================================================

    //! Set pointer to PDE(s) where the form is derived from 
    void SetPtPdes(SinglePDE * aPDE1, SinglePDE * aPDE2 );
    
    //! Set the result types and entities the bilinearform is working on
    void SetResults( shared_ptr<ResultDof> result1, 
                     shared_ptr<ResultDof> result2,
                     shared_ptr<EntityList> list1, 
                     shared_ptr<EntityList> list2 );
    
    //! Return first set of current entities
    shared_ptr<EntityList> GetFirstEntities() { return ent1_; }

    //! Return second set of current entities
    shared_ptr<EntityList> GetSecondEntities() { return ent2_; }

    //! Returns pointer to first pde
    SinglePDE * GetFirstPde ()  { return ptPde1_; }

    //! Returns pointer to second pde
    SinglePDE * GetSecondPde () { return ptPde2_; }

    //! Set function for SetCounterPart
    void SetCounterPart( bool setCounterPart ) 
    { setCounterPart_ = setCounterPart; }
    
    //! Check, if element matrix has to be assembled to 
    //! upper and lower part of global matrix
    bool IsSetCounterPart() const {return setCounterPart_; }

    //get the pointe rto damping layer object!
    DampLayer* getPtDamplayer() {
      return dampingLayer_;};

  protected:

    //! Pointer to bilinearform
    BaseForm * integrator_;

    //! Destination matrix type
    FEMatrixType destMat_;

    //! Secondary destination matrix
    FEMatrixType secDestMat_;

    //! Secondary matrix factor
    Double secMatFac_;

    //! Flag indicating assembling of counterpart
    bool setCounterPart_;
    
    //! Type of used materialData
    DataType matDataType_;

    //! for damping layer
    DampLayer* dampingLayer_;

    // ======================================================
    //  MAPPING DATA
    // ======================================================
    
    //! Pointer to first pde
    SinglePDE * ptPde1_;

    //! Pointer to second pde
    SinglePDE * ptPde2_;

    //! Pointer to first result type
    shared_ptr<ResultDof> result1_;

    //! Pointer to second result type
    shared_ptr<ResultDof> result2_;

    //! Pointer to first entity list
    shared_ptr<EntityList> ent1_;

    //! Pointer to second entity list
    shared_ptr<EntityList> ent2_;

    //! Pointer to first equation map
    shared_ptr<EqnMap> map1_;
    
    //! Pointer to second equation map
    shared_ptr<EqnMap> map2_;


  };

  // -------------------------------------------------------------------------

  //! Class for wrapping a linearform (assembled to the right hand side)
  class LinearFormContext {

  public:

    //! Constructor
    LinearFormContext( LinearForm* linearForm,
                       const std::string& dynamics = std::string() );

    //! Destructor
    virtual ~LinearFormContext();

    //! Return integrator
    LinearForm* GetIntegrator() { return integrator_; }

    //! Get dynamics
    const std::string& GetDynamics() { return dynamics_; }

    //! Returns true if a non-linear dependency (geometry,
    //! solution) is present for the wrapped linearform
    bool IsNonLin();
    
    // ======================================================
    //  MAPPING METHODS
    // ======================================================

    //! Map equations of linearform for a given entitylist
    void MapEqns( EntityIterator& it,
                  StdVector<Integer>& eqnVec,
                  PdeIdType& id );

    // ======================================================
    // ENTITIES / RESULTS
    // ======================================================

    //! Set pointer to pde where the form is defined from
    void SetPtPde(SinglePDE * ptPde );
    
    //! Set the result types and entities the linearform is defined on
    void SetResult( shared_ptr<ResultDof> result,
                    shared_ptr<EntityList> list );
    
    //! Return first set of current entities
    shared_ptr<EntityList> GetEntities() { return ent_; }

    //! returns pointer to first pde
    SinglePDE * GetPde () { return ptPde_; }

  protected:

    //! Pointer to bilinearform
    LinearForm * integrator_;
    
    // ======================================================
    //  LOAD DATA
    // ======================================================
    
    //! Load dynamics file
    std::string dynamics_;
    
    // ======================================================
    //  MAPPING DATA
    // ======================================================
    
    //! Pointer to pde
    SinglePDE * ptPde_;

    //! Pointer to result type
    shared_ptr<ResultDof> result_;

    //! Pointer to entity list
    shared_ptr<EntityList> ent_;

    //! Pointer to equation map
    shared_ptr<EqnMap> map_;
    
  };




  
  
 //  //! Class for performing the assembly of the system matrices
//   class Assemble {


//   public:

//   protected:

//   private:

//   };



 


} // end of namespace

#endif
