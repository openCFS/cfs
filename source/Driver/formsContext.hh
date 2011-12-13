// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_FORMS_CONTEXT_HH
#define CFS_FORMS_CONTEXT_HH


#include <string>

#include "General/defs.hh"
#include "General/environment.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {
class DampLayer;
template <class TYPE> class StdVector;
template <class TYPE> class Vector;
}  // namespace CoupledField

namespace CoupledField
{

  //! Forward class declarations
  class BaseForm;
  class EntityIterator;
  class EntityList;
  class EqnMap;
  class LinearForm;
  class SinglePDE;
  struct ResultInfo;



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
    	destMat_ = destMat;
      }

    //! Defines a secondary destination for the element matrix
    void SetSecDestMat( FEMatrixType aSecMat,
                        std::string aSecMatFac ); 

    //! initialize object for damping layer
    void SetDampLayer(std::string& dampingTypeFnc,
		      Vector<Double>& mPoint,
		      Double& dampFactor,
		      Double& dampFactorMax,
		      Double& startRadius,
		      Double& endRadius);

    //! Returns matrix type of the secondary matrix
    FEMatrixType GetSecDestMat() const { return secDestMat_; }

    //! Returns the factor the secondary matrix gets multiplied with (string representation)
    std::string GetSecMatFac() const;

    //! Returns the current value of the secondary matrix factor (evaluated, number representation)
    Double EvalSecMatFac() const;
    
    //! Returns the integrator
    BaseForm * GetIntegrator() {return integrator_; };

    //! Return entry type of matrix (real/imag part)
    Global::ComplexPart GetEntryType() {return entryType_;};

    //! Set entrytype for matrix (real/imag part)
    void SetEntryType( Global::ComplexPart &pEntryType ){
      entryType_ = pEntryType;};

    // ======================================================
    //  MAPPING METHODS
    // ======================================================

    //! Map equations for bilinear form for combination of two given entities
    virtual void MapEqns( EntityIterator& it1,
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
    void SetResults( shared_ptr<ResultInfo> result1,
                     shared_ptr<ResultInfo> result2,
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

    //! Returns information about first result info
    shared_ptr<ResultInfo> GetFirstResultInfo() { return result1_; }

    //! Returns information about second result info
    shared_ptr<ResultInfo> GetSecondResultInfo() { return result2_; }

    //get the pointer to damping layer object!
    DampLayer* getPtDamplayer() {
      return dampingLayer_;}

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

    //! create human readable debug output */
    std::string ToString();

  protected:

    //! Pointer to bilinearform
    BaseForm * integrator_;

    //! Destination matrix type
    FEMatrixType destMat_;

    //! Secondary destination matrix
    FEMatrixType secDestMat_;

    //! Handle for secondary matrix factor
    MathParser::HandleType secMatFacHandle_;
    
    //! Pointer to math parser instance
    MathParser* mathParser_;

    //! Entry type of matrix (real/imag part)
    Global::ComplexPart entryType_;

    //! for damping layer
    DampLayer* dampingLayer_;

    // Flag indicating assembling of the integrator
    // in the counterpart of the pde location
    bool setCounterPart_;

    // Flag indicating negating the entries of the element matrix
    bool negateEntries_;

    // ======================================================
    //  MAPPING DATA
    // ======================================================

    //! Pointer to first pde
    SinglePDE * ptPde1_;

    //! Pointer to second pde
    SinglePDE * ptPde2_;

    //! Pointer to first result type
    shared_ptr<ResultInfo> result1_;

    //! Pointer to second result type
    shared_ptr<ResultInfo> result2_;

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
    void SetResult( shared_ptr<ResultInfo> result,
                    shared_ptr<EntityList> list );

    //! Return first set of current entities
    shared_ptr<EntityList> GetEntities() { return ent_; }

    //! returns pointer to first pde
    SinglePDE * GetPde () { return ptPde_; }
    
    //! Returns information about result info
    shared_ptr<ResultInfo> GetResultInfo() { return result_; }

    std::string ToString() const;

  protected:

    //! Pointer to bilinearform
    LinearForm * integrator_;

    // ======================================================
    //  MAPPING DATA
    // ======================================================

    //! Pointer to pde
    SinglePDE * ptPde_;

    //! Pointer to result type
    shared_ptr<ResultInfo> result_;

    //! Pointer to entity list
    shared_ptr<EntityList> ent_;

    //! Pointer to equation map
    shared_ptr<EqnMap> map_;

  };

  //! Specialized context for non-conforming interfaces
  class NcBiLinFormContext : public BiLinFormContext  {

  public:

    //! Constructor
    //! \param biLinForm pointer to the bilinearform to be wrapped
    //! \param destMat destination Matrix (STIFFNESS, MASS, ...) of the
    //!                bilinearform
    NcBiLinFormContext( BaseForm* biLinForm,
                        FEMatrixType destMat );

    //! Destructor
    virtual ~NcBiLinFormContext();

    // ======================================================
    //  MAPPING METHODS
    // ======================================================

    //! Map equations for bilinear form for combination of two given entities
    void MapEqns( EntityIterator& it1,
                  EntityIterator& it2,
                  StdVector<Integer>& eqnVec1,
                  StdVector<Integer>& eqnVec2,
                  PdeIdType& id1, PdeIdType& id2 );
  };

} // end of namespace

#endif
