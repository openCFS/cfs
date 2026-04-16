#ifndef OPTIMIZATION_DESIGN_LOCALELEMENTCACHE_HH_
#define OPTIMIZATION_DESIGN_LOCALELEMENTCACHE_HH_

#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "General/Enum.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Matrix.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/ElemMapping/EntityLists.hh"

using std::complex;
using std::string;

namespace CoupledField
{

class BiLinearForm;

/** This class holds precomputed local element matrices in cache. This is to returned e.g. by
 * BDBInt::CalcElementMatrix(). It speeds up Assemble::AssembleMatrices*() and ErsatzMaterial::CalcUKU()
 * via OptimzationMaterial. It is closely related to OptimizationMaterial and handles regular and non-regular
 * meshed.
 * SIMP is easy, more critical are bi-material, multi-material, ParamMat (for the derivatives) and Bloch with the
 * wave vector dependence. Therefore the LocalElementCache can be disabled for not-yet implemented features.
 *
 * Trigger in the .xml via <designSpace local_element_cache="false"/> before design */
class LocalElementCache
{
public:
  LocalElementCache(DesignSpace* space);

  /** initialize the org data. FIXME: add context! */
  void InitOrg();

  /** initialize bimaterial for org and for shadow */
  void InitShadow(DesignSpace::DesignRegion* dr);

  /** Initialize mechanic material derivatives by checking if MECH_11 is in space */
  void InitMechMatDeriv(StdVector<RegionIdType>& reg);

  /** set a material derivative */
  void SetMatDeriv(const string& integrator, RegionIdType reg, DesignElement::Type dir, int design_id);

  /** deletes cached data */
  void ClearMatDeriv(const string& integrator, RegionIdType reg, DesignElement::Type dir, int design_id);

  /** This is essentially a replicate of CoefFunctionOpt::State which we cannot easily use
   * due to cyclic dependencies */
  typedef enum { NONE = -1, ORG = 0, SHADOW, DIRECTION } Type;

  bool IsActive() const { return active_; }

  /** set the element matrix if active. If active but not found, an exception is thrown */
  template<class T>
  bool CachedOrgElement(Matrix<T>& mat_out, BiLinearForm* form, const Elem* elem);

  /** variant for gradient evaluation without copying the stored data. */
  template<class T>
  const Matrix<T>& CachedElement(const string& integrator, Type type, const Elem* elem, DesignElement::Type dir = DesignElement::NO_DERIVATIVE, PtrCoefFct shadow_coef = PtrCoefFct(), int design_id = -1);

  bool HasCachedData(const string& integrator, Type type, DesignElement::Type dir = DesignElement::NO_DERIVATIVE, PtrCoefFct shadow_coef = PtrCoefFct(), int design_id = -1) {
    return active_ ? GetFormData(integrator, type, dir, shadow_coef, design_id) != NULL : false;
  }

  template<class T>
  const Matrix<T>& CachedOrgElement(const string& integrator, const Elem* elem) {
    return CachedElement<T>(integrator, ORG, elem);
  }

  template<class T>
  const Matrix<T>& CachedShadowElement(const string& integrator, const Elem* elem, PtrCoefFct shadow_coef) {
    return CachedElement<T>(integrator, SHADOW, elem, DesignElement::NO_DERIVATIVE, shadow_coef);
  }

  template<class T>
  const Matrix<T>& CachedMatDerivElement(const string& integrator, const Elem* elem, DesignElement::Type dir, int design_id = -1) {
    return CachedElement<T>(integrator, DIRECTION, elem, dir, PtrCoefFct(), design_id);
  }

  void ToInfo(PtrParamNode info);

private:

  /** Initialize (additional) material derivatives. This is always the tensor coefficients, even if the
   * design variables are different (e.g. Young's modulus)
   * @param dir only MECH_11, MECH_12, ... material derivative local element matrices can be cached!
   * @return false if dir was not feasible */
  bool InitMatDeriv(const string& integrator, RegionIdType reg, DesignElement::Type dir);

  Type GetType(bool lower_bimat,  DesignElement::Type direction) const;

  /** collects data for a form over all regions.  */
  struct FormData
  {
    /** the constructor which resizes the elements */
    void Init(const BiLinearForm* form, bool structured);

    /** Checks if context id and wave number do match with the current context */
    bool CheckContext() const
    {
      Context* curr = Optimization::context;
      if((int) curr->context_idx != this->contex)
        return false;

      assert((curr->DoBloch() && curr->num_bloch_wave_vectors > 0) || (!curr->DoBloch() && curr->num_bloch_wave_vectors == 0));
      if(curr->num_bloch_wave_vectors > 0)
        if(curr->GetEigenFrequencyDriver()->GetCurrentWaveVector() != this->wave)
          return false;

      return true;
    }

    /** index within data_real_ or data_cplx_ */
    int    idx = -1;

    /** the name of the bilinear form */
    string integrator;

    Type   type = NONE;

    /** the tensor coefficient MECH_11, ... for param mat derivatie */
    DesignElement::Type dir = DesignElement::NO_DERIVATIVE;

    /** the shadow coef for bimat and multimaterial */
    PtrCoefFct shadow;

    /** The wave vector for the bloch case. Taken from current context */
    Vector<double> wave;

    /** The context number. Always the current context is relevant, no option */
    int contex = -1;

    /** element data, indexed by elemNum. Only filled for regions which are not regular or when unregular is enforced  */
    StdVector<Matrix<double> > elem_real;
    StdVector<Matrix<complex<double> > >  elem_cplx;

    /** region data, indexed by RegionTypeId. If not set, the data shall be in elem_.
     * see LocalElementCache::regular_ */
    StdVector<Matrix<double> > region_real;
    StdVector<Matrix<complex<double> > >  region_cplx;

    /** which data to fill */
    bool isComplex = false;

    /** which iteration the data belongs to */
    int designID = -1;

    /** Helper to hide complex.
     * @param region write to region_ or elem_. region_ must be written only once! */
    void CalcElementMatrix(BiLinearForm* form, EntityIterator& entity, bool region);

    /** to be called by assert. Contains a lot of asserts */
    bool Validate(const BiLinearForm* form, Type type, DesignElement::Type dir);

    /** computes the memory in MBytes */
    double CalcMemory();

    /** for debugging purpos */
    std::string ToString();
  };

  /** searches the data and conditionally creates it */
  FormData& GetFormData(const BiLinearForm* form, Type type, DesignElement::Type dir, bool create = false);

  /** @return NULL if the data could not be found */
  LocalElementCache::FormData* GetFormData(const string& integrator, Type type,
      DesignElement::Type dir = DesignElement::NO_DERIVATIVE,
      PtrCoefFct shadow_coef = PtrCoefFct(),
      int design_id = -1)
  {
    // try to avoid as many string comparisons as possible!
    int res_idx = -1;
    switch(type)
    {
    case ORG:
      assert(dir == DesignElement::NO_DERIVATIVE);
      assert(!shadow_coef);
      assert(design_id == -1);
      for(unsigned int i = 0; res_idx == -1 && i < org_end_; i++)
        if(data_[i].CheckContext() && data_[i].integrator == integrator)
          res_idx = i;
      break;
    case SHADOW:
      assert(dir == DesignElement::NO_DERIVATIVE);
      assert(shadow_coef);
      assert(design_id == -1);
      for(unsigned int i = org_end_; res_idx == -1 && i < shadow_end_; i++)
        if(data_[i].CheckContext() && data_[i].integrator == integrator && data_[i].shadow == shadow_coef)
          res_idx = i;
      break;
    case DIRECTION:
      assert(dir != DesignElement::NO_DERIVATIVE);
      assert(!shadow_coef);
      for(unsigned int i= shadow_end_; res_idx == -1 && i < dir_end_; i++)
        if(data_[i].CheckContext() && data_[i].dir == dir && data_[i].integrator == integrator && data_[i].designID == design_id)
          res_idx = i;
      break;
    default:
      assert(false); // not yet implmented
    }
    if(res_idx == -1)
      return NULL;

    assert(data_[res_idx].integrator == integrator);
    assert(data_[res_idx].type == type);
    assert(data_[res_idx].contex == (int) Optimization::context->context_idx);
    assert(data_[res_idx].designID == design_id);
    return &data_[res_idx];
  }

  /** Helper for GetFormData() */
  FormData* AppendFormData(const BiLinearForm* form, Type type, DesignElement::Type dir = DesignElement::NO_DERIVATIVE, PtrCoefFct shadow_coef = PtrCoefFct(), int design_id = -1);

  /** fills either all elements for the region or all elements of the region.
   * Switches the Form to the type encoded in data!
   * @return false if either wrong form (SingleEntryInt, LBM) or nor region data */
  bool FillFormData(FormData& data, BiLinearForm* form, RegionIdType reg);

  /** checks of the form is not not org. Meant for debug assert() */
  bool CheckFormState(BiLinearForm* form, Type type);

  /** For InitShadow() and InitMatDeriv() */
  void Init(const string& integrator, RegionIdType reg, Type type, DesignElement::Type dir, PtrCoefFct coef);

  StdVector<FormData> data_;

  /** speedup access to data_. last org element (not included!) */
  unsigned int org_end_ = 0;
  /** shadow_start is org_end (not included!)*/
  unsigned int shadow_end_ = 0;
  /** index of last direction entry (not included!). dir_start is shadow_end.
   * Note that we can only cache MECH_11, MECH_12, ... material derivatives! */
  unsigned int dir_end_ = 0;

  /** if not active we fill the data and do no caching. */
  bool active_ = false;

  /** we have only a global regular. Either all regions are regular and "enforce_unstructured" is not set or
   * everything is non-regular. Has impact on where to store in FormData */
  bool regular_ = false;

  DesignSpace* space_ = NULL;

}; // end of LocalElementCache

} // end of namespace

#endif /* OPTIMIZATION_DESIGN_LOCALELEMENTCACHE_HH_ */
