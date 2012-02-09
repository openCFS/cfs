// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MECHSTRESSSTRAIN_04
#define FILE_MECHSTRESSSTRAIN_04

#include <complex>

#include "Forms/linElastInt.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"

namespace CoupledField {
class BaseMaterial;
class BaseResult;
class EntityIterator;
template <class TYPE> class Vector;
}  // namespace CoupledField

namespace CoupledField
{
  class MechPDE;

  //! class for calculation of mechanical stresses and strains
  template <class TYPE>
  class MechStressStrain : public linElastInt
  {
  public:

    /// Constructor
    MechStressStrain(BaseMaterial* matData, SubTensorType type);

    /// Destructor
    virtual ~MechStressStrain();  

    /** Calculates the results for stress/strain resulted stuff in MechPDE and PiezoCoupling
    * @param st either MECH_STRESS or MECH_STRAIN or VON_MISES_STRESS or VON_MISES_STRAIN
    * @param elec set to NULL when used for mech stress stuff only!
    * @param density divide by element volume? only for mechPDE */
    void CalcStressStrainResult(MechPDE* mech, shared_ptr<BaseResult> res, SolutionType st, bool density = false);

    /** calculates linear mech or piezoelectric stress */
    virtual void CalcStressVec(Vector<TYPE>& stressVec, UInt ip, EntityIterator& ent);

    // calculate linear part of Green Lagrangian Strain tensor
    void CalcStrainVec(Vector<TYPE>& strainVec, UInt ip, EntityIterator& ent);
    
    /// in stress calculations, the actual displacement of the element is needed
    /*!
      \param disp (input) Matrix with displacement d of all nodes of actual element
      \f[ \left( \begin{array}{ccc} 
      d_{x1} &  d_{x2} &  d_{x3} \\
      d_{y1} &  d_{y2} &  d_{y3} \\
      \end{array}\right) \f]         
    */
    virtual void SetActElemSol(Matrix<TYPE>& disp) { elemDisp_ = disp; }

  protected:  
  
    /// displacement of all nodes of actual element
    Matrix<TYPE> elemDisp_;
  };
  

} //end namespace

#endif // FILE_XXX


