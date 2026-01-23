// =====================================================================================
// 
//       Filename:  rhsBUInt.hh
// 
//    Description:  This class implements the general integrator for RHS integrators of 
//                  the form
//                  \int_K {\cal B} \cdot \vec{U} \ \text{d} K
//                  So we have a quantity U specified by the coefficient function
//                  passed to the constructor and some kind of BOperator 
// 
//        Version:  1.0
//        Created:  11/02/2011 10:09:14 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef FILE_RHSBUINTEGRATOR_
#define FILE_RHSBUINTEGRATOR_

#include "LinearForm.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Forms/Operators/BaseBOperator.hh"

namespace CoupledField{

template< class VEC_DATA_TYPE=Double,
    bool SURFACE = false>
class BUIntegrator : public LinearForm {
public:

  //! Constructor for volume integration
  BUIntegrator(BaseBOperator * bOp,
      VEC_DATA_TYPE factor,
      shared_ptr<CoefFunction > rhsCoef,
      bool coordUpdate = false,
      bool fullEvaluation = true,
      bool extractReal = false,
      const string& id = ""); // to save coil id

  //! Constructor for surface integration
  BUIntegrator(BaseBOperator * bOp,
      VEC_DATA_TYPE factor,
      shared_ptr<CoefFunction > rhsCoef,
      const std::set<RegionIdType>& volRegions,
      bool coordUpdate = false,
      bool fullEvaluation = true,
      bool extractReal = false,
      const string& id = ""); // to save coil id

  //! Copy constructor
  BUIntegrator(const BUIntegrator& right )
  :  LinearForm(right),
     fullEvaluation_(right.fullEvaluation_)
  {
    this->bOperator_ = right.bOperator_->Clone();
    this->factor_ = right.factor_;
    this->rhsCoefs_ = right.rhsCoefs_;
    this->volRegions_ = right.volRegions_;
    this->Bdim_ = right.Bdim_;
    this->id_ = right.id_;
  }

  //! \copydoc LinearForm::Clone
  virtual BUIntegrator* Clone(){
    return new BUIntegrator( *this );
  }

  virtual PtrCoefFct GetCoef() {
    return rhsCoefs_;
  }

  virtual string GetId() {
      return id_;
    }

  virtual ~BUIntegrator(){
    delete bOperator_;
  }

  void CalcElemVector(Vector<VEC_DATA_TYPE> & elemVec,EntityIterator& ent);

  bool IsComplex() const {
    return std::is_same<VEC_DATA_TYPE,Complex>::value;
  }

  virtual void SetFeSpace(shared_ptr<FeSpace> feSpace ){
    this->ptFeSpace_ = feSpace;
    UInt opDim = ptFeSpace_->GetNumDofs();
    intScheme_ = ptFeSpace_->GetIntScheme();
    Bdim_ = opDim;
  }

protected:

  //! Differential operator
  BaseBOperator * bOperator_;

  //! Additional factor for integrator
  VEC_DATA_TYPE factor_;

  //! String to store for example coil id
  string id_;

  //! Flag if full accuracy should be used for coefficient evaluation

  //! This flag denotes, if the coefficient function "u" should be evaluated
  //! at every integration point (true) or if only midpoint evaluation
  //! should be performed (false)
  const bool fullEvaluation_;

  //! Coefficient function "u"
  PtrCoefFct rhsCoefs_;

  //! set containing all volume regions for surface integrators
  std::set<RegionIdType> volRegions_;

  //! dimension of b-operator
  UInt Bdim_ = 0;

  // ---- Work buffers (reused across CalcElemVector calls to avoid allocations) ----
  // These are mutable because CalcElemVector is logically const but needs to reuse buffers
  // Since forms are cloned per-thread in parallel assembly, these are effectively thread-local
  mutable Matrix<Double> work_bMat_;
  mutable Vector<VEC_DATA_TYPE> work_cVec_;
  mutable StdVector<LocPoint> work_intPoints_;
  mutable StdVector<Double> work_weights_;
  mutable Vector<VEC_DATA_TYPE> work_pt1_;
  mutable Vector<VEC_DATA_TYPE> work_pt2_;
  mutable Matrix<Double> work_JacT_;
  mutable Matrix<Double> work_TF_;
  mutable Matrix<Double> work_TFinv_;

};

}
//Include template definition file
#include "BUInt.cc"
#endif
