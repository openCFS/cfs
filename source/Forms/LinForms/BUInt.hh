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
#include <boost/tr1/type_traits.hpp>
#include "Domain/CoefFunction/CoefFunction.hh"


namespace CoupledField{


template< class B_OP,
class VEC_DATA_TYPE=Double,
bool SURFACE = false>
class BUIntegrator : public LinearForm{
public:

  //! Constructor for volume integration
  BUIntegrator(VEC_DATA_TYPE factor,
               shared_ptr<CoefFunction > rhsCoef,
               bool coordUpdate = false );

  //! Constructor for surface integration
  BUIntegrator(VEC_DATA_TYPE factor,
               shared_ptr<CoefFunction > rhsCoef,
               const std::set<RegionIdType>& volRegions,
               bool coordUpdate = false );

  virtual ~BUIntegrator(){

  }

  void CalcElemVector(Vector<VEC_DATA_TYPE> & elemVec,EntityIterator& ent);

  bool IsComplex(){
    return std::tr1::is_same<VEC_DATA_TYPE,Complex>::value;
  }

  virtual void SetFeSpace(shared_ptr<FeSpace> feSpace ){
    this->ptFeSpace_ = feSpace;
    UInt opDim = feSpace->GetFeFunction()->GetResultInfo()->dofNames.GetSize();
    Bdim_ = opDim;
  }

  //! \copydoc LinearForm::IsSolDependent
  bool IsSolDependent() {
    return rhsCoefs_->GetDependency() == CoefFunction::SOLUTION;
  }

protected:
  B_OP operator_;

  VEC_DATA_TYPE factor_;

  shared_ptr<CoefFunction> rhsCoefs_;

  //! set containing all volume regions for surface integrators
  std::set<RegionIdType> volRegions_;

  //! dimension of b-operator
  UInt Bdim_;


};
}
//Include template definition file
#include "BUInt.cc"
#endif
