// =====================================================================================
// 
//       Filename:  rhsBUInt.hh
// 
//    Description:  This class implements the general integrator for RHS integrators of 
//                  the form
//                  \int_K {\cal B}  \[d\] \cdot \vec{U} \ \text{d} K
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

#ifndef FILE_RHS_BDUINTEGRATOR_
#define FILE_RHS_BDUINTEGRATOR_

#include "LinearForm.hh"
#include <boost/tr1/type_traits.hpp>
#include "Domain/CoefFunction/CoefFunction.hh"


namespace CoupledField{


template< class B_OP,
class VEC_DATA_TYPE=Double,
bool SURFACE = false>
class BDUIntegrator : public LinearForm{
public:

  //! Constructor for volume integration
  BDUIntegrator(VEC_DATA_TYPE factor,
                shared_ptr<CoefFunction > rhsCoef,
                shared_ptr<CoefFunction> dCoef,
                bool coordUpdate = false,
				bool extractReal = false);

  //! Constructor for surface integration
  BDUIntegrator(VEC_DATA_TYPE factor,
                shared_ptr<CoefFunction > rhsCoef,
                shared_ptr<CoefFunction> dCoef,
                const std::set<RegionIdType>& volRegions,
                bool coordUpdate = false,
				bool extractReal = false);

  //! Copy constructor
  BDUIntegrator(const BDUIntegrator& right )
    : LinearForm(right) {
    this->operator_ = right.operator_;
    this->factor_ = right.factor_;
    this->dCoef_ = right.dCoef_;
    this->rhsCoefs_ = right.rhsCoefs_;
    this->volRegions_ = right.volRegions_;
    this->Bdim_ = right.Bdim_;
  }

  //! \copydoc LinearForm::Clone
  virtual BDUIntegrator* Clone(){
    return new BDUIntegrator( *this );
  }

  virtual ~BDUIntegrator(){

  }

  void CalcElemVector(Vector<VEC_DATA_TYPE> & elemVec,EntityIterator& ent);

  bool IsComplex() const {
    return std::tr1::is_same<VEC_DATA_TYPE,Complex>::value;
  }

  virtual void SetFeSpace(shared_ptr<FeSpace> feSpace ){
    this->ptFeSpace_ = feSpace;
    UInt opDim = ptFeSpace_->GetNumDofs();
    intScheme_ = ptFeSpace_->GetIntScheme();
    Bdim_ = opDim;
  }

protected:
  B_OP operator_;

  VEC_DATA_TYPE factor_;

  //! Coefficient function for excitation
  PtrCoefFct rhsCoefs_;
  
  //! Coefficient function for material parameter (d-tensor)
  PtrCoefFct dCoef_;

  //! set containing all volume regions for surface integrators
  std::set<RegionIdType> volRegions_;

  //! dimension of b-operator
  UInt Bdim_;


};
}
//Include template definition file
#include "BDUInt.cc"
#endif
