#ifndef COEFFUNCTION_SURFACE_FORCE_BALANCE_HH
#define COEFFUNCTION_SURFACE_FORCE_BALANCE_HH

//#include <boost/tr1/type_traits.hpp>
#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"

// header for resultHandling
#include "DataInOut/ResultHandler.hh"
#include "Domain/Results/ResultFunctor.hh"
#include "PDE/SinglePDE.hh"


namespace CoupledField {

// ============================================================================
//  Coef Function Surface Force Balance
// ============================================================================
//! Coefficient function for calculating a scaled force that will always
//! satisfy \sum \bm{F} \cdot \bm{\eta} = 0.
//!
//! The coefFct gets passed a vector of surface pairs and the corresponding 
//! feFunctions. Furthermore, we pass the coefFct of a result as well as the
//! result functor containing the integrated value. We use the integrated value
//! to perform a scaling of the surface force density in order to ensure
//! \sum \bm{F} \cdot \bm{\eta} = 0. It has to be noted that each surface can
//! only be used once in one unique surface pair, otherwise, an exception is 
//! thrown.
//! Since we perform this calculation for each call and use the current location 
//! of the nodes, this function is applicable for moving mesh applications.
//! 
//! \note This class only works for real-valued data.

class CoefFunctionSurfaceForceBalance :
   public CoefFunction,
   public boost::enable_shared_from_this<CoefFunctionSurfaceForceBalance >
{
public:

  //! Constructor
  CoefFunctionSurfaceForceBalance( PtrCoefFct coef, 
                                    shared_ptr<ResultFunctor> fnct,
                                    shared_ptr<BaseFeFunction> feFct,
                                    StdVector<std::string> surfList1,
                                    StdVector<std::string> surfList2,
                                    StdVector<std::string> volumeList,
                                    StdVector<Vector<Double>> dirVec );

     //! Destructor
  virtual ~CoefFunctionSurfaceForceBalance(){
    ;
  }

  virtual string GetName() const { return "CoefFunctionSurfaceForceBalance"; }


  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat,
                 const LocPointMapped& lpm ) {
    EXCEPTION("Not implemented");
  }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Complex>& coefMat,
                 const LocPointMapped& lpm ) {
    EXCEPTION("Not implemented");
  }

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Double>& coefVec,
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const {
    assert(this->dimType_ == VECTOR );
    return size_;
  }

  //! \copydoc CoefFunction::IsZero
  bool IsZero() const;

  //! \copydoc CoefFunction::ToString
  std::string ToString() const;

  //! This function calculates a scaling factor based on averaging the two coefFunctions projected in a direction
  double CalcScalingFactor(Vector<Double> surfaceRes, Vector<Double> balanceRes, Vector<Double> curDirVec);

  //! This function defines the result functors based on the result handler.
  void DefineResFunctor(void);


protected:

  //! Size of square tensor
  UInt size_;

  //! CoefFunction, result functor, feFunction and surface pair lists
  PtrCoefFct coef_;
  SolutionType fnctResName_;
  shared_ptr<BaseFeFunction> feFct_;
  shared_ptr<ResultFunctor> fnct_;
  StdVector<std::string> volumeList_;
  StdVector<std::string> surfList1_;
  StdVector<std::string> surfList2_;
  StdVector<Vector<Double>> dirVec_;

  StdVector<shared_ptr<ResultFunctor>> fnctList1_;
  StdVector<shared_ptr<ResultFunctor>> fnctList2_;

  //! Set with all regionIdTypes
  std::set<RegionIdType> volRegions_;
  std::set<RegionIdType> surfRegions1_;
  std::set<RegionIdType> surfRegions2_;

  //! scaling factor
  double scalingFactor_;

  ResultHandler *resHandler_;

  //! bool to signal the first function call
  bool firstFuncCall_ = true;
};

}
#endif
