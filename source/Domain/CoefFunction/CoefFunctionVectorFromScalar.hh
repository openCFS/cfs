#ifndef COEFFUNCTION_VECTOR_FROM_SCALAR_HH
#define COEFFUNCTION_VECTOR_FROM_SCALAR_HH

//#include <boost/tr1/type_traits.hpp>
#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"


namespace CoupledField {

//! Provide a diagonal tensorial coefficient function defined by scalar ones

//! This class generates a tensorial coefficient function, where the diagonal
//! elements are defined in terms of scalar-valued coefficient functions.
class CoefFunctionVectorFromScalar :
   public CoefFunction,
   public boost::enable_shared_from_this<CoefFunctionVectorFromScalar >
{
public:

  //! Constructor
  CoefFunctionVectorFromScalar(shared_ptr<BaseFeFunction> feFct, 
                                StdVector<std::string> surfList, 
                                StdVector<PtrCoefFct>& scalVals, 
                                StdVector<UInt> dofInd, 
                                UInt size);

     //! Destructor
  virtual ~CoefFunctionVectorFromScalar(){
    ;
  }

  virtual string GetName() const { return "CoefFunctionVectorFromScalar"; }

  //! \copydoc CoefFunction::GetVector
  void GetVector( Vector<Double>& coefVec, const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const {
    assert(this->dimType_ == VECTOR );
    return size_;
  }

  //! \copydoc CoefFunction::IsZero
  bool IsZero() const;

  //! \copydoc CoefFunction::ToString
  std::string ToString() const;


protected:

  //! list of surfaces considered for evaluation
  StdVector<std::string> surfList_;

  //! Set with all regionIdTypes
  std::set<RegionIdType> surfRegions_;
  
  //! Size of the vectors
  UInt size_;

  //! Vector with coefFunction entries
  StdVector<PtrCoefFct> scalVals_;

  //! DoF name
  StdVector<UInt> dofInd_;

  //! Zero function
  PtrCoefFct constZero_;
  
  //! Mathparser
  MathParser* mp_ = nullptr;
};

}
#endif
