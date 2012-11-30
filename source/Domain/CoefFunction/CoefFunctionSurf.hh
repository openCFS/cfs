#ifndef FILE_COEFFUNCTION_SURF_HH
#define FILE_COEFFUNCTION_SURF_HH


#include "CoefFunction.hh"

namespace CoupledField {

//! Evaluates a coefficient function on a surface 

//! This class represents coefficient functions, which are defined just on a
//! surface. Thus, it internally holds pointers to coefficient functions for
//! each region. If queried for a scalar / vector / tensor entry of a given
//! surface element, it looks for the "correct" volume coefficient function
//! and returns it. In addition, it can also perform a normal-mapping of
//! the coefficient function.
template<typename T>
class CoefFunctionSurf : public CoefFunction {
public:

  //! Constructor
  CoefFunctionSurf( bool mapNormal );

  //! Destructor
  virtual ~CoefFunctionSurf();

  //! Pass volume coefficients
  void SetVolumeCoefs( std::map<RegionIdType, PtrCoefFct> coefs ); 

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<T>& coefMat, 
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<T>& coefVec, 
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(T& coefScalar, 
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVecSize
  UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;

  //! \copydoc CoefFunction::ToString
  std::string ToString() const;

private:

  //! Map with CoefFunctions for each region
  std::map<RegionIdType, PtrCoefFct> coefs_;
  
  //! Set with all regionIdTypes
  std::set<RegionIdType> regions_;

  //! Flag, if normal mapping should be performed
  bool mapNormal_;

  //! Uniform vector size
  UInt vecSize_;
};


} // end of namespace
#endif
