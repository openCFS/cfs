#ifndef FILE_COEFFUNCTION_SURF_HH
#define FILE_COEFFUNCTION_SURF_HH


#include "CoefFunction.hh"

namespace CoupledField {

// forward class declaration
class ResultInfo;

//! Evaluates a coefficient function on a surface 

//! This class represents coefficient functions, which are defined just on a
//! surface. Thus, it internally holds pointers to coefficient functions for
//! each region. If queried for a scalar / vector / tensor entry of a given
//! surface element, it looks for the "correct" volume coefficient function
//! and returns it. In addition, it can also perform a normal-mapping of
//! the coefficient function.
class CoefFunctionSurf : public CoefFunction {
public:

  //! Constructor
  
  //! Constructor for the class
  //! \param mapNormal If true, only the normal component w.r.t. to the
  //!                  surface element is taken into account.
  //! \param surfInfo Result info object for surface result
  CoefFunctionSurf( bool mapNormal, shared_ptr<ResultInfo> surfInfo =  shared_ptr<ResultInfo>());

  //! Destructor
  virtual ~CoefFunctionSurf();

  //! Set single volume coefficient function
  void AddVolumeCoef( RegionIdType, PtrCoefFct );
  
  //! Pass volume coefficients
  void SetVolumeCoefs( std::map<RegionIdType, PtrCoefFct> coefs ); 

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, 
                 const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetTensor
    void GetTensor(Matrix<Complex>& coefMat, 
                   const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Double>& coefVec, 
                 const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Complex>& coefVec, 
                 const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetScalar
    void GetScalar(Complex& coefScalar, 
                   const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVecSize
  UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;

  //! \copydoc CoefFunction::ToString
  std::string ToString() const;

  // ===========================
  //  NORMAL MAPPING OPERATIONS
  // ===========================
  //@{ \name Functions for normal mapping
  
  //! Mapping operation for vectors
  template<typename TYPE>
  static void MapVecNormal( TYPE& ret, const Vector<TYPE>& vec, 
                            const Vector<Double>& normal );
  
  //! Mapping operation for tensors in Voigt notation
  template<typename TYPE>
  static void MapTensorNormal( Vector<TYPE>& ret, const Vector<TYPE>& tensor,
                               const Vector<Double>& normal );
  
  //@}
  
private:

  //! Map with CoefFunctions for each region
  std::map<RegionIdType, PtrCoefFct> coefs_;
  
  //! Set with all regionIdTypes
  std::set<RegionIdType> regions_;

  //! Flag, if normal mapping should be performed
  bool mapNormal_;
};


} // end of namespace
#endif
