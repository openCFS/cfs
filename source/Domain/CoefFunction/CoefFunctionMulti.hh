// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEF_FUNCTION_MULTI_HH
#define COEF_FUNCTION_MULTI_HH

#include "CoefFunction.hh"

namespace CoupledField  {

//! This class is composed of other CoefFunctions for different regions

//! This class can be used to compose a CoefFunction spatially from other
//! feFunctions, e.g. for mixing analytical and interpolated results.
class CoefFunctionMulti : public CoefFunction {
public:
  //! Constructor
  
  //! Generic constructor.
  //! \param zeroEmptyRegions If true, the class returns a zero-CoefFunction
  //!                         for regions without assigned CoefFunction.
  //!                         Otherwise an exception is thrown.
  CoefFunctionMulti(bool zeroEmptyRegions = true);
  
  //! Destructor 
  virtual ~CoefFunctionMulti();
  
  //! Set coefficient function for a region
  void AddRegion( RegionIdType region, PtrCoefFct coef );

  //! Get all region definitions
  std::map<RegionIdType,PtrCoefFct > GetRegionCoefs() {
    return regionCoefs_;
  }
  
  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods
  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<Complex>& coefMat,
                         const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<Complex>& coeVec,
                         const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(Complex& coef,
                         const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<Double>& coefMat,
                         const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<Double>& coefVec,
                         const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(Double& coef,
                         const LocPointMapped& lpm );
  
  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;
  //@}

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;
  
private:
  
  //! Return coefficient function for a given region
  inline PtrCoefFct GetRegionCoef( RegionIdType region ) {
    std::map<RegionIdType,PtrCoefFct >::const_iterator it = 
        regionCoefs_.find(region);
    if(it == regionCoefs_.end()){
      if ( !zeroEmptyRegions_ )
        EXCEPTION("Region " << region << " is not contained in functor");
      return zeroCeof_;
    }
    return it->second;
  }
  
  //! Map storing coefFunction of the analytical regions
  std::map<RegionIdType,PtrCoefFct > regionCoefs_;
  
  //! "Zero" coefficient function, returned for empty regions
  PtrCoefFct zeroCeof_;  
  
  //! Flag, if zero coefficient function is return for non-set regions
  bool zeroEmptyRegions_;
};

} // end of namespace
#endif
