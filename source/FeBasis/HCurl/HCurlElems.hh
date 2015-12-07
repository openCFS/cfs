// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
#ifndef FILE_CFS_HCURL_ELEMENTS_HH
#define FILE_CFS_HCURL_ELEMENTS_HH

#include "FeBasis/BaseFE.hh"

namespace CoupledField {

  //! Base class for H(Curl)-conforming reference elements
  
  //! The HCurl element represents all elements with H(Curl) 
  //! conforming shape function
  class FeHCurl : public BaseFE {
  public:

    //! Constructor
    FeHCurl();

    //! Destructor
    virtual ~FeHCurl();

    //! Get (vectorial) shape functions
    
    //! Return HCurl shape functions 
    virtual void GetShFnc( Matrix<Double>& shape, 
                           const LocPointMapped& lp,
                           const Elem* elem, UInt comp = 1 );

    //! Return global curl of shape functions
    virtual void GetCurlShFnc( Matrix<Double>& curl, 
                               const LocPointMapped& lp,
                               const Elem* elem, UInt comp = 1 );

  protected:
    
    //! Calculate local basis functions (based on local gradient)
    virtual void CalcLocShFnc( Matrix<Double>& curl, 
                               const LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 ) {
      
    }
    
    //! Calculate local curl of basis function
    virtual void CalcLocCurlShFnc( Matrix<Double>& curl,
                                   const LocPointMapped& lp,
                                   const Elem* elem, UInt comp = 1 ) {
      
    }

    
  private:
  };
  
} // namespace CoupledField

#endif
