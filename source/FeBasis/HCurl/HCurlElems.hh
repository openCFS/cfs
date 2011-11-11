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
    virtual void GetShFnc( Matrix<Double>& shape, LocPointMapped& lp,
                           const Elem* elem, UInt comp = 1 );

    //! Return global curl of shape functions
    void GetCurlShFnc( Matrix<Double>& curl, LocPointMapped& lp,
                       const Elem* elem, UInt comp = 1 );

    //! Set the isotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) The desired order of the element
    virtual void SetIsoOrder(UInt order){
      // after debugging phase, this EXCEPTION can become a warning
      EXCEPTION("Trying to set the ISOTROPIC order of an element which does not support this opeeration.\
                The element will remain unchanged!");
    }

    //! Set the Anisotropic order of the Element. This methods gets overwritten 
    //! by the child classes to calculate the number of functions according to
    //! the given order
    //! \param order (input) vector of element orders for each space direction 
    virtual void SetAnisoOrder(StdVector<UInt> order){
      // after debugging phase, this EXCEPTION can become a warning
      EXCEPTION("Trying to set the ANISOTROPIC order of an element which does not support this operation.\
                The element will remain unchanged!");
    }

  protected:
    
    //! Calculate local basis functions (based on local gradient)
    virtual void CalcLocShFnc( Matrix<Double>& curl, LocPointMapped& lp,
                             const Elem* elem, UInt comp = 1 ) = 0;
    
    //! Calculate local curl of basis function
    virtual void CalcLocCurlShFnc( Matrix<Double>& curl, LocPointMapped& lp,
                                   const Elem* elem, UInt comp = 1 ) = 0;

    
  private:
  };
  
} // namespace CoupledField

#endif
