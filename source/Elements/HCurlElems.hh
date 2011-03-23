// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
#ifndef FILE_CFS_HCURL_ELEMENTS_HH
#define FILE_CFS_HCURL_ELEMENTS_HH

#include "basefe.hh"

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
    
    //! Calculate vectorial shape functions
    virtual void GetShFnc( Matrix<Double>& S, const LocPoint& lp,
                           const Elem* ptElem,  UInt comp = 1 )  = 0;
    
    //! Calculate Global Curl of shape functions
    void GetCurlShFcn( Matrix<Double>& s, const LocPoint& lp, 
                       const Elem* ptElem,  UInt comp = 1 );
     
    
    //! Return global derivative of shape functions
    void GetGlobDerivShFnc( Matrix<Double>& deriv, LocPointMapped& lp,
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
    
    //! Polynomial shape function
    //IntLegendrePol poly_;
    
  private:
  };

  
  // ========================================================================
  //  H1 Fe Lagrangian Elements of lowest order (1st / 2nd)
  // ========================================================================
  
} // namespace CoupledField

#endif
