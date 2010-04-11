// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
#ifndef FILE_CFS_H1_ELEMENTS_HI_HH
#define FILE_CFS_H1_ELEMENTS_HI_HH

#include "H1Elems.hh"

namespace CoupledField {

  //! Base class for hierarchical H1-conforming finite elements of arbitrary order
  
  //! This is the base class for all H1 elements of arbitrary order using 
  //! hierarchical polynomials. Currently we use the integrated Legendre polynomials,
  //! but in general we can plug in any arbitrary orthogonal polynomial (e.h. Jacobi,
  //! Gegenbauer, etc.) 
  class FeH1Hi : public FeH1 {
  public:

    //! Constructor
    FeH1Hi();

    //! Destructor
    virtual ~FeH1Hi();

    
    //! Calculate the shape functions / local derivatives for all integration points
    
    //! Not that this method is just a stub, as the current implementation of the
    //! hierarchical finite elements requires the global orientation of the element.
    //! Therewith we can not pre-calculate the shape functions.
    void SetFunctionsAtIp(const StdVector<LocPoint>& iPoints) {};
    
    
    //! Get number of shape functions for a given type (NODE/EDGE/FACE/ELEM)
    void GetNumFncs( StdVector<UInt>& numFcns,
                     EntityType fctEntityType,
                     UInt dof = 1 );

    //! This holds only for line,quad and hex,
    //! if other types are available the has to be reimplemented!
    //! Get the permutation Vector for a given Face or Edge
    //! e.g. If asked for a face, the element will check the flags
    //! of this face and return a vector of size NumberOfFncs on the Face
    //! holding the correct ordering 
    /*!
         \param fncPermutation (output) The Permuation Vector 
         \param ptElem (input) pointer to Grid Element to get grip of flags 
         \param fctEntityType (input) The Entity type, Node/Edge/Face where the 
                                      nodes are located at
         \param entNumber (input) The local entity number 
     */
    virtual void GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                      const Elem* ptElem,
                                      EntityType fctEntityType,
                                      UInt entNumber);

    
  protected:

    //@{
    //! Temporary helper functions until we use the final
    //! polynomial object 
    
    //! New version: Calculate all Legendre polynomials up to order p
    //void EvalPolynom( UInt p, Vector<Double>& vals );
    
    //! Evaluate polynom and its derivative using Honer's algorithm
    void EvalPolynom( Double& value, Double& deriv,
                 const UInt order, const Double* coeff,
                 const Double xVal );

    //! Coefficients of 1D Legendre coefficients up to order 8
    static Double lCoeff_[9][10];
    //@}
    
    //! Polynomial order per entity
    std::map<EntityType,StdVector<UInt> > entityOrder_;
  
    //! Number of shape functions per entity
    std::map<EntityType,StdVector<UInt> > entityFncs_;
  };

 
  
  //! H1 conforming hierarchical higher order line element
  class FeH1HiLine : public FeH1Hi {

  public:

    //! Constructor
    FeH1HiLine();

    //! Destructor
    virtual ~FeH1HiLine();
    
    //! Set isotropic polynomial order
    void SetIsoOrder(UInt order);

  protected:

    //! @see BaseFE::GetShFnc 
    void GetShFnc( Vector<Double>& shape, const LocPoint& lp,
                   const Elem* ptElem, UInt comp = 1 ); 

    //! @see BaseFE::GetDerivShFnc 
    void GetDerivShFnc( Matrix<Double> & deriv, const LocPoint& lp,
                        const Elem * elem,  UInt comp = 1);

  };

  //! H1 conforming hierarchical higher order quadrilateral element
  class  FeH1HiQuad : public FeH1Hi {

  public:

    //! Constructor
    FeH1HiQuad();
    

    //! Destructor
    virtual ~FeH1HiQuad();

    //! Set isotropic polynomial order
    void SetIsoOrder(UInt order);
    

    //! @see BaseFE::GetShFnc 
    void GetShFnc( Vector<Double>& shape, const LocPoint& lp,
                   const Elem * elem, UInt comp = 1 ); 

    //! @see BaseFE::GetDerivShFnc 
    void GetDerivShFnc( Matrix<Double> & deriv, const LocPoint& lp,
                        const Elem * elem,  UInt comp = 1);
  protected:
    };

  
  //! H1 conforming hierarchical higher order hexahedral element
   class  FeH1HiHex : public FeH1Hi {

   public:

     //! Constructor
     FeH1HiHex();
     
     //! Destructor
     virtual ~FeH1HiHex();
     
     //! Set isotropic polynomial order
     void SetIsoOrder(UInt order);

     protected:
     //! @see BaseFE::GetShFnc 
     void GetShFnc( Vector<Double>& shape, const LocPoint& lp,
                    const Elem * elem, UInt comp = 1 ); 

     //! @see BaseFE::GetDerivShFnc 
     void GetDerivShFnc( Matrix<Double> & deriv, const LocPoint& lp,
                         const Elem * elem,  UInt comp = 1);

   };

} // namespace CoupledField

#endif
