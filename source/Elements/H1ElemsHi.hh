// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode c
#ifndef FILE_CFS_H1_ELEMENTS_HI_HH
#define FILE_CFS_H1_ELEMENTS_HI_HH

#include "H1Elems.hh"
#include <boost/array.hpp>


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


    //! Set isotropic polynomial order
    void SetIsoOrder(UInt order);
    
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

    //! \copydoc BaseFE::GetIsoOrder
    virtual UInt GetIsoOrder() const;
    
    //! \copydoc BaseFE::GetMaxOrder
    virtual UInt GetMaxOrder() const;

    //! \copydoc BaseFE::GetMaxOrderLocDir
    virtual void GetMaxOrderLocDir(StdVector<UInt>& order );

  protected:
    
    //! Calculate number of unknowns
    
    //! This method calculates the number of unknowns functions
    //! for real tensor-product elements, i.e. line, quadrilateral
    //! and hexahedral elements.
    virtual void CalcNumUnknowns();
    
    //! Flag if re-calculation of number of unknowns is needed
    
    //! After changing the order of the element, a re-calculation of 
    //! the unumber of unknowns (actNumFncs_) is necessary.
    bool updateUnknowns_;
    
    //! Number of shape functions per entity
    std::map<EntityType,StdVector<UInt> > entityFncs_;
    
    //! Isotropic order. 0 if anisotropic
    UInt isoOrder_;
    
    // ========================================================================
    // DEFINITION OF (ANISOTROPIC) ORDER
    // ========================================================================
    //@{ \name Definition of (anisotropic) polynomial order
    
    //! Polynomial order of edges (#edges x 1 local direction)
    StdVector<UInt> orderEdge_;
    
    //! Polynomial order of faces (#faces x 2 local directions)
    StdVector<boost::array<UInt,2> > orderFace_;
    
    //! Polynomial order of inner (1 x 3 local directions)
    boost::array<UInt,3> orderInner_;
    //@}
    
    
    // ========================================================================
    // DEPRECATED / OLD METHODS
    // ========================================================================
    //@{ \name Deprecated section
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
  };

 
  
  // ======
  //  LINE
  // =======
  class FeH1HiLine : public FeH1Hi {

  public:

    //! Constructor
    FeH1HiLine();

    //! Destructor
    virtual ~FeH1HiLine();
    
  protected:

    //! @copydoc FeH1::CalcShFnc
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point,
                    const Elem* ptElem,
                    UInt comp = 1 );

    //! @copydoc FeH1::CalcLocDerivShFnc
    void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 );
  };
  
  // ==========
  //  TRIANGLE
  // ==========
  //! H1 conforming hierarchical higher order triangular element
  class  FeH1HiTria : public FeH1Hi {

  public:

    //! Constructor
    FeH1HiTria();

    //! Destructor
    virtual ~FeH1HiTria();

  protected:

    //! @copydoc FeH1::CalcShFnc
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point,
                    const Elem* ptElem,
                    UInt comp = 1 );

    //! @copydoc FeH1::CalcLocDerivShFnc
    void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 );


    //! Templatized version of calculation for shape function
    template<typename T_SCAL, typename T_VEC>
    void _CalcShFnc( const T_SCAL x, const T_SCAL y, 
                     const Elem * elem,
                     T_VEC& ret );

    //! @copydoc FeH1::GetNumFncs
    void CalcNumUnknowns();
  };

  // ===============
  //  QUADRILATERAL 
  // ===============
  
  //! H1 conforming hierarchical higher order quadrilateral element
  class  FeH1HiQuad : public FeH1Hi {

  public:

    //! Constructor
    FeH1HiQuad();

    //! Destructor
    virtual ~FeH1HiQuad();

  protected:

    //! @copydoc FeH1::CalcShFnc
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point,
                    const Elem* ptElem,
                    UInt comp = 1 );

    //! @copydoc FeH1::CalcLocDerivShFnc
    void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 );


    //! Templatized version of calculation for shape function
    template<typename T_SCAL, typename T_VEC>
    void _CalcShFnc( const T_SCAL x, const T_SCAL y, 
                     const Elem * elem,
                     T_VEC& ret );
  };

  // ============
  //  HEXAHEDRAL 
  // ============
  //! H1 conforming hierarchical higher order hexahedral element
  class  FeH1HiHex : public FeH1Hi {

  public:

    //! Constructor
    FeH1HiHex();

    //! Destructor
    virtual ~FeH1HiHex();

  protected:

    //! @copydoc FeH1::CalcShFnc
    void CalcShFnc( Vector<Double>& shape,
                    const Vector<Double>& point,
                    const Elem* ptElem,
                    UInt comp = 1 );

    //! @copydoc FeH1::CalcLocDerivShFnc
    void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 );

    //! Templatized version of calculation for shape function
    template<typename T_SCAL, typename T_VEC>
    void _CalcShFnc( const T_SCAL x, const T_SCAL y, const T_SCAL z,
                     const Elem * elem,
                     T_VEC& ret );

  };
  
  // =======
  //  WEDGE  
  // =======
  
  //! H1 conforming hierarchical higher order wedge element
   class  FeH1HiWedge : public FeH1Hi {

   public:

     //! Constructor
     FeH1HiWedge();

     //! Destructor
     virtual ~FeH1HiWedge();

   protected:

     //! @copydoc FeH1::CalcShFnc
     void CalcShFnc( Vector<Double>& shape,
                     const Vector<Double>& point,
                     const Elem* ptElem,
                     UInt comp = 1 );

     //! @copydoc FeH1::CalcLocDerivShFnc
     void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                             const Vector<Double>& point,
                             const Elem* ptElem,
                             UInt comp = 1 );

     //! Templatized version of calculation for shape function
     template<typename T_SCAL, typename T_VEC>
     void _CalcShFnc( const T_SCAL x, const T_SCAL y, const T_SCAL z,
                      const Elem * elem,
                      T_VEC& ret );
     
     //! @copydoc FeH1::GetNumFncs
     void CalcNumUnknowns();
   };

} // namespace CoupledField

#endif
