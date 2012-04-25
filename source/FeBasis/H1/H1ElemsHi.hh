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

    // ----------------------------------------------------------------------
    //  Polynomial Order
    // ----------------------------------------------------------------------
    
    //! Set isotropic polynomial order
    void SetIsoOrder( UInt order );
    
    //! Set anisotropic / directional dependent order
    
    //! Set the polynomial order direction dependent, i.e. different for each 
    //! direction. This can be used to prefer certain directions, i.e. for 
    //! thin structures or anisotropic field distribution.
    //! \param order Vector containing the polynomial order in each local 
    //!              direction (xi,eta,zeta).
    //! \param ptElem Pointer to geometric element
    //!
    //! \note Only simplex-based elements (quad, wedge, hex) support the 
    //!       prescription of anisotropic polynomial order.
    
    void SetAnisoOrder( const StdVector<UInt>& order,
                        const Elem* ptElem );

    //! Return set of basis functions, which exceeds a given polynomial order
    
    //! This method returns a set of indices of basis functions, which have
    //! a higher polynomial degree than the (spatially dependent) order
    //! provided as argument. The indices correspond to the positions of the
    //! basis functions / derivatives obtainable e.g. by the method
    //! FeH1::GetShFnc.
    //! \param nodes Contains all indices of basis functions, which have higher
    //!              polynomial degree than the one provided in the order
    //!              vector
    //! \param order Vector containing the polynomial order in each local 
    //!              direction (xi,eta,zeta) (see FeH1Hi::SetAnisoOrder)
    //! \param ptElem Pointer to geometric element
    void GetNodesExceedingOrder( std::set<UInt>& nodes, 
                                 const StdVector<UInt>& order,
                                 const Elem* ptElem );
    
    //! Set polynomial order for an edge
    
    //! Sets the order of a given edge number (0-based) to a given order. 
    void SetEdgeOrder( UInt edgeNum, UInt order  );
    
    //! Set polynomial order for a face
    
    //! Sets the order of a given face number (0-based) to a given order.
    //! The order is specified in terms of face-local directions.
    void SetFaceOrder( UInt faceNum, const boost::array<UInt,2>& order );
    
    //! Set polynomial order for element interior
    void SetInteriorOrder( const boost::array<UInt,3>& order ); 
    
    //! Return edge order
    const StdVector<UInt>& GetEdgeOrder( ) const {
      return orderEdge_; 
    }
    
    //! Return face order 
    const StdVector<boost::array<UInt,2> >& GetFaceOrder( ) const {
      return orderFace_;
    }
    
    //! Return interior order
    const boost::array<UInt,3>& GetInnerOrder( ) const {
      return orderInner_;
    }
    
    //! Calculate the shape functions / local derivatives for all integration points
    
    //! Not that this method is just a stub, as the current implementation of the
    //! hierarchical finite elements requires the global orientation of the element.
    //! Therewith we can not pre-calculate the shape functions.
    void SetFunctionsAtIp(const StdVector<LocPoint>& iPoints) {};
    
    //! Get total number of functions
    UInt GetNumFncs();
    
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
         \param fncPermutation (output) The Permutation Vector 
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

  protected:
    
    //! Calculate number of unknowns
    
    //! This method calculates the number of unknowns functions
    //! for real tensor-product elements, i.e. line, quadrilateral
    //! and hexahedral elements.
    virtual void CalcNumUnknowns();
    
    //! Return matrix with polynomial degrees for each shape function
    
    //! This method returns a matrix, which contains for every basis function
    //! and every local direction the polynomial degree. This can be used e.g.
    //! to determine the spatial order of a given node.
    virtual void GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                      const Elem* ptElem ) {
      EXCEPTION( "Not implemented" );
    }
    
    //! Flag if re-calculation of number of unknowns is needed
    
    //! After changing the order of the element, a re-calculation of 
    //! the number of unknowns (actNumFncs_) is necessary.
    bool updateUnknowns_;
    
    //! Number of shape functions per entity
    std::map<EntityType,StdVector<UInt> > entityFncs_;
    
    //! Isotropic order. 0 if anisotropic
    UInt isoOrder_;
    
    //! Maximum polynomial degree of element
    UInt maxOrder_;
    
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
    
    
    //! @copydoc FeH1Hi::GetPolyOrderOfNodes
    virtual void GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                      const Elem* ptElem );
    
    //! Templatized version of calculation for shape function
    template<typename T_SCAL, typename T_VEC>
    void _CalcShFnc( const T_SCAL x,  
                     const Elem * elem,
                     T_VEC& ret );
    
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

    //! @copydoc FeH1Hi::GetPolyOrderOfNodes
    virtual void GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                      const Elem* ptElem );

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

    //! @copydoc FeH1Hi::GetPolyOrderOfNodes
    virtual void GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                      const Elem* ptElem );
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

    //! @copydoc FeH1Hi::GetPolyOrderOfNodes
    virtual void GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                      const Elem* ptElem );
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

     //! @copydoc FeH1Hi::GetPolyOrderOfNodes
     virtual void GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                       const Elem* ptElem );
     
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
