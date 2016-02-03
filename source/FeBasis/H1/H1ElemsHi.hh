#ifndef FILE_CFS_H1_ELEMENTS_HI_HH
#define FILE_CFS_H1_ELEMENTS_HI_HH

#include "H1Elems.hh"
#include "FeBasis/FeHi.hh"

namespace CoupledField {
  //! Base class for hierarchical H1-conforming finite elements of arbitrary order
  
  //! This is the base class for all H1 elements of arbitrary order using 
  //! hierarchical polynomials. Currently we use the integrated Legendre polynomials,
  //! but in general we can plug in any arbitrary orthogonal polynomial (e.h. Jacobi,
  //! Gegenbauer, etc.) 
  class FeH1Hi : public FeH1, public FeHi {
  public:

    //! Constructor
    FeH1Hi(Elem::FEType feType  );

    //! Copy constructor
    FeH1Hi(const FeH1Hi & other);

    //! Destructor
    virtual ~FeH1Hi();

    //! Deep Copy
    virtual FeH1Hi* Clone() = 0;

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

    //! \copydoc BaseFE::IsIsotropic 
    virtual bool IsIsotropic() {
      return isIsotropic_;
    }
    
    //! \copydoc BaseFE::GetIsoOrder
    virtual UInt GetIsoOrder() const;
    
    //! \copydoc BaseFE::GetMaxOrder
    virtual UInt GetMaxOrder() const;
    
    //! \copydoc BaseFE::GetAnisoOrder
    virtual void GetAnisoOrder(StdVector<UInt>& order ) const;
    
    //! Compare two element for equality (= same shape and approximation);
    bool operator==( const FeH1Hi& comp) const;

  protected:
    
    //! Calculate number of unknowns

    //! This method calculates the number of unknowns functions
    //! for real tensor-product elements, i.e. line, quadrilateral
    //! and hexahedral elements.
    virtual void CalcNumUnknowns();
    
  };

 
  
  // ======
  //  LINE
  // =======
  class FeH1HiLine : public FeH1Hi {

  public:

    //! Constructor
    FeH1HiLine();

    FeH1HiLine(const FeH1HiLine & other)
      : FeH1Hi(other){

    }

    //! Destructor
    virtual ~FeH1HiLine();
    
    //! Dolly method
    virtual FeH1HiLine* Clone(){
      return new FeH1HiLine(*this);
    }

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

    //! Copy constructor
    FeH1HiTria(const FeH1HiTria & other)
      : FeH1Hi(other){
    }

    //! Destructor
    virtual ~FeH1HiTria();

    //! Dolly method
    virtual FeH1HiTria* Clone(){
      return new FeH1HiTria(*this);
    }
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

    //! Copy constructor
    FeH1HiQuad(const FeH1HiQuad & other)
      : FeH1Hi(other){
    }

    //! Destructor
    virtual ~FeH1HiQuad();

    //! Dolly method
    virtual FeH1HiQuad* Clone(){
      return new FeH1HiQuad(*this);
    }

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

    //! Copy constructor
    FeH1HiHex(const FeH1HiHex & other)
      : FeH1Hi(other){
    }

    //! Destructor
    virtual ~FeH1HiHex();

    //! Dolly method
    virtual FeH1HiHex* Clone(){
      return new FeH1HiHex(*this);
    }

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

     //! Copy constructor
     FeH1HiWedge(const FeH1HiWedge & other)
       : FeH1Hi(other){
     }

     //! Destructor
     virtual ~FeH1HiWedge();

     //! Dolly method
     virtual FeH1HiWedge* Clone(){
       return new FeH1HiWedge(*this);
     }

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
