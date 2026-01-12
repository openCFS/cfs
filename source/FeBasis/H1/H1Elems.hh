#ifndef FILE_CFS_H1_ELEMENTS_HH
#define FILE_CFS_H1_ELEMENTS_HH

#include "FeBasis/BaseFE.hh"
#include <unordered_map>
#include "def_use_openmp.hh"
#ifdef USE_OPENMP
  #include <omp.h>
#endif

namespace CoupledField {

  //! Base class for H1-conforming reference elements
  
  //! The H1 element represents all elements with H1 conforming shape function
  class FeH1 : public BaseFE {
  public:

    //! Enum defining the completenesss of the basis functions
    
    //! This enum defines (especially for higher order polynomials), if the
    //! set of polynomials spans a complete set up to a given order p (i.e.
    //! every single monomial occurring has maximum order p -> tensorial 
    //! approach / Lagrange approach) or if the only the sum of every 
    //! combination spans the polynomial degree p (Serendipity / trunk type)
    //! This has an influence on the possible choice of integration method,
    //! as 
    //@{
    typedef enum {
      SERENDIPITY_TYPE,      // Only side nodes (also referred as trunk space)
      TENSOR_TYPE            // Tensor approach (also referred as Lagrange space)
    } PolyCompleteType;
    static Enum<PolyCompleteType> polyCompleteEnum_;
    //@}
    
    
    //! Constructor
    FeH1();

    //! Copy Constructor
    FeH1(const FeH1 & other);

    //! Deep Copy
    virtual FeH1* Clone() = 0;

    //! Destructor
    virtual ~FeH1();

    // ======================================================================
    //  Standard Shape Functions
    // ======================================================================
    
    //@{ \name Standard bilinear shape function 
    //! Get value of all shape fnc at local poiont lp
    /*!
    \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param ip (input) Integration point
    */
    virtual void GetShFnc( Vector<Double>& S, const LocPoint& lp,
                           const Elem* ptElem,  UInt comp = 1 );
    
    virtual void GetShFnc( Matrix<Double>& S,  LocPoint& lp,
                           const Elem* ptElem,  UInt comp = 1 ){
      EXCEPTION("This GetShFnc is not implemented for H1 elements")
    }

    //! Return global derivative of shape functions
    void GetGlobDerivShFnc( Matrix<Double>& deriv, 
                            const LocPointMapped& lp,
                            const Elem* elem, UInt comp = 1 );

    //! Return local derivative of shape functions
    Matrix<Double>&  GetLocDerivShFnc( const LocPoint& lp,
                                       const Elem* elem,
                                       UInt comp = 1 );
    //@}

    // ======================================================================
    //  Incompatible modes for mechanical softening 
    // ======================================================================
    //@{ \name Incompatible modes  
    
    //! Query if element is capable of calculation of incompatible modes
    bool HasICModes() {
      return hasICModes_;
    }
    
    //! Return number of incompatible modes, if present
    virtual UInt GetNumICModes() {
      return 0;
    }
    
    //! Get value of all shape incompatible modes at specific integration point
    /*!
      \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param ip (input) Integration point
     */
    virtual void GetShFncICModes( Vector<Double>& S, const LocPoint& lp,
                                  const Elem* ptElem,  UInt comp = 1 );

    //! Return global derivative of incompatible modes
    void GetGlobDerivShFncICModes( Matrix<Double>& deriv, 
                                   const LocPointMapped& lp,
                                   const Elem* elem, UInt comp = 1 );

    //! Return local derivative of incompatible modes
    void GetLocDerivShFncICModes( Matrix<Double>& deriv, 
                                  const LocPoint& lp,
                                  const Elem* elem, UInt comp = 1 );

    virtual std::string GetFeSpaceName(){
      std::string r = "H1";
      return r;
    }

    //@}
    
  protected:

    //! Compute shape function at given position
    virtual void CalcShFnc( Vector<Double>& shape,
                            const Vector<Double>& point,
                            const Elem* ptElem,
                            UInt comp = 1 ) = 0;

    //! Get local derivatives of all shape fnc at arbitrary local point
    //! Local means here on the reference element
    /*! 
    \param S (output) Matrix with global derivatives of all shape functions
    \f [ \left( \begin{array}{ccc} N_{1,dx} & N_{1,dy} & \cdots \\
    N_{2,dx} & N_{2,dy} & \cdots \\
    \cdots     & \cdots      & \cdots \end{array}\right) \f ]
    \param LCoord (input) Local Coordinates of evaluation point
    \param CornerCoords (input) Coordinates of element corners
    \f [ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
    \cdots & \cdots & \cdots \end{array} \right) \f ]       
     */
    virtual void CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                    const Vector<Double>& point,
                                    const Elem* ptElem,
                                    UInt comp = 1 ) = 0;
    
    //! Flag defining the completeness of the polynomial
    PolyCompleteType completeType_;
    
    // =======================================================================
    //  PRE CALCULATION OF SHAPE FUNCTIONS AT INTEGRATION POINTS
    //  Changed here to map data structure as it is more versatile and
    //  seeing the fact, that we do not have millions of integration points,
    //  it should be ok
    // =======================================================================

    //! Stores Shape Functions for each integration point defined
    std::unordered_map<Integer, Vector<Double> > shapeFncsAtIp_;

    //! Stores shape function derivatives for each integration point
    std::unordered_map<Integer, Matrix<Double> > shapeFncDerivsAtIp_;
    
    // ======================================================================
    //  Incompatible modes for mechanical softening 
    // ======================================================================
    //@{ \name Incompatible modes for mechanical softening

    //! Flag, if incompatible modes are implemented for this element
    bool hasICModes_;

    //! Compute shape function of incompatible modes
    virtual void CalcShFncICModes( Vector<Double>& shape,
                                   const Vector<Double>& point,
                                   const Elem* ptElem,
                                   UInt comp = 1 ) {
      EXCEPTION( "Element does not support incompatible modes");
    }

    //! Compute local derivative of incompatible modes
    virtual void CalcLocDerivShFncICModes( Matrix<Double> & deriv, 
                                           const Vector<Double>& point,
                                           const Elem* ptElem,
                                           UInt comp = 1 ) {
      EXCEPTION( "Element does not support incompatible modes");
    }

    //! Stores Shape Functions for each integration point defined
    std::unordered_map<Integer, Vector<Double> > icModesAtIp_;

    //! Stores shape function derivatives for each integration point
    std::unordered_map<Integer, Matrix<Double> > icModesDerivsAtIp_;

    //@}
  private:
    //! Cached value for local derivative
    Matrix<Double> locDeriv_;
    
  };

} // namespace CoupledField

#endif
