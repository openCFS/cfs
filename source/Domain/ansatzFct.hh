// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_ANSATZFCT_HH
#define FILE_CFS_ANSATZFCT_HH

#include "General/defs.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {


  //! Base class representing an ansatz function space
  class AnsatzFct {

  public:
    typedef enum {CONST, LAGRANGE, LEGENDRE, NEDELEC, SPECTRAL} AnsatzFctType;

    typedef enum {NONE, ALL, NODE, EDGE, FACE, INTERIOR} FctEntityType;

    //! Constuctor
    AnsatzFct();

    //! Destructor
    virtual ~AnsatzFct();

    //! Return type of ansatz functions
    AnsatzFctType GetType() const { return type_; }

    //! Return true, if approximation is isotropic
    bool IsIsotropic() const { return isIsotropic_; }

    //! Return true, if Function is Discontinuous from Element to element
    bool IsDiscontinuous() const;

    //! Sets flag if the ansatz function is discontiuous from element to element
    void SetDiscontinuity(bool discontinuous);

  protected:

    //! Type of ansatz functions used
    AnsatzFctType type_;

    //! Flag indicating isotropy of ansatz functions
    bool isIsotropic_;

    //! Flag indicating if the AnsatzFnc is discontinous from element to element
    bool isDiscontinuous_;
  };

  //! Comparison operator for AnsatzFct
  bool operator==( const AnsatzFct& a, const AnsatzFct& b );

  //! Negative comparison operator for AnsatzFct
  bool operator!=( const AnsatzFct& a, const AnsatzFct& b );

  //! Class for constant ansatz functions
  class ConstFct : public AnsatzFct {

  public:

    //! Constructor
    ConstFct();

  };

  //! Class for standard 1st and 2nd orer Lagrange functions
  class LagrangeFct : public AnsatzFct {

  public:

    //! Constructor
    LagrangeFct();

  };

  //! Class for standard 1st and 2nd order Nedelec functions
  class NedelcFct : public AnsatzFct {

  public:

    //! Constructor
    NedelcFct();

  };

  //! Class for Spectral Element Ansatz functions
  class SpectralFct : public AnsatzFct {

  public:

    //! Constructor
    SpectralFct();

    //! Set order of lagrange polinomial
    void SetOrder( UInt order );

    //! Get order of lagrange polinomial
    UInt GetOrder( ) const;

    //! Calculate the derivatives for each supporting point for the given coordinate
    //@param deriv the vector for the resulting derivatives
    //@param coord is the given coordinate
    void EvaluateDerivPolynomial( Vector<Double> & deriv, Double coord );

    //! Calculate the value of the lagrange function for the given coordinates
    //@param shape vector of shape functions
    //@param coord the coordinate to be evaluated
    void EvaluatePolynomial( Vector<Double> & shape, Double coord );

  protected:

    //! Order of Polynomial
    UInt Order_;

    Vector< Double > supportingPoints_;

   	// Gauss Lobatto Integration Points Order: 1
    static Double l1[][1];

    //  Gauss-Lobatto  quadrature  points  and  weights  order  2
    static Double l2[][1];

    //  Gauss-Lobatto  quadrature  points  and  weights  order  3
    static Double l3[][1];

     //  Gauss-Lobatto  quadrature  points  and  weights  order 4
    static Double l4[][1];

    //  Gauss-Lobatto  quadrature  points  and  weights  order 5
    static Double l5[][1];

    //  Gauss-Lobatto  quadrature  points  and  weights  order 6
    static Double l6[][1];

    //  Gauss-Lobatto  quadrature  points  and  weights  order 7
    static Double l7[][1];
  };

  //! Class for standard 1st and 2nd orer Lagrange functions
   class NedelecFct : public AnsatzFct {

   public:

     //! Constructor
     NedelecFct();
     
   };
   
  //! Class for hierarchic Legendre ansatz functions
  class LegendreFct : public AnsatzFct{

  public:

    //! Define enum for different sub-spaces of Legendre polynomials
    typedef enum {TENSOR, PRODUCT} SubSpaceType;

    //! Constructor
    LegendreFct();

    //! Set isotropic order
    virtual void SetIsoOrder( UInt order );

    //! Return isotropic order
    virtual UInt GetIsoOrder() const;

    //! Set anisotropic order matrix
    void SetAnisoOrder( Matrix<UInt>& order );

    //! Return anisotropic order
    const Matrix<UInt>& GetAnisoOrder() const;

    //! Get maximum order w.r.t. to local coordinate direction
    UInt GetMaxOrderLocDir( UInt locDir ) const;

    //! Get maximum order w.r.t. to degree of freedom
    UInt GetMaxOrderDof( UInt dof ) const;

    //! Get maximum polynomial degree
    UInt GetMaxOrder() const;

    // =======================================================================
    // D A T A   M E M B E R S
    // =======================================================================

  protected:

    //! Isotropic order of method
    UInt isoOrder_;

    //! Anisotropic order (local dir x DOFs)
    Matrix<UInt> anOrder_;

    //! Maximum order w.r.t. to local coordinate direction
    StdVector<UInt> maxOrderLocDir_;

    //! Maximum order w.r.t. to degree of freedom
    StdVector<UInt> maxOrderDof_;

    //! Maximum order in complete array
    UInt maxOrder_;

    //! Type of subspace for approximation
    SubSpaceType subSpace_;

  };

}



#endif
