// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_ANSATZFCT_HH
#define FILE_CFS_ANSATZFCT_HH

#include "General/environment.hh"
#include "Matrix/matrix.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {


  //! Base class representing an ansatz function space
  class AnsatzFct {

  public:
    typedef enum {CONST, LAGRANGE, LEGENDRE, NEDELEC} AnsatzFctType;

    typedef enum {NONE, ALL, NODE, EDGE, FACE, INTERIOR} FctEntityType;

    //! Constuctor
    AnsatzFct();

    //! Destructor
    virtual ~AnsatzFct();
    
    //! Return type of ansatz functions
    AnsatzFctType GetType() const { return type_; }

    //! Return true, if approximation is isotropic
    bool IsIsotropic() const { return isIsotropic_; }

  protected:
    
    //! Type of ansatz functions used
    AnsatzFctType type_;

    //! Flag indicating isotropy of ansatz functions
    bool isIsotropic_;

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

  //! Class for hierarchic Legendre ansatz functions
  class LegendreFct : public AnsatzFct{ 

  public:

    //! Define enum for different sub-spaces of Legendre polynomials
    typedef enum {TENSOR, PRODUCT} SubSpaceType;
    
    //! Constructor
    LegendreFct();
    
    //! Set isotropic order
    void SetIsoOrder( UInt order );

    //! Return isotropic order
    UInt GetIsoOrder() const;

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
