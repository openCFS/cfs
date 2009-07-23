// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FESPACE_HH
#define FILE_CFS_FESPACE_HH

#include "General/environment.hh"
#include "Domain/entityList.hh"
#include "Domain/elem.hh"

namespace CoupledField {


// forward class declarations
class BaseFe;
class FeFunction;

//!  Base class for the Finite Element Space (FeSpace) 
/*!
  The Finite Element Space (FeSpace) representes the mathematic discrete 
  function spacein which a unknown function is defined.
  
  It has the following functionality:
  - Collect pointers to reference elements of related functional space
  - Perform equation numbering
  
  Each unknown function is associated with an appropriate FeSpace. 
  The spaces differ in the following senses
  
    - Continuity of derivatives (H1, HCurl, HDiv, L2, Constants)
    - Order of shape functions
      -* Lower Order: The approximation order is expliclty encoded in the 
                      element itself (i.e. Lagrange Elements)
      -* Higher Order: The functions of the elements are composed of 1D 
                       polynomials (e.g. Legendre, Jacobi, Gegenbauer).
  
  The FeSpaceConst plays a special role: It can be used in conjunction with
  discrete unknowns, e.g. additional network equations.
*/

class FeSpace {
public:

  //! Enumeration type for element space
  //@{
  typedef enum { CONST, H1_LO, H1_HI, HCURL_LO, HCURL_HI, HDIV_Lo, HDIV_HI, L2_LO, L2_HI } Type;
  static Enum<Type> TypeEnum;
  //@}

  //! Constructor
  FeSpace();

  //! Destructor
  ~FeSpace();

  //! Return type of FeSpace
  Type GetType() { return type_;}
  
  //! Generate instance of specific element space type
  shared_ptr<FeSpace> CreateInstance( Type type );
  
  // ========================================================================
  //  ELEMENT HANDLING
  // ========================================================================
  //@{ \name Element Handling

  //! Return pointer to reference element
  virtual shared_ptr<BaseFe> GetFe( const EntityIterator ent ) = 0;
  //@}
  
  // ========================================================================
  //  EQUATION HANDLING
  // ========================================================================
  //@{ \name Equation Handling
  
  //! Return equation numbers
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ) = 0; 
  
  //! Add result
  virtual void AddFeFunction( shared_ptr<FeFunction> fct ) = 0;
  
  //! Map equations i.e. intialize object
  virtual void Finalize() = 0;
  //@}
  
protected:

  //! Type of element space
  Type type_;
  
  //! Map for reference elmenets
  std::map<Elem::FEType, BaseFE*> refElems_;
};


///////////////////////////////////////////////////////////////////
// H1 - Lower Order / non hierarchical
///////////////////////////////////////////////////////////////////

class FeSpaceH1 : public FeSpace {


public:

  //! Type of basis used
  typedef enum {LAGRANGE, DUAL} BasisType; 

  //! Constructor
  FeSpaceH1();

  //! Destructor
  ~FeSpaceH1();

  //! Set type of basis
  void SetBasis( BasisType type );

  //! Determine order
  //! This method determines the order of the generated elements,
  //! as 
  void SetOrder( UInt order);

private:

  //! Type of basis
  BasisType basis_;

  //! Order of functions used
  UInt order_;
};

///////////////////////////////////////////////////////////////////
// H1 - Higher Order / hierarchical
///////////////////////////////////////////////////////////////////

class FeSpaceH1Hi : public FeSpace {

};

///////////////////////////////////////////////////////////////////
// HCurl - Lower Order (Nedelec) / not hierarchical
///////////////////////////////////////////////////////////////////

class FeSpaceHCur : public FeSpace {

};

///////////////////////////////////////////////////////////////////
// HCurl - Higher Order 
///////////////////////////////////////////////////////////////////

class FeSpaceHCurlHi : public FeSpace {

};


}


  // namespace CoupledField
#endif
