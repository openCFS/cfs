// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FEFUNCTION_HH
#define FILE_CFS_FEFUNCTION_HH

#include "General/environment.hh"
#include "MatVec/vector.hh"

namespace CoupledField {

//!  Base class for a function approximated by Finite Elements 
/*!
  The FeFunction class represents an function, approximated by in a finite
  element space (replaces the old Node-/ElemStoreSol class)
  Therefore it holds the following information:
    - Pointer to related Finite Element Space (which knows about the equation 
      numbers)
    - Description of the function (i.e. ResultInfo)
    - List of entities it is defined on (region, element list)
    - Boundary conditions (Dirichlet, Constraints, etc.)
    - Coefficient vector
  
  The class has methods to
  - Aquire the solution of an element
  - 
*/

class BaseFeFunction {
public:


  //! Constructor
  BaseFeFunction();
  
  //! Destructor
  ~BaseFeFunction();

  
protected:

  //! Identifier for the function
  FeFctIdType fctId_;

};


///////////////////////////////////////////////////////////////////
// Templatized Version
///////////////////////////////////////////////////////////////////

//! Templatized FeFunction (real/complex)
template<typename T>
class FeFunction : public BaseFeFunction {
public:

protected:

  //! Coefficient vector
  Vector<T> coeffs_;
};


}  // namespace CoupledField
#endif
