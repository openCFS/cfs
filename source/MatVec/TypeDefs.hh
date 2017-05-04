#ifndef OLAS_TYPEDEFS_HH
#define OLAS_TYPEDEFS_HH

//! \file TypeDefs.hh
//! This file collects several type and class definitions related to the use
//! of matrices and vectors with block entries.

#include <string>
#include <def_expl_templ_inst.hh>
#include "General/defs.hh"

namespace CoupledField {

  //! Class for determining matrix/vector block size

  //! This class is used to associate with a template that specifies a tiny
  //! matrix or vector the dimension of that matrix or vector. The class
  //! contains a single public attribute called <em>size</em> of type Integer
  //! that stores that dimension. Note that we only consider square matrices
  //! so one integer is enough to specify a matrix' dimension.
  template<class T>
  class BlockSize {
  public:
    static const Integer size = 1;
  };

} // namespace
#endif // OLAS_TYPEDEFS_HH
