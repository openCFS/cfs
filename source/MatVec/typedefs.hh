// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


//! \file typedefs.hh
//! This file collects several type and class definitions related to the use
//! of matrices and vectors with block entries.

#ifndef OLAS_TYPEDEFS_HH
#define OLAS_TYPEDEFS_HH

#include <string>
#include "def_expl_templ_inst.hh"
#include "General/defs.hh"

namespace CoupledField {

  //! Traits for associating Matrices with vectors and scalars
  template<class T>
  class AssocType{
  public:
    typedef T T_Mtype;        // matrix type 
    typedef T T_Vtype;        // vector type
    typedef T T_Stype;        // scalar type
    static std::string tagM;  // symbolic name for humans (matrix type)
    static std::string tagV;  // symbolic name for humans (vector type)
    static std::string tagS;  // symbolic name for humans (scalar type)
  };                                              

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
