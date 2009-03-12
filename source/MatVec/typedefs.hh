// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


//! \file typedefs.hh
//! This file collects several type and class definitions related to the use
//! of matrices and vectors with block entries.

#ifndef OLAS_TYPEDEFS_HH
#define OLAS_TYPEDEFS_HH

#include <string>

#include <def_expl_templ_inst.hh>

#include "General/defs.hh"
//#include "General/environment.hh"

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

#define DEFINE_ASSOC(M_TYPE,V_TYPE,S_TYPE)       \
template<>                                       \
class AssocType<M_TYPE>{                         \
public:                                          \
  typedef M_TYPE T_Mtype;                        \
  typedef V_TYPE T_Vtype;                        \
  typedef S_TYPE T_Stype;                        \
  static std::string tagM;                       \
  static std::string tagV;                       \
  static std::string tagS;                       \
};                                               \
template<>                                       \
class AssocType<V_TYPE>{                         \
public:                                          \
  typedef M_TYPE T_Mtype;                        \
  typedef V_TYPE T_Vtype;                        \
  typedef S_TYPE T_Stype;                        \
  static std::string tagM;                       \
  static std::string tagV;                       \
  static std::string tagS;                       \
}                                                

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

#define DEFINE_BLOCK_SIZE(M_TYPE,SIZE) \
template<>                             \
class BlockSize<M_TYPE>{               \
public:                                \
  static const Integer size = SIZE;    \
}


} // namespace

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "typedefs.cc"
#endif


#endif // OLAS_TYPEDEFS_HH
