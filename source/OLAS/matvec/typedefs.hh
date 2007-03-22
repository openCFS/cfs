// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


//! \file typedefs.hh
//! This file collects several type and class definitions related to the use
//! of matrices and vectors with block entries.

#include <string>

#ifndef OLAS_TYPEDEFS_HH
#define OLAS_TYPEDEFS_HH

#include "utils/defs.hh"
#include "utils/environment.hh"

namespace OLAS {

  //! Traits for associating Matrices with vectors and scalars
  template<class T>
  class assocType{
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
class assocType<M_TYPE>{                         \
public:                                          \
  typedef M_TYPE T_Mtype;                        \
  typedef V_TYPE T_Vtype;                        \
  typedef S_TYPE T_Stype;                        \
  static std::string tagM;                       \
  static std::string tagV;                       \
  static std::string tagS;                       \
};                                               \
template<>                                       \
class assocType<V_TYPE>{                         \
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

  // Function for determining matrix/vector entry type (i.e. Integer, Double,
  // Complex) for enum-type refer to environment.hh

  //! Class for determining the type of a matrix/vector entry on scalar level

  //! This class is used to associate with a template that specifies a tiny
  //! matrix or vector the type of the matrix' or vector's entry on the
  //! scalar level. The class contains a single public attribute called
  //! <em>M_EntryType</em> of type MatrixEntryType for storing this
  //! information.
  template<class T>
  class EntryType{
  public:
    //static const MatrixEntryType M_EntryType = NOENTRYTYPE;
  };

#define DEFINE_ENTRY_TYPE(M_TYPE,ENTRYTYPE)                     \
template<>                                                      \
class EntryType<M_TYPE>{                                        \
public:                                                         \
  static const MatrixEntryType M_EntryType = ENTRYTYPE;         \
}

  DEFINE_ENTRY_TYPE(float      , DOUBLE);
  DEFINE_ENTRY_TYPE(double     , DOUBLE);
  DEFINE_ENTRY_TYPE(Complex    , COMPLEX);

} // namespace


#endif // OLAS_TYPEDEFS_HH
