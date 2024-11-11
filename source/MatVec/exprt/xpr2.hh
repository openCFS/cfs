// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

//-*- C++ -*-
//
// Disambiguated Glommarble Expression Templates
// for Matrix Calculations
//
// ref. Computers in Physics, 11, 263 (1997)
//
//     Copyright (C) 2000 Masakatsu Ito
//     Institute for Molecular Science
//     Myodaiji  Okazaki  444-0867 Japan
//     E-mail mito@ims.ac.jp

//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 

#ifndef XPR2_H_
#define XPR2_H_

#include <stdexcept>
#include "xpr1.hh"

template <class P, class A, class F>
class Xpr2Func {
  A a;
public:
  inline Xpr2Func(const A& a_) : a(a_) {}
  inline P operator()(unsigned int i, unsigned int j) const { return F::apply( a(i,j) ); }
  inline unsigned int rows() const {return a.rows();}
  inline unsigned int cols() const {return a.cols();}

};

template <class P, class A, class F>
class Xpr2FuncTrans {
  A a;
public:
  inline Xpr2FuncTrans(const A& a_) : a(a_) {}
  inline P operator()(unsigned int i, unsigned int j) const { return F::apply( a(j,i) ); }
  inline unsigned int rows() const {return a.cols();}
  inline unsigned int cols() const {return a.rows();}

};

template <class P, class A, class F>
class Xpr2FuncHermitian {
  A a;
public:
  inline Xpr2FuncHermitian(const A& a_) : a(a_) {}
  inline P operator()(unsigned int i, unsigned int j) const { return F::apply( a(j,i) ); }
  inline unsigned int rows() const {return a.cols();}
  inline unsigned int cols() const {return a.rows();}

};

// Specialication for complex case
template <class A, class F>
class Xpr2FuncHermitian<Complex, A, F> {
  A a;
public:
  inline Xpr2FuncHermitian(const A& a_) : a(a_) {}
  inline Complex operator()(unsigned int i, unsigned int j) const { return F::apply( std::conj(a(j,i)) ); }
  inline unsigned int rows() const {return a.cols();}
  inline unsigned int cols() const {return a.rows();}

};

template <class P, class A, class B, class Op>
class Xpr2BinOp {
  A a;
  B b;
public:
  inline Xpr2BinOp(const A& a_, const B& b_) : a(a_), b(b_) {}
  inline P operator()(unsigned int i, unsigned int j) const { return Op::apply( a(i,j), b(i,j) ); }
  inline unsigned int rows() const {return a.rows();}
  inline unsigned int cols() const {return a.cols();}
};

// Specialziation for Scalar Addition and Substractin
// (Not Scalar Multiplication or Division)
template <class P, class A, class Op, class P2>
class Xpr2OpScalar {
  A a;
  P b;
public:
  inline Xpr2OpScalar(const A& a_, const P& b_) : a(a_), b(b_) {}
  inline PROMOTE(P,P2) operator()(unsigned int i, unsigned int j) const { return Op::apply( a(i,j), b ); }
  inline unsigned int rows() const {return a.rows();}
  inline unsigned int cols() const {return a.cols();}
};

template <class P, class B, class Op, class P2>
class Xpr2ScalarOp {
  P a;
  B b;
public:
  inline Xpr2ScalarOp(const P& a_, const B& b_) : a(a_), b(b_) {}
  inline PROMOTE(P,P2) operator()(unsigned int i, unsigned int j) const { return Op::apply( a, b(i,j) ); }
  inline unsigned int rows() const {return b.rows();}
  inline unsigned int cols() const {return b.cols();}
};


template <class P, class M>
class ConstRef2 {
  const M& m;
public:
  inline ConstRef2(const M& m_) : m(m_) {}
  inline P operator()(unsigned int i, unsigned int j) const { return m(i,j); }
  inline unsigned int rows() const {return m.rows();}
  inline unsigned int cols() const { return m.cols(); }
};

template <class P, class E>
class Xpr2 {
private:
  E e;
public:
  inline Xpr2(const E& e_) : e(e_) {}
  inline P operator()(unsigned int i, unsigned int j) const { return e(i,j); }
  inline unsigned int rows() const {return e.rows();}
  inline unsigned int cols() const {return e.cols();}
};

//template <unsigned int Nr,unsigned int Nc,class P,class I> class TDim2;

template <class P, class I>
class Dim2 {
private:
  void error(const char *msg) const {
    std::cerr << "Dim2 error: " << msg << std::endl;
    throw std::runtime_error("Dim3 error in Xpr2");
  }

public:
  explicit Dim2() {}
  virtual ~Dim2() {}
  unsigned int size() const { return static_cast<const I*>(this)->size(); }
  unsigned int rows() const { return static_cast<const I*>(this)->rows(); }
  unsigned int cols() const { return static_cast<const I*>(this)->cols(); }
  void Resize(unsigned int nRows, unsigned int nCols) {
    static_cast<const I*>(this)->Resize(nRows,nCols);
  }
  //P operator() (int n) const {
  //  return static_cast<const I*>(this)->operator()(n);
  //}
  P operator() (unsigned int i, unsigned int j) const {
    return static_cast<const I*>(this)->operator()(i,j);
  }
  template <class X> I& assignFrom(const Xpr2<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    me->Resize(rhs.rows(),rhs.cols());
    for (unsigned int i=0; i<me->rows(); i++) {
      for (unsigned int j=0; j<me->cols(); j++) me->operator()(i,j) = rhs(i,j);
    }
    return *me;
  }
  template <class M> I& assignFrom(const Dim2<P,M>& x) {
    I *me = static_cast<I*>(this);
    me->Resize(x.rows(),x.cols());
    //for (int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    for (unsigned  int i=0; i<me->rows(); i++) {
      for (unsigned int j=0; j<me->cols(); j++) me->operator()(i,j) = x(i,j);
    }
    return *me;
  }
  //template <int Nr,int Nc,class M> I& assignFrom(const TDim2<Nr,Nc,P,M>& x);
  I& assignFrom(P x) {  // for Square Matrix Only
    I *me = static_cast<I*>(this);
    unsigned int i,j, rsz = me->rows(), csz = me->cols(),
      //min = ( rsz < csz ) ? rsz : csz ;
      min = rsz;
    for (i=0; i<min; i++) {
      for (j=0; j<i; j++) me->operator()(i,j) = P(0);
      me->operator()(i,i) = x;
      for (j=i+1; j<csz; j++) me->operator()(i,j) = P(0);
    }
    //for (i=min; i<rsz; i++) {
    //  for (j=0; j<csz; j++) me->operator()(i,j) = P(0);
    //}
    return *me;
  }

  template <class X> Dim2<P,I>& operator+=(const Xpr2<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i<me->rows(); i++) {
      for (unsigned int j=0; j<me->cols(); j++) me->operator()(i,j) += rhs(i,j);
    }
    return *me;
  }
  template <class M> Dim2<P,I>& operator+=(const Dim2<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    //for (int i=0; i < me->size(); i++) me->operator()(i) += x(i);
    for (unsigned  int i=0; i<me->rows(); i++) {
      for (unsigned int j=0; j<me->cols(); j++) me->operator()(i,j) += rhs(i,j);
    }
    return *me;
  }
  Dim2<P,I>& operator+=(P x) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i<me->rows(); i++) 
      for (unsigned int j=0; j<me->cols(); j++)
        me->operator()(i,j) += x;
    return *me;
  }
 
  template <class X> Dim2<P,I>& operator-=(const Xpr2<P,X>& rhs) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i<me->rows(); i++) {
      for (unsigned int j=0; j<me->cols(); j++) me->operator()(i,j) -= rhs(i,j);
    }
    return *me;
  }
  template <class M> Dim2<P,I>& operator-=(const Dim2<P,M>& rhs) {
    I *me = static_cast<I*>(this);
    //for (int i=0; i < me->size(); i++) me->operator()(i) -= x(i);
    for (unsigned int i=0; i<me->rows(); i++) {
      for (unsigned int j=0; j<me->cols(); j++) me->operator()(i,j) -= rhs(i,j);
    }
    return *me;
  }
  Dim2<P,I>& operator-=(P x) { 
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i<me->rows(); i++) 
      for (unsigned int j=0; j<me->cols(); j++)
        me->operator()(i,j) -= x;
    return *me;
  }
 
  Dim2<P,I>& operator*=(P x) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i<me->rows(); i++) {
      for (unsigned int j=0; j<me->cols(); j++) me->operator()(i,j) *= x;
    }
    return *me;
  }

  Dim2<P,I>& operator/=(P x) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i<me->rows(); i++) {
      for (unsigned int j=0; j<me->cols(); j++) me->operator()(i,j) /= x;
    }
    return *me;
  }
};

template <class T,class A>
std::ostream& operator<<(std::ostream& s, const Dim2<T,A>& a) {
  for (unsigned int i=0; i< a.rows(); i++) {
    for (unsigned int j=0; j< a.cols(); j++)
      s << std::setw(6) << std::setprecision(3) << a(i,j);
    s << std::endl;
  }
  return s;
}


// Matrix Vector Multiplication
// Dim2 * Dim1
// Dim2 * Xpr1
// Xpr2 * Dim1
template <class P, class A, class B >
class Xpr1Reduct {
  A a;
  B b;
public:
  inline Xpr1Reduct(const A& a_, const B& b_) : a(a_), b(b_) {}
  inline P operator()(unsigned int i) const {
    P sum = a(i,0)*b(0);
    for (unsigned int j=1; j < a.cols(); j++) sum += a(i,j)*b(j);
    return sum;
  }
  inline unsigned int cols() const {return a.cols(); }
  inline unsigned int rows() const {return a.rows();}
  inline unsigned int size() const {return a.rows();}
};

// Matrix Matrix Multiplication
// Dim2 * Dim2
// Dim2 * Xpr2
// Xpr2 * Dim2
// Xpr2 * Xpr2
template <class P, class A, class B>
class Xpr2Reduct {
  A a;
  B b;
public:
  inline Xpr2Reduct(const A& a_, const B& b_) : a(a_), b(b_) {}
  inline P operator()(unsigned int i, unsigned int j) const {
    P sum = a(i,0)*b(0,j);
    for (unsigned int k=1; k < a.cols(); k++) sum += a(i,k)*b(k,j);
    return sum;
  }
  inline unsigned int cols() const {return b.cols(); }
  inline unsigned int rows() const {return a.rows();}
};



// Functions of Dim2 (NO TYPE PROMOTION)
#define XXX(f,ap) \
template <class P, class A> \
Xpr2<P, Xpr2Func<P, ConstRef2<P,Dim2<P,A> >, ap<P  > > > \
static inline f(const Dim2<P,A>& a) \
{\
   typedef Xpr2Func<P, ConstRef2<P,Dim2<P,A> >, ap<P  > > ExprT;\
   return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a))); \
}
//XXX(ident, Identity)
XXX(operator- , UnaryMinus)
XXX(Conj , Conjugate)
//XXX(exp, Exp)
#undef XXX

// Explicit Transpose of Dim2 (NO TYPE PROMOTION)
#define XXX(f,ap) \
template <class P, class A> \
Xpr2<P, Xpr2FuncTrans<P, ConstRef2<P,Dim2<P,A> >, ap<P  > > > \
static inline f(const Dim2<P,A>& a) \
{\
   typedef Xpr2FuncTrans<P, ConstRef2<P,Dim2<P,A> >, ap<P  > > ExprT;\
   return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a))); \
}
XXX(Transpose, Identity)
//XXX(exp, Exp)
#undef XXX

// Explicit Transpose of Xpr2 (NO TYPE PROMOTION)
#define XXX(f,ap) \
template <class P, class A> \
Xpr2<P, Xpr2FuncTrans<P, Xpr2<P,Dim2<P,A> >, ap<P  > > > \
static inline f(const Xpr2<P,Dim2<P,A> >& a) \
{\
   typedef Xpr2FuncTrans<P, Xpr2<P,Dim2<P,A> >, ap<P  > > ExprT;\
   return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a))); \
}
XXX(Transpose, Identity)
//XXX(exp, Exp)
#undef XXX

// Explicit Hermitian of Dim2 (NO TYPE PROMOTION)
#define XXX(f,ap) \
template <class P, class A> \
Xpr2<P, Xpr2FuncHermitian<P, ConstRef2<P,Dim2<P,A> >, ap<P  > > > \
static inline f(const Dim2<P,A>& a) \
{\
   typedef Xpr2FuncHermitian<P, ConstRef2<P,Dim2<P,A> >, ap<P  > > ExprT;\
   return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a))); \
}
XXX(Herm, Identity)
//XXX(exp, Exp)
#undef XXX


// Functions of Xpr2 (NO TYPE PROMOTION)
#define XXX(f,ap) \
template <class P, class E> \
Xpr2<P, Xpr2Func<P, Xpr2<P,E>, ap<P  > > > \
static inline f(const Xpr2<P,E>& a) \
{\
   typedef Xpr2Func<P, Xpr2<P,E>, ap<P> > ExprT;\
   return Xpr2<P,ExprT>(ExprT(a)); \
}
XXX(operator- , UnaryMinus)
XXX(exp, Exp)
#undef XXX

//Binary operations between Two Dim2s (type promotion)
#define XXX(op,ap) \
template <class P,class A,class B, class P2>                                   \
Xpr2<PROMOTE(P,P2), Xpr2BinOp<PROMOTE(P,P2), ConstRef2<P, Dim2<P,A> >, ConstRef2<P2,Dim2<P2,B> >, ap<P,P2> > > \
static inline op (const Dim2<P,A>& a, const Dim2<P2,B>& b) {\
typedef \
  Xpr2BinOp<PROMOTE(P,P2), ConstRef2<P,Dim2<P,A> >, ConstRef2<P2,Dim2<P2,B> >, ap<P,P2> > \
ExprT;\
return Xpr2<PROMOTE(P,P2),ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),    \
                                         ConstRef2<P2,Dim2<P2,B> >(b))); \
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

// Multiplication with Two Dim2s (type promotion)
#define XXX(op) \
  template <class P,class A,class B,class P2>                                   \
  Xpr2<PROMOTE(P,P2), Xpr2Reduct<PROMOTE(P,P2), ConstRef2<P, Dim2<P,A> >, ConstRef2<P2,Dim2<P2,B> > > > \
static inline op (const Dim2<P,A>& a, const Dim2<P2,B>& b) {\
  typedef \
    Xpr2Reduct<PROMOTE(P,P2), ConstRef2<P,Dim2<P,A> >, ConstRef2<P2,Dim2<P2,B> > > \
      ExprT;\
  return Xpr2<PROMOTE(P,P2),ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),    \
			      ConstRef2<P2,Dim2<P2,B> >(b)));\
}
XXX(operator*)
#undef XXX

// Multiplication between Xpr2 and Dim2 (type promotion)
#define XXX(op) \
  template <class P,class A,class B,class P2>                                   \
  Xpr2<PROMOTE(P,P2), Xpr2Reduct<PROMOTE(P,P2), Xpr2<P,A>, ConstRef2<P2,Dim2<P2,B> > > > \
static inline op (const Xpr2<P,A>& a, const Dim2<P2,B>& b) {\
  typedef \
    Xpr2Reduct<PROMOTE(P,P2), Xpr2<P,A >, ConstRef2<P2,Dim2<P2,B> > > \
      ExprT;\
  return Xpr2<PROMOTE(P,P2),ExprT>(ExprT(Xpr2<P,A>(a),    \
            ConstRef2<P2,Dim2<P2,B> >(b)));\
}
XXX(operator*)
#undef XXX

//Binary operations between Dim2 and Xpr2 (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P,class A,class B> \
Xpr2<P, Xpr2BinOp<P, ConstRef2<P, Dim2<P,A> >, Xpr2<P,B>, ap<P,P> > >\
static inline op (const Dim2<P,A>& a, const Xpr2<P,B>& b) {\
  typedef Xpr2BinOp<P, ConstRef2<P,Dim2<P,A> >, Xpr2<P,B>, ap<P,P> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

/* Not Efficient !! // Muitiplicatino with  Dim2 and Xpr2
#define XXX(op) \
template <class P,class A,class B> \
Xpr2<P, Xpr2Reduct<P, ConstRef2<P, Dim2<P,A> >, Xpr2<P,B> > >\
op (const Dim2<P,A>& a, const Xpr2<P,B>& b) {\
  typedef \
    Xpr2Reduct<P, ConstRef2<P, Dim2<P,A> >, Xpr2<P,B> > \
      ExprT;\
  return Xpr2<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),b));\
}
XXX(operator*)
#undef XXX */

//Binary operations between Dim2 and Scalar (type promotion)
#define XXX(op,ap) \
  template <class P,class A, class P2>                                   \
  Xpr2<PROMOTE(P,P2), Xpr2OpScalar<P2, ConstRef2<P, Dim2<P,A> >, ap<P,P2>, P > > \
static inline op (const Dim2<P,A>& a, const P2& b) {\
  typedef \
    Xpr2OpScalar<P2, ConstRef2<P,Dim2<P,A> >, ap<P,P2>, P >       \
      ExprT;\
  return Xpr2<PROMOTE(P,P2),ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),b)); \
}
  //XXX(operator+, OpAdd) //impossible!
  //XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

// impossible ! // Binary Operations between Dim2 and Dim1

// Multiplication with Dim2 and Dim1 (type promotion)
#define XXX(op)                                                         \
template <class P, class A, class B, class P2 >                       \
Xpr1<PROMOTE(P,P2), Xpr1Reduct<PROMOTE(P,P2), ConstRef2<P, Dim2<P,A> >, ConstRef1<P2,Dim1<P2,B> > > > \
static inline op (const Dim2<P,A>& a, const Dim1<P2,B>& b) {                          \
  typedef                                                               \
    Xpr1Reduct<PROMOTE(P,P2), ConstRef2<P,Dim2<P,A> >, ConstRef1<P2,Dim1<P2,B> > > \
    ExprT;                                                              \
  return Xpr1<PROMOTE(P,P2),ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),    \
                                         ConstRef1<P2,Dim1<P2,B> >(b))); \
}
XXX(operator*)
#undef XXX

// Multiplication with Xpr2 and Dim1 (type promotion)
#define XXX(op)                                                         \
template <class P, class A, class B, class P2 >                       \
Xpr1<PROMOTE(P,P2), Xpr1Reduct<PROMOTE(P,P2), Xpr2<P,A>, ConstRef1<P2,Dim1<P2,B> > > > \
static inline op (const Xpr2<P,A>& a, const Dim1<P2,B>& b) {                          \
  typedef                                                               \
    Xpr1Reduct<PROMOTE(P,P2), Xpr2<P,A>, ConstRef1<P2,Dim1<P2,B> > > \
    ExprT;                                                              \
  return Xpr1<PROMOTE(P,P2),ExprT>(ExprT(Xpr2<P,A>(a),    \
                                         ConstRef1<P2,Dim1<P2,B> >(b))); \
}
XXX(operator*)
#undef XXX


// // Multiplication with Dim2 and Dim1
// #define XXX(op)                                                         \
//   template <class P, class A, class B, class P2>                        \
// Xpr1<P, Xpr1Reduct<P, ConstRef2<P, Dim2<P,A> >, ConstRef1<P,Dim1<P,B> > > >\
// op (const Dim2<P,A>& a, const Dim1<P,B>& b) {\
//   typedef \
//     Xpr1Reduct<P, ConstRef2<P,Dim2<P,A> >, ConstRef1<P,Dim1<P,B> > > \
//       ExprT;\
//   return Xpr1<P,ExprT>(ExprT(ConstRef2<P,Dim2<P,A> >(a),\
// 			      ConstRef1<P,Dim1<P,B> >(b)));\
// }
// XXX(operator*)
// #undef XXX


//Binary operations between Two Xpr2s (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P, class A, class B> \
Xpr2<P, Xpr2BinOp<P, Xpr2<P,A>, Xpr2<P,B>, ap<P,P> > >\
static inline op (const Xpr2<P,A>& a, const Xpr2<P,B>& b) {\
  typedef Xpr2BinOp<P, Xpr2<P,A>, Xpr2<P,B>, ap<P,P> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(Xpr2<P,A>(a),Xpr2<P,B>(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX


//Binary operations between Xpr2 and Scalar (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P, class A> \
Xpr2<P, Xpr2OpScalar<P, Xpr2<P,A>, ap<P,P>, P > >         \
static inline op (const Xpr2<P,A>& a, const P& b) {\
  typedef Xpr2OpScalar<P, Xpr2<P,A>, ap<P,P>, P > ExprT;  \
  return Xpr2<P,ExprT>(ExprT(a,b));\
}
  //XXX(operator+, OpAdd) // impossible!
  //XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between Xpr2 and Dim2 (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P,class A,class B> \
Xpr2<P, Xpr2BinOp<P, Xpr2<P,A>, ConstRef2<P,Dim2<P,B> >, ap<P,P> > >\
static inline op (const Xpr2<P,A>& a, const Dim2<P,B>& b) {\
  typedef Xpr2BinOp<P, Xpr2<P,A>, ConstRef2<P,Dim2<P,B> >, ap<P,P> > ExprT;\
  return Xpr2<P,ExprT>(ExprT(a, ConstRef2<P,Dim2<P,B> >(b)));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

//Binary operations between Scalar and Dim2 (type promotion)
#define XXX(op,ap) \
template <class P,class B, class P2>                                   \
Xpr2<PROMOTE(P,P2), Xpr2ScalarOp<P2, ConstRef2<P,Dim2<P,B> >, ap<P2,P>, P > > \
static inline op (const P2& a, const Dim2<P,B>& b) {\
  typedef \
    Xpr2ScalarOp<P2, ConstRef2<P,Dim2<P,B> >, ap<P2,P>, P >       \
      ExprT;\
  return Xpr2<PROMOTE(P,P2),ExprT>(ExprT(a,ConstRef2<P,Dim2<P,B> >(b))); \
}
  //XXX(operator+, OpAdd) //impossible!
  //XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

//Binary operations between Scalar and Xpr2 (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P, class B> \
Xpr2<P, Xpr2ScalarOp<P, Xpr2<P,B>, ap<P,P>, P > >         \
static inline op (const P& a, const Xpr2<P,B>& b) {\
  typedef Xpr2ScalarOp<P, Xpr2<P,B>, ap<P,P>, P > ExprT;  \
  return Xpr2<P,ExprT>(ExprT(a,b));\
}
  //XXX(operator+, OpAdd) // impossible!
  //XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

// ALREADY DEFINED! //Binary operations between Scalars
// Already Defined! //Binary operations between Scalar and Dim1
// Already Defined! //Binary operations between Scalar and Xpr1
// Already Defined! //Binary operations between Scalar and DimD
// Already Defined! //Binary operations between Scalar and XprD

// impossible ! //Binary operations between Dim1 and Dim2
// impossible ! //Multiplication with Dim1 and Dim2
// impossible ! //Binary operations between Dim1 and Xpr2
// impossible ! //Multiplication with Dim1 and Xpr2
// Already Defined! //Binary operations between Dim1 and Scalar
// Already Defined! //Binary operations between Dim1 and Dim1
// Already Defined! //Binary operations between Dim1 and Xpr1
// Already Defined! //Binary operations between Dim1 and DimD
// Already Defined! //Binary operations between Dim1 and XprD

// impossible ! //Binary operations between Xpr1 and Dim2
// impossible ! //Multiplication with Xpr1 and Dim2
// impossible ! //Binary operations between Xpr1 and Xpr2
// impossible ! //Multiplication with Xpr1 and Xpr2
// Already Defined! //Binary operations between Xpr1 and Scalar
// Already Defined! //Binary operations between Xpr1 and Dim1
// Already Defined! //Binary operations between Xpr1 and Xpr1
// Already Defined! //Binary operations between Xpr1 and DimD
// Already Defined! //Binary operations between Xpr1 and XprD
// Already Defined! //Binary operations between DimD and Scalar

// impossible ! //Binary Operations between DimD and Dim1
// impossible ! //Binary Operations between DimD and Xpr1
// Already Defined! //Binary operations between DimD and DimD
// Already Defined! //Binary operations between DimD and XprD
// Already Defined! //Binary operations between XprD and Scalar

// impossible ! //Binary Operations between XprD and Dim1
// impossible ! //Binary Operations between XprD and Xpr1
// Already Defined! //Binary operations between XprD and DimD
// Already Defined! //Binary operations between XprD and XprD

#endif // XPR2_H_

