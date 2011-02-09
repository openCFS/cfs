// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

//-*- C++ -*-
//
// 1 Dimensional Expression Templates
// (Glommable and Disambiguated)
//
// ref. Computers in Physics, 11, 263 (1997)
//      Computers in Physics, 10, 552 (1996)
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

#ifndef XPR1_H_
#define XPR1_H_

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdexcept>

#include "MatVec/promote.hh"

template <class P>
class Identity {
public:
  static inline P apply(const P& a) { return a; }
  static inline unsigned int size() {return 1;}
};

template <class P>
class UnaryMinus {
public:
  static inline P apply(const P& a) { return -a; }
  static inline unsigned int size()  {return 1;}
};

template <class P>
class Exp {
public:
  static inline P apply(const P& a) { return P( exp(a) ); }
  static inline unsigned int size()  {return 1;}
};

template <class P, class A, class F>
class Xpr1Func {
  A a;
public:
  inline Xpr1Func(const A& a_) : a(a_) {}
  inline P operator()(unsigned int n) const { return F::apply( a(n) ); }
  inline unsigned int size() const {return a.size();}
};




#define XXX(ap,op) \
  template <class P, class P2>                           \
class ap {\
public:\
  static inline PROMOTE(P,P2) apply(const P& a, const P2& b) { return a op b; }       \
  static inline unsigned int size() {return 1;}\
};
XXX(OpAdd,+)
XXX(OpSub,-)
XXX(OpMul,*)
XXX(OpDiv,/)
#undef XXX

template <class P, class A, class B, class Op>
class Xpr1BinOp {
  A a;
  B b;
public:
  inline Xpr1BinOp(const A& a_, const B& b_) : a(a_), b(b_) {}
  inline P operator() (unsigned int n) const { return Op::apply( a(n), b(n) ); }
  inline unsigned int size() const {return a.size() ;}

};

// Specialziation for the all operations with Scalar
template <class P, class A, class Op, class P2>
class Xpr1OpScalar {
  A a;
  P b;
public:
  inline Xpr1OpScalar(const A& a_, const P& b_) : a(a_), b(b_) {}
  inline PROMOTE(P,P2) operator() (unsigned int n) const { return Op::apply( a(n), b ); }
  inline unsigned int size() const {return a.size();}

};

// Specialziation for the +,- and * operations with Scalar
template <class P, class B, class Op, class P2>
class Xpr1ScalarOp {
  P a;
  B b;
public:
  inline Xpr1ScalarOp(const P& a_, const B& b_) : a(a_), b(b_) {}
  inline PROMOTE(P,P2) operator() (unsigned int n) const { return Op::apply( a, b(n) ); }
  inline unsigned int size() const {return b.size();}
};

template <class P, class V>
class ConstRef1 {
  const V& v;
public:
  inline ConstRef1(const V& v_) : v(v_) {}
  inline P operator()(unsigned int n) const { return v(n); }
  inline unsigned int size() const {return v.size();}
};


template <class P, class E>
class Xpr1 {
private:
  E e;
public:
  inline Xpr1(const E& e_) : e(e_) {}
  inline P operator() (unsigned int n) const { return e(n); }
  inline unsigned int size() const {return e.size();}
};

//template <unsigned int N, class P, class I> class TDim1;

// 1 Dimensional Array Base Class
// for Glommable Expression Templates
template <class P, class I>
class Dim1 {
private:
  void error(const char *msg) const {
    std::cerr << "Dim1 error: " << msg << std::endl;
    throw std::runtime_error("Dim1 error in Xpr1");
  }

public:
  explicit Dim1() {}
  unsigned int size() const {
    return static_cast<const I*>(this)->size();
  }

  void resize( unsigned int n) const {
    static_cast<const I*>(this)->Resize(n);
  }

  P operator() (unsigned int n) const {
    return static_cast<const I*>(this)->operator()(n);
  }
  template <class E> I& 
  assignFrom(const Xpr1<P,E>& x) {
    I *me = static_cast<I*>(this);
    me->Resize(x.size());
    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    return *me;
  }
  template <class V> 
  I& assignFrom(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);
    me->Resize(x.size());

    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) = x(i);
    return *me;
  }
  I& assignFrom(P x) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) = x;
    return *me;
  }
  template <class E> Dim1<P,I>& operator+=(const Xpr1<P,E>& x) {
    I* me = static_cast<I*>(this);

    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) += x(i);
    return *me;
  }
  template <class V> Dim1<P,I>& operator+=(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);

    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) += x(i);
    return *me;
  }
  Dim1<P,I>& operator+=(P x) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) += x;
    return *me;
  }
  template <class E> Dim1<P,I>& operator-=(const Xpr1<P,E>& x) {
    I* me = static_cast<I*>(this);

    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) -= x(i);
    return *me;
  }
  template <class V> Dim1<P,I>& operator-=(const Dim1<P,V>& x) {
    I *me = static_cast<I*>(this);

    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) -= x(i);
    return *me;
  }
  Dim1<P,I>& operator-=(P x) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) -= x;
    return *me;
  }
  Dim1<P,I>& operator*=(P x) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) *= x;
    return *me;
  }
  Dim1<P,I>& operator/=(P x) {
    I *me = static_cast<I*>(this);
    for (unsigned int i=0; i < me->size(); i++) me->operator()(i) /= x;
    return *me;
  }
  
  template <class E> double in(const Xpr1<P,E>& x) const {
    double sum = 0.0;
    const I* me = static_cast<const I*>(this);

    for (unsigned int i=0; i < me->size(); i++) sum += me->operator()(i) * x(i);
    return sum;
  }
  template <class V> double in(const Dim1<P,V>& x) const {
    double sum = 0.0;
    const I* me = static_cast<const I*>(this);

    for (unsigned int i=0; i < me->size(); i++) sum += me->operator()(i) * x(i);
    return sum;
  }
};

template <class T,class A>
std::ostream& operator<<(std::ostream& s, const Dim1<T,A>& a) {
  for (unsigned int i=0; i< a.size(); i++)
    s << std::setw(6) << std::setprecision(3) << a(i) << std::endl;
  return s;
}



// Functions of Dim1 (NO TYPE PROMOTION)
#define XXX(f,ap) \
template <class P,class A> \
Xpr1<P, Xpr1Func<P, ConstRef1<P,Dim1<P,A> >, ap<P> > > \
static inline f(const Dim1<P,A>& a) \
{\
  typedef Xpr1Func<P, ConstRef1<P, Dim1<P,A> >, ap<P> > ExprT; \
   return Xpr1<P,ExprT>(ExprT(ConstRef1<P,Dim1<P,A> >(a))); \
}
XXX(ident, Identity)
XXX(operator- , UnaryMinus)
XXX(exp, Exp)
#undef XXX


// Functions of Xpr1 (NO TYPE PROMOTION)
#define XXX(f,ap) \
template <class P, class A> \
Xpr1<P, Xpr1Func<P, Xpr1<P,A>, ap<P> > > \
static inline f(const Xpr1<P,A>& a) \
{\
   typedef Xpr1Func<P, Xpr1<P,A> , ap<P> > ExprT;\
   return Xpr1<P,ExprT>(ExprT(a)); \
}
XXX(operator- , UnaryMinus)
XXX(exp, Exp)
#undef XXX

//Binary operations between Two Dim1s (type promotion)
#define XXX(op,ap) \
template <class P,class A,class B, class P2>                          \
Xpr1<PROMOTE(P,P2), Xpr1BinOp<PROMOTE(P,P2),ConstRef1<P, Dim1<P,A> >, ConstRef1<P2,Dim1<P2,B> >, ap<P,P2> > > \
static inline op (const Dim1<P,A>& a, const Dim1<P2,B>& b) {                        \
  typedef    Xpr1BinOp<PROMOTE(P,P2), ConstRef1<P, Dim1<P,A> >, ConstRef1<P2, Dim1<P2,B> >, ap<P,P2> > \
ExprT;\
return Xpr1<PROMOTE(P,P2),ExprT>(ExprT(ConstRef1<P,Dim1<P,A> >(a), \
                                       ConstRef1<P2,Dim1<P2,B> >(b)));  \
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

//Binary operations between Dim1 and Scalar (type promotion)
#define XXX(op,ap) \
  template <class P,class A, class P2>\
  Xpr1<PROMOTE(P,P2), Xpr1OpScalar<P2, ConstRef1<P, Dim1<P,A> >, ap<P,P2>, P > > \
static inline op (const Dim1<P,A>& a, const P2& b) {\
    typedef Xpr1OpScalar<P2, ConstRef1<P, Dim1<P,A> >, ap<P,P2>, P > ExprT; \
  return Xpr1<PROMOTE(P,P2),ExprT>(ExprT(ConstRef1<P,Dim1<P,A> >(a),b)); \
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between Scalar and Dim1 (type promotion)
#define XXX(op,ap) \
template <class P,class B,class P2>\
Xpr1<PROMOTE(P,P2), Xpr1ScalarOp<P2, ConstRef1<P, Dim1<P,B> >, ap<P2,P>, P > > \
static inline op (const P2& a, const Dim1<P,B>& b) {\
  typedef Xpr1ScalarOp<P2, ConstRef1<P, Dim1<P,B> >, ap<P2,P>, P > ExprT; \
return Xpr1<PROMOTE(P,P2),ExprT>(ExprT(a,ConstRef1<P,Dim1<P,B> >(b))); \
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

//Binary operations between Xpr1 and Scalar (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P,class A>\
Xpr1<P, Xpr1OpScalar<P, Xpr1<P,A>, ap<P,P>, P > >         \
static inline op (const Xpr1<P,A>& a, const P& b) {\
  typedef Xpr1OpScalar<P, Xpr1<P,A>, ap<P,P>, P > ExprT;  \
  return Xpr1<P,ExprT>(ExprT(a,b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
XXX(operator/, OpDiv)
#undef XXX

//Binary operations between Scalar and Xpr1 (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P,class B> \
Xpr1<P, Xpr1ScalarOp<P, Xpr1<P,B>, ap<P,P>, P > >         \
static inline op (const P& a, const Xpr1<P,B>& b) {\
  typedef Xpr1ScalarOp<P, Xpr1<P,B>, ap<P,P>, P > ExprT;  \
  return Xpr1<P,ExprT>(ExprT(a,b));\
}
XXX(operator+, OpAdd)
XXX(operator-, OpSub)
XXX(operator*, OpMul)
#undef XXX

//Binary operations between Dim1 and Xpr1 (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P,class A,class B>\
Xpr1<P, Xpr1BinOp<P, ConstRef1<P, Dim1<P,A> >, Xpr1<P,B>, ap<P,P> > > \
static inline op (const Dim1<P,A>& a, const Xpr1<P,B>& b) {\
  typedef Xpr1BinOp<P, ConstRef1<P, Dim1<P,A> >, Xpr1<P,B>, ap<P,P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(ConstRef1<P,Dim1<P,A> >(a), b));\
}

XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

//Binary operations between Xpr1 and Dim1 (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P,class A,class B>\
Xpr1<P, Xpr1BinOp<P, Xpr1<P,A>, ConstRef1<P, Dim1<P,B> >, ap<P,P> > > \
static inline op (const Xpr1<P,A>& a, const Dim1<P,B>& b) {\
  typedef Xpr1BinOp<P, Xpr1<P,A>, ConstRef1<P, Dim1<P,B> >, ap<P,P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a,ConstRef1<P,Dim1<P,B> >(b)));\
}

XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX

//Binary operations between Two Xpr1's (NO TYPE PROMOTION)
#define XXX(op,ap) \
template <class P, class A, class B>\
Xpr1<P, Xpr1BinOp<P, Xpr1<P,A>, Xpr1<P,B>, ap<P,P> > >\
static inline op (const Xpr1<P,A>& a, const Xpr1<P,B>& b) {\
  typedef Xpr1BinOp<P, Xpr1<P,A>, Xpr1<P,B>, ap<P,P> > ExprT;\
  return Xpr1<P,ExprT>(ExprT(a, b));\
}

XXX(operator+, OpAdd)
XXX(operator-, OpSub)
#undef XXX
#endif // XPR1_H_
