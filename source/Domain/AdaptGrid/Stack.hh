// - C++ -
/***************************************************************************
    File        : Stack.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Wed Dec 12 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

 // -*-
// General Description
//
//   a vector based stack implementation
//   for index-based random access


#ifndef STACK_H
#define STACK_H

#include <vector>


namespace grd {


template <class T, class Sequence = std::vector<T> >
class Stack {
public:
  typedef typename Sequence::value_type value_type;
  typedef typename Sequence::size_type size_type;
  typedef typename Sequence::reference reference;
  typedef typename Sequence::const_reference const_reference;

protected:
  Sequence c;

public:
  bool empty() const { return c.empty(); }
  //size_type size() const { return c.size(); }
  size_t size() const { return c.size(); }
  void resize(int ns) { c.resize(ns); }

  reference top() { return c.back(); }
  const_reference top() const { return c.back(); }
  void push(const value_type& x) { c.push_back(x); }
  void pop() { c.pop_back(); }

  // add vector subscripting
  reference operator[](size_type n) { return *(c.begin() + n); }
  const_reference operator[](size_type n) const { return *(c.begin() + n); }

  reference operator()(size_type n) { return *(c.begin() + n); }
  const_reference operator()(size_type n) const { return *(c.begin() + n); }

  // add memory pre-allocation for efficiency
  void reserve(size_type n) { c.reserve(n); }
};


} // namespace grd

#endif // STACK_H
