// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef _PLAIN_ALLOC_HH
#define _PLAIN_ALLOC_HH

// 
// The following code example is taken from the book
// The C++ Standard Library - A Tutorial and Reference
// by Nicolai M. Josuttis, Addison-Wesley, 1999
// © Copyright Nicolai M. Josuttis 1999

#include <limits>
#include <iostream>

namespace OLAS {


  //! Plain (de-)allocator class which does not use memory pools or caching

  //! This class implements a plain (de-)allocator, which directly uses the
  //! global new() and delete() methods without keeping a pool / cache of
  //! allocated memory for similar objects. 
  //! This can be used to prevent memory bottlenecks, which can occur for 
  //! example with big std::list objects, where the object is destroyed,
  //! but the memory is not available for other objects.
  //!
  //! A typical exmaple for the use with STL-containers could look like this:
  //! \code
  //! #include <list>
  //! #include "utils/plainAlloc.hh"
  //!
  //! std::list<UInt,PlainAlloc<UInt> > myList;
  //! \endcode
  //! \note Since the definition of an extra memory-allocator changes the 
  //! signature of the STL-container, it is convenient to use an own typedef:
  //! \code
  //! #include <list>
  //! #include "utils/plainAlloc.hh"
  //!
  //! typedef std::list<Integer,PlainAlloc<Integer> > SimpleIntList;
  //! SimpleIntList myList;
  //! SimpleIntList::iterator it = myList.begin();
  //! \endcode

  template <class T>
  class PlainAlloc {

  public:

    //@{ \name type definitions
    typedef T        value_type;
    typedef T*       pointer;
    typedef const T* const_pointer;
    typedef T&       reference;
    typedef const T& const_reference;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;
    //@}
    

    //! rebind allocator to type U
    template <class U>
    struct rebind {
      typedef PlainAlloc<U> other;
    };

    //@{ \name return address of values

    //! returns the address
    pointer address (reference value) const {
      return &value;
    }

    //! returns the const adress
    const_pointer address (const_reference value) const {
      return &value;
    }
    //@}
    
    //@{ \name constructors and destructor
    // nothing to do because the allocator has no state

    //! standard constructor
    PlainAlloc() throw() {
    }

    //! alternative constructor
    PlainAlloc(const PlainAlloc&) throw() {
    }

    //! alternative constructor
    template <class U>
    PlainAlloc (const PlainAlloc<U>&) throw() {
    }

    //! standard destructor
    ~PlainAlloc() throw() {
    }
    //@}

    //! return maximum number of elements that can be allocated
    size_type max_size () const throw() {
      return std::numeric_limits<std::size_t>::max() / sizeof(T);
    }

    //! allocate but don't initialize num elements of type T
    pointer allocate (size_type num, const void* = 0) {
      // allocate memory with global new
      pointer ret = (pointer)(::operator new(num*sizeof(T)));
      return ret;
    }

    //! initialize elements of allocated storage p with value value
    void construct (pointer p, const T& value) {
      // initialize memory with placement new
      new((void*)p)T(value);
    }

    //! destroy elements of initialized storage p
    void destroy (pointer p) {
      // destroy objects by calling their destructor
      p->~T();
    }

    //! deallocate storage p of deleted elements
    void deallocate (pointer p, size_type num) {
      // deallocate memory with global delete
      ::operator delete((void*)p);
    }
  };

  // return that all specializations of this allocator are interchangeable
  template <class T1, class T2>
  bool operator== (const PlainAlloc<T1>&,
                   const PlainAlloc<T2>&) throw() {
    return true;
  }
  
  template <class T1, class T2>
  bool operator!= (const PlainAlloc<T1>&,
                   const PlainAlloc<T2>&) throw() {
    return false;
  }
}
#endif
