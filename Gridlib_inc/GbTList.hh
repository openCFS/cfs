/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBTLIST_HH
#define  GBTLIST_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

// The TList class represents a singly linked list with header node.  The
// class T is intended to be native data (int, float, char*, etc.).  All
// list insertions occur at the front of the list.  If the list values were
// dynamically allocated by the application, it is the application's
// responsibility to delete those values.  This can be done by iterating
// over the list using GetFirst/GetNext and explicitly deleting the values.
//
// To iterate over a list to access values with an implicit iterator, use the
// following.  Because an implicit iterator is used, you cannot use this
// format for nested loops with the same list.
//
//     GbTList<T> kList;  // already built list
//     T tValue;
//     GbBool bFound = kList.getFirst(tValue);
//     while ( bFound )
//     {
//         // code to process tValue goes here...
//         bFound = kList.getNext(tValue);
//     }
//
// To iterate over a list to access values with an explicit iterator, use the
// following.  This format does support nested loops with the same list.
//
//     GbTList<T> kList;  // already built list
//     T tValue0, tValue1;
//     GbTList<T>::Iterator kIter0 = kList.getFirst();
//     for (kI0 = kList.getFirst(); kI0; kI0 = kList.getNext(kI0))
//     {
//         tValue0 = kList.getValue(kI0);
//         for (kI1 = kList.getFirst(); kI1; kI1 = kList.getNext(kI1))
//         {
//             tValue1 = kList.getValue(kI1);
//             // code to process pair <tValue0,tValue1> goes here
//         }
//     }

template <class T>
class GbTList
{
public:
  // construction and destruction
  GbTList ();
  ~GbTList ();

  // element access
  INLINE unsigned int getQuantity () const;
  INLINE void add (T tValue);
  INLINE GbBool addUnique (T tValue);  // adds only not already in list
  INLINE GbBool remove (T tValue);
  INLINE GbBool removeFront (T& rtValue);
  INLINE void removeAll ();

  // linear traversal of list (implicit iterator)
  INLINE GbBool getFirst (T& rtValue) const;
  INLINE GbBool getNext (T& rtValue) const;

  // linear traversal of list (explicit iterator)
  typedef const void* Iterator;
  INLINE Iterator getFirst () const;
  INLINE Iterator getNext (Iterator pkIterator) const;
  INLINE T& getValue (Iterator pkIterator) const;

  // list operations
  INLINE void reverseOrder ();
  INLINE GbBool contains (T tValue) const;

protected:
  class Node {
  public:
    Node (T tValue, Node* pkNext) {
      value_ = tValue;
      next_ = pkNext;
    }

    // GbTList is responsible for deleting the nodes properly

    T value_;
    Node* next_;
  };

  unsigned int quantity_;
  Node* front_;

  // iterator for traversal
  mutable Node* iterator_;
};


#ifdef __GNUG__

#include "GbTList.in"
#include "GbTList.T"

#else

// #pragma instantiate GbTList<float>
// #pragma instantiate GbTList<double>
// #pragma instantiate GbTList<GbBool>
#pragma instantiate GbTList<unsigned int>
class BoundaryList;
#pragma instantiate GbTList<BoundaryList*>
// #pragma instantiate GbTList<int>

#ifndef OUTLINE
#include "GbTList.in"
#endif

#endif  // g++

#endif  // GBTLIST_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/04/23 10:02:23  prkipfer
| introduced morphological image operations and TList class
|
|
+---------------------------------------------------------------------*/
