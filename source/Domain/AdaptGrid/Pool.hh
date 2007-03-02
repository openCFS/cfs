/***************************************************************************
    File        : Pool.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Thu Dec 6 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef POOL_H
#define POOL_H

#include <malloc.h>

namespace grd {

template<class T>
class Pool
{
private:
  struct Link { Link* next; };

  Link* head;
  const unsigned esize;

  // copy function
  void copy(const Pool &p) {
    if (this == &p) return;
    head  = p.head;
    esize = p.esize;
  }

  // dynamic memory allocation
  void grow()
  {
    // allocate memory for 1000 objects of type T
    char* start = ::new char[200*esize];
    char* last  = &start[(199)*esize];

    for (char* p = start; p < last; p += esize)
      ((Link*) p)->next = (Link*) (p+esize);

    ((Link*) last)->next = 0;
    head = (Link*) start;
  }

public:
  // Constructors
  Pool() : esize(sizeof(T) < sizeof(Link*) ? sizeof(Link*):sizeof(T))
  { head = 0; }

  Pool(const Pool &p) { copy(p); }
  // Destructor
  virtual ~Pool() { }

  Pool& operator=(const Pool& p) { copy(p); return *this; }

  // Member Functions

  void* alloc()
  {
    if (head == 0) grow();
    Link* p = head; head = p->next;
    return p;
  }

  void free(void* b)
  {
    Link* p = (Link*) b;
    p->next = head;
    head = p;
  }

};

} // namespace grd

#endif // POOL_H
