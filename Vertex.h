//  -*- C++ -*-
/***************************************************************************
    File        : Vertex.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Thu Dec 6 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef VERTEX_H
#define VERTEX_H


// project files
#include "Pool.h"

namespace grd {

class Vertex  {
public:
  // Constructors and Destructior
  Vertex() {
    tag = 0;
    id  = -1;
    pos[0] = 0.0;
    pos[1] = 0.0;
    pos[2] = 0.0;
  }
  Vertex(double x[3]) {
    tag = 0;
    id  = -1;
    pos[0] = x[0];
    pos[1] = x[1];
    pos[2] = x[2];
  }

  virtual ~Vertex() { }

  // Numbering
  inline int  getId();
  inline void setId(int Id);

  // Geometry
  inline double& operator [] (int);
  inline const double&  operator[] (int i) const;

  inline void    setPosition(double* p);
  inline double* getPosition() { return pos; }

  // Marks manipulation
  inline void ref();
  inline void unref();
  inline bool queryRef();

  // Memory management
  void * operator new(size_t) { return pool.alloc(); }
  void   operator delete(void *p) { pool.free(p); }

private:
  // numbering
  int  id;
  // database manipulation
  bool tag;
  // geometry
  double pos[3];

  // Memory management
  static Pool<Vertex> pool;
};


// ____________________________________________________________
//
//	
//	         inline functions
//
// ____________________________________________________________
//


// Description
//
inline int
Vertex::getId()
{
  return id;
}


// Description
//
inline void
Vertex::setId(int Id)
{
  id = Id;
}

// Description:
//
inline double&
Vertex::operator[](int i)
{
  return pos[i];
}


// Description:
//
inline const double&
Vertex::operator[](int i) const
{
  return pos[i];
}


// Description
//
inline void
Vertex::setPosition(double* p)
{
  pos[0] = p[0];
  pos[1] = p[1];
  pos[2] = p[2];
}



// Description
//
inline void
Vertex::ref()
{
  tag = true;
}


// Description
//
inline void
Vertex::unref()
{
  tag = false;
}


// Description
//
inline bool
Vertex::queryRef()
{
  return tag;
}


} // namespace grd

#endif // VERTEX_H