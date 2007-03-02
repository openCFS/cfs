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
#include "Pool.hh"

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

  // Numerics
  bool isBoundaryNode();
  void setBoundaryNode();
  void setInnerNode();

  void setColor(int bcol);
  int  color();

  // Memory management
  //void * operator new(size_t) { return pool.alloc(); }
  //void   operator delete(void *p) { pool.free(p); }

private:
  // numbering
  int  id;

  // Vertex marks
  unsigned char tag;

  // geometry
  double pos[3];

  // Memory management
  //static Pool<Vertex> pool;
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


//
// Vertex marks: one byte
//
//  1.- 6. bit: color number
//      7. bit: boundary mark:
//           1 is a boundary node
//           0 is an inner node
//      8. bit: reference or general mark
//           1 is referenced or marked
//           0 is not referenced or not marked
//
const int VERTEX_MARK_BOUNDARY = 64;
const int VERTEX_MARK_GENERAL  = 128;
const int VERTEX_MASK_LASTTWOBITS  = 192;
const int VERTEX_MASK_FIRSTSIXBIT  = 127;

// Description
// Mark last bit for general
// vertex mark
inline void
Vertex::ref()
{
  tag |= VERTEX_MARK_GENERAL;
}


// Description
// Unset last bit for general
// vertex mark
inline void
Vertex::unref()
{
  tag ^=VERTEX_MARK_GENERAL;
}


// Description
// Check last bit if
// it is set
inline bool
Vertex::queryRef()
{
  return (tag&VERTEX_MARK_GENERAL);
}


// Description
//  Boundary is the 7. bit
inline bool
Vertex::isBoundaryNode()
{
  return (tag&VERTEX_MARK_BOUNDARY);
}

// Description
inline void
Vertex::setBoundaryNode()
{
  tag |= VERTEX_MARK_BOUNDARY;
}

// Description
//
inline void
Vertex::setInnerNode()
{
  if (tag&VERTEX_MARK_BOUNDARY)
    tag ^= VERTEX_MARK_BOUNDARY;
}

// Description
//  Color is an int given by the first 6 bits
inline void 
Vertex::setColor(int bcol)
{
  // set first 6 bits to zero
  tag &= VERTEX_MASK_LASTTWOBITS;
  // set node color
  tag |= bcol;
}

// Description
//  read color from the first 6 bits
inline int  
Vertex::color()
{
  unsigned char res = tag;
  res &=VERTEX_MASK_FIRSTSIXBIT;
  return res;
}

} // namespace grd

#endif // VERTEX_H
