// - C++ -
/***************************************************************************
    File        : ConformingClosure.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Thu Jan 17 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef CONFORMINGCLOSURE_H
#define CONFORMINGCLOSURE_H

// System
#include <vector>
// Grid
#include "Element.h"
#include "Triangle.h"
#include "Quadrangle.h"
#include "Tetrahedron.h"
#include "Octahedron.h"

namespace grd {

#define NO_TRIANGLES   10
#define NO_QUADRANGLES 10
#define NO_TETRAHEDRA  40
#define NO_OCTAHEDRA   40
#define NO_HEXAHEDRA   10
#define NO_PRISMS      10
#define NO_PYRAMIDS    10

class ConformingClosure {
public:
  typedef Triangle**    triangleIterator;
  typedef Quadrangle**  quadrangleIterator;
  typedef Tetrahedron** tetrahedronIterator;
  typedef Octahedron**  octahedronIterator;
//  typedef Hexahedron**  hexahedronIterator;
//  typedef Prism**       prismIterator;
//  typedef Pyramid**     piramidIterator;
  typedef triangleIterator    TRI;
  typedef quadrangleIterator  QUA;
  typedef tetrahedronIterator TET;
  typedef octahedronIterator  OCT;
//  typedef hexahedronIterator  HEX;
//  typedef prismIterator       PRI;
//  typedef piramidIterator     PYR;
public:
  ConformingClosure() {
    int i;

    for (i = 0; i < NO_TRIANGLES; i++)
      triangle[i] = new Triangle;
    for (i = 0; i < NO_QUADRANGLES; i++)
      quadrangle[i] = new Quadrangle;
    for (i = 0; i < NO_TETRAHEDRA; i++)
      tetrahedron[i] = new Tetrahedron;
    for (i = 0; i < NO_OCTAHEDRA; i++)
      octahedron[i] = new Octahedron;
//    for (i = 0; i < NO_HEXAHEDRA; i++)
//      hexahedron[i] = new Hexahedron;
//    for (i = 0; i < NO_PRISMS; i++)
//      prism[i] = new Prism;
//    for (i = 0; i < NO_PYRAMIDS; i++)
//      pyramid[i] = new Pyramid;
    // set pointers
    init();
  }

  // Desctructor
  ~ConformingClosure() {
    int i;
    for (i = 0; i < NO_TRIANGLES; i++)
      delete triangle[i];
    for (i = 0; i < NO_QUADRANGLES; i++)
      delete quadrangle[i];
    for (i = 0; i < NO_TETRAHEDRA; i++)
      delete tetrahedron[i];
    for (i = 0; i < NO_OCTAHEDRA; i++)
      delete octahedron[i];
//    for (i = 0; i < NO_HEXAHEDRA; i++)
//      delete hexahedron[i];
//    for (i = 0; i < NO_PRIMS; i++)
//      delete prism[i];
//    for (i = 0; i < NO_PYRAMIDS; i++)
//      delete pyramid[i];
  }

  void init() {
    endTri = triangle;
    endQua = quadrangle;
    endTet = tetrahedron;
    endOct = octahedron;
  }



  void setNoOfTriangles(int noOfTri) {
    endTri = triangle + noOfTri;
  }
  void setNoOfQuadrangles(int noOfQua) {
    endQua = quadrangle + noOfQua;
  }
  void setNoOfTetrahedra(int noOfTet) {
    endTet = tetrahedron + noOfTet;
  }
  void setNoOfOctahedra(int noOfOct) {
    endOct = octahedron + noOfOct;
  }

  TRI beginTriangle()    { return triangle;    }
  QUA beginQuadrangle()  { return quadrangle;  }
  TET beginTetrahedron() { return tetrahedron; }
  OCT beginOctahedron()  { return octahedron;  }

  TRI endTriangle()    { return endTri; }
  QUA endQuadrangle()  { return endQua; }
  TET endTetrahedron() { return endTet; }
  OCT endOctahedron()  { return endOct; }


private:
  Triangle*    triangle[NO_TRIANGLES];
  Quadrangle*  quadrangle[NO_QUADRANGLES];
  Tetrahedron* tetrahedron[NO_TETRAHEDRA];
  Octahedron*  octahedron[NO_OCTAHEDRA];

  TRI endTri;
  QUA endQua;
  TET endTet;
  OCT endOct;
};


} // namespace grd
#endif // CONFORMINGCLOSURE_H