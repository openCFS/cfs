// - C++ -
/***************************************************************************
    File        : TetrahedronMeasure.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Tue Sep 10 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/


#ifndef TETRAHEDRONMEASURE_H
#define TETRAHEDRONMEASURE_H

// System
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

// Grid
#include "Vertex.h"
#include "Edge.h"
#include "Triangle.h"
#include "TetrahedronGeometry.h"
#include "Tetrahedron.h"

using namespace std;
namespace grd {

class TetrahedronMeasure {
public:
  TetrahedronMeasure() { }
  ~TetrahedronMeasure () { }

  // Method
  void operator() (Element* t)
  {
    int type = t->type();
    if(type == GRD_TETRAHEDRON)
    {
      getTetMeasure(t);
    }
  }

  const double& operator [] (int i) const { return measure[i]; }
  double& operator [] (int i) { return measure[i]; }

  void init() { measure.clear(); }

  inline double getArea(int face,Element* t);
  inline void   print(char* filename);

private:
  vector<double> measure;
  inline void getTetMeasure(Element* t);
};



//
//
inline void
TetrahedronMeasure::getTetMeasure(Element* t)
{
  double area;
  double volume;
  double edgeSize;
  Jacobian jac;
  Edge* eg = new Edge;
  TetrahedronGeometry tetGeo;

  // Get Area
  area = 0.0;
  for (int face = 0; face < 4; face++)
  {
    area += getArea(face,t);
  }

  // Get volume
  tetGeo.getJacobian(t,jac);
  volume = jac.getDeterminant();
  volume = fabs(volume);

  // get largest edge
  edgeSize = 0.0;
  for (int edge = 0; edge < 6; edge++)
  {
    eg->setVertex(0,t->getVertexAtEdge(edge,0));
    eg->setVertex(1,t->getVertexAtEdge(edge,1));
    double val = eg->length();
    edgeSize = (edgeSize > val) ? edgeSize : val;
  }

  double res = edgeSize*area / volume;
  measure.push_back(res);

}

//
//
inline double
TetrahedronMeasure::getArea(int face,Element* t)
{
  Vertex* v0 = t->getVertexAtFace(face,0);
  Vertex* v1 = t->getVertexAtFace(face,1);
  Vertex* v2 = t->getVertexAtFace(face,2);

  double x1[3];
  double x2[3];

  for (int i = 0; i < 3; i++)
  {
    x1[i] = (*v1)[i] - (*v0)[i];
    x2[i] = (*v2)[i] - (*v0)[i];
  }

  double normal[3];
  grd::crossProduct(normal,x1,x2);
  double area = fabs(grd::length(normal)) * 0.5;

  return area;
}

//
//
inline void
TetrahedronMeasure::print(char* filename)
{
  ofstream theFile;

  theFile.open(filename);

  // initialize
  for (uint i = 0; i < measure.size(); i++)
  {
    theFile << i << "   " << (measure[i]) << '\n';
  }

  // open parameter file
  theFile.close();
}

} // namespace grd
#endif // TETRAHEDRONMEASURE_H