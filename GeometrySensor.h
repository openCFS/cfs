// - C++ -
/***************************************************************************
    File        : GeometrySensor.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Wed Apr 24 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef GEOMETRYSENSOR_H
#define GEOMETRYSENSOR_H

#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Triangle.h"
#include "GridLevel.h"
#include "MultilevelGrid.h"
#include "TriSkeleton.h"

using namespace std;
namespace grd {

class GeometrySensor {
public:
  GeometrySensor() {}
  ~GeometrySensor() {}

  void markForRefinement(MultilevelGrid& grid);

};

} // namespace grd
#endif // GEOMETRYSENSOR_H