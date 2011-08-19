// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// - C++ -
/***************************************************************************
    File        : MeshWriter.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Apr 26 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef MESHWRITER_H
#define MESHWRITER_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

// Utilities
#include "Utility.hh"

// Grid
#include "Vertex.hh"
#include "Element.hh"
#include "Tetrahedron.hh"
#include "Octahedron.hh"
#include "MultilevelGrid.hh"

using namespace std;
namespace grd {

class MeshWriter {
public:
  MeshWriter () {}
  ~MeshWriter () {}

  inline void writeMesh(char* filename, MultilevelGrid& grid);

};



class Zaehler {
public:
  Zaehler () { counter = 0; }
  ~Zaehler () {}

  // Method
  void operator() (grd::Element* t)
  {
    counter++;
  }
  int  getCounter() { return counter; }
  void init() { counter = 0; }
private:
  int counter;
};

class WriteElement {
public:
  WriteElement(ofstream* fl) {
    theFile = fl;
    tets.resize(8);
    for (int l = 0; l < 8; l++)
    {
      tets[l] = (Element*) new Tetrahedron;
    }
  }
  ~WriteElement() {}

  void operator () (Element* t)
  {
    int type = t->type();
    if (type == GRD_TETRAHEDRON) {
      writeTet(t);
    }
    else if (type == GRD_OCTAHEDRON) {
      ((Octahedron*)t)->getTetras(tets);
      for (int l = 0; l < 8; l++)
      {
        writeTet(tets[l]);
      }
    }
  }

  void writeTet(Element* t)
  {
    int vt0,vt1,vt2,vt3;
    vt0 = (t->getVertex(0))->getId() + 1;
    vt1 = (t->getVertex(1))->getId() + 1;
    vt2 = (t->getVertex(2))->getId() + 1;
    vt3 = (t->getVertex(3))->getId() + 1;
    (*theFile) << "1 4"
               << " " << vt0
               << " " << vt1
               << " " << vt2
               << " " << vt3
               << '\n';
  }

private:
  ofstream* theFile;
  vector<Element*> tets;
};


class WriteVertex {
public:
  WriteVertex(ofstream* fl) {
    theFile = fl;
    counter = 0;
  }

  ~WriteVertex() {}

  //
  void operator() (Vertex* t)
  {
    int locId = t->getId();
    if (locId != counter)
      WARN("counting vertices");
    counter++;

    double x0,x1,x2;
    x0 = (*t)[0];
    x1 = (*t)[1];
    x2 = (*t)[2];
    (*theFile) << "    "
               << x0 << "   "
               << x1 << "   "
               << x2 << "   "
               << '\n';

  }
private:
  int counter;
  ofstream* theFile;
};

//
//
inline void
MeshWriter::writeMesh(char* filename,MultilevelGrid& grid)
{
  ofstream theFile;
  Zaehler zaehler;
  WriteElement writeElement(&theFile);
  WriteVertex writeVertex(&theFile);

  theFile.open(filename);
  // write header
  theFile << "mesh3d\n\n"
          << "surfaceelementsgi\n"
          << "620\n"
          << "1 1 1 0 3 1 2 3 1 1 1\n\n"
          << "volumeelements\n";

  grid.forEachElement(zaehler);
  int noOfElement = zaehler.getCounter();
  theFile << noOfElement << '\n';

  grid.forEachElement(writeElement);

  // write vertices
  theFile << "\npoints\n";
  int noOfVertex = grid.getNoOfVertices();
  theFile << noOfVertex << '\n';
  grid.forEachVertex(writeVertex);


}

} // namespace grd

#endif // MESHWRITER_H
