//  -*- C++ -*-
/***************************************************************************
    File        : main.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Tue Aug 13 17:11:32 CEST 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>


// Inventor
#include <Inventor/Xt/SoXt.h>
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoNormal.h>
#include <Inventor/nodes/SoNormalBinding.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/SoInput.h>
#include <Inventor/SoDB.h>


// Grid
#include "Chronometer.h"
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Tetrahedron.h"
#include "Octahedron.h"
#include "MultilevelGrid.h"
#include "GeometrySensor.h"
#include "MeshReader.h"
#include "TetrahedronMeasure.h"
#include "MeshWriter.h"

class MeshRefine {
public:
  MeshRefine() {}
  ~MeshRefine() {}
  void operator ()(grd::Element* t)
  {
    t->markForRefinement();
  }
};


class CheckSkeleton {
public:
  CheckSkeleton () {}
  ~CheckSkeleton () {}
  // Method
  void operator() (grd::Element* t)
  {
    grd::TriSkeleton* triSk = ((grd::Triangle*)t)->getSkeleton();
    if (triSk == 0)
      cerr << "ERROR skeleton not defined\n";
  }

};

class CountElement {
public:
  CountElement () { counter = 0; }
  ~CountElement () {}

  // Method
  void operator() (grd::Element* t)
  {
    int etype = t->type();

    if (etype == GRD_TRIANGLE){
      counter++;
    }
    else {
      int noOfFaces = t->getNoOfFaces();
      grd::Element* p = t->getParent();
      for (int i = 0; i < noOfFaces; i++)
      {
        grd::Element* n = t->getNeighbor(i);
        if (n)
        {
          int type = n->type();
          if (type == GRD_TRIANGLE)
            counter++;
        }
      }
    }
  }
  int getCounter() { return counter; }
  void init() { counter = 0; }
private:
  int counter;
};

int main(int argc, char *argv[])
{
  std::string meshFile = "/home/grosso/Projects/Grid/Data/head3.vol";
  std::string srgFile  = "/home/grosso/Projects/Grid/Data/head3.srg";
  grd::MeshReader mrBnd;
  grd::MeshReader mrVol;
  grd::MultilevelGrid grid;
  grd::GeometrySensor geomSensor;
  grd::Chronometer crono;
  grd::TetrahedronMeasure tetMeasure;
  CountElement cntElm;
  int noOfElem;

  // Read Corase mesh
  cerr << "Read meshes from file\n";
  mrBnd.readSRG(srgFile);
  mrVol.readGRUMMP(meshFile);

  // Set coarse mesh
  crono.start();
  grid.buildCoarseMesh(mrVol.getVertexList(),mrVol.getElementList());
  crono.stop();
  crono.print("Build coarse mesh");
  //grid.buildCoarseMesh(mrBnd.getVertexList(),mrBnd.getElementList());
  crono.start();
  grid.setBoundaryMesh(mrBnd.getElementList());
  crono.stop();
  crono.print("Set Boundary hierarchy mesh");

  //grid.buildMultilevelBoundary(mrBnd.getVertexStack(),mrBnd.getElementStack());
  CheckSkeleton chksk;
  grid.forEachBoundaryElement(chksk);
  grid.setVertexNumbers();
  cerr << "Mesh finished\n";


  // Check Mesh quality
  tetMeasure.init();
  grid.forEachVolumeElement(tetMeasure);
  tetMeasure.print("CoarseMesh.txt");

  // Count no. of elements
  cntElm.init();
  grid.forEachVolumeElement(cntElm);
  noOfElem = cntElm.getCounter();
  cerr << "No. of volume elements: " << noOfElem << '\n';
  cntElm.init();
  grid.forEachUnrefinedBoundaryElement(cntElm);
  noOfElem = cntElm.getCounter();
  cerr << "No. of surface elements: " << noOfElem << '\n';

  // Refine Mesh according to geometry
  // Refine 1.
  crono.start();
  geomSensor.markForRefinement(grid);
  crono.stop();
  crono.print("Mark for refinement");
  crono.start();
  grid.refine();
  crono.stop();
  crono.print("Mesh refinement");
  grid.setVertexNumbers();
  cntElm.init();
  grid.forEachUnrefinedVolumeElement(cntElm);
  noOfElem = cntElm.getCounter();
  cerr << "No. of volume elements: " << noOfElem << '\n';
  cntElm.init();
  grid.forEachUnrefinedBoundaryElement(cntElm);
  noOfElem = cntElm.getCounter();
  cerr << "No. of surface elements: " << noOfElem << '\n';

  // Check Mesh quality
  tetMeasure.init();
  grid.forEachVolumeElement(tetMeasure);
  tetMeasure.print("FirstMesh.txt");

  return 1;
}

