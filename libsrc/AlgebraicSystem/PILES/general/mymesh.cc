#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "mydefs.hh"
#include "output.hh"
#include "basemesh.hh"
#include "mymesh.hh"

namespace CoupledField
{

MyMesh :: MyMesh(char * ainputfile)
  : BaseMesh()
{
  filename = new char[100];

  strcpy(filename,ainputfile);

#ifdef TRACE
  (*trace) << "entering MyMesh::MyMesh" << endl;
#endif

#ifdef DEBUG
  (*debug) << "opening the file" << endl;
#endif
}
  
MyMesh :: ~MyMesh()
{
#ifdef TRACE
  (*trace) << "entering MyMesh::~MyMesh" << endl;
#endif

  delete filename;
}

void MyMesh :: Read()
{
#ifdef TRACE
  (*trace) << "entering MyMesh::Read" << endl;
#endif

  ifstream infile;

  infile.open(strcat(filename,".dat"));

  // sizes

  infile.ignore(100,'/');

  infile >> numnode;
  infile.ignore(100,'/'); 
  infile >> numelem;
  infile.ignore(100,'/');
  infile >> numdir;
  infile.ignore(100,'/');
  infile >> nummat;
  infile.ignore(100,'/');
  infile >> dim;
  infile.ignore(100,'/');
  infile >> elemsize;

  coord   = new Double[dim*numnode];
  bval    = new Double[numdir];
  connect = new Integer[elemsize*numelem];
  mat     = new Integer[numelem];
  bnode   = new Integer[numdir];

  // connectivity

  infile.ignore(100,'/');

  Integer i,j,elnum;
  
  for (i=0; i<numelem; i++)
    {
      infile >> elnum;
      infile >> mat[elnum-1];

      for (j=0; j<elemsize; j++)
	{
	  infile >> connect[(elnum-1)*elemsize+j];
	}
    }

  // coordinates

  infile.ignore(100,'/');

  for (i=0; i<numnode; i++)
    {
      for (j=0; j<dim; j++)
	{
	  infile >> coord[i*dim+j];
	}
    }

  // dirichlet nodes

  infile.ignore(100,'/');

  for (i=0; i<numdir; i++)
    {
      infile >> bnode[i];
      infile >> bval[i];
    }

#ifdef DEBUG
  Print();
#endif
}

}
