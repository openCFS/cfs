#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>

namespace CoupledField
{

BEMMesh :: BEMMesh(char * ainputfile)
  : BaseMesh()
{
  filename = new char[100];

  strcpy(filename,ainputfile);

#ifdef TRACE
  (*trace) << "entering BEMMesh::BEMMesh" << endl;
#endif

#ifdef DEBUG
  (*debug) << "opening the file" << endl;
#endif
}
  
BEMMesh :: ~BEMMesh()
{
#ifdef TRACE
  (*trace) << "entering BEMMesh::~BEMMesh" << endl;
#endif

  delete filename;

  delete [] coord;
  delete [] bval;
  delete [] connect;
  delete [] mat;
  delete [] bnode;

}

void BEMMesh :: Read()
{
#ifdef TRACE
  (*trace) << "entering BEMMesh::Read" << endl;
#endif

  ifstream infile;

  infile.open(strcat(filename,".dat"));

  // sizes

  Integer nx,ny,nd;

  infile.ignore(100,'/');
  infile >> nx;
  infile.ignore(100,'/'); 
  infile >> ny;

  numnode = 2*(2*nx + 2*ny - 4);
  numelem = numnode;

  infile.ignore(100,'/');
  infile >> nd;

  if (nd == 1) /// dirichlet
    {
      numdir = numnode;
      numneu = 0;
    }
  else if (nd == 2) ///neumann
    {
      numdir = 0;
      numneu = numnode;
    }

  infile.ignore(100,'/');
  infile >> nummat;
  infile.ignore(100,'/');
  infile >> dim;
  infile.ignore(100,'/');
  infile >> elemsize;

  coord   = new Double[dim*numnode];
  bval    = new Double[numdir+numneu];
  connect = new Integer[elemsize*numelem];
  mat     = new Integer[numelem];
  bnode   = new Integer[numdir+numneu];

  // connectivity

  Double hx,hy,x1,x2,y1,y2,y3,y4;
  Integer i,j,k,offset,p;
  
  offset = dim*(2*nx+2*ny-4);
  p      = 2*nx+2*ny-4;

  for (i=0; i<numelem; i++)
    {
      mat[i] = 1;
    }

  x1 = 0;
  x2 = 1;
  y1 = 0;
  y2 = 0.2;
  y3 = 0.3;
  y4 = 0.5;

  hx = (x2-x1)/(nx-1);
  hy = (y2-y1)/(ny-1);
  j  = 0;

  for (i=1; i<nx; i++)
    {
      connect[j]   = i;
      connect[j+1] = i+1;

      connect[offset+j]   = i+p;
      connect[offset+j+1] = i+p+1;

      j += dim;
    }

  k = nx-1;

  for (i=1; i<ny; i++)
    {
      connect[j]   = i+k;
      connect[j+1] = i+k+1;

      connect[offset+j]   = p+i+k;
      connect[offset+j+1] = p+i+k+1;

      j += dim;
    }

  k = nx+ny-2;
  
  for (i=1; i<nx; i++)
    {
      connect[j]   = i+k;
      connect[j+1] = i+k+1;

      connect[offset+j]   = p+i+k;
      connect[offset+j+1] = p+i+k+1;

      j += dim;
    }

  k = 2*nx+ny-3;
  
  for (i=1; i<ny-1; i++)
    {
      connect[j]   = i+k;
      connect[j+1] = i+k+1;

      connect[offset+j]   = p+i+k;
      connect[offset+j+1] = p+i+k+1;

      j += dim;
    }
  
  connect[j]   = p;
  connect[j+1] = 1;

  connect[offset+j]   = p+p;
  connect[offset+j+1] = p+1;

  /// coordinate
  
  j = 0;

  for (i=1; i<nx; i++)
    {
      coord[j]   = x1+(i-1)*hx;
      coord[j+1] = y1;

      coord[offset+j]   = x1+(i-1)*hx;
      coord[offset+j+1] = y3;

      j += dim;
    }

  for (i=1; i<ny; i++)
    {
      coord[j]   = x2;
      coord[j+1] = y1+(i-1)*hy;

      coord[offset+j]   = x2;
      coord[offset+j+1] = y3+(i-1)*hy;

      j += dim;
    }
  
  for (i=1; i<nx; i++)
    {
      coord[j]   = x2-(i-1)*hx;
      coord[j+1] = y2;

      coord[offset+j]   = x2-(i-1)*hx;
      coord[offset+j+1] = y4;

      j += dim;
    }
  
  for (i=1; i<ny; i++)
    {
      coord[j]   = x1;
      coord[j+1] = y2-(i-1)*hy;

      coord[offset+j]   = x1;
      coord[offset+j+1] = y4-(i-1)*hy;

      j += dim;
    }

  // dirichlet nodes

  j = numdir/2;

  for (i=0; i<j; i++)
    {
      bnode[i] = i+1;
      bval[i]  = 1;
    }

  for (i=j; i<numdir; i++)
    {
      bnode[i] = i+1;
      bval[i]  = 0;
    }

  //ConstructGraph();

#ifdef DEBUG
  (*debug) << "printing the mesh" << endl;
  Print();
#endif
}

}
