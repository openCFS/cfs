#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>

namespace CoupledField
{

SquareMesh :: SquareMesh(char * ainputfile)
  : BaseMesh()
{
  filename = new char[100];

  strcpy(filename,ainputfile);

#ifdef TRACE
  (*trace) << "entering SquareMesh::SquareMesh" << endl;
#endif

#ifdef DEBUG
  (*debug) << "opening the file" << endl;
#endif
}
  
SquareMesh :: ~SquareMesh()
{
#ifdef TRACE
  (*trace) << "entering SquareMesh::~SquareMesh" << endl;
#endif

  delete [] filename;

  delete [] coord;
  delete [] bval;
  delete [] connect;
  delete [] mat;
  delete [] bnode;

}

void SquareMesh :: Read()
{
#ifdef TRACE
  (*trace) << "entering SquareMesh::Read" << endl;
#endif

  ifstream infile;

  infile.open(strcat(filename,".dat"));

  // sizes

  Integer nx,ny,nd;

  infile.ignore(100,'/');
  infile >> nx;
  infile.ignore(100,'/'); 
  infile >> ny;

  numnode = nx*ny;
  numelem = (nx-1)*(ny-1);

  infile.ignore(100,'/');
  infile >> nd;

  numdir = nx;

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

  Integer i,j,k,elnum;
  
  for (i=0; i<numelem; i++)
    {
      mat[i] = 1;
    }

  k     = 0;
  elnum = 1;
  
  for (i=1; i<=ny-1; i++)
    {
      k++;

      for (j=1; j<=nx-1; j++)
	{
	  connect[(elnum-1)*elemsize+0] = k;
	  connect[(elnum-1)*elemsize+1] = k+1;
	  connect[(elnum-1)*elemsize+2] = nx+k;
	  connect[(elnum-1)*elemsize+3] = nx+k+1;
	  
	  k++;
	  elnum++;
	}
    }

  // coordinates

  Double hx,hy;

  hx = 1./(nx-1);
  hy = 1./(ny-1);
  k  = 0;

  for (i=1; i<=ny; i++)
    {
      for (j=1; j<=nx; j++)
	{
	  coord[k*dim+0] = (j-1)*hx;
	  coord[k*dim+1] = (i-1)*hy;
	  k++;
	}
    }

  // dirichlet nodes

  for (i=0; i<numdir; i++)
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
