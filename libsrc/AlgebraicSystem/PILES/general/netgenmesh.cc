#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>

namespace CoupledField
{

NETGENMesh :: NETGENMesh(char * ainputfile)
  : BaseMesh()
{
  filename = new char[100];
  pathname = new char[100];

  strcpy(pathname,ainputfile);

#ifdef TRACE
  (*trace) << "entering NETGENMesh::NETGENMesh" << endl;
#endif

#ifdef DEBUG
  (*debug) << "opening the file" << endl;
#endif
}
  
NETGENMesh :: ~NETGENMesh()
{
#ifdef TRACE
  (*trace) << "entering NETGENMesh::~NETGENMesh" << endl;
#endif

  delete filename;

  delete [] coord;
  delete [] bval;
  delete [] connect;
  delete [] mat;
  delete [] bnode;

}

void NETGENMesh :: Read()
{
#ifdef TRACE
  (*trace) << "entering NETGENMesh::Read" << endl;
#endif

  cout << "reading NETGEN ..." << endl;

  ifstream infile;

  strcpy(filename,pathname);
  infile.open(strcat(filename,".vol"));

  if (!infile.good())
    {
      cout << "shit" << endl;
    }

  Char * s = new Char[100];

  nummat = 1;

  infile >> s;
  infile >> s;
  infile >> dim;

  int i,j,k,num;

  for (i=0; i<10; i++)
    {
      infile >> s;
    }

  infile >> num;

  for (i=0; i<num; i++)
    {
      for (j=0; j<11; j++)
	{
	  infile >> k;
	}
    }

  numdir = 1;

  cout << "numdir " << numdir << endl;

  bval    = new Double[numdir];
  bnode   = new Integer[numdir]; 

  bval[0]  = 0;
  bnode[0] = 1;

  for (i=0; i<8; i++)
    {
      infile >> s;
    }

  infile >> numelem;

  cout << "numelem " << numelem << endl;

  elemsize = 4;

  connect = new Integer[elemsize*numelem];
  mat     = new Integer[numelem];

  for (i=0; i<numelem; i++)
    {
      infile >> mat[i]; 
      infile >> elemsize;

      for (j=0; j<elemsize; j++)
	{
	  infile >> k;
	  connect[i*elemsize+j] = k;
	}
    }

  for (i=0; i<5; i++)
    {
      infile >> s;
    }

  infile >> num;

  num *= 2;

  for (i=0; i<=num; i++)
    {
      infile.ignore(32000,'\n');
    }

  for (i=0; i<5; i++)
    {
      infile >> s;
    } 

  infile >> numnode;

  cout << "numnode " << numnode << endl;

  coord   = new Double[dim*numnode];

  for (i=0; i<numnode; i++)
    {
      for (j=0; j<dim; j++)
	{
	  infile >> coord[i*dim+j];
	}
    }

  infile.close();

  cout << "... finished" << endl;


#ifdef DEBUG
  //(*debug) << "printing the mesh" << endl;
  //Print();
#endif
}

}
