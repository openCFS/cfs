#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>

namespace CoupledField
{

NeuroFEMMesh :: NeuroFEMMesh(char * ainputfile)
  : BaseMesh()
{
  filename = new char[100];
  pathname = new char[100];

  strcpy(pathname,ainputfile);

#ifdef TRACE
  (*trace) << "entering NeuroFEMMesh::NeuroFEMMesh" << endl;
#endif

#ifdef DEBUG
  (*debug) << "opening the file" << endl;
#endif
}
  
NeuroFEMMesh :: ~NeuroFEMMesh()
{
#ifdef TRACE
  (*trace) << "entering NeuroFEMMesh::~NeuroFEMMesh" << endl;
#endif

  delete filename;

  delete [] coord;
  delete [] bval;
  delete [] connect;
  delete [] mat;
  delete [] bnode;

}

void NeuroFEMMesh :: Read()
{
#ifdef TRACE
  (*trace) << "entering NeuroFEMMesh::Read" << endl;
#endif

  cout << "reading NeuroFEM ..." << endl;

  ifstream infile;

  strcpy(filename,pathname);
  infile.open(strcat(filename,".geo"));

  if (!infile.good())
    {
      cout << "shit" << endl;
    }

  // sizes

  infile >> numnode;
  infile >> numelem;
  infile >> elemsize;
  infile >> dim;

  numdir = 1;
  nummat = 1;

  coord   = new Double[dim*numnode];
  bval    = new Double[numdir];
  connect = new Integer[elemsize*numelem];
  mat     = new Integer[numelem];
  bnode   = new Integer[numdir];

  // connectivity

  Integer i,j;
  Char s[20];
  
  // coordinates

  for (i=0; i< numnode; i++)
    {
      for (j=0; j<dim; j++)
	{
	  infile >> coord[i*dim+j];
	}
    }


  for (i=0; i<numelem; i++)
    {
      infile >> s;
      infile.get();

      for (j=0; j<elemsize; j++)
	{
	  infile.get(s,7);

	  connect[i*elemsize+j] = atoi(s);
	}
    }

  infile.close();

  strcpy(filename,pathname);
  infile.open(strcat(filename,".dir"));
  // dirichlet nodes

  for (i=0; i<numdir; i++)
    {
      infile >> bnode[i];
      infile >> bval[i];
    }

  infile.close();

  for (i=0;i<numelem; i++)
    {
      mat[i] = 1; 
    }

//   strcpy(filename,pathname);
//   infile.open(strcat(filename,".mat"));

//   Double h;

//   for (i=0; i<numelem; i++)
//     {
//       infile >> j;
//       infile >> s;
//       infile >> h;

//       if (h == 0.33)
// 	{
// 	  mat[i] = 1;
// 	}
//       else if (h == 0.0042)
// 	{
// 	  mat[i] = 2;
// 	}
//       else 
// 	{
// 	  cout << h << endl;
// 	  mat[i] = 3;
// 	}
//     }

//   infile.close();

  cout << "... finished" << endl;

  // calculate the graph

  //ConstructGraph();

#ifdef DEBUG
  //(*debug) << "printing the mesh" << endl;
  //Print();
#endif
}

}
