#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>

namespace CoupledField
{

BaseMesh :: BaseMesh()
{
#ifdef TRACE
  (*trace) << "entering BaseMesh::BaseMesh" << endl;
#endif

  offset = 10;
}
  
BaseMesh :: ~BaseMesh()
{
#ifdef TRACE
  (*trace) << "entering BaseMesh::~BaseMesh" << endl;
#endif

  delete [] rowsize;
  delete [] graph;
}

void BaseMesh :: Print() const
{
#ifdef TRACE
  (*trace) << "entering BaseMesh::Print" << endl;
#endif

  cout << "number of nodes      " << numnode << endl;
  cout << "number of elements   " << numelem << endl;
  cout << "number of dirichlet  " << numdir << endl;
  cout << "number of material   " << nummat << endl;
  cout << "number of nodes/elem " << elemsize << endl;
  cout << "dimension            " << dim << endl;

  Integer i,j;

  cout << "### connectivity" << endl;

  for (i=0; i<numelem; i++)
    {
      for (j=0; j<elemsize; j++)
	{
	  cout << connect[i*elemsize+j] << " ";
	}
      cout << endl;
    }

  cout << "### coordinate" << endl;

  for (i=0; i<numnode; i++)
    {
      for (j=0; j<dim; j++)
	{
	  cout << coord[i*dim+j] << " ";
	}
      cout << endl;
    }

  cout << "### dirichlet" << endl;

  for (i=0; i<numdir; i++)
    {
      cout << bnode[i] << " " << bval[i] << endl;
    }
  
  cout << "### end mesh" << endl;
}


void BaseMesh :: PrintMETIS() const
{
#ifdef TRACE
  (*trace) << "entering BaseMesh::PrintMETIS" << endl;
#endif

  
  Integer i,j;

  ofstream outfile;

  strcpy(filename,pathname);

  outfile.open(strcat(filename,".xyz"));

  cout << "+++ printing the coordinates for PMVIS: " << filename << endl;

  for (i=0; i<numnode; i++)
    {
      for (j=0; j<dim; j++)
	{
	  outfile << coord[i*dim+j] << " ";
	}
      outfile << endl;
    }

  outfile.close();

  /////////////////////////////////////////////////////

  strcpy(filename,pathname);

  outfile.open(strcat(filename,".cps"));

  cout << "+++ printing the connection for METIS: " << filename << endl;

  outfile << numelem << " " << 2 << endl;

  for (i=0; i<numelem; i++)
    {
      for (j=0; j<elemsize; j++)
	{
	  outfile << connect[i*elemsize+j] << " ";
	}
      outfile << endl;
    }

  outfile.close();

  /////////////////////////////////////////////////////

  strcpy(filename,pathname);

  outfile.open(strcat(filename,".con"));

  cout << "+++ printing the connection for PMVIS: " << filename << endl;

  for (i=0; i<numelem; i++)
    {
      for (j=0; j<elemsize; j++)
	{
	  outfile << connect[i*elemsize+j] << " ";
	}
      outfile << endl;
    }

  outfile.close();
}

void BaseMesh :: ConstructGraph()
{
  rowsize = new Integer[numnode];
  graph   = new Integer[numnode*offset];

  Integer i,j,k,l,p1,p2;
  Boolean set;

  for (i=0; i<numnode; i++)
    {
      for (j=0; j<offset; j++)
	{
	  graph[offset*i+j] = 0;
	}

      rowsize[i] = 0;
    }

  for (i=0; i<numelem; i++)
    {
      for (j=0; j<elemsize; j++)
	{
	  p1 = connect[i*elemsize+j];

	  for(k=0; k<elemsize; k++)
	    {
	      p2 = connect[i*elemsize+k];

	      set = TRUE;

	      for (l=0; l<rowsize[p1-1]; l++)
		{
		  if (graph[offset*(p1-1)+l] == p2)
		    {
		      set = FALSE;
		    }
		}

	      if (set == TRUE)
		{
		  rowsize[p1-1]++;

		  graph[offset*(p1-1)+rowsize[p1-1]-1] = p2;
		}
	    }
	}
    }

  // calc non-zero entries

  nne = 0;

  for (i=0; i<numnode; i++)
    {
      nne += rowsize[i];
    }

  // sort the graph: diagonal on first place, rest

  Integer help;

  for (i=1; i<=numnode; i++)
    {
      for (j=0; j<rowsize[i-1]; j++)
	{
	  if (graph[(i-1)*offset+j] == i)
	    {
	      help = j;
	    }
	}

      graph[(i-1)*offset+help] = graph[(i-1)*offset];
      graph[(i-1)*offset]      = i;

      for (j=1; j<rowsize[i-1]; j++)
	{
	  for (k=j+1; k<rowsize[i-1]; k++)
	    {
	      if (graph[(i-1)*offset+j] > graph[(i-1)*offset+k])
		{
		  help                  = graph[(i-1)*offset+j];
		  graph[(i-1)*offset+j] = graph[(i-1)*offset+k];
		  graph[(i-1)*offset+k] = help;
		}
	    }
	}
    }
  
#ifdef DEBUG

  for (i=1; i<=numnode; i++)
    {
      (*test) << "rowsize " << rowsize[i-1] << endl;

      for (j=0; j<rowsize[i-1]; j++)
	{
	  (*test) << graph[(i-1)*offset+j] << " ";
	}
      (*test) << endl;
    }

#endif
}

}
