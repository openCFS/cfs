#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include <general.hh>

namespace CoupledField
{

PILES_GMV_Writer::PILES_GMV_Writer (BaseMesh * amesh, const char * name, 
				    int abinary,const char * afromfile)
  : GMV_Writer (name, abinary, afromfile)
{
  mesh = amesh;

  WriteFileHeader();
}


PILES_GMV_Writer :: ~PILES_GMV_Writer ()
{
}


void PILES_GMV_Writer :: Write (double * sol)
{
  int i;

  // variable mode on ?
  if (! variablemode)
    {
      (*gmvfile) << "variable";
      if (! binary) 
	(*gmvfile) << endl << endl;

      variablemode = 1;
    }
     
  // nodal data - scalar : other stuff see wolfram's feppgmvwriter.*

  // transform name to 8 chars
  char funname[9];
  funname[8] = '\0';
  strncpy(funname, "scalar", 8);
  
  if (binary)
    for (i = 7; i >= 0; i--)
      if (funname[i] ==  '\0')
	funname[i] = ' ';
  
  (*gmvfile) << funname;
  if (! binary)
    (*gmvfile) << "   ";
  WriteLong (1);                // write node data
  if (! binary)
    (*gmvfile) << endl;
  
  
  // allocate temporary variables
  float * values = new float[nrnodes];

  for (i=0; i <= nrnodes; i++)
    {
      values[i] = sol[i];
    }

  // write values to file
  WriteFloat (nrnodes, values);

  delete [] values;
}

void PILES_GMV_Writer :: WriteNodes ()
{
  int i,j;

  // write nodes token
  (*gmvfile) << "nodes   ";
  WriteLong (nrnodes);
  if (! binary)
    (*gmvfile) << endl;

  // write coordinates of nodes
  float * coords = new float[nrnodes];

  int dim = mesh->GetDim();
  double * point;

  for (j=0; j<dim; j++)
    {
      for (i=1; i<=nrnodes; i++)
	{
	  point = mesh->GetCoord(i);
	  
	  coords[i-1] = point[j];
	}

      WriteFloat (nrnodes, coords);
      
      if (! binary)
	{
	  (*gmvfile) << endl;
	}
    }

  // free coords
  delete [] coords;
}



void PILES_GMV_Writer :: WriteCells ()
{
  int i, nrelem, elemsize, dim;

  // write cell token
  (*gmvfile) << "cells   ";
  
  nrelem   = mesh->GetNumElem();
  elemsize = mesh->GetElemSize();
  dim      = mesh->GetDim();

  WriteLong (nrcells);
  if (! binary)
    (*gmvfile) << endl;

  // write cells
  int * connect;

  for (i = 1; i <= nrelem; i++)
    {
      connect  = mesh->GetConnect(i);

      if (dim == 2)     // 2d elements
	{
	  if (elemsize == 3)     // triangle
	    {
	      if (binary)
		{
		  WriteString ("tri     ");
		  WriteLong(3);
		}
	      else
		(*gmvfile) << "tri 3" << endl;

	      long index[3];
	      index[0] = connect[0];
	      index[1] = connect[1];
	      index[2] = connect[2];
	      WriteLong (3, index);
	    }
	  else if (elemsize == 4)   // quadrilateral
	    {
	      if (binary)
		{
		  WriteString ("quad    ");
		  WriteLong (4);
		}
	      else
		(*gmvfile) << "quad 4" << endl;

	      long index[4];
	      index[0] = connect[0];
	      index[1] = connect[1];
	      index[2] = connect[3];
	      index[3] = connect[2];
	      WriteLong (4, index);
	    }
	}
      else if (dim == 3)
	{
	  if (elemsize == 4)     // tetrahedron
	    {
	      if (binary)
		{
		  WriteString ("tet     ");
		  WriteLong(4);
		}
	      else
		(*gmvfile) << "tet 4" << endl;

	      long index[4];
	      index[0] = connect[0];
	      index[1] = connect[1];
	      index[2] = connect[2];
	      index[3] = connect[3];
	      WriteLong (4, index);
	    }
	  else if (elemsize == 8)   // hexahedron
	    {
	      if (binary)
		{
		  WriteString ("thex     ");
		  WriteLong(8);
		}
	      else
		(*gmvfile) << "hex 8" << endl;

	      long index[8];
	      index[0] = connect[0];
	      index[1] = connect[1];
	      index[2] = connect[2];
	      index[3] = connect[3];
	      index[4] = connect[4];
	      index[5] = connect[5];
	      index[6] = connect[6];
	      index[7] = connect[7];
	      WriteLong (8, index);
	    }
	}  
    }

  if (! binary)
    (*gmvfile) << endl;

}

}
