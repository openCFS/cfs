#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "gmvwriter.hh"

GMV_Writer :: GMV_Writer (const char * name, int abinary, const char * afromfile)
{
  const char * gmv_default_name = "gmvoutput.dat";

  binary = abinary;

  // test file name
  if (name != NULL)
    {
      filename = new char[strlen(name) + 1];
      strcpy (filename, name);
    }
  else
    {
      filename = new char [strlen(gmv_default_name) + 1];
      strcpy (filename, gmv_default_name);
    }

  // test fromfile (file with grid information, useful for dynamic simulations)
  if (afromfile != NULL)
    {
      fromfile = new char[strlen(afromfile) + 1];
      strcpy (fromfile, afromfile);  
    }
  else
    {
      fromfile = NULL;
    }

  // open file
  if (binary)
    // gmvfile = new ofstream(filename, ios::out | ios::binary);
    gmvfile = new ofstream(filename, 0);
  else
    gmvfile = new ofstream(filename);

  if (! gmvfile->good())
    {
      cerr << "GMV_Writer : can not open file " << filename << endl;
      return;
    }

  variablemode = 0;
}


GMV_Writer :: ~GMV_Writer ()
{
  if (variablemode)
    {
      (*gmvfile) << "endvars ";
      if (! binary) 
	(*gmvfile) << endl << endl;
    }

  (*gmvfile) << "endgmv  ";
  if (! binary)
    (*gmvfile) << endl;

  // close file
  delete gmvfile;

  // delete filename
  delete [] filename;
  if (fromfile)
    delete [] fromfile;
}

/////////////////////////////////////////////////////////////////

void GMV_Writer :: WriteFileHeader () 
{
  // write header
  (*gmvfile) << "gmvinput";
  if (binary)
    (*gmvfile) << "ieee    ";
  else
    (*gmvfile) << " ascii" << endl << endl;


  // write nodes of the mesh (if necessary)
  nrnodes = NrNodes();

  if (fromfile)
    {
      if (binary)
	(*gmvfile) << "nodes   fromfile\"" << fromfile << "\"";
      else
	(*gmvfile) << "nodes fromfile \"" << fromfile << "\"" << endl << endl;
    }
  else
    {
      WriteNodes ();
    }

  // write cells of the mesh (if necessary)
  nrcells = NrCells ();

  if (fromfile)
    {
      if (binary)
	{
	  (*gmvfile) << "cells   fromfile\"" << fromfile << "\"";
	}
      else
	{
	  (*gmvfile) << "cells fromfile \"" << fromfile << "\"" << endl << endl;
	}
    }
  else
    {
      WriteCells ();
    }
}


void GMV_Writer :: WriteLong (long val)
{
  if (binary)
    gmvfile->write ((char *) & val, sizeof (long));
  else
    (*gmvfile) << val;
}


void GMV_Writer :: WriteLong (long len, const long * val)
{
  long i;

  if (binary)
    gmvfile->write ((char *) val, len * sizeof (long));
  else
    {
      for (i = 0; i < len; i++)
	{
	  (*gmvfile) << val[i] << "   ";
	  if ((i + 1) % 8 == 0)
	    (*gmvfile) << endl;
	}
      if (i % 8 != 0)
	(*gmvfile) << endl;
    }
}


void GMV_Writer :: WriteString (const char * str)
{
  (*gmvfile) << str;
}


void GMV_Writer :: WriteFloat (float val)
{
  if (binary)
    gmvfile->write ((char *) & val, sizeof (float));
  else
    (*gmvfile) << val;
}


void GMV_Writer :: WriteFloat (long len, const float * val)
{
  long i;

  if (binary)
    gmvfile->write ((char *) val, len * sizeof (float));
  else
    {
      for (i = 0; i < len; i++)
	{
	  (*gmvfile) << val[i] << "   ";
	  if ((i + 1) % 8 == 0)
	    (*gmvfile) << endl;
	}
      if (i % 8 != 0)
	(*gmvfile) << endl;
    }
}

