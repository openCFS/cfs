#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "mydefs.hh"
#include "output.hh"
#include "filesystem.hh"

namespace CoupledField
{

FileSystem :: FileSystem(char * ainputfile)
{
  filename = new char[100];

#ifdef TRACE
  strcpy(filename,ainputfile);
  trace = new ofstream(strcat(filename,".trace"));

  (*trace) << "entering FileSystem::FileSystem" << endl;
#endif

#ifdef DEBUG
  strcpy(filename,ainputfile);
  debug = new ofstream(strcat(filename,".deb"));

  strcpy(filename,ainputfile);
  test = new ofstream(strcat(filename,".test"));
  
#endif

#ifdef MEMTRACE
  strcpy(filename,ainputfile);
  memtrace = new ofstream(strcat(filename,".mem"));
#endif

  strcpy(filename,ainputfile);
  conv = new ofstream(strcat(filename,".con"));
}
  
FileSystem :: ~FileSystem()
{
#ifdef TRACE
  (*trace) << "entering FileSystem::~FileSystem" << endl;
#endif

  delete [] filename;
}

}
