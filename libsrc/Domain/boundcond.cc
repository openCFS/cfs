#include <iostream>
#include <fstream>
#include <string>

#include "bcs.hh"

namespace CoupledField
{

BCs :: BCs(FileType * const aInFile)
{
#ifdef TRACE
  (*trace) << "entering BCs::BCs" << std::endl;
#endif
 InFile     = aInFile; 
}

BCs :: ~BCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::~BCs" << std::endl;
#endif
  delete [] numDirichlet;
  delete [] numNeumann;
  delete [] numConstraints;
  
}

void BCs :: ReadBCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::ReadBCs" << std::endl;
#endif
  level = 0;

 //read restraints information
 InFile->ReadNumNodesforDirichletBC(restr[level].num);
 numDirichlet[level] = restr[level].num;
 restr[level].info = new Integer*[restr[level].num];
 for (Integer i=0;i<restr[level].num;i++)
   restr[level].info[i] = new Integer[3];
 restr[level].val = new Double[restr[level].num];
 InFile->ReadBoundRestr(restr[level].info, restr[level].num, restr[level].val);

}

}
