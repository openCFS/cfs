#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "bcs.hh"
#include "conffile.hh"

namespace CoupledField
{

BCs :: BCs(FileType * const aInFile)
{
#ifdef TRACE
  (*trace) << "entering BCs::BCs" << std::endl;
#endif
 InFile_     = aInFile; 

 Integer i;
 for (i=0; i<NUMLEVELGRID; i++) bcs_[i]=NULL;
  
 conf->getliststr("list_interfaces",levels_);
 bcs_[0]=new std::list<Integer>[levels_.size()];

}

BCs :: ~BCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::~BCs" << std::endl;
#endif

 Integer i;
 for (i=0; i<NUMLEVELGRID; i++)   
   if (bcs_[i]) delete [] bcs_[i];  

}

void BCs :: ReadBCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::ReadBCs" << std::endl;
#endif

 InFile_->ReadBCs(bcs_[0],levels_);

}

std::list<Integer> BCs::GetNodesLevel(const std::string level, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::GetNodesLevel" << std::endl;
#endif
 
 Integer i;
 for (i=0; i<levels_.size(); i++)
  if (level==levels_[i]) break;

 return bcs_[lev][i];
 
}

Integer BCs::GetNumNodesLevel(const std::string level, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::GetNumNodesLevel" << std::endl;
#endif

 Integer i;
 for (i=0; i<levels_.size(); i++)
  if (level==levels_[i]) break;

 return bcs_[lev][i].size();

}

void BCs :: Update(Grid * ptgrid)
{
   ;
}

}
