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
 toplevel_=0;

}

std::list<Integer> BCs::GetNodesLevel(const std::string color, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::GetNodesLevel" << std::endl;
#endif
 
  Integer level=lev;
 if (lev==-1) level=toplevel_;
 Integer i;
 for (i=0; i<levels_.size(); i++)
  if (color==levels_[i]) break;

 return bcs_[level][i]; 
}

Integer BCs::GetNumNodesLevel(const std::string color, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::GetNumNodesLevel" << std::endl;
#endif

  Integer level=lev;
  if (lev==-1) level=toplevel_; 
 Integer i;
 for (i=0; i<levels_.size(); i++)
  if (color==levels_[i]) break;

 return bcs_[level][i].size();

}

void BCs :: Update(Grid * ptgrid)
{   
    toplevel_++;
    bcs_[toplevel_]=new std::list<Integer>[levels_.size()]; 
    if (!bcs_[toplevel_]) Error(" Not enought memory",__FILE__,__LINE__); 

    ptgrid->UpdateBCs(bcs_[toplevel_]); 
}

void BCs :: printBCs(const Integer alevel)
{
  Integer level=alevel;
  if (level==-1) level=toplevel_;  

  Integer i,ilevel;
  for (i=0; i<levels_.size(); i++)
    {
      for (std::list<Integer>::const_iterator p=bcs_[ilevel][i].begin(); p!=bcs_[ilevel][i].end(); p++)
	{ std::cout << (*p) << std::endl;} 
      
    } 
}

}
