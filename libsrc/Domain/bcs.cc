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
 for (i=0; i<NUMLEVELGRID; i++) { bcs_[i]=NULL; bcsEdges_[i]=NULL;}
  
 conf->getliststr("list_nodes",levels_);
 bcs_[0]=new std::list<Integer>[levels_.size()];
 
 conf->getliststr("list_edges",color_edges_);
 bcsEdges_[0]=new std::vector<Elem*>[color_edges_.size()]; 

 conf->getliststr("list_neighelems",color_neighelems_);
 bcsNeighElems_[0]=new std::vector<Elem*>[color_neighelems_.size()];

}

BCs :: ~BCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::~BCs" << std::endl;
#endif

 Integer i,j,k;

 if (!levels_.size()) {
   for (i=0; i<NUMLEVELGRID; i++)   
    if (bcs_[i]) {
      for (j=0; j<levels_.size(); j++) {
	bcs_[i][j].clear();
      } // loop over colors of grid( over subdomains), index j
      delete [] bcs_[i];
    } // loop over levels of grid, index i
 } // end of if

  if (!color_edges_.size()) {
   for (i=0; i<NUMLEVELGRID; i++)   
    if (bcsEdges_[i]) {
      for (j=0; j<color_edges_.size(); j++) {
	for (k=0; k<bcsEdges_[i][j].size(); k++) {
	  Elem* tmpEl=bcsEdges_[i][j][k];
	  delete tmpEl;
	} // loop over elements of one color, index k
      } // loop over colors of grid( over subdomains), index j
      delete [] bcsEdges_[i];
    } // loop over levels of grid, index i
 } // end of if

  if (!color_neighelems_.size()) {
   for (i=0; i<NUMLEVELGRID; i++)
    if (bcsNeighElems_[i]) {
      for (j=0; j<color_neighelems_.size(); j++) {
        for (k=0; k<bcsNeighElems_[i][j].size(); k++) {
          Elem* tmpEl=bcsNeighElems_[i][j][k];
          delete tmpEl;
        } // loop over elements of one color, index k
      } // loop over colors of grid( over subdomains), index j
      delete [] bcsNeighElems_[i];
    } // loop over levels of grid, index i
 } // end of if

}

void BCs :: ReadBCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::ReadBCs" << std::endl;
#endif

 InFile_->ReadBCs(bcs_[0],levels_);
 toplevel_=0;

 if (color_edges_.size()) InFile_->ReadEl1d(bcsEdges_[0],color_edges_);

 if (color_neighelems_.size())  InFile_->ReadEl2d(bcsNeighElems_[0],color_neighelems_);
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

std::vector<Elem*> BCs::getEdgesBC(const std::string color, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::getEdgesBC" << std::endl;
#endif

  Integer level=lev;
  if (lev==-1) level=toplevel_;
  Integer i;
  for (i=0; i<color_edges_.size(); i++)
    if (color==color_edges_[i]) break;

  return bcsEdges_[level][i];
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

std::vector<Elem*> BCs::getNeighbors2DElemsFor1D(const std::string color, const Integer lev
)
{
#ifdef TRACE
  (*trace) << "entering BCs::getNeighElems" << std::endl;
#endif

  Integer level=lev;
  if (lev==-1) level=toplevel_;
  Integer i;
  for (i=0; i<color_neighelems_.size(); i++)
    if (color==color_neighelems_[i]) break;

  return bcsNeighElems_[level][i];
}

}
