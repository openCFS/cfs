#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "bcs.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

  // ***************************************************
  //   Old constructor version using conffile approach
  // ***************************************************
#ifndef XMLPARAMS
  BCs :: BCs(FileType * const aInFile)
  {
#ifdef TRACE
    (*trace) << "entering BCs::BCs" << std::endl;
#endif
    InFile_     = aInFile; 

 Integer i;
 for (i=0; i<NUMLEVELGRID; i++) { bcs_[i]=NULL; bcsEdges_[i]=NULL; bcsFaces_[0]=NULL; bcsNeighElems_[0]=NULL;}
    conf->ifgetliststr("list_nodes",levels_);
    if (levels_.size()) 
      bcs_[0]=new std::list<Integer>[levels_.size()];
 
    conf->ifgetliststr("list_edges",color_edges_);
    if (color_edges_.size()) 
      bcsEdges_[0]=new std::vector<Elem*>[color_edges_.size()]; 

    conf->ifgetliststr("list_faces",color_faces_);
    if (color_faces_.size()) 
      bcsFaces_[0]=new std::vector<Elem*>[color_faces_.size()]; 

    conf->ifgetliststr("list_neighelems",color_neighelems_);
    if (color_neighelems_.size())
      bcsNeighElems_[0]=new std::vector<Elem*>[color_neighelems_.size()];
  }

#else

  // ********************************************************
  //   New constructor version using XML parameter handling
  // ********************************************************
  BCs :: BCs(FileType * const aInFile) : InFile_(aInFile)
  {
    ENTER_FCN( "BCs::BCs" );

    // Initialise internal pointer arrays
    for( Integer i = 0; i < NUMLEVELGRID; i++ )
      {
	bcs_[i] = NULL;
	bcsEdges_[i] = NULL;
      }

    // Get list of special node sets
    params->GetList( "name", levels_, "domain", "nodes" );
    if( levels_.size() > 0 )
      {
	bcs_[0] = new std::list<Integer>[levels_.size()];
      }

    //
    // NOTE: This must still be converted !!!
    //
    Info->Warning( "BCs: list_edges an co. not supported by XML!?" );
//    conf->ifgetliststr("list_edges",color_edges_);
//    if (color_edges_.size()) 
//      bcsEdges_[0]=new std::vector<Elem*>[color_edges_.size()]; 
//
//    conf->ifgetliststr("list_faces",color_faces_);
//    if (color_faces_.size()) 
//      bcsFaces_[0]=new std::vector<Elem*>[color_faces_.size()]; 
//
//    conf->ifgetliststr("list_neighelems",color_neighelems_);
//    if (color_neighelems_.size())
//      bcsNeighElems_[0]=new std::vector<Elem*>[color_neighelems_.size()];
  }
#endif


BCs :: ~BCs()
{
#ifdef TRACE
  (*trace) << "entering BCs::~BCs" << std::endl;
#endif

 Integer i,j,k;

 if (levels_.size()) {
   for (i=0; i<NUMLEVELGRID; i++)   
    if (bcs_[i]) {
      for (j=0; j<levels_.size(); j++) {
	bcs_[i][j].clear();
      } // loop over colors of grid( over subdomains), index j
      delete [] bcs_[i];
    } // loop over levels of grid, index i
 } // end of if

  if (color_edges_.size()) {
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

 if (color_faces_.size()) {
   for (i=0; i<NUMLEVELGRID; i++)   
    if (bcsFaces_[i]) {
      for (j=0; j<color_faces_.size(); j++) {
	for (k=0; k<bcsFaces_[i][j].size(); k++) {
	  Elem* tmpEl=bcsFaces_[i][j][k];
	  delete tmpEl;
	} // loop over elements of one color, index k
      } // loop over colors of grid( over subdomains), index j
      delete [] bcsFaces_[i];
    } // loop over levels of grid, index i
 } // end of if

  if (color_neighelems_.size()) {
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

 if (levels_.size())
   InFile_->ReadBCs(bcs_[0],levels_);
 toplevel_=0;

 if (color_edges_.size())
   InFile_->ReadEl1d(bcsEdges_[0],color_edges_);

 if (color_faces_.size())
   InFile_->ReadEl2d(bcsFaces_[0],color_faces_);

 if (color_neighelems_.size()) {
    if (InFile_->ReadDim()==2)
       InFile_->ReadEl2d(bcsNeighElems_[0],color_neighelems_);
     else
       InFile_->ReadEl3d(bcsNeighElems_[0],color_neighelems_);   
   }
}

std::list<Integer>&  BCs::GetNodesLevel(const std::string color, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::GetNodesLevel" << std::endl;
#endif
  
  Boolean Found = FALSE;

  Integer level=lev;
  if (lev==-1) level=toplevel_;
  Integer i;
  for (i=0; i<levels_.size(); i++)
    if (color==levels_[i]) 
      {
	Found = TRUE;
	break;
      }

  if (!Found)
    {
      std::string ErrMsg = "Nodes for level \'" + color + "\' could not be found!";
      Error(ErrMsg.c_str(),__FILE__,__LINE__); 
    }

 return bcs_[level][i]; 
}

std::vector<Elem*>& BCs::getEdgesBC(const std::string color, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::getEdgesBC" << std::endl;
#endif

  Boolean Found = FALSE;

  Integer level=lev;
  if (lev==-1) level=toplevel_;
  Integer i;
  for (i=0; i<color_edges_.size(); i++)
    if (color==color_edges_[i]) 
      {
	Found = TRUE;
	break;
      }
  
  if (!Found)
    {
      std::string ErrMsg = "Edges for level \'" + color + "\' could not be found!";
      Error(ErrMsg.c_str(),__FILE__,__LINE__); 
    }

  return bcsEdges_[level][i];
}

std::vector<Elem*>& BCs::getFacesBC(const std::string color, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::getFacesBC" << std::endl;
#endif
  Boolean Found = FALSE;
 
  Integer level=lev;
  if (lev==-1) level=toplevel_;
  Integer i;

 if (color_faces_.size())
  for (i=0; i<color_faces_.size(); i++)
    if (color==color_faces_[i]) 
      {
	Found = TRUE;
	return bcsFaces_[level][i];
	break;
      }

 if (color_edges_.size())
  for (i=0; i<color_edges_.size(); i++)
    if (color==color_edges_[i]) 
      {
	Found = TRUE;
	return bcsEdges_[level][i];
	break;
      }

   if (!Found)
    {
      std::string ErrMsg = "Faces nor Edges for level \'" + color + "\' could not be found!";
      Error(ErrMsg.c_str(),__FILE__,__LINE__); 
    }
   

  return bcsFaces_[level][i];
}

Integer BCs::GetNumNodesLevel(const std::string color, const Integer lev)
{
#ifdef TRACE
  (*trace) << "entering BCs::GetNumNodesLevel" << std::endl;
#endif
   Boolean Found = FALSE;

   Integer level=lev;
   if (lev==-1) level=toplevel_; 
   Integer i;
   for (i=0; i<levels_.size(); i++)
     if (color==levels_[i]) 
       {
	 Found = TRUE;
	 break;
       }

   if (!Found)
     {
       std::string ErrMsg = "Nodes for level \'" + color
	 + "\' could not be found!";
       Error( ErrMsg.c_str(), __FILE__, __LINE__ ); 
     }

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

std::vector<Elem*> BCs::getNeighElemsForSurfaces(const std::string color, const Integer lev
)
{
#ifdef TRACE
  (*trace) << "entering BCs::getNeighElems" << std::endl;
#endif

  Boolean Found = FALSE;

  Integer level=lev;
  if (lev==-1) level=toplevel_;
  Integer i;
  for (i=0; i<color_neighelems_.size(); i++)
    if (color==color_neighelems_[i]) 
      {
	Found = TRUE;
	break;
      }

  if (!Found)
    {
      std::string ErrMsg = "Neighbor to surface Elements for level \'" + color + "\' could not be found!";
      Error(ErrMsg.c_str(),__FILE__,__LINE__); 
    }

  return bcsNeighElems_[level][i];
}

}
