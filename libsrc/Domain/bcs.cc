#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "bcs.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField {


  // ****************************************************
  //   Constructor version using XML parameter handling
  // ****************************************************
  BCs::BCs(FileType * const aInFile) {

    ENTER_FCN( "BCs::BCs" );

    InFile_     = aInFile; 
  
    // Initialise internal pointer arrays
    for( Integer i = 0; i < NUMLEVELGRID; i++ ) {
      bcs_[i] = NULL;
      bcsEdges_[i] = NULL;
    }

    // Get list of special node sets
    params->GetList( "name", levels_, "domain", "nodes" );
    if( levels_.GetSize() > 0 )
      bcs_[0] = new std::list<Integer>[levels_.GetSize()];
    
    //get geometric dimension
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    if ( probGeo == "3d" ) {
      //surface elements can occure
      params->GetList( "name", color_faces_, "domain", "elements" );
      if (color_faces_.GetSize()) {
        bcsFaces_[0] = new StdVector<Elem*>[color_faces_.GetSize()]; 
      }
    }    
    else {
      params->GetList( "name", color_edges_, "domain", "elements" );    
      if (color_edges_.GetSize()) {
        bcsEdges_[0] = new StdVector<Elem*>[color_edges_.GetSize()];
      }
    }
    
    //
    // NOTE: This must still be converted !!!
    //
#ifdef DEBUG
    //    Info->Warning( "BCs: list_edges an co. not supported by XML!?" );
#endif

  }


  // **************
  //   Destructor
  // **************
  BCs::~BCs() {

    ENTER_FCN( "BCs::~BCs" );

    Integer i, j, k;

    if ( levels_.GetSize() ) {

      // loop over levels of grid
      for ( i = 0; i < NUMLEVELGRID; i++ ) {
        if ( bcs_[i] ) {

          // loop over colors of grid / subdomains
          for ( j = 0; j < levels_.GetSize(); j++ ) {
            bcs_[i][j].clear();
          }
          delete [] bcs_[i];
        }
      }
    }

    if ( color_edges_.GetSize() ) {

      // loop over levels of grid
      for ( i = 0; i < NUMLEVELGRID; i++ ) {
        if (bcsEdges_[i]) {

          // loop over colors of grid / subdomains
          for ( j = 0; j < color_edges_.GetSize(); j++ ) {

            // loop over elements of one color
            for ( k = 0; k < bcsEdges_[i][j].GetSize(); k++ ) {
              delete (bcsEdges_[i][j][k]);
            }
          }
          delete [] bcsEdges_[i];
        }
      }
    }

    // NOTE: color_faces_ and color_edges_ are identical (see the code in
    //       the constructor), thus we should not try to do another deep
    //       deletion)
    delete [] bcsFaces_[i];
    // if ( color_faces_.GetSize() ) {

    // loop over levels of grid
    // for ( i = 0; i < NUMLEVELGRID; i++ ) {
    // if ( bcsFaces_[i] ) {

    // loop over colors of grid / subdomains
    // for ( j = 0; j < color_faces_.GetSize(); j++ ) {

    // loop over elements of one color
    // for ( k = 0; k < bcsFaces_[i][j].GetSize(); k++ ) {
    // delete (bcsFaces_[i][j][k]);
    // }
    // }
    // delete [] bcsFaces_[i];
    // }
    // }
    // }

    if ( color_neighelems_.GetSize() ) {

      // loop over levels of grid
      for ( i = 0; i < NUMLEVELGRID; i++ ) {
        if ( bcsNeighElems_[i] ) {

          // loop over colors of grid / subdomains
          for ( j = 0; j < color_neighelems_.GetSize(); j++ ) {

            // loop over elements of one color
            for ( k = 0; k < bcsNeighElems_[i][j].GetSize(); k++ ) {
              delete (bcsNeighElems_[i][j][k]);
            }
          }
          delete [] bcsNeighElems_[i];
        }
      }
    }
  }


  void BCs::ReadBCs() {

    ENTER_FCN( "BCs::ReadBCs" );

    // Create a dummy vector, since we do not need
    // an vector with pointers to the elements, sorted
    // by element numbers
    StdVector<Elem*> temp;

    if (levels_.GetSize())
      InFile_->ReadBCs(bcs_[0],levels_);
    toplevel_=0;

    if (color_edges_.GetSize())
      InFile_->ReadEl1d(bcsEdges_[0],temp,color_edges_);

    if (color_faces_.GetSize())
      {
        InFile_->ReadEl2d(bcsFaces_[0],temp,color_faces_);
      }

    if (color_neighelems_.GetSize()) {
      if (InFile_->ReadDim()==2)
        InFile_->ReadEl2d(bcsNeighElems_[0],temp,color_neighelems_);
      else
        InFile_->ReadEl3d(bcsNeighElems_[0],temp,color_neighelems_);   
    }
  }

  std::list<Integer>&  BCs::GetNodesLevel(const std::string color, const Integer lev)
  {
    ENTER_FCN( "BCs::GetNodesLevel" );
  
    Boolean Found = FALSE;

    Integer level=lev;
    if (lev==-1) level=toplevel_;
    Integer i;
    for (i=0; i<levels_.GetSize(); i++)
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

  StdVector<Elem*>& BCs::getEdgesBC(const std::string color, const Integer lev)
  {
    ENTER_FCN( "BCs::getEdgesBC" );

    Boolean Found = FALSE;

    Integer level=lev;
    if (lev==-1) level=toplevel_;
    Integer i;
    for (i=0; i<color_edges_.GetSize(); i++)
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

  StdVector<Elem*>& BCs::getFacesBC(const std::string color, const Integer lev)
  {
    ENTER_FCN( "BCs::getFacesBC" );

    Boolean Found = FALSE;
 
    Integer level=lev;
    if (lev==-1) level=toplevel_;
    Integer i;

    if (color_faces_.GetSize())
      for (i=0; i<color_faces_.GetSize(); i++)
        if (color==color_faces_[i]) 
          {
            Found = TRUE;
            return bcsFaces_[level][i];
            break;
          }

    if (color_edges_.GetSize())
      for (i=0; i<color_edges_.GetSize(); i++)
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
    ENTER_FCN( "BCs::GetNumNodesLevel" );

    Boolean Found = FALSE;

    Integer level=lev;
    if (lev==-1) level=toplevel_; 
    Integer i;
    for (i=0; i<levels_.GetSize(); i++)
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
    ENTER_FCN( "BCs::Update" );
    toplevel_++;
    bcs_[toplevel_]=new std::list<Integer>[levels_.GetSize()]; 
    if (!bcs_[toplevel_]) Error(" Not enought memory",__FILE__,__LINE__); 

    ptgrid->UpdateBCs(bcs_[toplevel_]); 
  }

  void BCs :: printBCs(const Integer alevel)
  {
    ENTER_FCN( "BCs::printBCs" );
    Integer level=alevel;
    if (level==-1) level=toplevel_;  

    Integer i,ilevel;
    for (i=0; i<levels_.GetSize(); i++)
      {
        for (std::list<Integer>::const_iterator p=bcs_[ilevel][i].begin(); p!=bcs_[ilevel][i].end(); p++)
          { std::cout << (*p) << std::endl;} 
      
      } 
  }

  StdVector<Elem*> BCs::getNeighElemsForSurfaces(const std::string color, const Integer lev
                                                 )
  {
    ENTER_FCN( "BCs::getNeighElems" );

    Boolean Found = FALSE;

    Integer level=lev;
    if (lev==-1) level=toplevel_;
    Integer i;
    for (i=0; i<color_neighelems_.GetSize(); i++)
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
