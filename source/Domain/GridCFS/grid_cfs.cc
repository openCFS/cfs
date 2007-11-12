// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <set>

#include "grid_cfs.hh"

#include "Elements/elements_header.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "General/exception.hh"
#include "Utils/coordSystem.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/domain.hh"

namespace CoupledField {

  // declare class specific logging stream
  DECLARE_LOG(gridcfs)
  DEFINE_LOG(gridcfs, "grid.cfs")

  GridCFS::GridCFS(UInt dim) : Grid( ) {


    isInitialized_ = false;
    isQuadratic_ = false;
    dim_ = dim;
    
    numNodes_ = 0;
    numElems_ = 0;
    numFaces_ = 0;
    numEdges_ = 0;
    edgesMapped_ = false;
    facesMapped_ = false;


  } 


  // **************
  //   Destructor
  // **************
  
  GridCFS::~GridCFS() {


    for ( UInt i = 0; i < numElems_; i++ ) {
      delete orderedElems_[i];
    }
    orderedElems_.Clear();
  } 

  /*
    void GridCFS::Read()
    {
 
    StdVector<StdVector<Elem*> >temp;

    // 1. Read Dimension
    dim_=ptFileType->GetDim();
    numNodes_ = ptFileType->GetNumNodes();
    numElems_ = ptFileType->GetNumElems();

    // 2. Read coordinatess
    ptFileType->GetCoordinates(coords_);

    // 3. Read in volume elements
    ptFileType->GetElements(volElems_, volRegionIds_, dim_);
     
    // 4. Create sorted list of volume elements
    // 5. If there are surfaces
    //    -  Iterate over all volume elements and make list
    //       of neighbouring elements for each node
    //    -   Read In Surface elements
    ptFileType->GetElements(temp, surfRegionIds_, dim_-1);
    
    // Obtain all region names
    ptFileType->GetAllRegionNames(regionNames_);
    // Create ordered list of elements
    Elem * ptVolElem;
    orderedElems_.Resize(numElems_);
    orderedElems_.Init(NULL);
    for ( UInt iRegion = 0; iRegion < volElems_.GetSize(); iRegion++ ) {
      
    for ( UInt iElem = 0; iElem < volElems_[iRegion].GetSize();
    iElem++ ) {
    ptVolElem =  volElems_[iRegion][iElem];
    BaseFE * ptFE = ptVolElem->ptElem;

    // Check, if quadratic elements are in the grid
    if  ( ptFE->GetNumCorners() < ptFE->GetNumNodes() ) {
    isQuadratic_ = true;
    }
        
    // Add type of FE to map
    numElemTypes_[ptFE->feType()]++;
        
    // Check, if element with same number is already contained
    // in the grid
    if ( orderedElems_[ptVolElem->elemNum-1]  != NULL) {
    EXCEPTION( "Element Nr. " << ptVolElem->elemNum 
    << " exists at least two times!\n"
    << "The first occurence is in region '"
    << regionNames_[orderedElems_[ptVolElem->elemNum-1]
    ->regionId] 
    << "', the second in region '"
    << regionNames_[volElems_[iRegion][iElem]->regionId]
    << "'.\nPlease check your mesh file!" );
    }
    orderedElems_[ptVolElem->elemNum-1] = ptVolElem;
    } // loop over elements
    } // loop over region
    
    if ( temp.GetSize() > 0 )
    CreateSurfaceElements(temp);

    // Perform final consistency check
    for( UInt i = 0; i<orderedElems_.GetSize(); i++ ) {
    if ( orderedElems_[i] == NULL ) {
    (*warning) << "Gap in numbering: No Element Nr. " << i+1 
    << " contained in the mesh! Errors can occur!";
    Warning( __FILE__, __LINE__ );
    }
    }

    
    // Get nodes of each volume and surface region
    ptFileType->GetNodesOfRegions( volElemNodes_, volRegionIds_ );
    ptFileType->GetNodesOfRegions( surfElemNodes_, surfRegionIds_ );
 
    
    // 8. Read in named nodes
    ptFileType->GetNamedNodes( namedNodes_, namedNodeNames_ );
    
    // 9. Read in named elements
    ptFileType->GetNamedElems( namedElems_, namedElemNames_ );

    // 10. Get nodes for directivity plots if defined in XML file
    GetNodesForDirectivity() ;

    // Print information about region mapping into
    PrintGridInfo();
    
    #ifdef ADAPTGRID
    FormNeighborsLists();
    #endif

    //     UInt inode = numNodes_+1;
    //     AddNodes(10);

    //     for(UInt i=0; i<10; i++)
    //     {
    //         Point p;
    //         p[0] = p[1] = p[2] = i;
    //         SetNodeCoordinate(inode++, p);
    //     }

    FinishInit();
    }
  */

  void GridCFS::CreateUserDefinedNodesElems() {
    
    // if no param object is present, just leave
    if (!param) return;
    
    Vector<Double> coord(3);
    std::string coordSys;
    std::string name;

    for( UInt iType = 0; iType < 2; iType++ ) {
      ParamNode * listNode = NULL;      
      StdVector<ParamNode*> nodes;
      bool isNode = true;
  
      if( iType == 0 ) {
        // iterate over nodes
        listNode = param->Get("domain")->Get("nodeList", false);
        isNode = true;
      } else {
        // iterate over nodes
        listNode = param->Get("domain")->Get("elemList", false);
        isNode = false;
      }

      if (listNode) {
        nodes = listNode->GetChildren();
        
        for( UInt i=0; i < nodes.GetSize(); i++ ) {
          
          // fetch name of nodes to be selected
          nodes[i]->Get("name", name );
          
          // check if node is defined by point coord
          ParamNode * coordNode = nodes[i]->Get("coord", false );
          if( coordNode ) {
            
            // ToDo: insert defintion for axisymmetric geometry
            coord.Init();
            coordNode->Get( "x", coord[0] );
            coordNode->Get( "y", coord[1] );
            coordNode->Get( "z", coord[2] );
            StdVector<UInt> entityNum(1);
            entityNum[0] = FindEntityMinDistance( isNode, coord );
            
            // add node / element 
            if( isNode ) {
              AddNamedNodes( name, entityNum );
            } else {
              AddNamedElems( name, entityNum );
            }
            
          }
          
          // check if node is defined by parametric description
          ParamNode * listNode = nodes[i]->Get("list", false );
          if( listNode ) {
            
            std::string coordSysId;

            // get coordinate system
            listNode->Get( "coordSysId", coordSysId );

            StdVector<PointSelection> selections;

            // get free component
            StdVector<ParamNode*>  freeNodes = 
              listNode->GetList("freeCoord");
            for( UInt i = 0; i < freeNodes.GetSize(); i++ ) {
              ParamNode * actNode = freeNodes[i];
              PointSelection actSel;
              actSel.isFree = true;
              actNode->Get( "comp", actSel.comp );
              actNode->Get( "start", actSel.start );
              actNode->Get( "stop", actSel.stop );
              actNode->Get( "inc", actSel.inc );
              selections.Push_back( actSel);
            }
            
            // get fixed component(s)
            StdVector<ParamNode*> fixedNodes = 
              listNode->GetList("fixedCoord");
            for( UInt i = 0; i < fixedNodes.GetSize(); i++ ) {
              ParamNode * actNode = fixedNodes[i];
              PointSelection actSel;
              actSel.isFree = false;
              actNode->Get( "comp", actSel.comp );
              actNode->Get( "value", actSel.value );
              selections.Push_back( actSel);
            }
            
            // Ensure, that we have as many entries in the 
            // vector as dimension
            if( dim_ != selections.GetSize() ) {
              EXCEPTION("There have been more coordinate components "
                        << "than there are space dimensions" );
            }
            AddEntityByParam( name, isNode, coordSysId,
                              selections );
          }
          
        }
      }
    }
  }
  
 

  void GridCFS::AddEntityByParam( const std::string& name, bool isNode, 
                                  const std::string& coordSysId,
                                  StdVector<PointSelection>& coords ) {

    // Check if entities with given name exist already
    if( nameTypeMap_.find( name) != nameTypeMap_.end() ) {
      EXCEPTION( "Entities with name " << name 
                 << " are already defined" );
    }
    
    // get coordinate system
    CoordSystem * cosy = NULL;
    try {
      cosy = domain->GetCoordSystem(coordSysId);
    } catch( Exception &e ) {
      RETHROW_EXCEPTION(e, "Could not select nodes within the"
                        << "coordinate system '"
                        << coordSysId << "'");
    }
    
    // map coordinate components to indices
    StdVector<UInt> dofs(dim_);

    for( UInt iDim = 0; iDim < dim_; iDim++ ) {
      dofs[iDim] = cosy->GetVecComponent( coords[iDim].comp )-1;
    }
    
    // require that the first entry is a free one 
    if( !coords[0].isFree ) {
      EXCEPTION( "First coordinate component must be free" );
    }

    // fetch math parser and register coordinate names
    MathParser * parser = domain->GetMathParser();
    StdVector<MathParser::HandleType> handles(dim_);
    StdVector<Vector<Double> > sampleVals(dim_);
    UInt totalNum = 1;
    for( UInt iDim = 0; iDim < dim_; iDim++ ) {
      if( !coords[iDim].isFree ) {
        handles[iDim] = parser->GetNewHandle();
        for( UInt iDim2 = 0; iDim2 < dim_; iDim2++ ) {
          if( !coords[iDim2].isFree )
            parser->SetValue( handles[iDim], coords[iDim2].comp, 0.0 );
        }
        parser->SetExpr( handles[iDim], coords[iDim].value );   
        sampleVals[iDim].Resize(1);
        sampleVals[iDim].Init(0.0);
      } else {
        UInt numSamples  = 
          UInt( floor ( (coords[iDim].stop-coords[iDim].start) 
                        / coords[iDim].inc ) )+1;
        sampleVals[iDim].Resize( numSamples );
        for( UInt iSample = 0; iSample < numSamples; iSample++ ) {
          sampleVals[iDim][iSample] = coords[iDim].start 
            + iSample * coords[iDim].inc;
        }
      }
      totalNum *= sampleVals[iDim].GetSize();
    }
    
    
    StdVector<UInt> entityNums(totalNum);
    Vector<Double> locCoord( dim_), globCoord( dim_ );
    Vector<Double> actEntCoord(dim_), temp;
    UInt pos = 0;
    
    // loop over all entries in the (first)free component vector
    // update first component, if it is free
    if( !coords[0].isFree ) {
      sampleVals[0][0] = parser->Eval(handles[0] );
    } 
    for( UInt iSample1 = 0; iSample1 < sampleVals[0].GetSize(); iSample1++ ) {
      locCoord[dofs[0]] = sampleVals[0][iSample1];      

      // update second component, if it is free
      if( !coords[1].isFree ) {
        parser->SetValue( handles[1], coords[0].comp, sampleVals[0][iSample1] );
        sampleVals[1][0] = parser->Eval(handles[1] );
      } 
      for( UInt iSample2 = 0; iSample2 < sampleVals[1].GetSize(); iSample2++ ) {

        locCoord[dofs[1]] = sampleVals[1][iSample2];        
        // update second component, if it is free
        if( dim_ == 3 ) {
          if( !coords[2].isFree ) {
            parser->SetValue( handles[2], coords[0].comp, sampleVals[0][iSample1] );
            parser->SetValue( handles[2], coords[1].comp, sampleVals[1][iSample1] );
            sampleVals[2][0] = parser->Eval(handles[2] );
          }
          for( UInt iSample3 = 0; iSample3 < sampleVals[2].GetSize(); iSample3++ ) {
            locCoord[dofs[2]] = sampleVals[2][iSample3];
            cosy->Local2GlobalCoord( globCoord, locCoord );
            entityNums[pos++] = FindEntityMinDistance( isNode, globCoord );
          }
          
        } // dim == 3 
        else {
          cosy->Local2GlobalCoord( globCoord, locCoord );
          entityNums[pos++] = FindEntityMinDistance( isNode, globCoord );
        }
      }
    }
    
    if( isNode) {
      
      // add named nodes to grid
      AddNamedNodes( name, entityNums );
    } else {
      AddNamedElems( name, entityNums );
    }
    
    // release handles of math parser
    for( UInt iDim = 0; iDim < dim_; iDim++ ) {
      if( !coords[iDim].isFree )
        parser->ReleaseHandle( handles[iDim] );
    }
  }

  UInt GridCFS::FindEntityMinDistance( bool isNode, Vector<Double>& coord ) {

    UInt entityNum;

    // iterate over all nodes/elements in the grid
    // vectors with node indices and distance
    std::vector<Double> entityDist;
    Vector<Double> actEntCoord, temp;
    
    if( isNode ) {
      // === loop over nodes ===
      entityDist.resize( numNodes_ );
      for( UInt iNode = 0; iNode < numNodes_; iNode++ ) {
        
        // calculate distance and store it in vector
        GetNodeCoordinate( actEntCoord, iNode+1, false );          
        temp = (actEntCoord-coord);
        entityDist[iNode] = temp.NormL2();
      } // nodes
    } else {
      
      // === loop over elements ===
      entityDist.resize(numElems_);
      
      Vector<Double> locMidPoint;
      Matrix<Double> connectCoord;
      
      // iterate over all elements
      for( UInt iElem = 0; iElem < numElems_; iElem++ ) {
        
        GetGlobalElemMidPoint( iElem+1, actEntCoord );
        temp = (actEntCoord-coord );
        entityDist[iElem] = temp.NormL2();
      } // elements
      
    }
    
    // find minimum entry in the vector
    std::vector<Double>::iterator it ;
    it = min_element(entityDist.begin(), entityDist.end());
    entityNum = std::distance(entityDist.begin(), it) + 1;
    return entityNum;


}

  void GridCFS::MapFaces() {

    LOG_TRACE(gridcfs) << "Starting to map faces ";

   
    // assert that any mesh was already read in
    assert( isInitialized_ == true );
   
    // If faces are already mapped ->leave
    if( facesMapped_ == true ) {
      return;
    }

    // create counter for number of faces
    UInt actFaceNum = 1;

    // vectors for local / global nodes
    StdVector<UInt> faceIndices, faceNodes;
   
    // iterate over all elements
    for( UInt iElem = 0; iElem < orderedElems_.GetSize(); iElem++ ) {

      // remember current element
      Elem & actElem = *orderedElems_[iElem];

      // if element is of wrong dimension (surface element )
      // ->leave
      //if ( actElem.ptElem->GetDim() < dim_ ) { continue; }

      // get number of element faces
      UInt numFaces = actElem.ptElem->GetNumFaces();

      // adapt size of faces and orientation array of element
      actElem.faces.Resize( numFaces );
      actElem.faceFlags.Resize( numFaces );

      // iterate over all faces of this element
      for( UInt iFace = 0; iFace < numFaces; iFace++ ) {
       
        // get local nodal indices of current face
        actElem.ptElem->GetFaceIndices( faceIndices, iFace );

        // create new Face object
        Face actFace;
        actFace.nodes.Resize( faceIndices.GetSize() );

        // insert node numbers into current face definition
        for( UInt iNode = 0; iNode < faceIndices.GetSize(); iNode++ ) {
          actFace.nodes[iNode] = actElem.connect[faceIndices[iNode]-1];
        }
       
        LOG_DBG3(gridcfs) << "Cecking face with nodes : " 
                          << actFace.nodes.Serialize();
       
        // Re-orientate face to match global orientation and
        // obtain the orientation flags
        std::bitset<3> orientation;
        actFace.Normalize( orientation );
       
        LOG_DBG3(gridcfs) << "Normalized : " 
                          << actFace.nodes.Serialize();
        LOG_DBG3(gridcfs) << "Orientation: " << orientation;
       

        // check if face was already numbered
        if( faceNums_.find( actFace ) == faceNums_.end() ) {
          LOG_DBG2(gridcfs) << "Adding face number " 
                            << actFaceNum << std::endl;
          faceNums_[actFace] = actFaceNum;
          faces_.Push_back( actFace );
          actElem.faces[iFace] = actFaceNum;
          actFaceNum++;
        } else {
          actElem.faces[iFace] = faceNums_[actFace];
          LOG_DBG2(gridcfs) << "--> already defined\n";
        }

        // Set also orientation flags for face
        actElem.faceFlags[iFace] = orientation;
      } // loop over faces
     
      // Print information about connectivity and faces
      LOG_DBG2(gridcfs) << "Elem Nr. " << actElem.elemNum;
      LOG_DBG2(gridcfs) << "===================";
      LOG_DBG2(gridcfs) << "Connectivity: " << actElem.connect.Serialize();
      LOG_DBG2(gridcfs) << "Faces: " << actElem.faces.Serialize();

      LOG_TRACE(gridcfs) << "Finished to map faces\n";

    } // loop over elements

    // Set flag for mapping of sub-entities
    numFaces_ = actFaceNum-1;
    LOG_DBG2(gridcfs) << "Total number of faces: " << numFaces_;

    // Set flag
    facesMapped_ = true;

    LOG_TRACE(gridcfs) << "Finished mapping faces\n";

  }
  


  void GridCFS::MapEdges() {

    LOG_TRACE(gridcfs) << "Starting to map edges";

    // assert that any mesh was already read in
    assert( isInitialized_ == true );

    // If edges/surfaces were already mapped ->leave
    if( edgesMapped_ == true ) {
      return;
    }

    // create counter for number of edges
    UInt actEdgeNum = 1;
    StdVector<UInt> locEdge(2), globEdge(2);

    // iterate over all elements
    for( UInt iElem = 0; iElem < orderedElems_.GetSize(); iElem++ ) {

      // remember current element
      Elem & actElem = *orderedElems_[iElem];

      // get number of edges
      UInt numEdges= actElem.ptElem->GetNumEdges();

      // adapt size of edge number array of element
      actElem.edges.Resize( numEdges );

      // iterate over all edges of this element
      for( UInt iEdge = 0; iEdge < numEdges; iEdge++ ) {
       
        // get local edge indices 
        actElem.ptElem->GetEdgeIndices( locEdge, iEdge );

        // create new edge
        Edge actEdge;
        actEdge.nodes[0] = actElem.connect[locEdge[0]-1];
        actEdge.nodes[1] = actElem.connect[locEdge[1]-1];

        // check if ordering is correct
        Integer orientation = 1;
       
        if( actEdge.nodes[1] < actEdge.nodes[0] ) {
          UInt secNode = actEdge.nodes[1];
          actEdge.nodes[1] = actEdge.nodes[0];
          actEdge.nodes[0] = secNode;

          // swap factor for orientation
          orientation = -1;
        }

        // check if edge was already numbered
        if( edgeNums_.find( actEdge ) == edgeNums_.end() ) {
          LOG_DBG2(gridcfs) << "Adding edge number " << actEdgeNum;
          LOG_DBG3(gridcfs) << "with nodes: " << actEdge.nodes[0] << ","
                            << actEdge.nodes[1] << std::endl;
          edgeNums_[actEdge] = actEdgeNum;
          edges_.Push_back( actEdge );
          actElem.edges[iEdge] = actEdgeNum*orientation;
          actEdgeNum++;
        } else {
          actElem.edges[iEdge] = edgeNums_[actEdge]*orientation;       
        }

      }
     
      // Print information about connectivity and edges
      LOG_DBG2(gridcfs) << "Elem Nr. " << actElem.elemNum;
      LOG_DBG2(gridcfs) << "===================";
      LOG_DBG2(gridcfs) << "Connectivity: " << actElem.connect.Serialize();
      LOG_DBG2(gridcfs) << "Edges: " << actElem.edges.Serialize();

      LOG_TRACE(gridcfs) << "Finished to map edges\n";
    }

    // Set flag for mapping of sub-entities
    numEdges_ = actEdgeNum-1;
    LOG_DBG2(gridcfs) << "Total number of edges: " << numEdges_ << std::endl;

    // Set flag
    edgesMapped_ = true;
  }

  UInt GridCFS::GetNumEdges() {
    return numEdges_;
    
  }

  UInt GridCFS::GetNumFaces() {
    return numFaces_;
    
  }


  const Edge&  GridCFS::GetEdge( UInt edgeNr ) {
    if( !edgesMapped_ ) {
      EXCEPTION( "Edges are not mapped yet!");
    }
  
    Edge const & ret = edges_[edgeNr-1];
    
    return ret;
  }
  
  const Face&  GridCFS::GetFace( UInt faceNr ) {
    if( !facesMapped_ ) {
      EXCEPTION( "Surfaces are not mapped yet!" );
    }
    
    Face const & ret = faces_[faceNr-1];
    
    return ret;
  }

  void GridCFS::FinishInit()
  {

    volElemNodes_.Resize(0);
    volRegionIds_.Resize(0);
    surfElemNodes_.Resize(0);
    surfRegionIds_.Resize(0);
    volElems_.Resize(0);
    surfElems_.Resize(0);
    maxNumElemNodes_ = 0;
        
        
    UInt e;
    UInt numElems = orderedElems_.GetSize();
    UInt numNodes;
    UInt numVolRegions, numSurfRegions;

    std::map<RegionIdType, StdVector<Elem*> > volRegionElems, surfRegionElems;
    std::map<RegionIdType, std::set<UInt> > volRegionNodes, surfRegionNodes;

    // set of elements, which get surface-mapped
    std::set<Elem*> surfElems;
    for(e=0; e<numElems; e++)
    {
      bool isSurfElem = false;
      Elem* el = orderedElems_[e];
      FEType type = el->ptElem->feType();
      numNodes = NUM_ELEM_NODES[type];

      maxNumElemNodes_ = maxNumElemNodes_ < numNodes ?
                         numNodes : maxNumElemNodes_;

      switch(type)
      {
      case ET_LINE2:
      case ET_LINE3:
        if(dim_ == 2) {
          isSurfElem = true;
        }
        break;
      case ET_TRIA3:
      case ET_TRIA6:
      case ET_QUAD4:
      case ET_QUAD8:
      case ET_QUAD9:
        if(dim_ == 3) {
          isSurfElem = true;
        }
        break;
      case ET_TET4:
      case ET_TET10:
      case ET_HEXA8:
      case ET_HEXA20:
      case ET_HEXA27:
      case ET_PYRA5:
      case ET_PYRA13:
      case ET_WEDGE6:
      case ET_WEDGE15:
        break;
      default:
        break;
      }
      
      // decide, what to do with the element
      if( isSurfElem ) {
        surfElems.insert( el );
      } else  {
        volRegionElems[el->regionId].Push_back(el);

        volRegionNodes[el->regionId].insert(&el->connect[0],
                                            (&el->connect[0])+numNodes);
      }
    }


    std::map<RegionIdType, StdVector<Elem*> >::iterator regionElemIt, regionElemEnd;
    std::map<RegionIdType, std::set<UInt> >::iterator regionNodeIt;
    std::set<UInt>::iterator setIt, setEnd;
    UInt region = 0;

    numVolRegions = volRegionElems.size();
    volElemNodes_.Resize(numVolRegions);
        
    regionElemIt = volRegionElems.begin();
    regionElemEnd = volRegionElems.end();
    regionNodeIt = volRegionNodes.begin();

    for( ; regionElemIt != regionElemEnd; regionElemIt++, regionNodeIt++, region++)
    {
      volElems_.Push_back(regionElemIt->second);
      volRegionIds_.Push_back(regionElemIt->first);
            
      setIt = regionNodeIt->second.begin();
      setEnd = regionNodeIt->second.end();
            
      for( ; setIt != setEnd; setIt++)
      {
        volElemNodes_[region].insert(*setIt);
      }
            
    }

    // Call creation of surface elements
    std::map<UInt, SurfElem* > mappedSurfElems;
    CreateSurfaceElements( surfElems, mappedSurfElems );
    
    // Iterate over all surface elements and put their nodes and elements
    // according to their region id
    std::map<UInt, SurfElem*>::iterator surfElemIt;
    for( surfElemIt = mappedSurfElems.begin();
         surfElemIt != mappedSurfElems.end();
         surfElemIt++ ) {
      SurfElem * surfEl = surfElemIt->second;
      UInt elemNum = surfElemIt->first;
      orderedElems_[elemNum-1] = surfEl;
      LOG_DBG3(gridcfs) << "Adding element #" << elemNum
                         << " to list of surface elements"; 
      if( surfEl->regionId != NO_REGION_ID ) {
        surfRegionElems[surfEl->regionId].Push_back( surfEl );
        surfRegionNodes[surfEl->regionId].insert( &surfEl->connect[0],
                                                  (&surfEl->connect[0]) + numNodes );
      }
    }
    numSurfRegions = surfRegionElems.size();
    surfElemNodes_.Resize(numSurfRegions );
                                  
    regionElemIt = surfRegionElems.begin();
    regionElemEnd = surfRegionElems.end();
    regionNodeIt = surfRegionNodes.begin();
    region = 0;

    for( ; regionElemIt != regionElemEnd; regionElemIt++, regionNodeIt++, region++)
    {
      surfElems_.Push_back(regionElemIt->second);
      surfRegionIds_.Push_back(regionElemIt->first);
            
      setIt = regionNodeIt->second.begin();
      setEnd = regionNodeIt->second.end();
            
      for( ; setIt != setEnd; setIt++)
      {
        surfElemNodes_[region].insert(*setIt);
      }
            
    }

    isInitialized_ = true;

    // Select nodes / elements according to the users specification in the
    // parameter file
    CreateUserDefinedNodesElems();

    // print information to file
    PrintGridInfo();
  }
    

  // ======================================================
  // GENERAL GRID INFORMATION
  // ======================================================

  
  UInt GridCFS::GetNumElemOfType( FEType type ) {
    return numElemTypes_[type];
  }

  void GridCFS::AddNodes(const UInt numNodes)
  { 
    coords_.Resize(this->numNodes_ + numNodes);
    numNodes_ += numNodes;
  }
    

  void GridCFS::SetNodeCoordinate(const UInt inode, const Point & rfPoint)
  {
    if ( inode > numNodes_ || inode < 0 ) {
      EXCEPTION( "GridCFS: There are only " << numNodes_
                 << " nodes in the grid. You wanted to set coordinates for "
                 << "node number " << inode );
    }
        
    if ( (dim_ == 2) && (rfPoint[2] != 0) ) {
      EXCEPTION(  "GridCFS: Dimension of grid is 2D. "
                  << "But you wanted to set a 3D coordinate for"
                  << "node number " << inode );
    }

    coords_[inode-1] = rfPoint;
  }  

  void GridCFS::SetNodeCoordinate(const UInt inode, const Vector<Double> & rfPoint)
  {
    if ( inode > numNodes_ || inode < 0 ) {
      EXCEPTION( "GridCFS: There are only " << numNodes_
                 << " nodes in the grid. You wanted to set coordinates for "
                 << "node number " << inode );
    }
        
    if ( (dim_ == 2) && (rfPoint[2] != 0) ) {
      EXCEPTION( "GridCFS: Dimension of grid is 2D. "
                 << "But you wanted to set a 3D coordinate for"
                 << "node number " << inode );
    }
        
    UInt idx = inode-1;
    coords_[idx][0] = rfPoint[0];
    coords_[idx][1] = rfPoint[1];
    coords_[idx][2] = rfPoint[2];
  }  

  
  UInt GridCFS::GetDim() {
    return dim_;
  }

  
  UInt GridCFS::GetNumNodes(){
    return numNodes_;
  }

  
  UInt GridCFS::GetNumNodes( const StdVector<RegionIdType> & regions ) {
    
    UInt numNodes = 0;
    Integer index = 0;

    for ( UInt i=0; i<regions.GetSize(); i++ ) {
      
      // look in volume regions
      index = volRegionIds_.Find(regions[i]);
      if ( index != -1 ) {
        numNodes += volElemNodes_[index].size();
      } else {
        
        // look in surface regions
        index = surfRegionIds_.Find(regions[i]);
        if ( index != -1 ) {
          numNodes += surfElemNodes_[index].size();
        } else {
          EXCEPTION("GridCFS: The region with id '" << regions[i]
                    << "' was not found in the grid!" );
        }
      }
    }
    
    return numNodes;
    
  }

  
  UInt GridCFS::GetNumNodes( const std::string & nodesName ) {

    UInt numNodes = 0;
    Integer index = namedNodeNames_.Find(nodesName);

    if ( index != -1 ) {
      numNodes = namedNodes_[index].GetSize();
    } else {
      EXCEPTION( "GridCFS: The Nodes with name '" << nodesName
                 << "' were not found in the grid!" );
    }

    return numNodes;
  }
  
  
  UInt GridCFS::GetNumElems() {
    return numElems_;
  }

  
  UInt GridCFS::GetNumSurfElems() {
    
    UInt numSurfElems = 0;
    
    for (UInt i=0; i<surfElems_.GetSize(); i++) {
      numSurfElems += surfElems_[i].GetSize();
    }
    
    return numSurfElems;
  }
  
  
  UInt GridCFS::GetNumVolElems() {
    
    UInt numVolElems = 0;
    
    for (UInt i=0; i<volElems_.GetSize(); i++) {
      numVolElems += volElems_[i].GetSize();
    }
    
    return numVolElems;
  }
  
  
  UInt GridCFS::GetNumElems( const StdVector<RegionIdType> & regions ){

    
    UInt numElems = 0;
    Integer index = 0;
    
    for ( UInt i=0; i<regions.GetSize(); i++ ) {
      
      
      // check if region Id is ALL_REGIONS
      if ( regions[i] == ALL_REGIONS ) {

        // iterate over all volume elements
        for ( UInt i = 0; i < volElems_.GetSize(); i++) {
          numElems += volElems_[index].GetSize();
        }
        
        // iterate over all surface elements
        for ( UInt i = 0; i < surfElems_.GetSize(); i++) {
          numElems += surfElems_[index].GetSize();
        }
        
      } else {
        
        
        // look in volume regions
        index = volRegionIds_.Find(regions[i]);
        if ( index != -1 ) {
          numElems += volElems_[index].GetSize();
        } else {
          
          // look in surface regions
          index = surfRegionIds_.Find(regions[i]);
          if ( index != -1 ) {
            numElems += surfElems_[index].GetSize();
          } else {
            EXCEPTION( "GridCFS: The region with id '" << regions[i]
                       << "' was not found in the grid!" );
          }
        }
      }
    }
    return numElems;
    
  }
  
  
  void GridCFS::AddNamedNodes( std::string name, StdVector<UInt> & nodeNums)
  {
    // Check if entities with given name exist already
    if( nameTypeMap_.find( name) != nameTypeMap_.end() ) {
      EXCEPTION( "Entities with name " << name 
                 << " are already defined" );
    }
    namedNodeNames_.Push_back(name);
    namedNodes_.Push_back(nodeNums);
    nameTypeMap_[name] = EntityList::NAMED_NODES;
  }

  void GridCFS::AddNamedElems( std::string name, StdVector<UInt> & elemNums)
  {
    // Check if entities with given name exist already
    if( nameTypeMap_.find( name) != nameTypeMap_.end() ) {
      EXCEPTION( "Entities with name " << name 
                 << " are already defined" );
    }
    namedElemNames_.Push_back(name);
    namedElems_.Push_back(elemNums);
    nameTypeMap_[name] = EntityList::NAMED_ELEMS;
    
    // get unique node number of elements
    UInt size = elemNums.GetSize();
    std::set<UInt> nodes;
    StdVector<UInt> nodeVec;
    for( UInt i = 0; i < size; i++ ) {
      const StdVector<UInt> & connect = orderedElems_[elemNums[i]-1]->connect;
      nodes.insert( connect.Begin(), connect.End() );
    }
    nodeVec.Resize( nodes.size() );
    std::copy( nodes.begin(), nodes.end(), nodeVec.Begin() );
    namedElemNodes_.Push_back( nodeVec );
  }

  void GridCFS::GetListNodeNames( StdVector<std::string> & nodeNames) {
    nodeNames = namedNodeNames_;
  }

  
  void GridCFS::GetListElemNames( StdVector<std::string> & elemNames) {
    elemNames = namedElemNames_;
  }
  
  // ======================================================
  // NODE ACCESS FUNCTIONS
  // ======================================================
  
  void GridCFS::GetNodesByName( StdVector<UInt> & nodeList,
                                const std::string & name ) {

    // Check if entities with given name exists at all
    if( nameTypeMap_.find( name) == nameTypeMap_.end() ) {
      EXCEPTION( "Entities with name " << name 
                 << " are already defined" );
    }
    
    // check, which entity type the name belongs to
    EntityList::DefineType defType = nameTypeMap_[name];
    Integer index = -1;
        
    switch( defType ) {
    
    case EntityList::REGION:
      GetNodesByRegion( nodeList, RegionNameToId( name ) );
      break;
    
    case EntityList::NAMED_NODES:
      index = namedNodeNames_.Find(name);
      if ( index != -1 ) {
        nodeList = namedNodes_[index];
      } else {
        EXCEPTION( "GridCFS: There are no nodes with name '" << name
                   << "' in the grid!" );
      }
      break;
    
    case EntityList::NAMED_ELEMS:
      index = namedElemNames_.Find(name);
      if ( index != -1 ) {
        nodeList = namedElemNodes_[index];
      } else {
        EXCEPTION( "GridCFS: There are no nodes with name '" << name
                   << "' in the grid!" );
      }
      
      break;
      
    default:
      EXCEPTION( "Can obtain nodes only for one region, named elements and "
                 << "named nodes" );
    }
  }

  void GridCFS::GetNodesByRegion( StdVector<UInt> & nodeList,
                                  const RegionIdType regionId ) {

    Integer index = 0;
    
    // look in volume regions
    index = volRegionIds_.Find(regionId);
    if ( index != -1 ) {
      nodeList.Resize(volElemNodes_[index].size());
      std::copy(volElemNodes_[index].begin(),
                volElemNodes_[index].end(),
                nodeList.Begin());
    } else {
      
      // look in surface regions
      index = surfRegionIds_.Find(regionId);
      if ( index != -1 ) {
        nodeList.Resize(surfElemNodes_[index].size());
        std::copy(surfElemNodes_[index].begin(),
                  surfElemNodes_[index].end(),
                  nodeList.Begin());
      } else {
        EXCEPTION( "GridCFS: The region with id '" << regionId
                   << "' was not found in the grid!" );
      }
    }
    
  }
  
  
  void GridCFS::GetNodeCoordinate( Point & rfPoint,
                                   const UInt inode, 
                                   bool updated ) {
    
    if ( inode > numNodes_ || inode < 0 ) {
      EXCEPTION( "GridCFS: There are only " << numNodes_
                 << " nodes in the grid. You requested coordinates for "
                 << "node number " << inode <<". Go check your mesh file!" );
    }

    rfPoint = coords_[inode-1];

    if( updated == true && deltCoords_.GetSize() != 0 ) {
      rfPoint += deltCoords_[inode-1];
    }


  }
  
  
  void GridCFS::GetNodeCoordinate( Vector<Double> & rfPoint,
                                   const UInt inode, 
                                   bool updated ) {

    if ( inode > numNodes_ || inode < 0 ) {
      EXCEPTION( "GridCFS: There are only " << numNodes_
                 << " nodes in the grid. You requested coordinates for "
                 << "node number " << inode <<". Go check your mesh file!" );
    }

    UInt idx = inode-1;
    rfPoint.Resize(dim_);
    rfPoint[0] = coords_[idx][0];
    rfPoint[1] = coords_[idx][1];
    if( dim_ == 3 ) {
      rfPoint[2] = coords_[idx][2];
    }

  }
    
  // ======================================================
  // ELEMENT ACCESS FUNCTIONS
  // ======================================================
  
  void GridCFS::AddElems(UInt nElems)
  {
    orderedElems_.Resize(numElems_ + nElems);

    UInt i=0;
    UInt idx=numElems_;

    for(; i<nElems; i++, idx++)
    {
      orderedElems_[idx] = new Elem();
      orderedElems_[idx]->elemNum = idx+1;
    }

    numElems_ = idx;
  }
    
      
  void GridCFS::SetElemData(UInt ielem,
                            FEType type,
                            RegionIdType region,
                            const UInt* connect)
  {
    UInt idx=ielem-1;
    Elem* el = orderedElems_[idx];
    UInt d = 2;
    UInt numNodes = NUM_ELEM_NODES[type];

    numElemTypes_[type]++;
        
    switch(type)
    {
    case ET_LINE2:
      el->ptElem = ptL1;
      break;
    case ET_LINE3:
      el->ptElem = ptL2;
      isQuadratic_ = true;
      break;
    case ET_TRIA3:
      el->ptElem = ptTr1;
      break;
    case ET_TRIA6:
      el->ptElem = ptTr2;
      isQuadratic_ = true;
      break;
    case ET_QUAD4:
      el->ptElem = ptQ1;
      break;            
    case ET_QUAD8:
      el->ptElem = ptQ2;
      isQuadratic_ = true;
      break;
    case ET_QUAD9:
      el->ptElem = NULL;
      break;
    case ET_TET4:
      d=3;
      el->ptElem = ptTet1;
      break;
    case ET_TET10:
      d=3;
      el->ptElem = ptTet2;
      isQuadratic_ = true;
      break;
    case ET_HEXA8:
      d=3;
      el->ptElem = ptHexa1;
      break;
    case ET_HEXA20:
      d=3;
      el->ptElem = ptHexa2;
      isQuadratic_ = true;
      break;
    case ET_HEXA27:
      d=3;
      el->ptElem = NULL;
      break;
    case ET_PYRA5:
      d=3;
      el->ptElem = ptPyra1;
      break;
    case ET_PYRA13:
      d=3;
      el->ptElem = ptPyra2;
      isQuadratic_ = true;
      break;
    case ET_WEDGE6:
      d=3;
      el->ptElem = ptWedge1;
      break;
    case ET_WEDGE15:
      d=3;
      el->ptElem = ptWedge2;
      isQuadratic_ = true;
      break;
    default:
      break;
    }

    if((dim_ == 2) && (d == 3))
    {
      EXCEPTION( "GridCFS: Cannot add 3D element type "
                 << ELEM_TYPE_NAMES[type]
                 << " to 2D grid." );
    }

    el->regionId = region;
    el->connect.Resize(numNodes);
    memcpy(&el->connect[0], connect, numNodes*sizeof(UInt));
  }
    
  void GridCFS::GetElemData(const UInt ielem,
                            FEType & type,
                            RegionIdType & region,
                            UInt* connect) const
  {
 #ifdef DEBUG
    if ( ielem > numElems_ || ielem < 0) {  
      EXCEPTION( "GridCFS: There are only " << numElems_ 
                 << " elements in the grid! You requested element number " 
                 << ielem << ". Go check your mesh file!" );
    }
    if ( orderedElems_[ielem-1] == NULL ) {
      EXCEPTION( "Element with Nr. " << ielem << " is not contained in mesh!" );
    }
 #endif

    UInt numNodes;

    type = orderedElems_[ielem-1]->ptElem->feType();
    region = orderedElems_[ielem-1]->regionId;
    numNodes = NUM_ELEM_NODES[type];
    memcpy(connect, &orderedElems_[ielem-1]->connect[0], numNodes*sizeof(UInt));
    
  }
    
  void GridCFS::SetElemRegion(UInt ielem, RegionIdType region)
  {
 #ifdef DEBUG
    if ( ielem > numElems_ || ielem < 0) {  
      EXCEPTION( "GridCFS: There are only " << numElems_ 
                 << " elements in the grid! You requested element number " 
                 << ielem << ". Go check your mesh file!" );
    }
    if ( orderedElems_[ielem-1] == NULL ) {
      EXCEPTION( "Element with Nr. " << ielem << " is not contained in mesh!" );
          
    }
 #endif

    orderedElems_[ielem-1]->regionId = region;
  }
    


  const Elem * GridCFS::GetElem( UInt elemNr ) {
    
 #ifdef DEBUG
    if ( elemNr > numElems_ || elemNr < 0) {  
      EXCEPTION( "GridCFS: There are only " << numElems_ 
                 << " elements in the grid! You requested element number " 
                 << elemNr << ". Go check your mesh file!" );
    }
    if ( orderedElems_[elemNr-1] == NULL ) {
      EXCEPTION( "Element with Nr. " << elemNr << " is not contained in mesh!" );
    }
 #endif
   
    return orderedElems_[elemNr-1];

  }


  
  void GridCFS::GetElems( StdVector<Elem*> & elems, 
                          const RegionIdType regionId ) {
    elems.Clear();
    
    // check if region Id is ALL_REGIONS
    if ( regionId == ALL_REGIONS ) {
      elems.Reserve( numElems_ );

      // iterate over all elements
      for ( UInt i = 0; i < numElems_; i++) {
        elems.Push_back(orderedElems_[i]);
      }
    } else {
      // look in volume regions
      Integer index = volRegionIds_.Find(regionId);
      if ( index != -1 ) {
        elems = volElems_[index];
      } else {    
        // look in surface regions
        index = surfRegionIds_.Find(regionId);
        if ( index != -1 ) {
          elems.Reserve( surfElems_[index].GetSize());
          for (UInt iElem=0; iElem<surfElems_[index].GetSize(); iElem++ ) {
            elems.Push_back(surfElems_[index][iElem]);
          }
        } else {    
          EXCEPTION( "GridCFS: The region with id '" << regionId
                     << "' was not found in the grid!" );
        }
      }
    }
  }
  

  
  void GridCFS::GetVolElems( StdVector<Elem*> & elems, 
                             const RegionIdType regionId ) {
    
    // check if region Id is ALL_REGIONS
    if ( regionId == ALL_REGIONS ) {
      elems.Reserve( GetNumVolElems() );
      for ( UInt i = 0; i < volElems_.GetSize(); i++) {
        for (UInt iElem = 0; iElem < volElems_[i].GetSize(); iElem++ ) {
          elems.Push_back(volElems_[i][iElem]);
        }
      }
    } else {
      Integer index = volRegionIds_.Find(regionId);
      if ( index != -1 ) {
        elems = volElems_[index];
      } else {
        EXCEPTION( "GridCFS: The volume region with id '" << regionId
                   << "' was not found in the grid!" );
      }
    }
  }
  
  
  void GridCFS::GetSurfElems( StdVector<SurfElem*> & elems, 
                              const RegionIdType regionId ) {
    
    Integer index = surfRegionIds_.Find(regionId);
    if ( index != -1 ) {
      UInt numElems = surfElems_[index].GetSize();
      
      for(UInt i=0; i<numElems; i++) 
      {    
        elems.Push_back(dynamic_cast<SurfElem*>(surfElems_[index][i]));
      }
      
    } else {
      EXCEPTION( "GridCFS: The surface region with id '" << regionId
                 << "' was not found in the grid!" );
    }
  }

  
  void GridCFS::GetElemsByName( StdVector<Elem*> & elems,
                                const std::string & elemsName ) {

    StdVector<UInt> elemNumbers;
    Integer index = namedElemNames_.Find(elemsName);
    

    if ( index != -1 ) {
      elemNumbers = namedElems_[index];
      elems.Resize( elemNumbers.GetSize() ); 
      for ( UInt i = 0; i < elemNumbers.GetSize(); i++ ) {
        elems[i] = orderedElems_[elemNumbers[i]-1 ];
      }
    } else {
      EXCEPTION( "GridCFS: There are no named elements with name '" 
                 << elemsName << "' in the grid!" );
    }
    
  }

  void GridCFS::GetElemNodes( StdVector<UInt> & connect, 
                              const UInt iElem ) {
    
    if ( iElem > numElems_ || iElem < 0) {  
      EXCEPTION( "GridCFS: There are only " << numElems_ 
                 << " elements in the grid! You requested element number " 
                 << iElem << ". Go check your mesh file!" );
    }
    
 #ifdef DEBUG
    if ( orderedElems_[iElem-1] == NULL ) {
      EXCEPTION( "Element with Nr. " << iElem << " is not contained in mesh!" );
    }
 #endif
    
    connect = orderedElems_[iElem-1]->connect;
    
  }

  
  void GridCFS::GetElemNodesCoord( Matrix<Double> & coordMat,  
                                   const StdVector<UInt> & connect,
                                   bool updated ) {

    coordMat.Resize(dim_, connect.GetSize());

    if( updated == true && deltCoords_.GetSize() != 0 ) {
      for (UInt k=0; k < connect.GetSize(); k++) {
        for (UInt actDim=0; actDim < dim_; actDim++) {
          coordMat[actDim][k] = coords_[connect[k]-1][actDim];
          //std::cerr << "\n coordMat before:\n" << coordMat << "\n\n";
          coordMat[actDim][k] += deltCoords_[connect[k]-1][actDim];
          //std::cerr << "\n coordMat after:\n" << coordMat << "\n\n";
        }
      }
    } else {
      for (UInt k=0; k < connect.GetSize(); k++)
      {
        for (UInt actDim=0; actDim < dim_; actDim++)
          coordMat[actDim][k] = coords_[connect[k]-1][actDim];
      }
        
    }


  }
  
  
  void GridCFS::GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                                     const StdVector<UInt> & nodeList,
                                     const StdVector<RegionIdType> 
                                     & regionIds ) {
    bool belongs2Interface;

    StdVector<UInt> map;
    Integer index = 0;
    
    // loop over all regionIDs
    for (UInt isd=0; isd<regionIds.GetSize(); isd++)
    {

      // get index for id in regionlist
      index = volRegionIds_.Find(regionIds[isd]);
      if ( index == -1 ) {
        EXCEPTION( "GetElemsNextToNodes: A region with id '" 
                   << regionIds[isd] << "' was not found in the list of "
                   << "of volume regions." );
      }
      // Get element of given region
      StdVector<Elem*> const & elems = volElems_[index];
        
      // loop over all elements in subdomain
      for (UInt iNS=0; iNS < elems.GetSize(); iNS++)
      {
        Elem *aux = elems[iNS];
        StdVector<UInt>  const & aux_connect = aux->connect;
        
        belongs2Interface = false;
        
        // check if any node is common in Interface
        for (UInt inode=0; inode<aux_connect.GetSize(); inode++) {

          for (UInt nnode=0; nnode<nodeList.GetSize(); nnode++) {

            if (nodeList[nnode] == aux_connect[inode]) {
              belongs2Interface = true;
              break;
            }
          }
        }
        
        if (belongs2Interface)
          elemList.Push_back(elems[iNS]);
      }
    }
  }
  
  
  void GridCFS::GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                       const StdVector<Elem*> & surfElems,
                                       const StdVector<RegionIdType> 
                                       &neighRegions ) {
    EXCEPTION( "Not implemented" );
  }
    



  // ======================================================
  // MISCELLANEOUS
  // ======================================================
  
  void GridCFS::GetNodesOfElemList( StdVector<UInt> & nodeList,
                                    const StdVector<Elem*> & elemList,
                                    bool onlyLinNodes) {

    std::set<UInt> elemNodes;
    std::set<UInt>::iterator it;
    UInt iElem, iNode, numElemCorners;

    // First, create a set with node numbers of elements
    for ( iElem = 0; iElem < elemList.GetSize(); iElem++ ) {
      StdVector<UInt> const & connecth = elemList[iElem]->connect;
      
      if (onlyLinNodes == true)
        numElemCorners = elemList[iElem]->ptElem->GetNumCorners();
      else
        numElemCorners = connecth.GetSize();
      
      for ( iNode = 0; iNode < numElemCorners; iNode++ ) {
        elemNodes.insert(connecth[iNode]);
      }
    }
    
    // Then copy this set into the nodeList vector
    nodeList.Resize(elemNodes.size());
    iNode = 0;
    for ( it = elemNodes.begin(); it != elemNodes.end(); it++) {
      nodeList[iNode++] = *it;
    }
  }
  
  
  
  void GridCFS::GetRegionNames( StdVector<std::string> & regionNames ) {
        
    regionNames = regionNames_;
        
  }

  
  void GridCFS::SetNodeOffset( const StdVector<UInt>& nodes, 
                               const Vector<Double>& offsets ) {

    // Check if node offsets were already set
    if( deltCoords_.GetSize() == 0 ) {
      deltCoords_.Resize( coords_.GetSize() );
      deltCoords_.Init();
    }

    // Set delta coordinates
    for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
      Point actOffset;
      for( UInt iDim = 0; iDim < dim_; iDim++ ) {
        actOffset[iDim] = offsets[iNode*dim_ + iDim];
      }
      deltCoords_[nodes[iNode]-1] = actOffset;
    }
  }

  
  bool GridCFS::HasNodalOffset() {

    if( deltCoords_.GetSize() != 0 ) {
      return true;
    } else {
      return false;
    }
     
  }

  // =======================================================================
  // Helper Methods
  // =======================================================================
  
  void GridCFS::CreateSurfaceElements( std::set<Elem*>& elems,
                                       std::map<UInt, SurfElem*>& surfElems ) {

    LOG_TRACE(gridcfs) << "Starting to map surface elements";
    
    // 1.) Create vector of vector of elems
    StdVector<StdVector<UInt> > elemNrPerNode;
    UInt nrNodes, iRegion, iElem;
    Elem * ptVolElem = NULL;
    elemNrPerNode.Resize(numNodes_);
    elemNrPerNode.Init();

    // 2.) Iterate over all volume elements and add for each
    //     element node the element number 

    for ( iRegion = 0; iRegion < volElems_.GetSize(); iRegion++ ) {
      for ( iElem = 0; iElem < volElems_[iRegion].GetSize(); iElem++ ) {
        ptVolElem = volElems_[iRegion][iElem];
        
        nrNodes = ptVolElem->connect.GetSize();
        
        for (UInt iNode = 0; iNode < nrNodes; iNode++ ) {
          elemNrPerNode[ptVolElem->connect[iNode]-1].
            Push_back(ptVolElem->elemNum);
          
        } // loop over nodes
      } // loop over elements
    } // loop over regions


    // iterate over all temporary elements and convert them into
    // surface elements
    SurfElem *myElem;
    Elem* oldElem;
    
    std::set<Elem*>::iterator elIt;
    
    LOG_DBG(gridcfs) << "There are " << elems.size() << " surface element to be mapped";
    for( elIt = elems.begin(); elIt != elems.end(); elIt++ ) {
      
      oldElem = *elIt;
     // create new surface element
      myElem = new SurfElem();
      myElem->connect = oldElem->connect;
      myElem->regionId = oldElem->regionId;
      myElem->elemNum = oldElem->elemNum;
      myElem->ptElem = oldElem->ptElem;
      surfElems[myElem->elemNum] = myElem;

      // delete old volume element
      delete oldElem;
    }

    // 3.) Iterate over all surface elements and look for each
    //     element, if all of its nodes can be assigned to one or
    //     two neighbours

    UInt surfNodeNr = 0;
    UInt elemsFound = 0;
    UInt elemsAssigned = 0;
    
    std::map<UInt,SurfElem*>::iterator surfElIt;
    for( surfElIt = surfElems.begin(); 
         surfElIt != surfElems.end(); 
         surfElIt++ ) {

      elemsAssigned = 0;

      myElem = surfElIt->second;

      // get number of nodes of surface element
      nrNodes = myElem->connect.GetSize();
      StdVector<UInt> const & connect = 
        myElem->connect;

      // get first node of surface element
      surfNodeNr = myElem->connect[0];

      // make loop over all elements, which have first node
      // of surface element in common
      for (UInt iVolElem = 0; 
      iVolElem < elemNrPerNode[surfNodeNr-1].GetSize(); iVolElem++ ) {
        elemsFound = 1;         

        // look if this element is also defined by the other nodes
        // of the surface element
        for (UInt iNode = 1; iNode < nrNodes; iNode++ ) {

          UInt index = connect[iNode]-1;
          for (UInt iElem2 = 0 ; iElem2 < elemNrPerNode[index].GetSize(); 
          iElem2++ ) {

            if ( elemNrPerNode[index][iElem2] == 
              elemNrPerNode[surfNodeNr-1][iVolElem] ) {
              elemsFound++;
              break;
            }

          } // loop over all elements of other nodes
        } // loop over all other nodes

        if ( elemsFound == nrNodes ) {

          ptVolElem = orderedElems_[elemNrPerNode[surfNodeNr-1][iVolElem]-1];

          if ( elemsAssigned == 0 ) {
            myElem->ptVolElem1 = ptVolElem;
          }
          else {
            myElem->ptVolElem2 = ptVolElem;
          }

          elemsAssigned++;
        }
      } // loop over element numbers of first node

      // sanity check (avoid the impossible ;-)
      if ( elemsAssigned > 2 ) {
        EXCEPTION( "Found " << elemsAssigned
                   << " volume elements for surface element no. "
                   << myElem->elemNum );
      }

    } // loop over surface elements


    // 4.) Iterate over all surface elements and calculate surface
    //     flag by comparing the directed and the undirected surface
    //     normal. If both differ, the surfaceNormalSign = -1, otherwise 1.
    Vector<Double> normalUndefSign, normalDefSign;
    Double sign;

    for( surfElIt = surfElems.begin(); 
         surfElIt != surfElems.end(); 
         surfElIt++ ) {

      myElem = surfElIt->second;

      // check, if each surface element has at least one volume neighbour
      if ( myElem->ptVolElem1 == NULL ) {
        //  EXCEPTION( "Pointer to first volume element is NULL for surface"
        //                    << " element no. "  
        //                    << surfElems_[iRegion][iSurfElem]->elemNum << ".\n"
        //                    << "Please check your mesh-file!" );
        //         }
        myElem->normalSign = 0;
      } else {

        CalcSurfNormal( normalUndefSign, *myElem, false );
        CalcSurfNormalOutOfVol( normalDefSign, 
                                *myElem,
                                *myElem->ptVolElem1,
                                false );


        // Check if all entries have the same sign by calulating
        // a scalar product between both vectors.
        // If it is positive, they point in the smae direction,
        // otherwise an angle of 180 lies in between.
        sign = normalUndefSign * normalDefSign;

        if ( sign > 0.0 ) {
          myElem->normalSign = 1;
        } else {
          myElem->normalSign = -1;
        }
      }
    }
    
  }

  
  void GridCFS::PrintGridInfo() const {
    
    std::string help, empty;
    if( !Info) return;
    Info->PrintF("", "\n--- GridCFS: Region Map ---\n\n");
    Info->PrintF("", "ID\t| Name\n");
    Info->PrintF("", "-------------------\n");
    
    for( UInt i = 0; i < regionNames_.GetSize(); i++ ) {
      help = GenStr(i);
      help += "\t| ";
      help += regionNames_[i];
      help += '\n';
      Info->PrintF("",help.c_str());
    }
  }
  
  
  void GridCFS::CalcSurfNormal( Vector<Double> & n, 
                                const Elem & surfElem,
                                bool updated ) {
    
    //compute normal vector
    Matrix<Double>  ptCoord;

    GetElemNodesCoord(ptCoord, surfElem.connect, updated );
    UInt surfCorners = surfElem.ptElem->GetNumCorners();
 
    // Check for dimension:
    if (surfElem.ptElem->GetDim() == 1) {
  
      // 1. step: compute vector perpendicular to line element
      // but without defined sign
      Double dx  = ptCoord[0][1] - ptCoord[0][0];
      Double dy  = ptCoord[1][1] - ptCoord[1][0];
      Double len = sqrt(dx*dx + dy*dy);
      if (len <= 0.0) {
        EXCEPTION( "length of normal vector is zero!" );
      }
      n.Resize(2);
      n[0] = dy/len;
      n[1] = -dx/len;
    }
    else {
      // 1. step: compute vector perpendicular to surface element
      // but without defined sign
      
      //compute the two vectors in the plane
      Vector<Double> vec1(3), vec2(3);
      for (UInt i=0; i<3; i++) {
        vec1[i] = ptCoord[i][1]             - ptCoord[i][0];
        vec2[i] = ptCoord[i][surfCorners-1] - ptCoord[i][0];
      }
      //compute cross product
      n.Resize(3);
      n[0] = vec1[1] * vec2[2] - vec1[2]*vec2[1];
      n[1] = vec1[2] * vec2[0] - vec1[0]*vec2[2];
      n[2] = vec1[0] * vec2[1] - vec1[1]*vec2[0];
      //normalize the length to 1
      Double length = n.NormL2();
      n /= length;
      
    }
  }

  
  void GridCFS::CalcSurfNormalOutOfVol(Vector<Double> & n,
                                       const Elem & surfElem,
                                       const Elem & volElem,
                                       bool updated )
  {

    //compute normal vector
    Matrix<Double>  ptVolCoord, ptSurfCoord;

    // First, calculate undefined normal
    CalcSurfNormal(n, surfElem, updated );

    GetElemNodesCoord(ptSurfCoord, surfElem.connect, updated );
    GetElemNodesCoord(ptVolCoord, volElem.connect, updated );
    

    UInt volCorners = volElem.ptElem->GetNumCorners();
 
    // Check for dimension:
    // A 2D volume element has only one face
    // -> we are in 2D
    if ( n.GetSize() == 2 ) {
  
      // compute direction
      
      Integer indexNode1=-1;
      Integer indexNode2=-1;
      
      for(UInt actNode=0; actNode < volCorners; actNode++)
      {
        if (volElem.connect[actNode] == surfElem.connect[0])
          indexNode1 = actNode;
        if (volElem.connect[actNode] == surfElem.connect[1])
          indexNode2 = actNode;
      }
      // if not clockwise orientation of nodes (difference of node indizes is -1)    
      if (indexNode1==-1 || indexNode2==-1)
        EXCEPTION("Nodes of neighbouring element not found!" );
    
    
      // counterclockwise orientation of nodes (difference of node indizes is +1)
      if ( ( indexNode2-indexNode1  == -1 || 
             (indexNode2-indexNode1)- (Integer) volCorners == -1 ) ) {
        n *= -1;
      }
    
      else
        // counterclockwise orientation of nodes (difference of node indizes is +1)

        if (! (indexNode2-indexNode1 == 1 || 
               (indexNode2-indexNode1)+volCorners == 1) )
          EXCEPTION("Nodes of interface don't lie beneath each other in neighbouring element!" );
    }
  
    else {
    
      // compute direction
    
      // find first common vertex index
      Integer firstCommonIndex = -1;
      for (UInt i=0; i<volCorners; i++)
        if (volElem.connect[i] == surfElem.connect[0]){
          firstCommonIndex = i;
          break;
        }
    
      // calculate barycenter of volume element
      Vector<Double> barycenter(3);
      for (UInt i=0; i<volCorners; i++){
        barycenter[0] += ptVolCoord[0][i];
        barycenter[1] += ptVolCoord[1][i];
        barycenter[2] += ptVolCoord[2][i];
      }

      barycenter /= volCorners;

      // check, if scalar product with vector (going from barycenter to
      // common edge) and perpendicular vector  are pointing in same direction
      Vector<Double> innerVec(3);
      Double product = 0;
      innerVec[0] = ptVolCoord[0][firstCommonIndex] - barycenter[0];
      innerVec[1] = ptVolCoord[1][firstCommonIndex] - barycenter[1];
      innerVec[2] = ptVolCoord[2][firstCommonIndex] - barycenter[2];
    
      product = innerVec * n;
      if (product < 0) {
        n *= -1;
      }

    }
  }
  
  Double GridCFS::CalcVolumeOfRegion( const RegionIdType regionId,
                                      bool isaxi,
                                      bool updated ) {

    StdVector<Elem*> elems;
    Matrix<Double> cornerCoords;
    Double volume = 0.0;

    GetElems(elems,regionId);

    for( UInt i = 0; i < elems.GetSize(); i++ ) {
      GetElemNodesCoord(cornerCoords, elems[i]->connect, updated );
      volume += elems[i]->ptElem->CalcVolume(cornerCoords, isaxi);
    }

    return volume;
  }
  
  void GridCFS::GetGlobalElemMidPoint( UInt elemNum, Vector<Double>& coord ) {

    if( elemNum > numElems_ ) {
      EXCEPTION("Eleement number " << elemNum << " is bigger than total "
                << "number of elements within the grid" );
    }
    Vector<Double> locMidPoint;
    Matrix<Double> connectCoord;
    
    Elem * actElem = orderedElems_[elemNum-1];
    BaseFE * ptFE = actElem->ptElem;
    
    GetElemNodesCoord( connectCoord, actElem->connect, false );
    ptFE->GetCoordMidPoint(locMidPoint);
    ptFE->Local2GlobalCoord( coord, locMidPoint, 
                             connectCoord, actElem );
    
  }


#ifdef ADAPTGRID
  
  void GridCFS::putNodesFromGrid_RG(grd::MultilevelGrid * grid,
                                    const UInt level)
  {

    Integer maxnumnodes = (*grid).getNoOfVertices();
    maxnumnodes_=maxnumnodes;
    ptCoordinate_=new Point<dim>[maxnumnodes];

    typedef std::list<grd::Vertex*>::iterator VerI;  
  
    Double * ps;
    Integer ilev, i=0;
    std::cout << "\t\033[32m no. of vertices: \033[0m "
              << (*grid).getNoOfVertices() << std::endl;

    Integer topLevel = grid->getTopLevel();
    for (ilev=0; ilev<=topLevel; ilev++) {
      std::list<grd::Vertex*>*le=(*grid).getGridLevel(ilev)->getVertexList();
      for (VerI p=le->begin(); p!=le->end(); ++p) {
        Integer index = (*p)->getId();
        index--;
        if ( index >= maxnumnodes ) {
          std::cerr << " ERROR: catastrophic error, index overflow\n "
                    << " The index: " << index << "  the maxnumnodes: "
                    << maxnumnodes << '\n';
        }
        ps=(*p)->getPosition();
        ptCoordinate_[index][0]=ps[0];
        ptCoordinate_[index][1]=ps[1];
        if (dim==3)
          ptCoordinate_[index][2]=ps[2];
        i++;
      }
    }
  } // end of function putNodesFromGrid_RG


  template<>
  void GridCFS<2>::putElemsFromGrid_RG(grd::MultilevelGrid * grid,
                                       const Integer level)
  {

    typedef std::list<grd::Element*>::iterator ElmI;
    ElmI p;

    Integer i, j, nnodes,type;

    std::list<grd::Element*> * le;
    std::list<grd::Element*> ** lt;

    // Element Map
    if (!elemMap_.empty()) {
      for (i = 0; i < elemMap_.GetSize(); i++) {
        ElementMap* tmp = elemMap_[i];
        delete tmp;
      }
      elemMap_.Clear();
    } 
 
    Integer noOfLevels=grid->getNoOfLevels();
    for (j=0; j< noOfLevels; j++)
    {
      lt=grid->getGridLevel(j)->getElementList();
      type=0;

      while(lt[type]) {
        for ( p=lt[type]->begin() ; p!= lt[type]->end(); ++p)
        {
          if ((*p)->isRegular())
          {
            Elem * el=new Elem();
            // Element maping
            ElementMap* tmpMap = new ElementMap;
             
            Integer nnodes=(*p)->getNoOfVertices();
            Integer etype = (*p)->type();
             
            switch(etype)
            {
            case GRD_TRIANGLE:
              el->ptElem=ptTr1;
              break;
            case GRD_QUADRANGLE:
              el->ptElem=ptQ;
              break;
            default:
              EXCEPTION("Unknown type of element in GridRG" );
              break;
            }
             
            (*el).connect.Resize(nnodes);
             
            for (i=0; i<nnodes; i++)
            {   
              (*el).connect[i]=((*p)->getVertex(i))->getId();    
            }
             
            Integer sd=(*p)->getValue();
            if (sd >= sd_.GetSize()) {
              EXCEPTION(" Value in element from Grid_RG is incorrect" );
            }
            (*el).namesd=sd_[sd];

            elems_[sd].Push_back(el);
             
            // put info in elemMap
            Integer position = elems_[sd].GetSize()-1;
            tmpMap->map.Resize(1);
            tmpMap->map[0] = position;
            tmpMap->sd = sd;
            elemMap_.Push_back(tmpMap);
             
          } // if isRegular
          else if ((*p)->isIrregular()) {
            ElementMap* tmpMap = new ElementMap;
            tmpMap->sd = (*p)->getValue();
            grd::ConformingClosure closure;
            typedef grd::ConformingClosure::triangleIterator TriI;
            typedef grd::ConformingClosure::quadrangleIterator QuadI;
            (*p)->close(closure);
           
            // Process closing triangles
            for (TriI tri = closure.beginTriangle();
                 tri != closure.endTriangle(); ++tri) {
              Elem * el=new Elem();
              Integer nnodes=(*tri)->getNoOfVertices();
              el->ptElem=ptTr1;
              (*el).connect.Resize(nnodes);
              for (i=0; i<nnodes; i++)
              {
                (*el).connect[i]=((*tri)->getVertex(i))->getId();
              }
              Integer sd=(*tri)->getValue();
              if (sd >= sd_.GetSize()) {
                EXCEPTION(" Value in element from Grid_RG is incorrect" );
              }
              (*el).namesd=sd_[sd];
              elems_[sd].Push_back(el);
             
              // maping
              Integer position = elems_[sd].GetSize()-1;
              tmpMap->map.Push_back(position);
            } // for tri
           
            // Process now the quads
            for (QuadI quad = closure.beginQuadrangle();
                 quad != closure.endQuadrangle(); ++quad) {
              Elem * el=new Elem();
              Integer nnodes=(*quad)->getNoOfVertices();
              el->ptElem=ptQ;
              (*el).connect.Resize(nnodes);
              for (i=0; i<nnodes; i++)
              {
                (*el).connect[i]=((*quad)->getVertex(i))->getId();
              }
              Integer sd=(*quad)->getValue();
              if (sd >= sd_.GetSize()) {
                EXCEPTION(" Value in element from Grid_RG is incorrect" );
              }
              (*el).namesd=sd_[sd];
              elems_[sd].Push_back(el);

              // maping
              Integer position = elems_[sd].GetSize() - 1;
              tmpMap->map.Push_back(position);
            } // for quad
            // update element map list
            elemMap_.Push_back(tmpMap);
          } // else if isIrregular
        } // for element
        type++;
      } // end while(); list of elements types
    } // for level
 
    // Info->PrintF("Total number of elements (only for first subdomain): %i", elems_[0].GetSize());
    
    FormNeighborsLists();
     
  } // end of function 


  template<>
  void GridCFS<3>::putElemsFromGrid_RG(grd::MultilevelGrid * grid,
                                       const Integer level)
  {

    typedef std::list<grd::Element*>::iterator ElmI;
    ElmI p;

    StdVector<grd::Element*> tets;

    Integer i, j, nnodes,type;

    std::list<grd::Element*> * le;
    std::list<grd::Element*> ** lt;

    // Init Tets buffer
    tets.Resize(4);
    for (i = 0; i < 4; i++)
    {
      tets[i] = new grd::Tetrahedron;
    }

    // Element Map
    if (!elemMap_.empty()) {
      for (i = 0; i < elemMap_.GetSize(); i++) {
        ElementMap* tmp = elemMap_[i];
        delete tmp;
      }
      elemMap_.Clear();
    } 
 
    Integer i1,i2,lev;
    Integer noOfLevels=grid->getNoOfLevels();
    for (j=0; j< noOfLevels; j++)
    {
      lt=grid->getGridLevel(j)->getElementList();
      type=2; // take only volume elementes

      while(lt[type])
      {
        for ( p=lt[type]->begin() ; p!= lt[type]->end(); ++p)
        {
          if ((*p)->isRegular())
          {
            Elem * el=new Elem();
            // Element maping
            ElementMap* tmpMap = new ElementMap;

            Integer nnodes=(*p)->getNoOfVertices();
            Integer etype = (*p)->type();

            Integer sd,position,nnds;

            switch(etype)
            {
            case GRD_TRIANGLE:
              break;
            case GRD_QUADRANGLE:
              break;
            case GRD_TETRAHEDRON:
              el->ptElem=ptTet;
              (*el).connect.Resize(nnodes);

              for (i=0; i<nnodes; i++) {
                (*el).connect[i]=((*p)->getVertex(i))->getId();
              }

              sd=(*p)->getValue();

              if (sd >= sd_.GetSize())
                EXCEPTION("Value in element from Grid_RG is incorrect");

              (*el).namesd=sd_[sd];

              elems_[sd].Push_back(el);

              // put info in elemMap
              position = elems_[sd].GetSize()-1;
              tmpMap->map.Resize(1);
              tmpMap->map[0] = position;
              tmpMap->sd = sd;

              elemMap_.Push_back(tmpMap);
              break;
            case GRD_OCTAHEDRON:
              (*p)->getTetras(tets);
              Elem * elT[4];

              Integer it;
              nnds=4;

              // loop over tetrahedrals of octahedron
              for (it=0; it<4; it++)
              {
                elT[it]=new Elem();
                elT[it]->ptElem=ptTet;
                elT[it]->connect.Resize(nnds);
                // copy of connection array
                for (i=0; i<nnds; i++)
                {
                  elT[it]->connect[i]=(tets[it]->getVertex(i))->getId();
                }

                sd=(*p)->getValue();

                if (sd >= sd_.GetSize()) 
                  EXCEPTION(" Value in element from Grid_RG is incorrect" );
                elT[it]->namesd=sd_[sd];

                elems_[sd].Push_back(elT[it]);
                // mapping
                Integer position = elems_[sd].GetSize()-1;
                tmpMap->map.Push_back(position);
                tmpMap->sd=sd;
              } // end: loop over tetrahedrals of octahedron

              // update element map list
              elemMap_.Push_back(tmpMap);
              break;
            case GRD_HEXAHEDRON:
              el->ptElem=ptHexa;
              break;
            default:
              EXCEPTION("Unknown type of element in GridRG" );
              break;
            } // end of switch
          }
          else if ((*p)->isIrregular()) {
            ElementMap* tmpMap = new ElementMap;
            tmpMap->sd = (*p)->getValue();
            grd::ConformingClosure closure;
            typedef grd::ConformingClosure::tetrahedronIterator TetI;
            (*p)->close(closure);
       
            // Process closing tetrahedrons
            for (TetI tet = closure.beginTetrahedron(); tet != closure.endTetrahedron(); ++tet)
            {
              Elem * el=new Elem();
              Integer nnodes=(*tet)->getNoOfVertices();
              el->ptElem=ptTet;
              (*el).connect.Resize(nnodes);
              for (i=0; i<nnodes; i++)
              {
                (*el).connect[i]=((*tet)->getVertex(i))->getId();
              }
              Integer sd=(*tet)->getValue();
              if (sd >= sd_.GetSize()) 
                EXCEPTION(" Value in element from Grid_RG is incorrect" );
              (*el).namesd=sd_[sd];
              elems_[sd].Push_back(el);

              // maping
              Integer position = elems_[sd].GetSize()-1;
              tmpMap->map.Push_back(position);
            } // for tet
            // update element map list
            elemMap_.Push_back(tmpMap);
          } // else if isIrregular
        } // for element
        type++;
      } // end while(); list of elements types
    } // for level

    // clean buffer of tets
    if (!tets.empty()) 
    {
      for (i = 0; i < 4; i++)
      {
        delete tets[i];
      }
      tets.Clear();
    }


    //     Info-PrintF("", "Total number of elements: %i", elems_[0].GetSize());

    FormNeighborsLists();
  }

  void GridCFS::Refine(grd::MultilevelGrid& grid)
  {

    Integer i,j;
    Integer counter = 0;

    // Mesh refinement
    Integer noOfLevels = grid.getNoOfLevels();
    typedef std::list<grd::Element*>::iterator ElmI;

    Integer k;
    for (Integer j = 0; j < noOfLevels; j++)
    {
      grd::GridLevel* gridlv = grid.getGridLevel(j);
      list<grd::Element*>** lt = gridlv->getElementList();
      Integer type;
      if (dim==3) type=2;
      else type=0;
      while (lt[type])
      {
        for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
        {
          if (!(*p)->isRefined())
          {
            Integer  sde = (*p)->getValue();
            Integer flag = 0;
            ElementMap* map = elemMap_[counter];
                     
            Integer sdm = map->sd;
            if (sde != sdm) 
              EXCEPTION("Wrong number of subdomain");
                     
            for (i = 0; i < map->map.GetSize(); i++)
            {
              Integer elmId = map->map[i];
              flag = elems_[sdm][elmId]->refinementFlag;
              Integer numRefs=elems_[sdm][elmId]->refinementNumber;
              if (flag == 1) 
              {
                            
                (*p)->markForRefinement(numRefs-1);
                    
                break;
              }

              if (flag == -1)
              {
                (*p)->markForCoarsening();
                break;
              }
                         
            }
                     
            // update counter
            counter++;
          }
        } // for loop elems
        // next element type
        type++;
      } // type loop
    } // level loop

    grid.refine();      
  }

  void GridCFS::ReRefine(grd::MultilevelGrid& grid)
  {
    
    Integer i,j;
    Integer counter = 0;

    // Mesh refinement
    Integer topLevel = grid.getTopLevel();
    typedef std::list<grd::Element*>::iterator ElmI;

    Integer k;
    grd::GridLevel* gridlv = grid.getGridLevel(topLevel);
    list<grd::Element*>** lt = gridlv->getElementList();
    Integer type;
    if (dim==3) type=2;
    else type=0;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if (!(*p)->isRefined())
        {
          Integer  sde = (*p)->getValue();
          grd::Element *parent = (*p)->getParent();
                 
          Integer numRefs = parent->getNumOfRefinements(); 
          if (numRefs)
            (*p)->markForRefinement(numRefs-1);
        } 
      } // for loop elems
      // next element type
      type++;
    } // type loop
    grid.refine();     
  }

  void GridCFS::RefineUniform(grd::MultilevelGrid& grid)
  {
    Integer i,j;

    // Mesh refinement
    Integer noOfLevels = grid.getNoOfLevels();
    typedef std::list<grd::Element*>::iterator ElmI;

    for (Integer j = 0; j < noOfLevels; j++) {
      grd::GridLevel* gridlv = grid.getGridLevel(j);
      list<grd::Element*>** lt = gridlv->getElementList();
      Integer type = 0;
      if (dim == 3)
        type = 2;
      while (lt[type]) {
        for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p) {
          if (!(*p)->isRefined())
            (*p)->markForRefinement();
        } // for loop elems
        // next element type
        type++;
      } // type loop
    } // level loop

    grid.refine();
  }

  void GridCFS::SetRefinementFlag()
  {
    Integer i,j;
    for (i=0; i<sd_.GetSize(); i++) {
      for (j=0; j<elems_[i].GetSize(); j++) {
        elems_[i][j]->refinementFlag=true;
      }
    }
  }

#endif // end of ADAPTGRID

  void GridCFS::AddNode( const Point & coord, UInt & inode)
  {
    if(!isInitialized_)
      EXCEPTION("Cannot add node to uninitialized grid!");

    coords_.Push_back(coord);
    inode = ++numNodes_;
  }
  
  void GridCFS::AddNode( const Vector<Double> & coord, UInt & inode )
  {
    if(!isInitialized_)
      EXCEPTION("Cannot add node to uninitialized grid!");
    
    if(coord.GetSize() != 3)
      EXCEPTION("Node to be added has wrong dimension!");
    
    Point p;
        
    for(UInt i=0; i<3; i++)
      p[i] = coord[i];
        
    coords_.Push_back(p);
    inode = ++numNodes_;
  }
  
  void GridCFS::AddNodes( const StdVector< Point > & coords,
                          StdVector< UInt > & inodes)
  {
    if(!isInitialized_)
      EXCEPTION("Cannot add nodes to uninitialized grid!");

    Point p;
    UInt i, n;

    n=coords.GetSize();
    inodes.Resize(n);

    for(i=0; i<n; i++)
    {
      coords_.Push_back(coords[i]);
      numNodes_++;
      inodes[i] = numNodes_;
    }
  }

  void GridCFS::AddSurfaceElems( const RegionIdType regionid,
                                 const StdVector< SurfElem* > & surfelems,
                                 StdVector< UInt > & elemids)
  {
    UInt i, n;    
    Integer index = 0;
    UInt *ptConn;
    UInt numNodes;
    
    if(!isInitialized_)
      EXCEPTION("Cannot add surface elements to uninitialized grid!");

    UInt regionIdx = surfRegionIds_.Find(regionid);
      
    if(regionIdx == -1)
      EXCEPTION("Surface regionid not found!");

    n=surfelems.GetSize();
    elemids.Resize(n);

    for(i=0; i<n; i++)
    {
      // a check should be added to avoid insertions
      // of already existing elements
      surfelems[i]->regionId = regionid;
      numElems_++;
      surfelems[i]->elemNum = numElems_;
        
      orderedElems_.Push_back(surfelems[i]);
      surfElems_[regionIdx].Push_back(surfelems[i]);
      elemids[i] = numElems_;

      ptConn = &surfelems[i]->connect[0];
      numNodes = surfelems[i]->connect.GetSize();
      surfElemNodes_[regionIdx].insert(ptConn,
                                       ptConn+numNodes);
    }
  }

  void GridCFS::AddVolumeElems( const RegionIdType regionid,
                                const StdVector< Elem* > & volelems,
                                StdVector< UInt > & elemids)
  {
    UInt i, n;
    UInt *ptConn;
    UInt numNodes;

    if(!isInitialized_)
      EXCEPTION("Cannot add volume elements to uninitialized grid!");

    UInt regionIdx = volRegionIds_.Find(regionid);
      
    if(regionIdx == -1)
      EXCEPTION("Volume regionid not found!");


    n=volelems.GetSize();
    elemids.Resize(n);

    for(i=0; i<n; i++)
    {
      // a check should be added to avoid insertions
      // of already existing elements
      volelems[i]->regionId = regionid;
      numElems_++;
      volelems[i]->elemNum = numElems_;

      orderedElems_.Push_back(volelems[i]);
      volElems_[regionIdx].Push_back(volelems[i]);
      elemids[i] = numElems_;

      ptConn = &volelems[i]->connect[0];
      numNodes = volelems[i]->connect.GetSize();
      volElemNodes_[regionIdx].insert(ptConn,
                                      ptConn+numNodes);
    }
  }

  void GridCFS::AddSurfaceRegion( const std::string name,
                                  RegionIdType& regionid)
  {
    
    // Check if entities with given name exist already
    if( nameTypeMap_.find( name) != nameTypeMap_.end() ) {
      EXCEPTION( "Entities with name " << name 
                 << " are already defined" );
    }
    
    if(!isInitialized_)
      EXCEPTION("Cannot add a surface region to an uninitialized grid!");

    if(regionNames_.Find(name) != -1)
    {
      EXCEPTION("A surface region with the same name already exists!");
      regionid = -1;
    }

    regionNames_.Push_back(name);
    regionid = regionNames_.GetSize()-1;
    surfRegionIds_.Push_back(regionid);
    
    StdVector<Elem*> dummy_elems;
    surfElems_.Push_back(dummy_elems);
    
    std::set<UInt> dummy_nodes;
    surfElemNodes_.Push_back(dummy_nodes);                                      
  }
  
  void GridCFS::AddVolumeRegion( const std::string name,
                                 RegionIdType& regionid)
  {
    // Check if entities with given name exist already
    if( nameTypeMap_.find( name) != nameTypeMap_.end() ) {
      EXCEPTION( "Entities with name " << name 
                 << " are already defined" );
    }
    
    if(!isInitialized_)
      EXCEPTION("Cannot add a volume region to an uninitialized grid!");

    if(regionNames_.Find(name) != -1)
    {
      EXCEPTION("A volume region with the same name already exists!");
      regionid = -1;
    }

    regionNames_.Push_back(name);
    regionid = regionNames_.GetSize()-1;
    volRegionIds_.Push_back(regionid);
    
    StdVector<Elem*> dummy_elems;
    volElems_.Push_back(dummy_elems);
    
    std::set<UInt> dummy_nodes;
    volElemNodes_.Push_back(dummy_nodes);                                      
  }

  void GridCFS::ClearRegion( const RegionIdType regionid )
  {
    StdVector<Elem*> newOrderedElems;
    StdVector<Elem*> elems;
    UInt numElems;
    UInt i, n;
    
    // look in volume regions
    UInt index = volRegionIds_.Find(regionid);
    if ( index != -1 ) {
      n = volElems_[index].GetSize();

      for(i=0; i<n; i++)
      {
        orderedElems_[volElems_[index][i]->elemNum-1] = NULL;
        delete volElems_[index][i];
      }

      volElems_[index].Clear();
    } else {    
      // look in surface regions
      index = surfRegionIds_.Find(regionid);
      if ( index != -1 ) {
        n = surfElems_[index].GetSize();

        for(i=0; i<n; i++)
        {
          orderedElems_[surfElems_[index][i]->elemNum-1] = NULL;
          delete surfElems_[index][i];
        }

        surfElems_[index].Clear();
      } else {
        EXCEPTION("GridCFS: The region with id '" << regionid
                  << "' was not found in the grid!");
      }
    }


    numElems = 0;
    for(i=0; i<numElems_; i++)
    {
      if(orderedElems_[i] != NULL)
      {
        newOrderedElems.Push_back(orderedElems_[i]);
        //        std::cout << "Clear Region: " << orderedElems_[i]->elemNum << " -> " << numElems << std::endl;
        orderedElems_[i]->elemNum = ++numElems;
      }
    }
    
    orderedElems_ = newOrderedElems;
    numElems_ = numElems;
  }


} // end namespace
