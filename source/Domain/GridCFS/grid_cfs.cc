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
#include "Domain/domain.hh"

namespace CoupledField {

  // declare class specific logging stream
  DECLARE_LOG(gridcfs)
  DEFINE_LOG(gridcfs, "grid.cfs")

  GridCFS::GridCFS() : Grid() {

    ENTER_FCN( "GridCFS::GridCFS" );

    isInitialized_ = false;
    isQuadratic_ = false;
    dim_ = 3;
    
    std::string probGeo;
    param->Get("domain")->Get( "geometryType", probGeo );
    
    if ( probGeo == "3d")
      dim_ = 3;
    else if (probGeo == "axi")
      dim_ = 2;
    else if (probGeo == "plane")
      dim_ = 2;

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

    ENTER_FCN( "GridCFS::GridCFS" );

    for ( UInt i = 0; i < numElems_; i++ ) {
      delete orderedElems_[i];
    }
    orderedElems_.Clear();
  } 

  /*
    void GridCFS::Read()
    {
    ENTER_FCN( "GridCFS::Read" );
 
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

  void  GridCFS::GetNodesForDirectivity() {
    ENTER_FCN( "GridCFS::GetNodesForDirectivity" );
    
    LOG_TRACE(gridcfs) << "Reading nodes for directivity pattern";
    Vector<Double> temp;

    // get parameter node of directivity definition
    ParamNode * dirNode = param->Get("domain")->Get("directivityNodes", false);
    if( !dirNode) return;

    // get reference coordinate system
    std::string coordName = dirNode->Get("coordSys")->Get("name")->AsString();
    CoordSystem * coordSys = domain->GetCoordSystem( coordName );

    // get new nodes parameter nodes
    StdVector<ParamNode*> nodesList = dirNode->GetList("nodes");

    // iterate over all nodeLists
    Double dist, startAng, stopAng, incAng;
    std::string name;
    UInt numAngles;
    StdVector<Double> angleList;
    for( UInt iList = 0; iList < nodesList.GetSize(); iList++ ) {

      // get name, distance, start/stop/incAngle
      ParamNode * actNode = nodesList[iList];
      actNode->Get( "name", name );
      actNode->Get( "distance", dist );
      actNode->Get( "startAngle", startAng );
      actNode->Get( "stopAngle", stopAng );
      actNode->Get( "incAngle", incAng );

      LOG_DBG(gridcfs) << "Creating nodes '" << name 
                       << "' with distance " << dist;

      // fill angleList
      numAngles = UInt( floor ( (stopAng-startAng) / incAng ) );
      if ( incAng <= (360.0-(stopAng-startAng)) 
           && (stopAng-startAng)>0.0) {
        numAngles++ ;
      }
      angleList.Resize( numAngles );
      for( UInt iAng = 0; iAng < numAngles; iAng++ ) {
        angleList[iAng] = startAng + iAng * incAng;
      }

      // create new node vector
      StdVector<UInt> nodes (numAngles);

      // iterate over all angles
      Vector<Double> locCoord(dim_), globCoord(dim_);
      Vector<Double> actNodeCoord(dim_);
      for( UInt iDeg = 0; iDeg < numAngles; iDeg++ ) {
        
        // map next local coordinate to global one
        locCoord[0] = dist;
        locCoord[1] = angleList[iDeg];
        coordSys->Local2GlobalCoord( globCoord, locCoord );
        LOG_DBG2(gridcfs) << iDeg+1 
                          << ". node with " << angleList[iDeg] 
                          << "° has coordinates " << globCoord.Serialize();

        // vectors with node indices and distance
        std::vector<Double> nodeDist( numNodes_ );

        // iterate over all nodes
        for( UInt iNode = 0; iNode < numNodes_; iNode++ ) {
        
          // calculate distance and store it in vector
          GetNodeCoordinate( actNodeCoord, iNode+1, false );          
          temp = (actNodeCoord-globCoord);
          nodeDist[iNode] = temp.NormL2();
        } // nodes
        
        // find minimum entry in the vector
        std::vector<Double>::iterator it ;
        it = min_element(nodeDist.begin(), nodeDist.end());
        LOG_DBG2(gridcfs) << "Found node " 
                          << std::distance(nodeDist.begin(), it) + 1;
        nodes[iDeg] = std::distance(nodeDist.begin(), it) + 1;
        
      } // angles
      AddNamedNodes( name, nodes );
      LOG_DBG2(gridcfs) << "Adding nodelist '" << name 
                        << "' with nodes " << nodes.Serialize();
    } // nodeLists

    LOG_TRACE(gridcfs) << "Finished reading nodes for directivity pattern";
  }

//   void  GridCFS::GetNodesForDirectivity() {
//     ENTER_FCN( "GridCFS::GetNodesForDirectivity" );


//     StdVector<std::string> keyVec, attrVec, valVec;
//     StdVector<Double> radiiVec;
//     attrVec = "","","";
//     valVec = "","","";
//     UInt j,k,l;
    
//     // get parameter node of directivity definition
//     ParamNode * dirNode = param->Get("domain")->Get("directivityNodes", false);
//     if( !dirNode) return;


//     // get "distance nodes" and "radius"
//     StdVector<ParamNode*> distNodes = dirNode->GetList("distance");
//     for( UInt i = 0; i < distNodes.GetSize(); i++ ) {
//       radiiVec.Push_back( distNodes[i]->Get("radius")->AsDouble() );
//     }

//     StdVector<Double> val;

//     //Read the coordinates of center of circle
//     StdVector<Double> center;
//     center.Resize(3);
//     center.Init();
//     ParamNode * centerNode = dirNode->Get("center");
//     centerNode->Get("x", center[0] );
//     centerNode->Get("y", center[1] );
//     if( dim_ == 3 ) {
//       centerNode->Get("z", center[2] );
//     }
     
//     //Read planes for saving nodes
//     StdVector<UInt> nodePlanes;
//     if ( dim_ == 3 )
//     {
//       nodePlanes.Resize(dim_);
//     }
//     else
//     {
//       nodePlanes.Resize(1);
//     }


//     // get parameter nodes defining planes
//     ParamNode* planeNode = dirNode->Get("planes");


//     nodePlanes.Init();
//     StdVector<Double> Begin, End;
//     Begin.Resize(3);
//     End.Resize(3);
//     Begin.Init();
//     End.Init();

//     StdVector<std::string> planeNames;
//     planeNames.Resize(3);
//     planeNames.Init();
    
//     planeNode->Get( "xyName", planeNames[0] );
//     planeNode->Get( "xyBegin", Begin[0] ); // flag for XY plane
//     planeNode->Get( "xyEnd", End[0] );
   
//     if ( dim_ == 3 ) {
      
//       planeNode->Get( "xzName", planeNames[1] );
//       planeNode->Get( "yzName", planeNames[2] );
      
//       planeNode->Get( "xzBegin", Begin[1] );
//       planeNode->Get( "xzEnd", End[1] );
//       planeNode->Get( "yzBegin", Begin[2] );
//       planeNode->Get( "yzEnd", End[2] );
//     }

//     //Read angle interval size for saving nodes
//     Double angleStep = 1.0;
//     dirNode->Get("saveIncr", angleStep, false );
   
//     StdVector<Integer> numDiv;
//     numDiv.Resize(3);
//     numDiv.Init();

//     nodePlanes[0]=1;
//     if (((End[0]-Begin[0])/angleStep)<1.0)
//     {
//       (*warning) << "Not saving directivity nodes on plane: " << planeNames[0]
//                  << " since sweep angle (End-Begin) "<<(End[0]-Begin[0])
//                  <<" is smaller than step angle " << angleStep << "";
//       Warning( __FILE__, __LINE__ );
//       End[0]=0;
//       Begin[0]=0;
//     }
//     numDiv[0]=Integer(floor((End[0]-Begin[0])/angleStep));
//     if (angleStep <= (360.0-(End[0]-Begin[0])) && (End[0]-Begin[0])>0.0)
//     {
//       numDiv[0]++;
//     }
    
            
//     if(dim_==3)
//     {
//       for (UInt actPlane=1; actPlane<=2; actPlane++)
//       {
//         nodePlanes[actPlane]=1;
//         if (((End[actPlane]-Begin[actPlane])/angleStep)<1.0)
//         {
//           (*warning) << "Not saving directivity nodes on plane: " << planeNames[actPlane]
//                      << " plane since sweep angle (End-Begin) "<<(End[actPlane]-Begin[actPlane])
//                      <<" is smaller than step angle " << angleStep << "";
//           Warning( __FILE__, __LINE__ );
//           End[actPlane]=0;
//           Begin[actPlane]=0;
//         }
//         numDiv[actPlane]=Integer(floor((End[actPlane]-Begin[actPlane])/angleStep));
//         if (angleStep <= (360.0-(End[actPlane]-Begin[actPlane])) && (End[actPlane]-Begin[actPlane])>0.0)
//         {
//           numDiv[actPlane]++;
//         }
//       }
//     }
    
//     //Here finally we start the process for searching nodes for directivity
//     //Create vector with angles to save solutions
//     Matrix<Double> angleList;
//     UInt totalNumDiv=numDiv[0]+numDiv[1]+numDiv[2];
//     UInt plane=3;

//     angleList.Resize(plane,totalNumDiv);
//     angleList.Init();
    
//     for (UInt actPlane=0; actPlane<numDiv.GetSize() ; actPlane++)
//     {
//       for (int i=0; i<numDiv[actPlane] ; i++)
//       {
//         angleList[actPlane][i] = (Begin[actPlane] + i*angleStep);
//       }
//     }
  
//     StdVector<UInt> directNodes;

//     if (dim_==3)
//       directNodes.Resize(radiiVec.GetSize()*(numDiv[0]*nodePlanes[0]+numDiv[1]*nodePlanes[1]+numDiv[2]*nodePlanes[2]));
//     else
//       directNodes.Resize(radiiVec.GetSize()*numDiv[0]);
    
//     UInt ctr=0;
//     for (UInt actPlane=0; actPlane<nodePlanes.GetSize(); actPlane++)
//     {
//       if (nodePlanes[actPlane]==1 && numDiv[actPlane]>0)
//       {
//         namedNodeNames_.Push_back( planeNames[actPlane]);
//         namedNodes_.Push_back( StdVector<UInt>() );
//         Info->PrintF("", "Saved directivity nodes on plane: %s\n", planeNames[actPlane].c_str()); 
//         Info->PrintF("", "Angle list: \nangles = [");  
//         for (int i=0; i<numDiv[actPlane] ; i++)
//         {
//           Info->PrintF( "", " %.1f;", angleList[actPlane][i] );
//         }
//         Info->PrintF( "", "]\n" );
            
//         if (actPlane==0)//XY
//         {
//           j=0;
//           k=1;
//           l=2;
//         }
//         if (dim_==3)
//         {
//           if (actPlane==1)//XZ
//           {
//             j=0;
//             k=2;
//             l=1;
//           }             
//           if (actPlane==2)//YZ
//           {
//             j=1;
//             k=2;
//             l=0;
//           }             
//         }  
//         for (UInt actRadIndex=0; actRadIndex<radiiVec.GetSize(); actRadIndex++)
//         {
//           std::vector<Double>::iterator it;
//           Matrix<Double> save_point;
//           save_point.Resize(numDiv[actPlane],dim_);
//           std::string nodename; 
//           Info->PrintF("", "Radius R= %.4f\nnodes = [", radiiVec[actRadIndex]);  

//           for (int i=0; i<numDiv[actPlane] ; i++)
//           {
//             save_point[i][j] = center[j] + cos(angleList[actPlane][i]/ 180 * PI)*radiiVec[actRadIndex];
//             save_point[i][k] = center[k] + sin(angleList[actPlane][i]/ 180 * PI)*radiiVec[actRadIndex];  
//             if (dim_==3)
//             {
//               save_point[i][l] = center[l] + 0.0; 
//             }
                    
//             //Computing distance between the savedpoint and the mesh nodes
//             std::vector<Double> dist2;
//             dist2.resize(numNodes_);
//             for (UInt nodeIndex=0; nodeIndex<numNodes_; nodeIndex++ )
//             {
//               Point p= coords_[nodeIndex];
                          
//               dist2[nodeIndex] = pow((p[j] - save_point[i][j]),2) +
//                                  pow((p[k] -  save_point[i][k]),2);
//               if (dim_==3)
//                 dist2[nodeIndex] +=  pow((p[l] -  save_point[i][l]),2);
//             }
//             it=min_element(dist2.begin(), dist2.end());
//             Integer numNodeMin = -1;
//             for (UInt nodeIndex=0; nodeIndex<numNodes_; nodeIndex++ )
//             {
//               if (*it==dist2[nodeIndex])
//               {
//                 numNodeMin=nodeIndex+1;
//               }
//             } 
//             directNodes[ctr*radiiVec.GetSize()+actRadIndex*numDiv[actPlane] + i] = numNodeMin;
                   
//             //                     nodename = "plane";
//             //                     nodename.append( planeNames[actPlane] );
//             //                     nodename.append( "_r" );
//             //                     nodename.append( GenStr(actRadIndex+1) );
//             //                     nodename.append( "_a" );
//             //                     nodename.append( GenStr(angleList[actPlane][i]) );
//             //namedNodeNames_.Push_back(nodename);

//             UInt lastIndex = namedNodes_.GetSize()-1; 
//             namedNodes_[lastIndex].Push_back(directNodes[ctr*radiiVec.GetSize()+actRadIndex*numDiv[actPlane] + i]);
//             Info->PrintF( "", " %i;", directNodes[ctr*radiiVec.GetSize()+actRadIndex*numDiv[actPlane] + i] ); 
//           }
//           Info->PrintF( "", "]\n" );
//         }
//         Info->PrintF( "", "\n" );
//         ctr+=numDiv[actPlane];
//       }
//     }
//   }

  void GridCFS::MapFaces() {
    ENTER_FCN( "GridCFS::MapFaces" );

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

      // offset for face index
      UInt offset=0;

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
    ENTER_FCN( "GridCFS::MapEdges" );

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

        //std::cerr << "EdgeIndices = \n" << locEdge << std::endl;
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

    for(e=0; e<numElems; e++)
    {
      Elem* el = orderedElems_[e];
      FEType type = el->ptElem->feType();
      numNodes = NUM_ELEM_NODES[type];

      maxNumElemNodes_ = maxNumElemNodes_ < numNodes ?
                         numNodes : maxNumElemNodes_;
            
      switch(type)
      {
      case ET_LINE2:
      case ET_LINE3:
        if(dim_ == 2)
        {
          surfRegionElems[el->regionId].Push_back(el);

          surfRegionNodes[el->regionId].insert(&el->connect[0],
                                               (&el->connect[0])+numNodes);
        }
        break;
      case ET_TRIA3:
      case ET_TRIA6:
      case ET_QUAD4:
      case ET_QUAD8:
      case ET_QUAD9:
        if(dim_ == 2)
        {
          volRegionElems[el->regionId].Push_back(el);

          volRegionNodes[el->regionId].insert(&el->connect[0],
                                              (&el->connect[0])+numNodes);
        }
        else
        {
          surfRegionElems[el->regionId].Push_back(el);

          surfRegionNodes[el->regionId].insert(&el->connect[0],
                                               (&el->connect[0])+numNodes);
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
        {
          volRegionElems[el->regionId].Push_back(el);

          volRegionNodes[el->regionId].insert(&el->connect[0],
                                              (&el->connect[0])+numNodes);
        }
      default:
        break;
      }
    }


    std::map<RegionIdType, StdVector<Elem*> >::iterator regionElemIt, regionElemEnd;
    std::map<RegionIdType, std::set<UInt> >::iterator regionNodeIt;
    std::set<UInt>::iterator setIt, setEnd;
    UInt region = 0;
    UInt regionId;

    numVolRegions = volRegionElems.size();
    numSurfRegions = surfRegionElems.size();
    volElemNodes_.Resize(numVolRegions);
    surfElemNodes_.Resize(numSurfRegions);
        
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
        
    CreateSurfaceElements(surfElems_);

    isInitialized_ = true;

    // Read in directivity node from parameter file
    GetNodesForDirectivity();

    // print information to file
    PrintGridInfo();
  }
    

  // ======================================================
  // GENERAL GRID INFORMATION
  // ======================================================

  
  UInt GridCFS::GetNumElemOfType( FEType type ) {
    ENTER_FCN( "GridCFS::GetNumElemOfType" );
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
    ENTER_FCN( "GridCFS::GetDim" );
    return dim_;
  }

  
  UInt GridCFS::GetNumNodes(){
    ENTER_FCN( "GridCFS::GetNumNodes" );
    return numNodes_;
  }

  
  UInt GridCFS::GetNumNodes( const StdVector<RegionIdType> & regions ) {
    ENTER_FCN( "GridCFS::GetNumNodes" );
    
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
    ENTER_FCN( "GridCFS::GetNumNodes" );

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
    ENTER_IFCN( "GridCFS::GetNumElems" );
    return numElems_;
  }

  
  UInt GridCFS::GetNumSurfElems() {
    ENTER_FCN( "GridCFS::GetNumSurfElems" );
    
    UInt numSurfElems = 0;
    
    for (UInt i=0; i<surfElems_.GetSize(); i++) {
      numSurfElems += surfElems_[i].GetSize();
    }
    
    return numSurfElems;
  }
  
  
  UInt GridCFS::GetNumVolElems() {
    ENTER_FCN( "GridCFS::GetNumVolElems" );
    
    UInt numVolElems = 0;
    
    for (UInt i=0; i<volElems_.GetSize(); i++) {
      numVolElems += volElems_[i].GetSize();
    }
    
    return numVolElems;
  }
  
  
  UInt GridCFS::GetNumElems( const StdVector<RegionIdType> & regions ){
    ENTER_FCN( "GridCFS::GetNumElems" );

    
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
    namedNodeNames_.Push_back(name);
    namedNodes_.Push_back(nodeNums);
  }

  void GridCFS::AddNamedNodes( std::string name, std::vector<UInt> & nodeNums)
  {
    UInt numNodes = nodeNums.size();
    StdVector<UInt> nn;

    namedNodeNames_.Push_back(name);
    nn.Resize(numNodes);
    for(UInt i=0; i<numNodes; i++)
    {
      nn[i] = nodeNums[i];
    }
        
    namedNodes_.Push_back(nn);
  }
    

  void GridCFS::AddNamedElems( std::string name, StdVector<UInt> & elemNums)
  {
    namedElemNames_.Push_back(name);
    namedElems_.Push_back(elemNums);
  }

  void GridCFS::AddNamedElems( std::string name, std::vector<UInt> & elemNums)
  {
    UInt numElems = elemNums.size();
    StdVector<UInt> en;

    namedElemNames_.Push_back(name);
    en.Resize(numElems);
    for(UInt i=0; i<numElems; i++)
    {
      en[i] = elemNums[i];
    }
        
    namedElems_.Push_back(en);
  }
    
  
  void GridCFS::GetListNodeNames( StdVector<std::string> & nodeNames) {
    ENTER_FCN( "GridCFS::GetListNodeNames" );
    nodeNames = namedNodeNames_;
  }

  
  void GridCFS::GetListElemNames( StdVector<std::string> & elemNames) {
    ENTER_FCN( "GridCFS::GetListElemNames" );
    elemNames = namedElemNames_;
  }
  
  void GridCFS::GetListNodeNames( std::vector<std::string> & nodeNames)
  {
    UInt numNamedNodes = namedNodeNames_.GetSize();
    nodeNames.resize(numNamedNodes);
        
    for(UInt i=0; i<numNamedNodes; i++)
    {
      nodeNames[i] = namedNodeNames_[i];
    }
  }
    
  void GridCFS::GetListElemNames( std::vector<std::string> & elemNames)
  {
    UInt numNamedElems = namedElemNames_.GetSize();
    elemNames.resize(numNamedElems);
        
    for(UInt i=0; i<numNamedElems; i++)
    {
      elemNames[i] = namedElemNames_[i];
    }
  }

  
  // ======================================================
  // NODE ACCESS FUNCTIONS
  // ======================================================
  
  void GridCFS::GetNodesByName( StdVector<UInt> & nodeList,
                                const std::string & name ) {
    ENTER_FCN( "GridCFS::GetNodesByName" );

    Integer index = namedNodeNames_.Find(name);
    if ( index != -1 ) {
      nodeList = namedNodes_[index];
    } else {
      EXCEPTION( "GridCFS: There are no nodes with name '" << name
                 << "' in the grid!" );
    }
    
  }

  void GridCFS::GetNodesByName( std::vector<UInt> & nodeList,
                                const std::string & name ) {
    ENTER_FCN( "GridCFS::GetNodesByName" );

    Integer index = namedNodeNames_.Find(name);
    if ( index != -1 ) {
      UInt numNodes = namedNodes_[index].GetSize();
      nodeList.resize(numNodes);
        
      for(UInt i=0; i<numNodes; i++)
        nodeList[i] = namedNodes_[index][i];
    } else {
      EXCEPTION( "GridCFS: There are no nodes with name '" << name
                 << "' in the grid!" );
    }
    
  }
  
  
  void GridCFS::GetNodesByRegion( StdVector<UInt> & nodeList,
                                  const RegionIdType regionId ) {
    ENTER_FCN( "GridCFS::GetNodesByRegion" );

    Integer index = 0;
    
    // look in volume regions
    index = volRegionIds_.Find(regionId);
    if ( index != -1 ) {
      nodeList.Resize(volElemNodes_[index].size());
      std::copy(volElemNodes_[index].begin(),
                volElemNodes_[index].end(),
                &nodeList[0]);
    } else {
      
      // look in surface regions
      index = surfRegionIds_.Find(regionId);
      if ( index != -1 ) {
        nodeList.Resize(surfElemNodes_[index].size());
        std::copy(surfElemNodes_[index].begin(),
                  surfElemNodes_[index].end(),
                  &nodeList[0]);
      } else {
        EXCEPTION( "GridCFS: The region with id '" << regionId
                   << "' was not found in the grid!" );
      }
    }
    
  }
  
  
  void GridCFS::GetNodeCoordinate( Point & rfPoint,
                                   const UInt inode, 
                                   bool updated ) {
    ENTER_FCN( "GridCFS::GetNodeCoordinate" );
    
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
    ENTER_FCN( "GridCFS::GetNodeCoordinate" );

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
    ENTER_FCN( "GridCFS::GetElem" );
    
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
    ENTER_FCN( "GridCFS::GetElems" );
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
    ENTER_FCN( "GridCFS::GetVolElems" );
    
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
    ENTER_FCN( "GridCFS::GetSurfElems" );
    
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
    ENTER_FCN( "GridCFS::GetElemsByName" );

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

  void GridCFS::GetElemsByName( std::vector<Elem*> & elems,
                                const std::string & elemsName ) {
    ENTER_FCN( "GridCFS::GetElemsByName" );

    Integer index = namedElemNames_.Find(elemsName);
    

    if ( index != -1 ) {
      UInt numElems = namedElems_[index].GetSize();
      elems.resize(numElems);
        
      for(UInt i=0; i<numElems; i++)
        elems[i] = orderedElems_[namedElems_[index][i]-1];
    } else {
      EXCEPTION( "GridCFS: There are no named elements with name '" 
                 << elemsName << "' in the grid!" );
    }
    
  }


  void GridCFS::GetElemNumsByName( std::vector<UInt> & elemNums,
                                   const std::string & elemsName )
  {
    ENTER_FCN( "GridCFS::GetElemNumsByName" );

    Integer index = namedElemNames_.Find(elemsName);
    

    if ( index != -1 ) {
      UInt numElems = namedElems_[index].GetSize();
      elemNums.resize(numElems);
        
      for(UInt i=0; i<numElems; i++)
        elemNums[i] = namedElems_[index][i];
    } else {
      EXCEPTION( "GridCFS: There are no named elements with name '" 
                 << elemsName << "' in the grid!" );
    }
    
  }
  
  void GridCFS::GetElemNodes( StdVector<UInt> & connect, 
                              const UInt iElem ) {
    ENTER_FCN( "GridCFS::GetElemNodes" );
    
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
    ENTER_FCN( "GridCFS::GetElemNodesCoord" );

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
    ENTER_FCN( "GridCFS::GetElemsNextToNodes" );
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
    ENTER_FCN( "GridCFS::GetElemsNextToSurface" );
    EXCEPTION( "Not implemented" );
  }
    



  // ======================================================
  // MISCELLANEOUS
  // ======================================================
  
  void GridCFS::GetNodesOfElemList( StdVector<UInt> & nodeList,
                                    const StdVector<Elem*> & elemList,
                                    bool onlyLinNodes) {
    ENTER_FCN( "GridCFS::GetNodesOfElemList" );

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
    ENTER_FCN( "GridCFS::GetAllRegionNames");
        
    regionNames = regionNames_;
        
  }

  void GridCFS::GetRegionNames( std::vector<std::string> 
                                & regionNames )
  {
    UInt numRegions = regionNames_.GetSize();

    regionNames.resize(numRegions);
        

    for(UInt i=0; i< numRegions; i++)
    {
      regionNames[i] = regionNames_[i];
    }
  }
    


  
  void GridCFS::SetNodeOffset( const StdVector<UInt>& nodes, 
                               const Vector<Double>& offsets ) {
    ENTER_FCN( "GridCFS::SetNodeOffset" );

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
    ENTER_FCN( "GridCFS::HasNodalOffset()" );

    if( deltCoords_.GetSize() != 0 ) {
      return true;
    } else {
      return false;
    }
     
  }

  // =======================================================================
  // Helper Methods
  // =======================================================================
  
  void GridCFS::CreateSurfaceElements( StdVector<StdVector<Elem*> > & elems) {
    ENTER_FCN( "GridCFS::CreateSurfaceElements" );
   

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
    surfElems_.Resize(elems.GetSize());
    //    surfElems_.Init();
    
    for ( UInt iRegion = 0; iRegion < elems.GetSize(); iRegion++ ) {

      surfElems_[iRegion].Resize(elems[iRegion].GetSize());
      //      surfElems_[iRegion].Init();

      for (UInt iSurfElem = 0; iSurfElem < elems[iRegion].GetSize();
           iSurfElem++ ) {

        oldElem = elems[iRegion][iSurfElem];
          

        // create new surface element
        myElem = new SurfElem();
        myElem->connect = oldElem->connect;
        myElem->regionId = oldElem->regionId;
        myElem->elemNum = oldElem->elemNum;
        myElem->ptElem = oldElem->ptElem;
        surfElems_[iRegion][iSurfElem] = myElem;
        
        // add surface element into vector with ordered Elements
        orderedElems_[myElem->elemNum-1] = myElem;

        // delete old volume element
        delete oldElem;
      }
    }

    // 3.) Iterate over all surface elements and look for each
    //     element, if all of its nodes can be assigned to one or
    //     two neighbours

    UInt surfNodeNr = 0;
    UInt elemsFound = 0;
    UInt elemsAssigned = 0;
    
    for ( iRegion = 0; iRegion < surfElems_.GetSize(); iRegion++ ) {

      
      for (UInt iSurfElem = 0; iSurfElem < surfElems_[iRegion].GetSize();
           iSurfElem++ ) {
        elemsAssigned = 0;

        myElem = dynamic_cast<SurfElem*>(surfElems_[iRegion][iSurfElem]);

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
                     << surfElems_[iRegion][iSurfElem]->elemNum );
        }

      } // loop over surface elements
    } // loop over regions




    // 4.) Iterate over all surface elements and calculate surface
    //     flag by comparing the directed and the undirected surface
    //     normal. If both differ, the surfaceNormalSign = -1, otherwise 1.
    Vector<Double> normalUndefSign, normalDefSign;
    Double sign;

    for ( UInt iRegion = 0; iRegion < surfElems_.GetSize(); iRegion++ ) {
      for (UInt iSurfElem = 0; iSurfElem < surfElems_[iRegion].GetSize();
           iSurfElem++ ) {
        
        myElem = dynamic_cast<SurfElem*>(surfElems_[iRegion][iSurfElem]);

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
    
  }

  
  void GridCFS::PrintGridInfo() const {
    ENTER_FCN( "GRIDCFS::PrintGridInfo()" );
    
    std::string help, empty;

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
    ENTER_FCN( "GridCFS::CalcSurfNormal" );
    
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
    ENTER_IFCN( "GridCFS::CalcSurfNormalOutOfVol" );

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
    ENTER_FCN( "GridCFS::GetVolumeOfRegion" );

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
  

#ifdef ADAPTGRID
  
  void GridCFS::putNodesFromGrid_RG(grd::MultilevelGrid * grid,
                                    const UInt level)
  {
    ENTER_FCN( "GridCFS::putNodesFromGrid_RG" );

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
    ENTER_FCN( "GridCFS::putElemsFromGrid_RG" );

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
    ENTER_FCN( "GridCFS::putElemsNodesFromGrid_RG" );

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
    ENTER_FCN( "GridCFS::Refine" );

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
    ENTER_FCN( "GridCFS<Dim>::ReRefine" );
    
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
