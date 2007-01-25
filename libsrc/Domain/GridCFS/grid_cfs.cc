#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <set>

#include "grid_cfs.hh"

#include "Domain/grid.hh"
#include "Elements/elements_header.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField {

// declare class specific logging stream
  DECLARE_LOG(gridcfs)
  DEFINE_LOG(gridcfs, "grid.cfs")

  template<UInt DIM>
  GridCFS<DIM>::GridCFS( FileType * const aptFileType) {

    ENTER_FCN( "GridCFS::GridCFS" );

    inFile_ = aptFileType;

    isInitialized_ = false;
    isQuadratic_ = false;
    dim_ = 0;
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
  template<UInt DIM>
  GridCFS<DIM>::~GridCFS() {

    ENTER_FCN( "GridCFS::GridCFS" );

    // Delete volume elements
    for ( UInt i = 0; i < volElems_.GetSize(); i++ ) {
      for ( UInt j = 0; j < volElems_[i].GetSize(); j++ ) {
        delete volElems_[i][j];
      }
    }
    volElems_.Clear();

    // Delete surface elements
    for ( UInt i = 0; i < surfElems_.GetSize(); i++ ) {
      for ( UInt j = 0; j < surfElems_[i].GetSize(); j++ ) {
        delete surfElems_[i][j];
      }
    }
    surfElems_.Clear();

  } 


  template<UInt DIM>
  void GridCFS<DIM>::Read()
  {
    ENTER_FCN( "GridCFS::Read" );
 
    StdVector<StdVector<Elem*> >temp;

    // 1. Read Dimension
    dim_=inFile_->GetDim();
    numNodes_ = inFile_->GetNumNodes();
    numElems_ = inFile_->GetNumElems();

    // 2. Read coordinatess
    inFile_->GetCoordinates(coords_);

    // 3. Read in volume elements
    inFile_->GetElements(volElems_, volRegionIds_, DIM);
     
    // 4. Create sorted list of volume elements
    // 5. If there are surfac
    //    -  Iterate over all volume elements and make list
    //       of neighbouring elements for each node
    //    -   Read In Surface elements
    inFile_->GetElements(temp, surfRegionIds_, DIM-1);
    
    // Obtain all region names
    inFile_->GetAllRegionNames(regionNames_);
    
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
          (*error) << "Element Nr. " << ptVolElem->elemNum 
                   << " exists at least two times!\n"
                   << "The first occurence is in region '"
                   << regionNames_[orderedElems_[ptVolElem->elemNum-1]
                                     ->regionId] 
                   << "', the second in region '"
                   << regionNames_[volElems_[iRegion][iElem]->regionId]
                   << "'.\nPlease check your mesh file!";
          Error( __FILE__, __LINE__ );
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
    inFile_->GetNodesOfRegions( volElemNodes_, volRegionIds_ );
    inFile_->GetNodesOfRegions( surfElemNodes_, surfRegionIds_ );
 
    
    // 8. Read in named nodes
    inFile_->GetNamedNodes( namedNodes_, namedNodeNames_ );
    
    // 9. Read in named elements
    inFile_->GetNamedElems( namedElems_, namedElemNames_ );

    // 10. Get nodes for directivity plots if defined in XML file
    GetNodesForDirectivity() ;


    // Print information about region mapping into
    PrintGridInfo();
    
#ifdef ADAPTGRID
    FormNeighborsLists();
#endif

    isInitialized_ = true;
  }

  template<UInt DIM>
  void  GridCFS<DIM>::GetNodesForDirectivity() {
    ENTER_FCN( "GridCFS::GetNodesForDirectivity" );


    StdVector<std::string> keyVec, attrVec, valVec;
    StdVector<Double> radiiVec;
    attrVec = "","","";
    valVec = "","","";
    UInt j,k,l;

    keyVec = "domain", "directivityNodes", "distance", "radius";
    params->GetList(keyVec, attrVec, valVec, radiiVec);
    // Check if there are any radii for saving nodes,
    // otherwise, jump out
    if ( radiiVec.GetSize() == 0 ) 
      {
        return;
      }

    StdVector<Double> val;

    //Read the coordinates of center of circle
    StdVector<Double> center;
    center.Resize(3);
    center.Init();
    //x
    keyVec = "domain" , "directivityNodes" , "center" , "x";
    params->GetList( keyVec, attrVec, valVec, val);
    center[0] =  val[0];
    
    //y
    keyVec = "domain" , "directivityNodes" , "center" , "y";
    params->GetList( keyVec, attrVec, valVec, val);
    center[1] =  val[0];
    
    if ( dim_ == 3 ) {
    //z
    keyVec = "domain" , "directivityNodes" , "center" , "z";
    params->GetList( keyVec, attrVec, valVec, val); 
    center[2] =  val[0];
    }
    else{
      center[2] = 0;
    }
     
    //Read planes for saving nodes
    Double planeFlag;
    
    StdVector<UInt> nodePlanes;
    if ( dim_ == 3 )
      {
        nodePlanes.Resize(dim_);
      }
    else
      {
        nodePlanes.Resize(1);
      }
    nodePlanes.Init();
    
    StdVector<Double> Begin, End;
    Begin.Resize(3);
    End.Resize(3);
    Begin.Init();
    End.Init();

    StdVector<std::string> planeNames;
    planeNames.Resize(3);
    planeNames.Init();
    
    keyVec = "domain" , "directivityNodes" , "planes" , "xyName";
    params->Get( keyVec, attrVec, valVec, planeNames[0]);
    
    keyVec = "domain" , "directivityNodes" , "planes" , "xyBegin";
    params->Get( keyVec, attrVec, valVec, planeFlag);
    Begin[0]=planeFlag; // flag for XY plane

    keyVec = "domain" , "directivityNodes" , "planes" , "xyEnd";
    params->Get( keyVec, attrVec, valVec, planeFlag);
    End[0]=planeFlag; // flag for XY plane
   
    if ( dim_ == 3 ) {
      
      keyVec = "domain" , "directivityNodes" , "planes" , "xzName";
      params->Get( keyVec, attrVec, valVec, planeNames[1]);

      keyVec = "domain" , "directivityNodes" , "planes" , "yzName";
      params->Get( keyVec, attrVec, valVec, planeNames[2]);
    
      keyVec = "domain" , "directivityNodes" , "planes" , "xzBegin";    
      params->Get( keyVec, attrVec, valVec, planeFlag);
      Begin[1]=planeFlag; // flag for XZ plane
      
      keyVec = "domain" , "directivityNodes" , "planes" , "xzEnd";    
      params->Get( keyVec, attrVec, valVec, planeFlag);
      End[1]=planeFlag; // flag for XZ plane

      keyVec = "domain" , "directivityNodes" , "planes" , "yzBegin";    
      params->Get( keyVec, attrVec, valVec, planeFlag);
      Begin[2]=planeFlag; // flag for YZ plane

      keyVec = "domain" , "directivityNodes" , "planes" , "yzEnd";    
      params->Get( keyVec, attrVec, valVec, planeFlag);
      End[2]=planeFlag; // flag for YZ plane
    }

    //Read angle interval size for saving nodes
    Double angleStep;
    params->Get( "saveIncr", angleStep, "directivityNodes");
   
    StdVector<Integer> numDiv;
    numDiv.Resize(3);
    numDiv.Init();

    nodePlanes[0]=1;
    if (((End[0]-Begin[0])/angleStep)<1.0)
      {
        (*warning) << "Not saving directivity nodes on plane: " << planeNames[0]
                   << " since sweep angle (End-Begin) "<<(End[0]-Begin[0])
                   <<" is smaller than step angle " << angleStep << "";
        Warning( __FILE__, __LINE__ );
        End[0]=0;
        Begin[0]=0;
      }
    numDiv[0]=Integer(floor((End[0]-Begin[0])/angleStep));
    if (angleStep <= (360.0-(End[0]-Begin[0])) && (End[0]-Begin[0])>0.0)
      {
        numDiv[0]++;
      }
    
            
    if(dim_==3)
      {
        for (UInt actPlane=1; actPlane<=2; actPlane++)
          {
            nodePlanes[actPlane]=1;
            if (((End[actPlane]-Begin[actPlane])/angleStep)<1.0)
              {
                (*warning) << "Not saving directivity nodes on plane: " << planeNames[actPlane]
                           << " plane since sweep angle (End-Begin) "<<(End[actPlane]-Begin[actPlane])
                           <<" is smaller than step angle " << angleStep << "";
                Warning( __FILE__, __LINE__ );
                End[actPlane]=0;
                Begin[actPlane]=0;
              }
            numDiv[actPlane]=Integer(floor((End[actPlane]-Begin[actPlane])/angleStep));
            if (angleStep <= (360.0-(End[actPlane]-Begin[actPlane])) && (End[actPlane]-Begin[actPlane])>0.0)
              {
                numDiv[actPlane]++;
              }
          }
      }
    
    //Here finally we start the process for searching nodes for directivity
    //Create vector with angles to save solutions
    Matrix<Double> angleList;
    UInt totalNumDiv=numDiv[0]+numDiv[1]+numDiv[2];
    UInt plane=3;

    angleList.Resize(plane,totalNumDiv);
    angleList.Init();
    
    for (UInt actPlane=0; actPlane<numDiv.GetSize() ; actPlane++)
      {
        for (int i=0; i<numDiv[actPlane] ; i++)
          {
            angleList[actPlane][i] = (Begin[actPlane] + i*angleStep);
          }
      }
  
    StdVector<UInt> directNodes;

    if (dim_==3)
      directNodes.Resize(radiiVec.GetSize()*(numDiv[0]*nodePlanes[0]+numDiv[1]*nodePlanes[1]+numDiv[2]*nodePlanes[2]));
    else
      directNodes.Resize(radiiVec.GetSize()*numDiv[0]);
    
    UInt ctr=0;
    for (UInt actPlane=0; actPlane<nodePlanes.GetSize(); actPlane++)
      {
        if (nodePlanes[actPlane]==1 && numDiv[actPlane]>0)
          {
            namedNodeNames_.Push_back( planeNames[actPlane]);
            namedNodes_.Push_back( StdVector<UInt>() );
            Info->PrintF("", "Saved directivity nodes on plane: %s\n", planeNames[actPlane].c_str()); 
            Info->PrintF("", "Angle list: \nangles = [");  
            for (int i=0; i<numDiv[actPlane] ; i++)
              {
                Info->PrintF( "", " %.1f;", angleList[actPlane][i] );
              }
            Info->PrintF( "", "]\n" );
            
            if (actPlane==0)//XY
              {
                j=0;
                k=1;
                l=2;
              }
            if (dim_==3)
              {
                if (actPlane==1)//XZ
                  {
                    j=0;
                    k=2;
                    l=1;
                  }             
                if (actPlane==2)//YZ
                  {
                    j=1;
                    k=2;
                    l=0;
                  }             
              }  
            for (UInt actRadIndex=0; actRadIndex<radiiVec.GetSize(); actRadIndex++)
              {
                std::vector<Double>::iterator it;
                Matrix<Double> save_point;
                save_point.Resize(numDiv[actPlane],dim_);
                std::string nodename; 
                Info->PrintF("", "Radius R= %.4f\nnodes = [", radiiVec[actRadIndex]);  

                for (int i=0; i<numDiv[actPlane] ; i++)
                  {
                    save_point[i][j] = center[j] + cos(angleList[actPlane][i]/ 180 * PI)*radiiVec[actRadIndex];
                    save_point[i][k] = center[k] + sin(angleList[actPlane][i]/ 180 * PI)*radiiVec[actRadIndex];  
                    if (dim_==3)
                      {
                        save_point[i][l] = center[l] + 0.0; 
                      }
                    
                    //Computing distance between the savedpoint and the mesh nodes
                    std::vector<Double> dist2;
                    dist2.resize(numNodes_);
                    for (UInt nodeIndex=0; nodeIndex<numNodes_; nodeIndex++ )
                      {
                        dist2[nodeIndex] = pow((coords_[nodeIndex][j] - save_point[i][j]),2) +
                          pow((coords_[nodeIndex][k] -  save_point[i][k]),2);
                        if (dim_==3)
                          dist2[nodeIndex] +=  pow((coords_[nodeIndex][l] -  save_point[i][l]),2);
                        //std::cout<<"dist2: "<<dist2[nodeIndex]<<std::endl;
                      }
                    it=min_element(dist2.begin(), dist2.end());
                    Integer numNodeMin = -1;
                    for (UInt nodeIndex=0; nodeIndex<numNodes_; nodeIndex++ )
                      {
                        if (*it==dist2[nodeIndex])
                          {
                            numNodeMin=nodeIndex+1;
                          }
                      } 
                    directNodes[ctr*radiiVec.GetSize()+actRadIndex*numDiv[actPlane] + i] = numNodeMin;
                   
//                     nodename = "plane";
//                     nodename.append( planeNames[actPlane] );
//                     nodename.append( "_r" );
//                     nodename.append( GenStr(actRadIndex+1) );
//                     nodename.append( "_a" );
//                     nodename.append( GenStr(angleList[actPlane][i]) );
                    //namedNodeNames_.Push_back(nodename);

                    UInt lastIndex = namedNodes_.GetSize()-1; 
                    namedNodes_[lastIndex].Push_back(directNodes[ctr*radiiVec.GetSize()+actRadIndex*numDiv[actPlane] + i]);
                    Info->PrintF( "", " %i;", directNodes[ctr*radiiVec.GetSize()+actRadIndex*numDiv[actPlane] + i] ); 
                  }
                Info->PrintF( "", "]\n" );
              }
                Info->PrintF( "", "\n" );
            ctr+=numDiv[actPlane];
          }
      }
  }

  template<UInt DIM>
  void GridCFS<DIM>::MapFaces() {
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
     if ( actElem.ptElem->GetDim() < DIM ) { continue; }

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
  



  template<UInt DIM>
  void GridCFS<DIM>::MapEdges() {
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

  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumEdges() {
    return numEdges_;
    
  }

  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumFaces() {
    return numFaces_;
    
  }


  template<UInt DIM>
  const Edge&  GridCFS<DIM>::GetEdge( UInt edgeNr ) {
    if( !edgesMapped_ ) {
      Error( "Edges are not mapped yet!",
             __FILE__, __LINE__ );
    }
  
    Edge const & ret = edges_[edgeNr-1];
    
    return ret;
  }
  
  template<UInt DIM>
  const Face&  GridCFS<DIM>::GetFace( UInt faceNr ) {
    if( !facesMapped_ ) {
      Error( "Surfaces are not mapped yet!",
             __FILE__, __LINE__ );
    }
  
    Face const & ret = faces_[faceNr-1];
    
    return ret;
  }

  // ======================================================
  // GENERAL GRID INFORMATION
  // ======================================================

  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumElemOfType( FEType type ) {
    ENTER_FCN( "GridCFS::GetNumElemOfType" );
    return numElemTypes_[type];
  }

  template<UInt DIM>
  UInt GridCFS<DIM>::GetDim() {
    ENTER_FCN( "GridCFS::GetDim" );
    return dim_;
  }

  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumNodes(){
    ENTER_FCN( "GridCFS::GetNumNodes" );
    return numNodes_;
  }

  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumNodes( const StdVector<RegionIdType> & regions ) {
    ENTER_FCN( "GridCFS::GetNumNodes" );
    
    UInt numNodes = 0;
    Integer index = 0;

    for ( UInt i=0; i<regions.GetSize(); i++ ) {
      
      // look in volume regions
      index = volRegionIds_.Find(regions[i]);
      if ( index != -1 ) {
        numNodes += volElemNodes_[index].GetSize();
      } else {
        
        // look in surface regions
        index = surfRegionIds_.Find(regions[i]);
        if ( index != -1 ) {
          numNodes += surfElemNodes_[index].GetSize();
        } else {
          (*error) << "GridCFS: The region with id '" << regions[i]
                   << "' was not found in the grid!";
          Error( __FILE__, __LINE__ );
        }
      }
    }
    
    return numNodes;
    
  }

  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumNodes( const std::string & nodesName ) {
    ENTER_FCN( "GridCFS::GetNumNodes" );

    UInt numNodes = 0;
    Integer index = namedNodeNames_.Find(nodesName);

    if ( index != -1 ) {
      numNodes = namedNodes_[index].GetSize();
    } else {
      (*error) << "GridCFS: The Nodes with name '" << nodesName
               << "' were not found in the grid!";
      Error( __FILE__, __LINE__ );
    }

    return numNodes;
  }
  
  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumElems() {
    ENTER_IFCN( "GridCFS::GetNumElems" );
    return numElems_;
  }

  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumSurfElems() {
    ENTER_FCN( "GridCFS::GetNumSurfElems" );
    
    UInt numSurfElems = 0;
    
    for (UInt i=0; i<surfElems_.GetSize(); i++) {
      numSurfElems += surfElems_[i].GetSize();
    }
    
    return numSurfElems;
  }
  
  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumVolElems() {
    ENTER_FCN( "GridCFS::GetNumVolElems" );
    
    UInt numVolElems = 0;
    
    for (UInt i=0; i<volElems_.GetSize(); i++) {
      numVolElems += volElems_[i].GetSize();
    }
    
    return numVolElems;
  }
  
  template<UInt DIM>
  UInt GridCFS<DIM>::GetNumElems( const StdVector<RegionIdType> & regions ){
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
            (*error) << "GridCFS: The region with id '" << regions[i]
                     << "' was not found in the grid!";
            Error( __FILE__, __LINE__ );
          }
        }
      }
    }
    return numElems;
    
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetRegionIds( StdVector<RegionIdType> & regions ) {
    ENTER_FCN( "GridCFS::GetRegionIds" );
    
    regions.Clear();
    regions.Reserve(volRegionIds_.GetSize() + surfRegionIds_.GetSize());
    
    for( UInt i = 0; i < volRegionIds_.GetSize(); i++) {
      regions.Push_back(volRegionIds_[i]);
    }

    for( UInt i = 0; i < surfRegionIds_.GetSize(); i++) {
      regions.Push_back(surfRegionIds_[i]);
    }
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetVolRegionIds( StdVector<RegionIdType> & volRegions ) {
    ENTER_FCN( "GridCFS::GetVolRegionIds" );
    
    volRegions = volRegionIds_;
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetSurfRegionIds( StdVector<RegionIdType> & surfRegions ) {
    ENTER_FCN( "GridCFS::GetSurfRegionIds" );
    
    surfRegions = surfRegionIds_;
  }

  
  template<UInt DIM>
  void GridCFS<DIM>::GetListNodeNames( StdVector<std::string> & nodeNames) {
    ENTER_FCN( "GridCFS<DIM>::GetListNodeNames" );
    nodeNames = namedNodeNames_;
  }

  template<UInt DIM>
  void GridCFS<DIM>::GetListElemNames( StdVector<std::string> & elemNames) {
    ENTER_FCN( "GridCFS<DIM>::GetListElemNames" );
    elemNames = namedElemNames_;
  }
  
  
  // ======================================================
  // NODE ACCESS FUNCTIONS
  // ======================================================
  template<UInt DIM>
  void GridCFS<DIM>::GetNodesByName( StdVector<UInt> & nodeList,
                                     const std::string & name ) {
    ENTER_FCN( "GridCFS::GetNodesByName" );

    Integer index = namedNodeNames_.Find(name);
    if ( index != -1 ) {
      nodeList = namedNodes_[index];
    } else {
      (*error) << "GridCFS: There are no nodes with name'" << name
               << "' in the grid!";
      Error( __FILE__, __LINE__ );
    }
    
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetNodesByRegion( StdVector<UInt> & nodeList,
                                       const RegionIdType regionId ) {
    ENTER_FCN( "GridCFS::GetNodesByRegion" );

    Integer index = 0;
    
    // look in volume regions
    index = volRegionIds_.Find(regionId);
    if ( index != -1 ) {
      nodeList = volElemNodes_[index];
    } else {
      
      // look in surface regions
      index = surfRegionIds_.Find(regionId);
      if ( index != -1 ) {
        nodeList  = surfElemNodes_[index];
      } else {
        (*error) << "GridCFS: The region with id '" << regionId
                 << "' was not found in the grid!";
        Error( __FILE__, __LINE__ );
      }
    }
    
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetNodeCoordinate( Point<DIM> & rfPoint,
                                        const UInt inode, 
                                        bool updated ) {
    ENTER_FCN( "GridCFS::GetNodeCoordinate" );
    
    if ( inode > numNodes_ || inode < 0 ) {
      (*error) << "GridCFS: There are only " << numNodes_
               << " nodes in the grid. You requested coordinates for "
               << "node number " << inode <<". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }

    rfPoint = coords_[inode-1];
    
    if( updated == true && deltCoords_.GetSize() != 0 ) {
      rfPoint += deltCoords_[inode-1];
    }


  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetNodeCoordinate( Vector<Double> & rfPoint,
                                        const UInt inode, 
                                        bool updated ) {
    ENTER_FCN( "GridCFS::GetNodeCoordinate" );

    if ( inode > numNodes_ || inode < 0 ) {
      (*error) << "GridCFS: There are only " << numNodes_
               << " nodes in the grid. You requested coordinates for "
               << "node number " << inode <<". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    rfPoint = coords_[inode-1];
    if( updated == true && deltCoords_.GetSize() != 0 ) {
      for( UInt iDim = 0; iDim < DIM; iDim++ ) {
        rfPoint[iDim] += deltCoords_[inode-1][iDim];
      }
    }
  }
    
  // ======================================================
  // ELEMENT ACCESS FUNCTIONS
  // ======================================================
  template<UInt DIM>
  const Elem * GridCFS<DIM>::GetElem( UInt elemNr ) {
    ENTER_FCN( "GridCFS::GetElem" );
    
#ifdef DEBUG
    if ( elemNr > numElems_ || elemNr < 0) {  
      (*error) << "GridCFS: There are only " << numElems_ 
               << " elements in the grid! You requested element number " 
               << elemNr << ". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    if ( orderedElems_[elemNr-1] == NULL ) {
      (*error) << "Element with Nr. " << elemNr << " is not contained in mesh!";
      Error( __FILE__, __LINE__ );
    }
#endif
   
    return orderedElems_[elemNr-1];

  }


  template<UInt DIM>
  void GridCFS<DIM>::GetElems( StdVector<Elem*> & elems, 
                                  const RegionIdType regionId ) {
    ENTER_FCN( "GridCFS::GetElems" );
    elems.Clear();
    
    // check if region Id is ALL_REGIONS
    if ( regionId == ALL_REGIONS ) {
      elems.Reserve( GetNumVolElems() + GetNumSurfElems() );

      // iterate over all volume elements
      for ( UInt i = 0; i < volElems_.GetSize(); i++) {
        for (UInt iElem = 0; iElem < volElems_[i].GetSize(); iElem++ ) {
          elems.Push_back(volElems_[i][iElem]);
        }
      }
      // iterate over all surface elements
      for ( UInt i = 0; i < surfElems_.GetSize(); i++) {
        for (UInt iElem = 0; iElem < surfElems_[i].GetSize(); iElem++ ) {
          elems.Push_back(surfElems_[i][iElem]);
        }
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
          (*error) << "GridCFS: The region with id '" << regionId
                   << "' was not found in the grid!";
          Error( __FILE__, __LINE__ );
        }
      }
    }
  }
  

  template<UInt DIM>
  void GridCFS<DIM>::GetVolElems( StdVector<Elem*> & elems, 
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
        (*error) << "GridCFS: The volume region with id '" << regionId
                 << "' was not found in the grid!";
        Error( __FILE__, __LINE__ );
      }
    }
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetSurfElems( StdVector<SurfElem*> & elems, 
                                   const RegionIdType regionId ) {
    ENTER_FCN( "GridCFS::GetSurfElems" );
    
    Integer index = surfRegionIds_.Find(regionId);
    if ( index != -1 ) {
      elems = surfElems_[index];
    } else {    
      (*error) << "GridCFS: The surface region with id '" << regionId
               << "' was not found in the grid!";
      Error( __FILE__, __LINE__ );
    }
  }

  template<UInt DIM>
  void GridCFS<DIM>::GetElemsByName( StdVector<Elem*> & elems,
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
      (*error) << "GridCFS: There are no named elements with name '" 
               << elemsName << "' in the grid!";
      Error( __FILE__, __LINE__ );
    }
    
  }

  template<UInt DIM>
  void GridCFS<DIM>::GetElemNodes( StdVector<UInt> & connect, 
                                   const UInt iElem ) {
    ENTER_FCN( "GridCFS::GetElemNodes" );
    
    if ( iElem > numElems_ || iElem < 0) {  
      (*error) << "GridCFS: There are only " << numElems_ 
               << " elements in the grid! You requested element number " 
               << iElem << ". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    
#ifdef DEBUG
    if ( orderedElems_[iElem-1] == NULL ) {
      (*error) << "Element with Nr. " << iElem << " is not contained in mesh!";
      Error( __FILE__, __LINE__ );
    }
#endif
    
    connect = orderedElems_[iElem-1]->connect;
    
  }

  template<UInt DIM>
  void GridCFS<DIM>::GetElemNodesCoord( Matrix<Double> & coordMat,  
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
        for (UInt actDim=0; actDim < dim_; actDim++)
          coordMat[actDim][k] = coords_[connect[k]-1][actDim];
    }


  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetElemsNextToNodes( StdVector<Elem*> & elemList, 
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
          (*error) << "GetElemsNextToNodes: A region with id '" 
                   << regionIds[isd] << "' was not found in the list of "
                   << "of volume regions.";
          Error( __FILE__, __LINE__ );
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
  
  template<UInt DIM>
  void GridCFS<DIM>::GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                            const StdVector<Elem*> & surfElems,
                                            const StdVector<RegionIdType> 
                                            &neighRegions ) {
    ENTER_FCN( "GridCFS::GetElemsNextToSurface" );
    Error( "Not implemented", __FILE__, __LINE__ );
  }
    



  // ======================================================
  // MISCELLANEOUS
  // ======================================================
  template<UInt DIM>
  void GridCFS<DIM>::GetNodesOfElemList( StdVector<UInt> & nodeList,
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
  
  
  template<UInt DIM>
  void GridCFS<DIM>::GetAllRegionNames( StdVector<std::string> & regionNames ) {
    ENTER_FCN( "GridCFS::GetAllRegionNames");
    
    regionNames = regionNames_;
    
  }

  template<UInt DIM>
  void GridCFS<DIM>::SetNodeOffset( const StdVector<UInt>& nodes, 
                                    const Vector<Double>& offsets ) {
    ENTER_FCN( "GridCFS::SetNodeOffset" );

    // Check if node offsets were already set
    if( deltCoords_.GetSize() == 0 ) {
      deltCoords_.Resize( coords_.GetSize() );
      deltCoords_.Init();
    }

    // Set delta coordinates
    for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
      Point<DIM> actOffset;
      for( UInt iDim = 0; iDim < DIM; iDim++ ) {
        actOffset[iDim] = offsets[iNode*DIM + iDim];
      }
      deltCoords_[nodes[iNode]-1] = actOffset;
    }
  }

  template<UInt DIM>
  bool GridCFS<DIM>::HasNodalOffset() {
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
  template<UInt DIM>
  void GridCFS<DIM>::CreateSurfaceElements( StdVector<StdVector<Elem*> > & elems) {
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
    surfElems_.Resize(elems.GetSize());
    surfElems_.Init();
    
    for ( UInt iRegion = 0; iRegion < elems.GetSize(); iRegion++ ) {

      surfElems_[iRegion].Resize(elems[iRegion].GetSize());
      surfElems_[iRegion].Init();

      for (UInt iSurfElem = 0; iSurfElem < elems[iRegion].GetSize();
           iSurfElem++ ) {

        // create new surface element
        myElem = new SurfElem();
        myElem->connect = elems[iRegion][iSurfElem]->connect;
        myElem->regionId = elems[iRegion][iSurfElem]->regionId;
        myElem->elemNum = elems[iRegion][iSurfElem]->elemNum;
        myElem->ptElem = elems[iRegion][iSurfElem]->ptElem;
        surfElems_[iRegion][iSurfElem] = myElem;
        
        // add surface element into vector with ordered Elements
        orderedElems_[myElem->elemNum-1] = myElem;

        // delete old volume element
        delete elems[iRegion][iSurfElem];
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

        // get number of nodes of surface element
        nrNodes = surfElems_[iRegion][iSurfElem]->connect.GetSize();
        StdVector<UInt> const & connect = 
          surfElems_[iRegion][iSurfElem]->connect;

        // get first node of surface element
        surfNodeNr = surfElems_[iRegion][iSurfElem]->connect[0];
        
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
              surfElems_[iRegion][iSurfElem]->ptVolElem1 = ptVolElem;
            }
            else {
              surfElems_[iRegion][iSurfElem]->ptVolElem2 = ptVolElem;
            }

            elemsAssigned++;
          }
        } // loop over element numbers of first node

        // sanity check (avoid the impossible ;-)
        if ( elemsAssigned > 2 ) {
          (*error) << "Found " << elemsAssigned
                   << " volume elements for surface element no. "
                   << surfElems_[iRegion][iSurfElem]->elemNum;
            Error( __FILE__, __LINE__ );
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
        

        // check, if each surface element has at least one volume neighbour
        if ( surfElems_[iRegion][iSurfElem]->ptVolElem1 == NULL ) {
         //  (*error) << "Pointer to first volume element is NULL for surface"
//                    << " element no. "  
//                    << surfElems_[iRegion][iSurfElem]->elemNum << ".\n"
//                    << "Please check your mesh-file!";
//           Error( __FILE__, __LINE__ );
//         }
          surfElems_[iRegion][iSurfElem]->normalSign = 0;
        } else {
          
          CalcSurfNormal( normalUndefSign, *surfElems_[iRegion][iSurfElem], false );
          CalcSurfNormalOutOfVol( normalDefSign, 
                                  *surfElems_[iRegion][iSurfElem],
                                  *surfElems_[iRegion][iSurfElem]->ptVolElem1,
                                  false );
          
          
          // Check if all entries have the same sign by calulating
          // a scalar product between both vectors.
          // If it is positive, they point in the smae direction,
          // otherwise an angle of 180 lies in between.
          sign = normalUndefSign * normalDefSign;
          
          if ( sign > 0.0 ) {
            surfElems_[iRegion][iSurfElem]->normalSign = 1;
          } else {
            surfElems_[iRegion][iSurfElem]->normalSign = -1;
          }
        }
      }
    }
    
  }

  template<UInt DIM>
  void GridCFS<DIM>::PrintGridInfo() const {
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
  
  template<UInt DIM>
  void GridCFS<DIM>::CalcSurfNormal( Vector<Double> & n, 
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
        Error("length of normal vector is zero!", __FILE__, __LINE__ );
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

  template<UInt DIM>
  void GridCFS<DIM>::CalcSurfNormalOutOfVol(Vector<Double> & n,
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
        Error("Nodes of neighbouring element not found!", __FILE__, __LINE__);
    
    
      // counterclockwise orientation of nodes (difference of node indizes is +1)
      if ( ( indexNode2-indexNode1  == -1 || 
            (indexNode2-indexNode1)- (Integer) volCorners == -1 ) ) {
        n *= -1;
      }
    
      else
        // counterclockwise orientation of nodes (difference of node indizes is +1)

        if (! (indexNode2-indexNode1 == 1 || 
               (indexNode2-indexNode1)+volCorners == 1) )
          Error("Nodes of interface don't lie beneath each other in neighbouring element!", __FILE__, __LINE__);   
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
  template<UInt DIM>
  Double GridCFS<DIM>::CalcVolumeOfRegion( const RegionIdType regionId,
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
  template<UInt DIM>
  void GridCFS<DIM>::putNodesFromGrid_RG(grd::MultilevelGrid * grid,
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
                      Error("Unknown type of element in GridRG", __FILE__,
                            __LINE__);
                      break;
                    }
             
                  (*el).connect.Resize(nnodes);
             
                  for (i=0; i<nnodes; i++)
                    {   
                      (*el).connect[i]=((*p)->getVertex(i))->getId();    
                    }
             
                  Integer sd=(*p)->getValue();
                  if (sd >= sd_.GetSize()) {
                    Error(" Value in element from Grid_RG is incorrect",
                          __FILE__,__LINE__);
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
                    Error(" Value in element from Grid_RG is incorrect",
                          __FILE__,__LINE__);
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
                    Error(" Value in element from Grid_RG is incorrect",
                          __FILE__,__LINE__);
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
                          Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);

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

                            if (sd >= sd_.GetSize()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);
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
                        Error("Unknown type of element in GridRG", __FILE__,__LINE__);
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
                      if (sd >= sd_.GetSize()) Error(" Value in element from Grid_RG is incorrect",__FILE__,__LINE__);
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

  template<Integer DIM>
  void GridCFS<DIM>::Refine(grd::MultilevelGrid& grid)
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
                      Error("Wrong number of subdomain",__FILE__,__LINE__);
                     
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

  template<Integer DIM>
  void GridCFS<DIM>::ReRefine(grd::MultilevelGrid& grid)
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

  template<Integer DIM>
  void GridCFS<DIM>::RefineUniform(grd::MultilevelGrid& grid)
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

  template<Integer DIM>
  void GridCFS<DIM>::SetRefinementFlag()
  {
    Integer i,j;
    for (i=0; i<sd_.GetSize(); i++) {
      for (j=0; j<elems_[i].GetSize(); j++) {
        elems_[i][j]->refinementFlag=true;
      }
    }
  }

#endif // end of ADAPTGRID



  // ecplicit template instantiation for gcc
#ifdef __GNUC__
  template class GridCFS<3>;
  template class GridCFS<2>;
#endif

  // explicit template instantitation for sgi
#ifdef __sgi
#pragma instantiate GridCFS<3>
#pragma instantiate GridCFS<2>
#endif
} // end namespace


