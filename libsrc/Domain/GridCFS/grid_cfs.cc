#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <set>

#include "grid_cfs.hh"

#include <Domain/grid.hh>
#include <Elements/elements_header.hh>
#include <DataInOut/ParamHandling/BaseParamHandler.hh>

namespace CoupledField
{

  template<UInt DIM>
  GridCFS<DIM>::GridCFS(FileType * const aptFileType)
  {
    ENTER_FCN( "GridCFS::GridCFS" );

    inFile_ = aptFileType;

    isInitialized_ = FALSE;
    dim_ = 0;
    numNodes_ = 0;
    numElems_ = 0;
    

  } 

  template<UInt DIM>
  GridCFS<DIM>::~GridCFS()
  {
    ENTER_FCN( "GridCFS::GridCFS" );
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
    //inFile_->GetNamedElems( namedElems_, namedElemNames_ );

    
    // Print information about region mapping into
    PrintGridInfo();
    
#ifdef ADAPTGRID
    FormNeighborsLists();
#endif

    isInitialized_ = TRUE;
  }


  // ======================================================
  // GENERAL GRID INFORMATION
  // ======================================================
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
                                        const UInt inode ) {
    ENTER_FCN( "GridCFS::GetNodeCoordinate" );
    
    if ( inode > numNodes_ || inode < 0 ) {
      (*error) << "GridCFS: There are only " << numNodes_
               << " nodes in the grid. You requested coordinates for "
               << "node number " << inode <<". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    
    rfPoint = coords_[inode-1];
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetNodeCoordinate( Vector<Double> & rfPoint,
                                        const UInt inode ) {
    ENTER_FCN( "GridCFS::GetNodeCoordinate" );

    if ( inode > numNodes_ || inode < 0 ) {
      (*error) << "GridCFS: There are only " << numNodes_
               << " nodes in the grid. You requested coordinates for "
               << "node number " << inode <<". Go check your mesh file!";
      Error( __FILE__, __LINE__ );
    }
    rfPoint = coords_[inode-1];
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
    Error( "Not implemented", __FILE__, __LINE__ );
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
                                        const StdVector<UInt> & connect ) {
    ENTER_FCN( "GridCFS::GetElemNodesCoord" );

    coordMat.Resize(dim_, connect.GetSize());
    for (UInt k=0; k < connect.GetSize(); k++)    
      for (UInt actDim=0; actDim < dim_; actDim++)
        coordMat[actDim][k] = coords_[connect[k]-1][actDim];
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                                          const StdVector<UInt> & nodeList,
                                          const StdVector<RegionIdType> 
                                          & regionIds ) {
    ENTER_FCN( "GridCFS::GetElemsNextToNodes" );
    Boolean belongs2Interface;

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
					 Boolean onlyLinNodes) {
    ENTER_FCN( "GridCFS::GetNodesOfElemList" );

    std::set<UInt> elemNodes;
    std::set<UInt>::iterator it;
    UInt iElem, iNode, numElemCorners;

    // First, create a set with node numbers of elements
    for ( iElem = 0; iElem < elemList.GetSize(); iElem++ ) {
      StdVector<UInt> const & connecth = elemList[iElem]->connect;
      
      if (onlyLinNodes == TRUE)
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
    SurfElem * myElem;
    surfElems_.Resize(elems.GetSize());
    
    for ( UInt iRegion = 0; iRegion < elems.GetSize(); iRegion++ ) {

      surfElems_[iRegion].Resize(elems[iRegion].GetSize());

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
          (*error) << "Pointer to first volume element is NULL for surface"
                   << " element no. "  
                   << surfElems_[iRegion][iSurfElem]->elemNum << ".\n"
                   << "Please check your mesh-file!";
          Error( __FILE__, __LINE__ );
        }

        CalcSurfNormal( normalUndefSign, *surfElems_[iRegion][iSurfElem] );
        CalcSurfNormalOutOfVol( normalDefSign, 
                                *surfElems_[iRegion][iSurfElem],
                                *surfElems_[iRegion][iSurfElem]->ptVolElem1 );


        // Check if all entries have the same sign by calulating
        // a scalar product between both vectors.
        // If it is positive, they point in the smae direction,
        // otherwise an angle of 180° lies in between.
        sign = normalUndefSign * normalDefSign;
        
        if ( sign > 0.0 ) {
          surfElems_[iRegion][iSurfElem]->normalSign = 1;
        } else {
          surfElems_[iRegion][iSurfElem]->normalSign = -1;
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
      help = Info->GenStr(i);
      help += "\t| ";
      help += regionNames_[i];
      help += '\n';
      Info->PrintF("",help.c_str());
    }
  }
  
  template<UInt DIM>
  void GridCFS<DIM>::CalcSurfNormal( Vector<Double> & n, 
                                     const Elem & surfElem ) {
    ENTER_FCN( "GridCFS::CalcSurfNormal" );
    
    //compute normal vector
    Matrix<Double>  ptCoord;

    GetElemNodesCoord(ptCoord, surfElem.connect);
    UInt surfCorners = surfElem.ptElem->GetNumCorners();
 
    // Check for dimension:
    if (surfElem.ptElem->GetDim() == 1) {
  
      // 1. step: compute vector perpendicular to line element
      // but without defined sign
      Double dx  = ptCoord[0][1] - ptCoord[0][0];
      Double dy  = ptCoord[1][1] - ptCoord[1][0];
      Double len = sqrt(dx*dx + dy*dy);
      if (len <= 0.0)
        Error("length of normal vector is zero!");
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
                                            const Elem & volElem)
  {
    ENTER_IFCN( "GridCFS::CalcSurfNormalOutOfVol" );

    //compute normal vector
    Matrix<Double>  ptVolCoord, ptSurfCoord;

    // First, calculate undefined normal
    CalcSurfNormal(n, surfElem);

    GetElemNodesCoord(ptSurfCoord, surfElem.connect);
    GetElemNodesCoord(ptVolCoord, volElem.connect);
  

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
      if ( (indexNode2-indexNode1 == -1 || 
            (indexNode2-indexNode1)- volCorners == -1 ) ) {
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
                                           Boolean isaxi) {
    ENTER_FCN( "GridCFS::GetVolumeOfRegion" );

    StdVector<Elem*> elems;
    Matrix<Double> cornerCoords;
    Double volume = 0.0;

    GetElems(elems,regionId);

    for( UInt i = 0; i < elems.GetSize(); i++ ) {
      GetElemNodesCoord(cornerCoords, elems[i]->connect);
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
        if (index >= maxnumnodes) {
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
        elems_[i][j]->refinementFlag=TRUE;
      }
    }
  }

#endif // end of ADAPTGRID


 

  //   template<Integer DIM>
  //   void GridCFS<DIM>::FormNeighborsLists()
  //   {
  //     ENTER_FCN( "GridCFS::FormNeighborsLists" );

  //     Integer i,j,k,n,m;

  //     Integer size=sd_.GetSize(); // number of subdomains in domain
  //     elNeighbors_ = new  StdVector<StdVector<Elem*> >*[size];
  //     // initialize                                           
  //     // for each subdomain of grid we create vector( noOfElems) with vector of
  //     // neighbors of each element
  //     for (i = 0; i < sd_.GetSize(); i++) {
  //       size=elems_[i].GetSize();
  //       elNeighbors_[i] = new StdVector<StdVector<Elem*> >(size);    
  //     }

  //     // vector with vector of neighbors-elements for each node
  //     vtNeighbors_.Resize(maxnumnodes_);

  //     // first main loop: form list with element-neighbors for each node
  //     for (i = 0; i < sd_.GetSize(); i++) {
  //       for (j = 0; j < elems_[i].GetSize(); j++) {
  //    Integer noOfVertices = elems_[i][j]->connect.GetSize();
  //    for (k = 0; k < noOfVertices; k++) {
  //      Integer id = elems_[i][j]->connect[k];
  //      vtNeighbors_[id-1].Push_back(elems_[i][j]);
  //    } // for loop over vertices, index k
  //       } // for loop over elements of the subdomain i, index j
  //     } // for loop over sd, index i

  //     //   // check of nodes list
  //     //   for (i=0; i<maxnumnodes_; i++) {
  //     //     for (j=0; j<vtNeighbors_[i].GetSize(); j++) {
  //     //       cout << " size " << vtNeighbors_[i].GetSize() << endl;
  //     //       Elem* ptE=vtNeighbors_[i][j];
  //     //       cout << " no of nodes " << i+1 << endl;
  //     //       cout << ptE->connect << endl;
  //     //     }
  //     //   }
      
  //     // second main loop: form list with element-neighbors for each element
  //     for (i = 0; i < sd_.GetSize(); i++) {   // do loop over subdomains
  //       for (j = 0; j < elems_[i].GetSize(); j++) {   // do loop over elements of subdomain
  //    Integer noOfVertices = elems_[i][j]->connect.GetSize(); 
  //    for (k = 0; k < noOfVertices; k++) {   // do loop over vertices of elements of the subdomain
  //      Integer id = elems_[i][j]->connect[k]; 
  //      for (n= 0; n < vtNeighbors_[id-1].GetSize(); n++) {  // do loop over list of neighbors for node
  //        Elem* ptel = vtNeighbors_[id-1][n];           

  //        if (ptel != elems_[i][j]) {           // check that this element is not the same element for which we are looking for neighbors
  //          Boolean flag = false;
  //          for (m = 0; m < (*elNeighbors_[i])[j].GetSize(); m++) { // check that this element is new in list of neighbors
  //            Elem* ptel_tmp = (*elNeighbors_[i])[j][m];
  //            if (ptel_tmp == ptel)
  //              flag = true;
  //          } // for loop over neighbors, index m
  //          if (!flag)
  //            (*elNeighbors_[i])[j].Push_back(ptel);

  //        } // if (ptle != elem)
  //      } // for loop over vtNeighbors, index n
  //    } // for loop over vertices, index k
  //       } // for loop over elements of the subdomain, index j
  //     } // for loop over sd, index i

  //     //check of element-neighbors
  //     //   for (i=0; i<sd_.GetSize(); i++) {
  //     //       for (j = 0; j < elems_[i].GetSize(); j++) {  
  //     //     cout << " element: " << j << " with connect " << elems_[i][j]->connect << endl; 

  //     //     for (m = 0; m < (*elNeighbors_[i])[j].GetSize(); m++) { // loop over neigbors for elem j
  //     //       cout << " neighbors: " << m << endl;   
  //     //       Elem * ptElm=(*elNeighbors_[i])[j][m];
  //     //       cout << ptElm->connect;
  //     //     }
  //     //       }
  //     //   }

  //   } // end of function FormNeighborsLists

  //   //! return pointer to vector of element-neighbors for the element with number noOfElem
  //   template<Integer DIM>
  //   StdVector<Elem*> * GridCFS<DIM>::GetptNeighboursOfElem(const Integer noOfElem, const std::string color)
  //   {
  
  //     StdVector<Elem*> * result;

  //     if (elNeighbors_==NULL) Error("You can't use function GetptNeighboursOfElem, since list of neighbors is not formed. You should use function FormNeighborsList before.");

  //     Integer i;
  //     for (i=0; i<sd_.GetSize(); i++) {
  //       if (sd_[i]==color) {          
  //         result=&(*elNeighbors_[i])[noOfElem];
  //       }
  //     }

  //     return result;
  //   }

  

  //   template<Integer DIM>
  //   void GridCFS<DIM>::FormNeighbors4NodesOfElements(const StdVector<Elem*>& elems, 
  //                                               StdVector<StdVector<Elem*> > &nodeNeighbors, 
  //                                               StdVector<Integer> & map)
  //   {
  //     ENTER_FCN( "GridCFS::FormNeighbors4NodeOfElements" );

  //     Integer iel,k;
  //     CalcNumberOfNodesInPatch(elems,map);
  //     Integer maxnumnodes=map.GetSize();
  //     // calculation number of nodes in patch of elements

  //     // vector with vector of neighbors-elements for each node
  //     nodeNeighbors.Resize(maxnumnodes);    

  //     // first main loop: form list with element-neighbors for each node
  //     for (iel = 0; iel < elems.GetSize(); iel++) {
  //       Integer noOfVertices = elems[iel]->connect.GetSize();
  //       for (k = 0; k < noOfVertices; k++) {
  //    Integer id = elems[iel]->connect[k];
  //    Integer imp;
  //    for (imp=0;imp<map.GetSize();imp++) 
  //      {
  //        if (id==map[imp]) 
  //          {
  //            break;
  //          }     
  //      }     
  //    nodeNeighbors[imp].Push_back(elems[iel]);
  //       } // for loop over vertices, index k
  //     } // for loop over elements iel

  //   } // end of function FormNeighbors4NodesOfElements




  //  template<Integer DIM> void
  //   GridCFS<DIM>::DefineBelonging4Elems(const StdVector<Elem*>& elemsSurf, 
  //                                  const StdVector<Elem*>&elems, 
  //                                  StdVector<Elem*> & belongingSE)
  //   {
  //     ENTER_FCN( "GridCFS:DefineBelonging4Elems" );

  //     Integer noOfSurfElems=elemsSurf.GetSize();
  //     belongingSE.Resize(noOfSurfElems);


  //     // form list with neighbors for each nodes in patch of boundary elements
  //     StdVector<Integer> map;
  //     StdVector<StdVector<Elem*> > listNeighbors;
  //     FormNeighbors4NodesOfElements(elems,listNeighbors,map);

  //     Integer ise,je;
  //     for (ise=0; ise<noOfSurfElems; ise++) 
  //       { // loop over surface elements
      
  //    Boolean FoundNd=FALSE;
  //    Elem * ptAuxElem;
      
  //    StdVector<Integer> &connectSE=elemsSurf[ise]->connect;
      
  //    // get list of neighbors for first node of the surface element
  //    Integer imp;   // get local number for this node
  //    for (imp=0; imp<map.GetSize(); imp++)
  //      if (connectSE[0]==map[imp]) 
  //        break;
      
    
  //    StdVector<Elem*> &listNeigh4Elem=listNeighbors[imp];
      
  //    // loop over list of neighbors
  //    Integer ine;
  //    for (ine=0;ine<listNeigh4Elem.GetSize();ine++)
  //      {
  //        ptAuxElem=listNeigh4Elem[ine];
          
  //        // check is there other vertices of element
  //        // loop over other nodes of surf element
  //        for (je=1;je<connectSE.GetSize();je++) {
  //          Integer verSE=connectSE[je];
            
  //          //loop over vertices of the element
  //          FoundNd=FALSE;
  //          StdVector<Integer> &vertices=ptAuxElem->connect;
  //          Integer ivt;
  //          for(ivt=0;ivt<vertices.GetSize();ivt++) {
  //            if (verSE==vertices[ivt]) {
  //              FoundNd=TRUE;
  //              break;
  //            }
  //          } // end of loop over vertices of neigh-element
            
  //          if (!FoundNd) break;
  //        } // end of loop over nodes of surf element
          
  //        if (FoundNd) {
  //          belongingSE[ise]=ptAuxElem;
  //          break;
  //        }
  //      } // end loop over neighbors 
      
  //       } // loop over Surf element 
  //   }



//   template<UInt DIM> void
//   GridCFS<DIM>::CalcNumberOfNodesInPatch(const StdVector<Elem*> & patch,
// 					 StdVector<Integer> & map, 
// 					 Boolean OnlyLinNodes)
//   {
//     ENTER_FCN( "GridCFS::CalcNumberOfNodesInPatch" );

//     Integer iels,ivc,imp;
//     StdVector<Integer> vec_connect;
//     Boolean NewNode;

//     for (iels=0; iels<patch.GetSize(); iels++) // loop over elements in patch
//       {
     
// 	vec_connect=patch[iels]->connect;
// 	Integer numElemCorners;
// 	if (OnlyLinNodes == TRUE)
// 	  numElemCorners = patch[iels]->ptElem->GetNumCorners();
// 	else
// 	  numElemCorners = vec_connect.GetSize();
	
// 	for (ivc=0; ivc<numElemCorners; ivc++) {
// 	  NewNode=TRUE;
// 	  // loop over vector with global nodes for previous elements
// 	  for (imp=0; imp<map.GetSize(); imp++) {
// 	    // check that this node is not new
// 	    if (map[imp] == vec_connect[ivc]) { 
// 	      NewNode=FALSE;
// 	    }	 
// 	  }

// 	  if (NewNode) {
// 	    map.Push_back(vec_connect[ivc]);
// 	  }

// 	} // end of loop over nodes in element

//       } // end of loop over elements in patch   
//   }




  //   template<Integer DIM> void
  //   GridCFS<DIM>::GetVolNeighboursForSurf(const StdVector<Elem*> & surfElems,
  //                                    const StdVector<std::string> & neighRegions,
  //                                    StdVector<Elem*> & volElems,
  //                                    const Integer level)
  //   {
  //     ENTER_FCN( "GridCFS::GetVolNeighboursForSurf" );

  //     std::string errMsg;
  //     StdVector<Elem*> auxElems, neighElems;
  //     Integer noOfSurfElems=surfElems.GetSize();
  //     volElems.Resize(noOfSurfElems);


  //     // form list with neighbors for each nodes in patch of boundary elements
  //     StdVector<Integer> map;
  //     StdVector<StdVector<Elem*> > listNeighbors;
    
  //     // get all elements in neighborin region, from where later
  //     // the volume elems are picked
  //     for (Integer iSD=0; iSD<neighRegions.GetSize(); iSD++)
  //       {
  //    GetElemSD(auxElems, neighRegions[iSD], level);
  //    for (Integer iEl=0; iEl<auxElems.GetSize(); iEl++)
  //      neighElems.Push_back(auxElems[iEl]);
  //       }

  //     // create a list of elements in 'neighRegion', which have at least
  //     // one node in common with the element lying in neighRegions
  //     FormNeighbors4NodesOfElements(neighElems,listNeighbors,map);
  //     Integer ise,je;
  //     for (ise=0; ise<noOfSurfElems; ise++) 
  //       { // loop over surface elements
        
  //    Boolean FoundNd=FALSE;
  //    Elem * ptAuxElem;
        
  //    StdVector<Integer> &connectSE=surfElems[ise]->connect;
        
  //    // get list of neighbors for first node of the surface element
  //    Integer imp;   // get local number for this node
  //    for (imp=0; imp<map.GetSize(); imp++)
  //      if (connectSE[0]==map[imp]) 
  //        break;
        
        
  //    StdVector<Elem*> &listNeigh4Elem=listNeighbors[imp];
        
  //    // loop over list of neighbors
  //    Integer ine;
  //    for (ine=0;ine<listNeigh4Elem.GetSize();ine++)
  //      {
  //        ptAuxElem=listNeigh4Elem[ine];
            
  //        // check are there other vertices of element
  //        // loop over other nodes of surf element
  //        for (je=1;je<connectSE.GetSize();je++) {
  //          Integer verSE=connectSE[je];
              
  //          //loop over vertices of the element
  //          FoundNd=FALSE;
  //          StdVector<Integer> &vertices=ptAuxElem->connect;
  //          Integer ivt;
  //          for(ivt=0;ivt<vertices.GetSize();ivt++) {
  //            if (verSE==vertices[ivt]) {
  //              FoundNd=TRUE;
  //              break;
  //            }
  //          } // end of loop over vertices of neigh-element
              
  //          if (!FoundNd) break;
  //        } // end of loop over nodes of surf element
            
  //        if (FoundNd) {
  //          volElems[ise]=ptAuxElem;
  //          break;
  //        }
  //      } // end loop over all elements in neighbouring region
  //    if (!FoundNd)
  //      {
  //        errMsg  = "GridCFS::GetVolNeighboursForSurf: For the surface element with Nr. ";
  //        errMsg += Info->GenStr(surfElems[ise]->elemNum);
  //        errMsg += " an according volume element was not ";
  //        errMsg += "found in the regions '";
  //        for (Integer j=0; j<neighRegions.GetSize()-1; j++) 
  //          errMsg += neighRegions[j] + "', '";
  //        errMsg += neighRegions[neighRegions.GetSize()-1] + "'.\n";
  //        errMsg += "Please make sure, that for each ";
  //        errMsg += "surface element there is exactly ONE volume element ";
  //        errMsg += "in the speciefied neighbouring region.";
            
            
  //        Error(errMsg.c_str(), __FILE__, __LINE__);
  //      }
        
        
  //       } // loop over all Surface elements 
  //   }
  
  



  //  template<>
  //   Double GridCFS<2>::CalcAreaElem(const Elem* elem)
  //   {
  //     ENTER_FCN( "GridCFS<Dim>::CalcAreaElem" );

  //     Double res;
  //     Double a,b,c,s;
  //     Point<2> A,B,C,D;

  //     Integertype=elem->ptElem->feType();
  //     switch(type) {
  //     case LINE: // line
  //       A=ptCoordinate_[elem->connect[0]-1];
  //       B=ptCoordinate_[elem->connect[1]-1];
  //       res=dist(A,B);
  //       break;
  //     case TRIA:     // triangle
  //       A=ptCoordinate_[elem->connect[0]-1];
  //       B=ptCoordinate_[elem->connect[1]-1];
  //       C=ptCoordinate_[elem->connect[2]-1];
  
  //       a=dist(A,B);
  //       b=dist(B,C);
  //       c=dist(A,C);
  
  //       s=(a+b+c)/2;
  //       res=sqrt(s*(s-a)*(s-b)*(s-c));

  //       break; 
  //     case QUAD:     // quadrilateral
  //       Double res1, res2;
  //       A=ptCoordinate_[elem->connect[0]-1];
  //       B=ptCoordinate_[elem->connect[1]-1];
  //       C=ptCoordinate_[elem->connect[2]-1];
  //       D=ptCoordinate_[elem->connect[3]-1];
   
  //       a=dist(A,B);
  //       b=dist(B,C);
  //       c=dist(A,C);
  //       s=(a+b+c)/2;
  //       res1=sqrt(s*(s-a)*(s-b)*(s-c));
  //       a=dist(C,D);
  //       b=dist(A,D);
  //       s=(a+b+c)/2;
  //       res2=sqrt(s*(s-a)*(s-b)*(s-c));
  //       res=res1+res2;
  //       break;

  //     default:
  //       Error("There isn't such kind of element's type",__FILE__,__LINE__);
  //     }

  //     return res;
  //   }




  //   template<>
  //   Double GridCFS<3>::CalcAreaElem(const Elem* elem)
  //   {
  //     Error("Not implemented yet",__FILE__,__LINE__);

  //     //just for SUN compiler
  //     Double dummy =0.0;
  //     return dummy;
  //   }


#ifdef __sgi
#pragma instantiate GridCFS<3>
#pragma instantiate GridCFS<2>
#endif
} // end namespace


