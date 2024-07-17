#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <set>
#include <time.h>

// icl on Windows does not know not/and/.. as in boosts return not af.test_and_set()
#if defined(WIN32) && defined(__INTEL_COMPILER)
  #include <ciso646>
#endif

// includes needed for writing geometry file
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include <def_use_cgal.hh>
#include <def_use_eigen.hh>

#include <boost/sort/sort.hpp>
#include "GridCFS.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "General/Exception.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"

#ifdef WIN32
#include <direct.h>
#endif

namespace CoupledField {

  // declare class specific logging stream
  DEFINE_LOG(gridcfs, "grid.cfs")

  GridCFS::GridCFS(UInt dim, PtrParamNode param, PtrParamNode info,
      const std::string &id, bool buildExtend)
  : Grid( param, info ) {
    gridId_ = id;
    isQuadratic_ = false;
    dim_ = dim;
    assert(dim > 0);

    numNodes_ = 0;
    numElems_ = 0;
    numFaces_ = 0;
    numEdges_ = 0;
    numNcSurfElems_ = 0;
    edgesMapped_ = false;
    facesMapped_ = false;
    maxNumElemNodes_ = 0;
    mapNodeToElems_.Resize(0);
    buildExtendedElemInfo_ = buildExtend;
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

  void GridCFS::CreateUserDefinedNodesElems() {

    // if no param object is present, just leave
    if (!param_) return;

    Vector<Double> locCoord(dim_), globCoord(dim_);
    Double compoValue;
    std::string coordSysId, compoName;
    std::string name;

    for( UInt iType = 0; iType < 2; iType++ ) {
      PtrParamNode listNode;
      ParamNodeList nodes;
      bool isNode = true;

      if( iType == 0 ) {
        // iterate over nodes
        listNode = param_->Get("domain")->Get("nodeList", ParamNode::PASS);
        isNode = true;
      } else {
        // iterate over elements
        listNode = param_->Get("domain")->Get("elemList", ParamNode::PASS);
        isNode = false;
      }

      if (listNode) {
        nodes = listNode->GetChildren();

        for( UInt i=0; i < nodes.GetSize(); i++ ) {

          // fetch name of nodes to be selected
          nodes[i]->GetValue("name", name);

          //check if named nodes are defined by specifying a region
          if (nodes[i]->Has("allNodesInRegion")) {
            std::string regName;
            nodes[i]->Get("allNodesInRegion")->GetValue("regName",regName);
            RegionIdType regId = region_.Parse(regName);

            StdVector<UInt> nodesInRegion;
            GetNodesByRegion(nodesInRegion,regId);
            AddNamedNodes(name,nodesInRegion);
          }

          // check if node is defined by point coord
          PtrParamNode coordNode = nodes[i]->Get("coord", ParamNode::PASS );
          if( coordNode ) {
            UInt compoIndex;
            locCoord.Init(0.0);
            globCoord.Init(0.0);
            
            coordSysId = "default";
            coordNode->GetValue("coordSysId", coordSysId, ParamNode::PASS);
            const CoordSystem* cosy = domain->GetCoordSystem(coordSysId);
            
            ParamNodeList compoList = coordNode->GetChildren();
            ParamNodeList::iterator compoIt = compoList.Begin(),
                                    endIt = compoList.End();
            
            for ( ; compoIt != endIt; ++compoIt ) {
              compoName = (*compoIt)->GetName();
              if ( compoName == "coordSysId" ) continue;
              
              compoValue = (*compoIt)->MathParse<Double>();
              if ( compoValue == 0.0) continue;
              
              if ( isAxi_ && (compoName == "x" || compoName == "y") ) {
                EXCEPTION(listNode->GetName() << " '" << name
                          << "': Coordinate components must be 'r', 'z' "
                          << "in an axisymmetric simulation.");
              }
              if ( isAxi_ && compoName == "r" ) compoName = "x";
              if ( isAxi_ && compoName == "z" ) compoName = "y";
              
              try {
                compoIndex = cosy->GetVecComponent(compoName)-1;
              } catch (Exception &ex) {
                RETHROW_EXCEPTION(ex, "Unable to create "
                                  << listNode->GetName()
                                  << " '" << name << "'.");
              }
              
              locCoord[compoIndex] = compoValue;
            }

            cosy->Local2GlobalCoord(globCoord, locCoord);
            
            StdVector<UInt> entityNum(1);
            entityNum[0] = FindEntityMinDistance( isNode, globCoord );

            // add node / element
            if( isNode ) {
              AddNamedNodes( name, entityNum );
            } else {
              AddNamedElems( name, entityNum );
            }

          }

          // check if node is defined by parametric description
          PtrParamNode listNode = nodes[i]->Get("list", ParamNode::PASS );
          if( listNode ) {

            std::string gridId, coordSysId;

            // make sure, that this is the correct grid
            if ( listNode->Get("gridId", ParamNode::INSERT)->As<std::string>()
                != gridId_ ) return;

            // get coordinate system
            listNode->GetValue( "coordSysId", coordSysId );

            StdVector<PointSelection> selections;

            // get free component
            ParamNodeList  freeNodes =
              listNode->GetList("freeCoord");
            for( UInt i = 0; i < freeNodes.GetSize(); i++ ) {
              PtrParamNode actNode = freeNodes[i];
              PointSelection actSel;
              actSel.isFree = true;
              actSel.comp = actNode->Get( "comp" )->As<std::string>();
              actSel.start = actNode->Get( "start" )->MathParse<Double>();
              actSel.stop = actNode->Get( "stop" )->MathParse<Double>();
              actSel.inc = actNode->Get( "inc" )->MathParse<Double>();
              selections.Push_back( actSel);
            }

            // get fixed component(s)
            ParamNodeList fixedNodes =
                listNode->GetList("fixedCoord");
            for( UInt i = 0; i < fixedNodes.GetSize(); i++ ) {
              PtrParamNode actNode = fixedNodes[i];
              PointSelection actSel;
              actSel.isFree = false;
              actSel.comp = actNode->Get( "comp" )->As<std::string>();
              actSel.value = actNode->Get( "value" )->As<std::string>();
              selections.Push_back( actSel);
            }

            // Ensure, that we have as many entries in the
            // vector as dimension
            if( dim_ != selections.GetSize() )
              EXCEPTION("The node parametrization for '" << name  << "' does not match the mesh dimensions");

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
    //if( nameTypeMap_.find( name) != nameTypeMap_.end() ) {
    //  EXCEPTION( "Entities with name " << name
    //             << " are already defined" );
    //}

    // get coordinate system
    CoordSystem * cosy = nullptr;
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
    StdVector<unsigned int> handles(dim_);
    StdVector<Vector<Double> > sampleVals(dim_);
    for( UInt iDim = 0; iDim < dim_; iDim++ ) {
      if( !coords[iDim].isFree ) {
        handles[iDim] = parser->GetNewHandle();
        for( UInt iDim2 = 0; iDim2 < dim_; iDim2++ ) {
          if( !coords[iDim2].isFree )
            parser->SetValue( handles[iDim], coords[iDim2].comp, 0.0 );
        }
        parser->SetExpr( handles[iDim], coords[iDim].value );
        sampleVals[iDim].Resize(1);
        sampleVals[iDim].Init();
      } else {
        UInt numSamples  =
          UInt( floor ( abs(coords[iDim].stop-coords[iDim].start)
                        / coords[iDim].inc ) )+1;
        sampleVals[iDim].Resize( numSamples );
        for( UInt iSample = 0; iSample < numSamples; iSample++ ) {
          sampleVals[iDim][iSample] = coords[iDim].start
            + iSample * coords[iDim].inc;
        }
      }
    }


    std::set<UInt> entityNums;
    Vector<Double> locCoord( dim_), globCoord( dim_ );
    Vector<Double> actEntCoord(dim_), temp;

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
            parser->SetValue( handles[2], coords[1].comp, sampleVals[1][iSample2] );
            sampleVals[2][0] = parser->Eval(handles[2] );
          }
          for( UInt iSample3 = 0; iSample3 < sampleVals[2].GetSize(); iSample3++ ) {
            locCoord[dofs[2]] = sampleVals[2][iSample3];
            cosy->Local2GlobalCoord( globCoord, locCoord );
            entityNums.insert(FindEntityMinDistance( isNode, globCoord ));
          }

        } // dim == 3
        else {
          cosy->Local2GlobalCoord( globCoord, locCoord );
          entityNums.insert(FindEntityMinDistance( isNode, globCoord ));
        }
      }
    }

    // copy unique set of nodes / elements to vector
    StdVector<UInt> entityNumVec(entityNums.size());
    std::set<UInt>::const_iterator setIt = entityNums.begin();
    for(UInt pos = 0 ; setIt != entityNums.end(); ++setIt, ++pos ) {
      entityNumVec[pos] = *setIt;
    }
    
    if( isNode) {

      // add named nodes to grid
      AddNamedNodes( name, entityNumVec );
    } else {
      AddNamedElems( name, entityNumVec );
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

        // Check, if element has same dimension as grid
        // -> We want to find only volume elements
        if( Elem::shapes[orderedElems_[iElem]->type].dim == dim_ ) {
          shared_ptr<ElemShapeMap> esm = GetElemShapeMap(orderedElems_[iElem]);
          esm->GetGlobMidPoint(actEntCoord);
          temp = (actEntCoord-coord );
          entityDist[iElem] = temp.NormL2();
        } else {
          entityDist[iElem] = 1e16;
        }
      } // elements

    }

    // find minimum entry in the vector
    std::vector<Double>::iterator it ;
    it = min_element(entityDist.begin(), entityDist.end());
    entityNum = std::distance(entityDist.begin(), it) + 1;
    return entityNum;

}

  void GridCFS::MapFaces() {

    LOG_DBG(gridcfs) << "Starting to map faces ";


    // assert that any mesh was already read in
    assert( isInitialized_ == true );

    // If faces are already mapped ->leave
    if( facesMapped_ == true ) {
      return;
    }

    // map containing the nodes for each face
    std::map<std::set<UInt>, UInt > faceNums;
    
    // create counter for number of faces
    UInt actFaceNum = 1;

    // vectors for local / global nodes
    StdVector<UInt> faceIndices, faceNodes;

    // iterate over all elements
    for( UInt iElem = 0; iElem < orderedElems_.GetSize(); iElem++ ) {

      // remember current element
      Elem & actElem = *orderedElems_[iElem];
      ElemShape & actShape = Elem::shapes[actElem.type];

      // if element is of wrong dimension (surface element )
      // ->leave
      //if ( actElem.ptElem->GetDim() < dim_ ) { continue; }

      // get number of element faces
      UInt numFaces = actShape.numFaces;

      // adapt size of faces and orientation array of element
      actElem.extended->faces.Resize( numFaces );
      actElem.extended->faceFlags.Resize( numFaces );

      // iterate over all faces of this element
      for( UInt iFace = 0; iFace < numFaces; iFace++ ) {

        // get local nodal indices of current face
        faceIndices = actShape.faceVertices[iFace];

        // create new Face object
        Face actFace;
        StdVector<UInt> faceNodes;
        UInt numFaceIndices = faceIndices.GetSize();
        faceNodes.Resize( numFaceIndices );

        // insert node numbers into current face definition
        for( UInt iNode = 0; iNode < numFaceIndices; iNode++ ) {
          faceNodes[iNode] = actElem.connect[faceIndices[iNode]-1];
        }

        LOG_DBG3(gridcfs) << "Checking face with nodes : "
                          << faceNodes.Serialize();

        // Re-orientate face to match global orientation and
        // obtain the orientation flags
        std::bitset<5> orientation;
        actFace.Normalize( orientation, faceNodes );

        LOG_DBG3(gridcfs) << "Normalized : "
                          << faceNodes.Serialize();
        LOG_DBG3(gridcfs) << "Orientation: " << orientation;

        // copy all nodes in a set
        std::set<UInt> nodeSet;
        nodeSet.insert(faceNodes.Begin(), faceNodes.End());

        // check if face was already numbered
        if( faceNums.find(nodeSet) == faceNums.end() ) {
        
        //if( faceNums_.find( actFace ) == faceNums_.end() ) {
          LOG_DBG2(gridcfs) << "Adding face number "
                            << actFaceNum << std::endl;
          faceNums[nodeSet] = actFaceNum;
          actFace.neighbors.Push_back( orderedElems_[iElem]);
          actElem.extended->faces[iFace] = actFaceNum;
          faces_.Push_back( actFace );
          actFaceNum++;
        } else {
          UInt faceNum = faceNums[nodeSet];
          actElem.extended->faces[iFace] = faceNum;
          faces_[faceNum-1].neighbors.Push_back( orderedElems_[iElem]);
          LOG_DBG2(gridcfs) << "--> already defined\n";
        }

        // Set also orientation flags for face
        actElem.extended->faceFlags[iFace] = orientation;
      } // loop over faces

      // Print information about connectivity and faces
      LOG_DBG2(gridcfs) << "Elem Nr. " << actElem.elemNum;
      LOG_DBG2(gridcfs) << "===================";
      LOG_DBG2(gridcfs) << "Connectivity: " << actElem.connect.Serialize();
      LOG_DBG2(gridcfs) << "Faces: " << actElem.extended->faces.Serialize();

      LOG_DBG(gridcfs) << "Finished to map faces\n";

    } // loop over elements


    // Trim vector containing faces
    faces_.Trim();
    
    // Set flag for mapping of sub-entities
    numFaces_ = actFaceNum-1;
    LOG_DBG2(gridcfs) << "Total number of faces: " << numFaces_;

    // Set flag
    facesMapped_ = true;

    LOG_DBG(gridcfs) << "Finished mapping faces\n";

  }



  void GridCFS::MapEdges() {

    LOG_DBG(gridcfs) << "Starting to map edges";

    // assert that any mesh was already read in
    assert( isInitialized_ == true );

    // If edges/surfaces were already mapped ->leave
    if( edgesMapped_ == true ) {
      return;
    }

    // map with edge numbers
    std::map<std::set<UInt>, UInt > edgeNums;
    
    // create counter for number of edges
    UInt actEdgeNum = 1;
    StdVector<UInt> locEdge(2), globEdge(2);

    // iterate over all elements
    for( UInt iElem = 0; iElem < orderedElems_.GetSize(); iElem++ ) {

      // remember current element
      Elem & actElem = *orderedElems_[iElem];
      ElemShape & actShape = Elem::shapes[actElem.type];
      
      // get number of edges
      UInt numEdges= actShape.numEdges;

      //in case the extended element info is not yet created we do it here
      if(!actElem.extended){
        actElem.extended = new ExtendedElementInfo();
      }

      // adapt size of edge number array of element
      actElem.extended->edges.Resize( numEdges );

      // iterate over all edges of this element
      for( UInt iEdge = 0; iEdge < numEdges; iEdge++ ) {

        // get local edge indices
        locEdge = actShape.edgeVertices[iEdge];
        // create new edge
        Edge actEdge;
        StdVector<UInt> edgeNodes(2);
        edgeNodes[0] = actElem.connect[locEdge[0]-1];
        edgeNodes[1] = actElem.connect[locEdge[1]-1];

        // check if ordering is correct
        Integer orientation = 1;

        if( edgeNodes[1] < edgeNodes[0] ) {
          std::swap(edgeNodes[1], edgeNodes[0]);
          // swap factor for orientation
          orientation = -1;
        }
        
        std::set<UInt> nodeSet;
        nodeSet.insert(edgeNodes.Begin(), edgeNodes.End());

        // check if edge was already numbered
        if( edgeNums.find( nodeSet) == edgeNums.end() ) {
          LOG_DBG2(gridcfs) << "Adding edge number " << actEdgeNum;
          LOG_DBG3(gridcfs) << "with nodes: " << edgeNodes.Serialize() << std::endl;
          edgeNums[nodeSet] = actEdgeNum;
          actEdge.neighbors.Push_back(orderedElems_[iElem]);
          actElem.extended->edges[iEdge] = actEdgeNum*orientation;
          edges_.Push_back( actEdge );
          actEdgeNum++;
        } else {
          UInt edgeNum = edgeNums[nodeSet];
          LOG_DBG3(gridcfs) << "Edge was already mapped with number #" << edgeNum << std::endl;
          actElem.extended->edges[iEdge] = edgeNum*orientation;
          edges_[edgeNum-1].neighbors.Push_back(orderedElems_[iElem]);
        }
        

      }

      // Print information about connectivity and edges
      LOG_DBG2(gridcfs) << "Elem Nr. " << actElem.elemNum;
      LOG_DBG2(gridcfs) << "===================";
      LOG_DBG2(gridcfs) << "Connectivity: " << actElem.connect.Serialize();
      LOG_DBG2(gridcfs) << "Edges: " << actElem.extended->edges.Serialize();

      LOG_DBG(gridcfs) << "Finished to map edges\n";
    }

    // Trim vector containing faces
    edges_.Trim();
        
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

  RegionIdType GridCFS::CheckPatternRegion(StdVector<bool>& replace)
  {
    // check the input for a pattern region
    std::string name = "";
    std::string file = "";

    // be sensitive to the cfstool case
    if(!param_ || !param_->Has("domain/regionList")) return NO_REGION_ID;

    ParamNodeList list = param_->Get("domain/regionList")->GetList("region");
    for(UInt i = 0; i < list.GetSize(); i++)
    {
      if(list[i]->Has("pattern"))
      {
        if(name != "") EXCEPTION("Only a single region with pattern attribute allowed");
        name = list[i]->Get("name")->As<std::string>();
        file = list[i]->Get("pattern")->As<std::string>();
      }
    }

    // quick exit
    if(name == "") return NO_REGION_ID;

    // read the pattern file
    PtrParamNode xml = XmlReader::ParseFile(file);

    // check this file
    if (xml->Count("set") == 0)
      throw Exception("There are no design sets in the pattern file " + file);

    ParamNodeList elems = xml->GetList("set").Last()->GetList("element");

    if(elems.GetSize() != orderedElems_.GetSize())
      EXCEPTION("There are " << elems.GetSize() << " elements in pattern file " << file << " but " << orderedElems_.GetSize() << " in the mesh.");

    // blow up replace to the size of all elements
    replace.Resize(orderedElems_.GetSize() + 1, false); // numbers are 1-based

    for(UInt i = 0, n = elems.GetSize(); i < n; i++)
    {
      double val = elems[i]->Get("design")->As<double>();
      if(val >= 1.0)
      {
        UInt nr = elems[i]->Get("nr")->As<UInt>();
        replace[nr] = true;
      }
    }

    // finally add the region
    return AddRegion(name, false);
  }

  void GridCFS::FinishInit()
  {
    
    LOG_TRACE(gridcfs) << "Finalizing GridCFS (FinishInit)";
    volRegionIds_.Clear();
    surfRegionIds_.Clear();
    volElems_.Clear();
    surfElems_.Clear();
    maxNumElemNodes_ = 0;

    UInt e;
    UInt numElems = orderedElems_.GetSize();
    UInt numNodes = 0;

    // the pattern regions are defined by density files in domain/regionList/region with pattern set
    StdVector<bool> pattern;
    RegionIdType pattern_reg = CheckPatternRegion(pattern); // assumed to be volume region only!

    std::map<RegionIdType, StdVector<Elem*> > volRegionElems, surfRegionElems;
    std::map<RegionIdType, StdVector<UInt> > volRegionNodes, surfRegionNodes;
    std::map<RegionIdType, UInt > regionDims;

    LOG_DBG(gridcfs) << "Determine list of surface elements: " << numElems;
    
    // set of elements, which get surface-mapped
    StdVector<Elem*> surfElems;

    // loop over all elements
    for(e=0; e<numElems; e++)
    {
      Elem* el = orderedElems_[e];
      
      Elem::FEType type = el->type;
      numNodes = Elem::shapes[type].numNodes;

      maxNumElemNodes_ = maxNumElemNodes_ < numNodes ?
                         numNodes : maxNumElemNodes_;

      // in pattern case replace element region
      if(pattern_reg != NO_REGION_ID && pattern[el->elemNum])
        el->regionId = pattern_reg;

      // Insert dimension of first element in region into regionDims map
      // If elements with different dimension are encountered issue an exception
      if(!regionDims[el->regionId]) 
      {
        regionDims[el->regionId] = Elem::shapes[type].dim;
        LOG_DBG3(gridcfs) << "\tRegion '" 
                          << region_.ToString(el->regionId)
                          << "' has dimension " << regionDims[el->regionId];
        
      }
      else
      {
        // Elements in the region with id NO_REGION_ID may have arbitrary
        // dimensions.
        if( el->regionId != NO_REGION_ID &&
            regionDims[el->regionId] != Elem::shapes[type].dim )
        {
          EXCEPTION("Elements with different dimensions have been "
                    << "encountered in region '" << region_.ToString(el->regionId)
                    << "'!\nThe error occurred while examining element "
                    << el->elemNum << ".\nPlease check your mesh file!");
        }    
      }


      // get dimension of element
      UInt elemDim = Elem::shapes[type].dim;
      bool isSurfElem = (dim_ - elemDim) == 1 ? true : false; 

      // decide, what to do with the element
      if( isSurfElem ) {
        surfElems.Push_back( el );
        LOG_DBG3(gridcfs) << "\tadding elem #" << el->elemNum 
                          << " to list of surface elements.";
      } else  {
        if( el->regionId != NO_REGION_ID ) {
          volRegionElems[el->regionId].Push_back(el);
          LOG_DBG3(gridcfs) << "\tadding elem #" << el->elemNum 
                                    << " to list of volume elements.";
          const StdVector<UInt> & connect = el->connect;
          UInt numNodes = connect.GetSize();
          for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
              volRegionNodes[el->regionId].Push_back(connect[iNode]);
          }
        }
      } // if surfElem
    } // loop ordered elems

    
    std::map<RegionIdType, StdVector<Elem*> >::iterator regionElemIt, regionElemEnd;
    std::map<RegionIdType, StdVector<UInt> >::iterator regionNodeIt;
    std::set<UInt>::iterator setIt, setEnd;
    UInt region = 0;

    
    regionElemIt = volRegionElems.begin();
    regionElemEnd = volRegionElems.end();
    regionNodeIt = volRegionNodes.begin();

    for( ; regionElemIt != regionElemEnd; 
        regionElemIt++, regionNodeIt++, region++) {
      StdVector<UInt> elemNodes;
      regionElemIt->second.Trim();
      volElems_.Push_back(regionElemIt->second);
      regionData[regionElemIt->first].type = VOLUME_REGION;
      regionData[regionElemIt->first].type_idx = volRegionIds_.GetSize();
      volRegionIds_.Push_back(regionElemIt->first);
      // make number of region nodes unique
      StdVector<UInt> & regionNodes = regionNodeIt->second;

      // spreadsort is about 20 times faster than std::sort (tested for 1e7 nodes)
      boost::sort::spreadsort::spreadsort( regionNodes.Begin(), regionNodes.End() );
      //std::sort( regionNodes.Begin(), regionNodes.End() );

      StdVector<UInt>::iterator uIt;
      uIt = std::unique( regionNodes.Begin(), regionNodes.End() );
      UInt numRegionNodes = std::distance(regionNodes.Begin(), uIt);
      numVolElemNodes_.Push_back(numRegionNodes);
    }
    volRegionNodes.clear();

    // Call creation of surface elements
    StdVector<SurfElem*> mappedSurfElems;
    CreateSurfaceElements( surfElems, mappedSurfElems );

    // Iterate over all surface elements and put their nodes and elements
    // according to their region id
    UInt numSurfElems = mappedSurfElems.GetSize();
    for( UInt iSurfEl = 0; iSurfEl < numSurfElems; ++iSurfEl ) {
      SurfElem * surfEl = mappedSurfElems[iSurfEl];
      UInt elemNum = surfEl->elemNum;
      orderedElems_[elemNum-1] = surfEl;
      LOG_DBG3(gridcfs) << "Adding element #" << elemNum
          << " to list of surface elements";
      if( surfEl->regionId != NO_REGION_ID ) {
        surfRegionElems[surfEl->regionId].Push_back( surfEl );
        const StdVector<UInt> & connect = surfEl->connect;
        UInt numNodes = connect.GetSize();
        for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
          surfRegionNodes[surfEl->regionId].Push_back(connect[iNode]);
        }
      }
    }

    regionElemIt = surfRegionElems.begin();
    regionElemEnd = surfRegionElems.end();
    regionNodeIt = surfRegionNodes.begin();

    for( ; regionElemIt != regionElemEnd; regionElemIt++, regionNodeIt++, region++)
    {
      RegionIdType region = regionElemIt->first;
      regionElemIt->second.Trim();
      surfElems_.Push_back(regionElemIt->second);
      regionData[region].type = SURFACE_REGION;
      regionData[region].type_idx = surfRegionIds_.GetSize();
      surfRegionIds_.Push_back(region);
      
      // make number of region nodes unique
      StdVector<UInt> & regionNodes = regionNodeIt->second;
      std::sort( regionNodes.Begin(), regionNodes.End() );
      StdVector<UInt>::iterator uIt;
      uIt = std::unique( regionNodes.Begin(), regionNodes.End() );
      UInt numRegionNodes = std::distance(regionNodes.Begin(), uIt);
      numSurfElemNodes_.Push_back( numRegionNodes );
    }

    if (regionData.GetSize() != (volElems_.GetSize() + surfElems_.GetSize()))
    {
      for(UInt i = 0; i < regionData.GetSize(); i++)
      {
        if (volRegionIds_.Find(regionData[i].id) == -1
            && surfRegionIds_.Find(regionData[i].id) == -1)
        {
          EXCEPTION("Region has neither volume nor surface elements: " \
              << region_.ToString(regionData[i].id));
        }
      }
    }

    // in case of internalMesh the region is already marked as regular
    // so we can skip the test here
    for(unsigned int i = 0; i < regionData.GetSize(); i++) {
      regionData[i].regular = (regionData[i].regular || CheckForRegularRegion(i));
      LOG_DBG(gridcfs) << "FI: region " << regionData[i].name << " regular? " << regionData[i].regular;
    }

    isInitialized_ = true;

    // Try to fix problems due to negative Jacobian determinants
    LOG_DBG(gridcfs) << "Trying to correct negative Jacobians. -> CoupledField::LagrangeElemShapeMap::CalcJDet -> CoupledField::FeH1::GetLocDerivShFnc";
    CorrectElementConnectivities();

    // make named nodes from lines
    makeNameNodesFromLines();
    
    // Select nodes / elements according to the users specification in the
    // parameter file
    CreateUserDefinedNodesElems();
    
    // In the end, trim all vectors, i.e. delete any non-used memory from its
    // capacity.
    coords_.Trim();
    deltCoords_.Trim();
    orderedElems_.Trim();

    // print information to file - checks for exportGrid
    if(domain && info_) // we have no domain in the cfstool case
    {
      ToInfo(info_->Get(ParamNode::HEADER)->Get("domain"));

      StdVector<unsigned int> reg = CalcRegulardGridDiscretization();
      if(!reg.IsEmpty()) {
        MathParser* mp = domain->GetMathParser();
        mp->SetValue(MathParser::GLOB_HANDLER, "nx", reg[0]);
        mp->SetValue(MathParser::GLOB_HANDLER, "ny", reg[1]);
        mp->SetValue(MathParser::GLOB_HANDLER, "nz", reg[2]);
      }
    }
  }
  
  void GridCFS::GenerateDGSurfaceElemes(std::set<RegionIdType> regionList,
                                        StdVector<shared_ptr<NcSurfElem> > & interiorSurfElems,
                                        StdVector<shared_ptr<NcSurfElem> > & exteriorSurfElems){

    assert(edgesMapped_ == true);
    assert(facesMapped_ == true);
    std::cout << "   ++ Generating Surface Elements for DG-Methods ...";
    std::cout.flush();
    //ok now we generate an entitylist conaining all volume elements
    //ElemList volElemList;
    std::set<RegionIdType>::iterator regIt = regionList.begin();
    std::vector<Elem*> elemList;
    for(;regIt != regionList.end();++regIt){
      //volElemList.
      StdVector<Elem*> tmpEl;
      GetElems(tmpEl,*regIt);
      elemList.insert(elemList.end(),tmpEl.Begin(),tmpEl.End());
    }
    //reserve some memory such that pushbacks are less expensive
    //so we assume four subelements per element in 2d and 6 in 3d
    UInt memReserve = 0;
    if(dim_==2){
      memReserve = elemList.size() * 4;
    }else{
      memReserve = elemList.size() * 6;
    }
    StdVector<shared_ptr<NcSurfElem> > allelems;
    allelems.Reserve(memReserve);
    //UInt curNumelems = numElems_;
    Vector<Double> normalUndefSign, normalDefSign;
    for(UInt i =0;i<elemList.size();++i){
      Elem* curE = elemList[i];
      ElemShape curShape = Elem::shapes[curE->type];
      StdVector<Elem::FEType> subelems = curShape.surfElemTypes;
      UInt numSubElems = curShape.numSurfElems;
      for(UInt aSub =0; aSub<numSubElems ; ++aSub,++numNcSurfElems_){
        //create the subelement
        shared_ptr<NcSurfElem> curSurf(new NcSurfElem());
        //give it a number
        curSurf->elemNum = numElems_+numNcSurfElems_+1;
        curSurf->type = subelems[aSub];
        curSurf->regionId = curE->regionId;
        if(curShape.dim == 3){

          Face curFace = faces_[curE->extended->faces[aSub]-1];
          //for the connect array we have to extract the connectivity
          // right from the volume element as the normalized node array has no longer the expected
          // ordering vertexNodes-edgeNodes-FaceNodes

          //obtain the faceNodeArray
          StdVector<UInt> fNodes = curShape.faceNodes[aSub];
          UInt numNodes = curShape.faceNodes[aSub].GetSize();
          curSurf->connect.Resize(numNodes);
          //now we loop over the face nodes and assign the corresponding nodes from the volume element
          for(UInt i=0;i<curShape.faceNodes[aSub].GetSize();++i){
        	  curSurf->connect[i] = curE->connect[fNodes[i]-1];
          }
          curSurf->extended->faces.Push_back(curE->extended->faces[aSub]);
          curSurf->extended->faceFlags.Push_back(curE->extended->faceFlags[aSub]);
          //ok we know hat the face vertex array contains the edges of the face
          // in order to obtain the correct ordering of edges we loop over
          //the array and search for each pair of face vertices the correct
          // edge e.g. {1,2} {2,3} {3,4} {4,1}
          UInt numFVertices = curShape.faceVertices[aSub].GetSize();
          for(UInt i=0;i<=numFVertices;++i){
        	  UInt idx1 = (i%numFVertices);
        	  UInt idx2 = ((i+1)%numFVertices);
        	  UInt node1 = curShape.faceVertices[aSub][idx1];
        	  UInt node2 = curShape.faceVertices[aSub][idx2];
        	  for(UInt i=0;i<curShape.edgeVertices.GetSize();++i){
        		  if(curShape.edgeVertices[i].Contains(node1) &&
                     curShape.edgeVertices[i].Contains(node2)){
        			  curSurf->extended->edges.Push_back(curE->extended->edges[i]);
        		  }
        	  }
          }
        }else if(curShape.dim==2){
          Edge curEdge = edges_[abs(curE->extended->edges[aSub])-1];
          curSurf->extended->edges.Push_back(abs(curE->extended->edges[aSub]));


          if(curE->extended->edges[aSub]<0){
            //first the corner nodes
            for( Integer i = 1; i >= 0; i-- ) {
              curSurf->connect.Push_back(curE->connect[curShape.edgeVertices[aSub][(UInt)i]-1]);
              curSurf->localCoords.Push_back(curShape.nodeCoords[curShape.edgeVertices[aSub][(UInt)i]-1]);
            }
          }else{
            //first the corner nodes
            for( UInt i = 0; i < 2; i++ ) {
              curSurf->connect.Push_back(curE->connect[curShape.edgeVertices[aSub][i]-1]);
              curSurf->localCoords.Push_back(curShape.nodeCoords[curShape.edgeVertices[aSub][i]-1]);
            }
          }


          //now we check for more nodes on the edge
          //and assume the the corners were the first two entries in the edgeNodes array
          UInt numNodes = curShape.edgeNodes[aSub].GetSize();
          if(curE->extended->edges[aSub]>0){
            for( UInt i = 2; i < numNodes; i++ ) {
              curSurf->connect.Push_back(curE->connect[curShape.edgeNodes[aSub][i]-1]);
              curSurf->localCoords.Push_back(curShape.nodeCoords[curShape.edgeNodes[aSub][i]-1]);
            }
          }else{
            for( UInt i = 0; i < (numNodes-2); i++ ) {
              curSurf->connect.Push_back(curE->connect[curShape.edgeNodes[aSub][numNodes-i-1]-1]);
              curSurf->localCoords.Push_back(curShape.nodeCoords[curShape.edgeNodes[aSub][numNodes-i-1]-1]);
            }
          }
        }else{
          EXCEPTION("surfelem creation for one dimensional shapes not implemented");
        }

        curSurf->ptVolElems[0] = curE;

        // Mapping of normal sign is not needed anymore
//        shared_ptr<ElemShapeMap> es = GetElemShapeMap(curSurf.get(), false);
//        LocPoint lp = Elem::shapes[curSurf->type].midPointCoord;
//        es->CalcNormal( normalUndefSign, lp );
//        es->CalcNormalOutOfVol( normalDefSign, lp, *curSurf->ptVolElems[0] );
//
//        sign = normalUndefSign * normalDefSign;
//
//        if ( sign > 0.0 ) {
//          curSurf->normalSign = 1;
//        } else {
//          curSurf->normalSign = -1;
//        }

        allelems.Push_back(curSurf);

      }//foreach surface
    }//foreach volumen element
    //now we have created the vector with surface elements so we gonna find the neighbors
    //TODO: perhaps we have checked before if the regions are conforming.....
    ComputeSurfElemNeighbors(allelems,interiorSurfElems,exteriorSurfElems,true);

    LOG_DBG3(gridcfs) << std::endl << "==================================================" << std::endl;
    LOG_DBG3(gridcfs) << "Interior surface elements for DG-Methods" << std::endl;
    LOG_DBG3(gridcfs) << "==================================================" << std::endl;
    for(UInt iElem = 0; iElem < interiorSurfElems.GetSize();++iElem){
      LOG_DBG3(gridcfs) << "Added surfElem # " << interiorSurfElems[iElem]->elemNum;
      LOG_DBG3(gridcfs) << " with Edge # " << interiorSurfElems[iElem]->extended->edges[0];
      LOG_DBG3(gridcfs) << " and volElem # " << interiorSurfElems[iElem]->ptVolElems[0]->elemNum;
      LOG_DBG3(gridcfs) << " connect ";
      for(UInt i =0;i<interiorSurfElems[iElem]->connect.GetSize();++i){
        LOG_DBG3(gridcfs) << interiorSurfElems[iElem]->connect[i] << " ";
      }
      LOG_DBG3(gridcfs) << std::endl << "\t Neighbors:" << std::endl;
      StdVector<shared_ptr<NcSurfElem> > neigh = interiorSurfElems[iElem]->neighbors;
      for(UInt i=0;i < neigh.GetSize();++i){
        LOG_DBG3(gridcfs) << "\t " << neigh[i]->elemNum << " ";
      }
      LOG_DBG3(gridcfs) << std::endl;
      LOG_DBG3(gridcfs) << "--------------------------------------------------" << std::endl;
    }

    LOG_DBG3(gridcfs) << "==================================================" << std::endl;
    LOG_DBG3(gridcfs) << "exterior surface elements for DG-Methods" << std::endl;
    LOG_DBG3(gridcfs) << "==================================================" << std::endl;
    for(UInt iElem = 0; iElem < exteriorSurfElems.GetSize();++iElem){
      LOG_DBG3(gridcfs) << "Added surfElem # " << exteriorSurfElems[iElem]->elemNum;
      LOG_DBG3(gridcfs) << " with Edge # " << exteriorSurfElems[iElem]->extended->edges[0];
      LOG_DBG3(gridcfs) << " and volElem # " << exteriorSurfElems[iElem]->ptVolElems[0]->elemNum;
      LOG_DBG3(gridcfs) << " connect ";
      for(UInt i =0;i<exteriorSurfElems[iElem]->connect.GetSize();++i){
        LOG_DBG3(gridcfs) << exteriorSurfElems[iElem]->connect[i] << " ";
      }
      LOG_DBG3(gridcfs) << std::endl << "\t Neighbors:" << std::endl;
      StdVector<shared_ptr<NcSurfElem> > neigh = exteriorSurfElems[iElem]->neighbors;
      for(UInt i=0;i < neigh.GetSize();++i){
        LOG_DBG3(gridcfs) << "\t " << neigh[i]->elemNum << " ";
      }
      LOG_DBG3(gridcfs) << std::endl;
      LOG_DBG3(gridcfs) << "--------------------------------------------------" << std::endl;
    }


    std::cout << "done" << std::endl;
    std::cout.flush();
  }

  void GridCFS::ComputeSurfElemNeighbors(StdVector<shared_ptr<NcSurfElem> > surfElemList,
                                         StdVector<shared_ptr<NcSurfElem> > & interiorSurfElems,
                                         StdVector<shared_ptr<NcSurfElem> > & exteriorSurfElems,
                                         bool conforming){

    // TODO: perhaps we can do some resize commands to avoid the dynamic memory allocation...

    if(conforming){
      std::map<Integer,StdVector<UInt> > EdgeToSurfElemMap;
      std::map<Integer,StdVector<UInt> > FaceToSurfElemMap;
      //first we loop over each element and associate the elementindex to an edge/face number
      for(UInt surfIdx =0;surfIdx < surfElemList.GetSize(); ++surfIdx){
        shared_ptr<NcSurfElem> curSurf = surfElemList[surfIdx];
        ElemShape curShape = Elem::shapes[curSurf->type];
        if(curShape.dim == 2){
          //so we check for faces
          FaceToSurfElemMap[abs(curSurf->extended->faces[0])].Push_back(surfIdx);
        }else if(curShape.dim == 1){
          EdgeToSurfElemMap[abs(curSurf->extended->edges[0])].Push_back(surfIdx);
        }else{
          EXCEPTION("Unsupported dimension");
        }
      }
      //now lets check if we have consitent information
      // only one of the two maps may have elements
      if(EdgeToSurfElemMap.size()>0 && FaceToSurfElemMap.size() > 0)
        EXCEPTION("Supplied an invlid SufelemList: Edge and facemaps have entries. This is not allowed rigth now");

      shared_ptr<NcSurfElem> cS = surfElemList[0];
      ElemShape curShape = Elem::shapes[cS->type];
      //now loop ofer the correct map and assign neighbors to elements
      std::map<Integer,StdVector<UInt> >::iterator fiter;
      std::map<Integer,StdVector<UInt> >::iterator end;
      if(curShape.dim == 2){
        fiter = FaceToSurfElemMap.begin();
        end = FaceToSurfElemMap.end();
      }else if(curShape.dim == 1){
        fiter = EdgeToSurfElemMap.begin();
        end = EdgeToSurfElemMap.end();
      }

      //well we assume, that about 5% of the surfaces are external so reserve them
      interiorSurfElems.Reserve(UInt(std::ceil(surfElemList.GetSize() * 0.95)));
      exteriorSurfElems.Reserve(UInt(std::ceil(surfElemList.GetSize() * 0.05)));

      for(;fiter != end;++fiter){
        //get the index vector
        if(fiter->second.GetSize() > 1){
          // so we have an interior element so we let the guys know about each other
          StdVector<UInt> fIdx = fiter->second;
          for(UInt curE = 0; curE < fIdx.GetSize(); ++curE){
            //and loop again but only assign if cuE != j
            for(UInt j = 0; j < fIdx.GetSize(); ++j){
              if(curE != j){
                surfElemList[fIdx[curE]]->neighbors.Push_back(surfElemList[fIdx[j]]);
              }
            }
            interiorSurfElems.Push_back(surfElemList[fIdx[curE]]);
          }
        }else if (fiter->second.GetSize() == 1){
          // so we have an exterior element we just push it to the vector
          exteriorSurfElems.Push_back(surfElemList[fiter->second[0]]);
        }else if(fiter->second.GetSize() == 0 ){
          EXCEPTION("No association for face to elements. Can this be true????");
        }
      }

    }else{
      EXCEPTION("The nonconforming case is not implemented for the calculation of surf elem neighbors");
    }

  }

  void GridCFS::GetElemCentroid(Vector<Double>& center, UInt eNum, bool isupdated)
  {
    Elem* el = this->orderedElems_[eNum-1];

    if(!el->extended)
    {
      Point cP;
      shared_ptr<ElemShapeMap> esm = this->GetElemShapeMap(el,isupdated);
      esm->CalcBarycenter(cP);
      center = cP.data;
    }
    else
    {
      SetElementBarycenters(el->regionId,isupdated);
      center = el->extended->barycenter.data;
    }
  }

  void GridCFS::CreateGridInformation(ResultHandler* ptRes, std::map<std::string, CoordSystem*>& coordSysMap)
  {
    // This method crates a "dummy" multisequence step, in
    // which some grid-information results are created:
    // - Local directions (xi,eta,zeta) of elements
    // - Jacobian determinant
    // - surface element normals
    //
    // In addition we create for every local coordinate system
    // a results for the local directions.
    
    // Register results
    shared_ptr<BaseResult> sol;
    shared_ptr<EntityList> ent;
    
    // create result descriptions
    Vector<Double> locDir_xi(dim_), locDir_eta(dim_), locDir_zeta(dim_);
    shared_ptr<ResultInfo> locDir1( new ResultInfo );
    shared_ptr<ResultInfo> locDir2( new ResultInfo );
    shared_ptr<ResultInfo> locDir3( new ResultInfo );
    shared_ptr<ResultInfo> jacRes( new ResultInfo );
    shared_ptr<ResultInfo> surfNormal( new ResultInfo );
    shared_ptr<ResultInfo> aspectRatio( new ResultInfo );
    shared_ptr<ResultInfo> area( new ResultInfo );
    shared_ptr<ResultInfo> volume( new ResultInfo );

    StdVector<std::string> dirVec (dim_);
    if( dim_ == 3 ) {
      dirVec = "x", "y", "z";
    } else {
      dirVec = "x", "y";
    }
    // 1) Local directions
    locDir1->resultType = ELEM_LOC_DIR;
    locDir1->resultName = "xi";
    locDir1->definedOn = ResultInfo::ELEMENT;
    locDir1->entryType = ResultInfo::VECTOR;
    locDir1->dofNames = dirVec;
    locDir_xi[0] = 1.0;
    locDir1->unit = "";
    
    locDir2->resultType = ELEM_LOC_DIR;
    locDir2->resultName = "eta";
    locDir2->definedOn = ResultInfo::ELEMENT;
    locDir2->entryType = ResultInfo::VECTOR;
    locDir2->dofNames = dirVec;
    locDir_eta[1] = 1.0;
    locDir2->unit = "";
    
    if (dim_ == 3 )  {
      locDir3->resultType = ELEM_LOC_DIR;
      locDir3->resultName = "zeta";
      locDir3->definedOn = ResultInfo::ELEMENT;
      locDir3->entryType = ResultInfo::VECTOR;
      locDir3->dofNames = dirVec;
      locDir_zeta[2] = 1.0;
      locDir3->unit = "";
    }
    // 2) Jacobian Determinant
    jacRes->resultType = JACOBIAN;
    jacRes->resultName = "Jacobian";
    jacRes->definedOn = ResultInfo::ELEMENT;
    jacRes->entryType = ResultInfo::SCALAR;
    jacRes->dofNames = "";
    jacRes->unit = "";
    
    // 3) Surface element normals
    surfNormal->resultType = ELEM_LOC_DIR;
    surfNormal->resultName = "SurfaceNormal";
    surfNormal->definedOn = ResultInfo::ELEMENT;
    surfNormal->entryType = ResultInfo::VECTOR;
    surfNormal->dofNames = dirVec;
    surfNormal->unit = "";
    
    // 4) Aspect Ratio
    aspectRatio->resultType = ASPECT_RATIO;
    aspectRatio->resultName = "aspectRatio";
    aspectRatio->definedOn = ResultInfo::ELEMENT;
    aspectRatio->entryType = ResultInfo::SCALAR;
    aspectRatio->dofNames = "";
    aspectRatio->unit = "";
    
    // 5) Area
    area->resultType = VOLUME;
    area->resultName = "area";
    area->definedOn = ResultInfo::ELEMENT;
    area->entryType = ResultInfo::SCALAR;
    area->dofNames = "";
    area->unit = "m^2";
    
    // 5) Volume / Area
    volume->resultType = VOLUME;
    volume->resultName = "volume";
    volume->definedOn = ResultInfo::ELEMENT;
    volume->entryType = ResultInfo::SCALAR;
    volume->dofNames = "";
    volume->unit = "m^3";
    
    // === create results for all coordinate systems available ===
    std::map<std::string, StdVector<shared_ptr<ResultInfo> > > coordDirs;
    std::map<std::string, CoordSystem*>::const_iterator coordIt;
    
    // Empty functor
    shared_ptr<ResultFunctor> fnc;
    
    // loop over all coordinate systems
    for( coordIt = coordSysMap.begin(); 
         coordIt != coordSysMap.end(); ++coordIt ) {
       
      std::string id = coordIt->first;
      CoordSystem * actCosy = coordIt->second;
      
      // Create ResultInfo objects for every local direction
      StdVector<shared_ptr<ResultInfo> > dirInfo(dim_);
      for( UInt iDim = 0; iDim < dim_; ++iDim ) {
        shared_ptr<ResultInfo> locDir( new ResultInfo );
        locDir->resultType = ELEM_LOC_DIR;
        locDir->resultName = "CooSy-"+id+"-"+actCosy->GetDofName(iDim+1);
        locDir->definedOn = ResultInfo::ELEMENT;
        locDir->entryType = ResultInfo::VECTOR;
        locDir->dofNames = dirVec;
        locDir->unit = "";
        
        dirInfo[iDim] = locDir;
      }
      coordDirs[id] = dirInfo;
    }

    
    // loop over all volume regions
    StdVector<std::string> outDest;
    outDest.Push_back("");
    StdVector<shared_ptr<BaseResult> > resultList;
    std::map<std::string, StdVector<shared_ptr<BaseResult> > > coordResultList;
    for ( UInt i = 0, numRegions = volRegionIds_.GetSize();
        i < numRegions; i++) {
      
      // get elements
      ent = GetEntityList( EntityList::ELEM_LIST, region_.ToString(volRegionIds_[i]) );
      
    
      
      // create result objects
      sol = shared_ptr<BaseResult> (new Result<Double>());
      sol->SetEntityList(ent);
      sol->SetResultInfo(jacRes);
      ptRes->RegisterResult( sol,fnc,0,0,1,1,outDest,"",true,false);
      resultList.Push_back(sol);
      
      sol = shared_ptr<BaseResult> (new Result<Double>());
      sol->SetEntityList(ent);
      sol->SetResultInfo(aspectRatio);
      ptRes->RegisterResult( sol,fnc,0,0,1,1,outDest,"",true,false);
      resultList.Push_back(sol);

      sol = shared_ptr<BaseResult> (new Result<Double>());
      sol->SetEntityList(ent);
      sol->SetResultInfo(locDir1);
      ptRes->RegisterResult( sol,fnc,0,0,1,1,outDest,"",true,false);
      resultList.Push_back(sol);

      sol = shared_ptr<BaseResult> (new Result<Double>());
      sol->SetEntityList(ent);
      sol->SetResultInfo(locDir2);
      ptRes->RegisterResult( sol,fnc,0,0,1,1,outDest,"",true,false);
      resultList.Push_back(sol);


      if( dim_ == 3 ) {
        sol = shared_ptr<BaseResult> (new Result<Double>());
        sol->SetEntityList(ent);
        sol->SetResultInfo(locDir3);
        ptRes->RegisterResult( sol,fnc,0, 0,1,1,outDest,"",true,false);
        resultList.Push_back(sol);
        
        sol = shared_ptr<BaseResult> (new Result<Double>());
        sol->SetEntityList(ent);
        sol->SetResultInfo(volume);
        ptRes->RegisterResult(sol, fnc, 0, 0, 1, 1, outDest, "", true, false);
        resultList.Push_back(sol);
      }
      else {
        sol = shared_ptr<BaseResult> (new Result<Double>());
        sol->SetEntityList(ent);
        sol->SetResultInfo(area);
        ptRes->RegisterResult(sol, fnc, 0, 0, 1, 1, outDest, "", true, false);
        resultList.Push_back(sol);
      }
      
      // loop over all coordinate systems
      std::map<std::string, StdVector<shared_ptr<ResultInfo> > >
      ::iterator coordResIt;
      for( coordResIt = coordDirs.begin(); 
           coordResIt != coordDirs.end();
           ++coordResIt )   {
        
      // loop over all local directions
        StdVector<shared_ptr<ResultInfo> >& dirResults = coordResIt->second;
        for( UInt iDir = 0; iDir < dirResults.GetSize(); ++iDir ) {
          sol = shared_ptr<BaseResult> (new Result<Double>());
          sol->SetEntityList(ent);
          sol->SetResultInfo(dirResults[iDir]);
          ptRes->RegisterResult( sol,fnc,0, 0,1,1,outDest,"",true,false);
          coordResultList[coordResIt->first].Push_back(sol);
        }
      }
    }
    
    // loop over all surface regions
    StdVector<shared_ptr<BaseResult> > surfResultList;
    for ( UInt i = 0, numSurfaces = surfRegionIds_.GetSize();
        i < numSurfaces; i++) {

      // get elements
      std::string surfRegName = region_.ToString(surfRegionIds_[i]);
      ent = GetEntityList( EntityList::SURF_ELEM_LIST, surfRegName);

      // check whether the surface region doesn't belong to an NC-interface:
      // we have to make sure that the elements of the former have adjacent volume elements
      // required for the calculation of surface normals
      if (nciNameMap_.find(surfRegName) == nciNameMap_.end())
      {
        // create result objects
        sol = shared_ptr<BaseResult> (new Result<Double>());
        sol->SetEntityList(ent);
        sol->SetResultInfo(surfNormal);
        ptRes->RegisterResult( sol,fnc,0, 0,1,1,outDest,"",true,false);
        surfResultList.Push_back(sol);
      }
      
      if ( dim_ == 3 ) {
        sol = shared_ptr<BaseResult> (new Result<Double>());
        sol->SetEntityList(ent);
        sol->SetResultInfo(area);
        ptRes->RegisterResult(sol, fnc, 0, 0, 1, 1, outDest, "", true, false);
        resultList.Push_back(sol);
      }
    }
    
    // begin writing of results
    ptRes->BeginMultiSequenceStep( 0, BasePDE::STATIC, 1);
    ptRes->BeginStep(0,0);
    
    Matrix<Double> actCoord, jac, jacInv, jacInvT;
    Vector<Double> midPoint, globVec, actLocDir;
    
    // loop over all volume region results
    for( UInt i = 0; i < resultList.GetSize(); ++i ) {
      
      // fetch result object
      Result<Double> &  actSol = 
          dynamic_cast<Result<Double>&>(*resultList[i]);
      EntityIterator it = actSol.GetEntityList()->GetIterator();
      Vector<Double> & actVal = actSol.GetVector();

      // check, which result is required
      switch (actSol.GetResultInfo()->resultType) {
      case ELEM_LOC_DIR:
        if (actSol.GetResultInfo()->resultName == "xi") {
          actLocDir = locDir_xi;
        } else if (actSol.GetResultInfo()->resultName == "eta") {
          actLocDir = locDir_eta;
        } else if (actSol.GetResultInfo()->resultName == "zeta") {
          actLocDir = locDir_zeta;
        }
        actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );

        // loop over elements
        for ( it.Begin(); !it.IsEnd(); it++ ) {
          const Elem * ptElem = it.GetElem();
          shared_ptr<ElemShapeMap>  esm = GetElemShapeMap(ptElem);
          Vector<Double> midPoint = 
              Elem::shapes[ptElem->type].midPointCoord;
          esm->CalcJ( jac, midPoint );
          // calculate for every element jacobian and map "local" direction vector
          if (jac.GetNumCols() == 1) {
            globVec.Resize(dim_); // do nothing for the case of surface/line elements in 3D
          }
          else {
            jac.Invert( jacInv );
            globVec = Transpose(jacInv) * actLocDir;
            globVec /= globVec.NormL2();
          }
          for(UInt iDim = 0; iDim < dim_; iDim++ ) {
            actVal[it.GetPos()*dim_ + iDim] = globVec[iDim];
          }
        }
        break;

      case ASPECT_RATIO:
        actVal.Resize( actSol.GetEntityList()->GetSize() );
        // loop over elements
        for ( it.Begin(); !it.IsEnd(); it++ ) {

          const Elem * ptElem = it.GetElem();
          shared_ptr<ElemShapeMap>  esm = GetElemShapeMap(ptElem);
          Double minEdge, maxEdge;
          esm->GetMaxMinEdgeLength(maxEdge,minEdge);
          actVal[it.GetPos()] = maxEdge / minEdge;
        } // loop over elements
        break;
        
      case JACOBIAN:
        actVal.Resize( actSol.GetEntityList()->GetSize() );
        // loop over elements
        for ( it.Begin(); !it.IsEnd(); it++ ) {

          const Elem * ptElem = it.GetElem();
          shared_ptr<ElemShapeMap>  esm = GetElemShapeMap(ptElem);
          Vector<Double> midPoint = 
              Elem::shapes[ptElem->type].midPointCoord;
          esm->CalcJ( jac, midPoint );
          actVal[it.GetPos()] = esm->CalcJDet(jac,midPoint);
        } // loop over elements
        break;
        
      case VOLUME:
        actVal.Resize( actSol.GetEntityList()->GetSize() );
        // loop over elements
        for ( it.Begin(); !it.IsEnd(); it++ ) {

          const Elem * ptElem = it.GetElem();
          shared_ptr<ElemShapeMap>  esm = GetElemShapeMap(ptElem);
          actVal[it.GetPos()] = esm->CalcVolume();
        } // loop over elements
        break;
        
      default:
        WARN("Grid info '" << actSol.GetResultInfo()->resultName
             << "is enabled, but calculation is not implemented.");
      }
      ptRes->UpdateResult(resultList[i]);
    } // loop over volume results
    
    // loop over all coordinate systems
    for( coordIt = coordSysMap.begin(); coordIt != coordSysMap.end(); ++coordIt){

      std::string id = coordIt->first;
      CoordSystem * actCosy = coordIt->second;

      // obtain result vector
      StdVector<shared_ptr<BaseResult> >  & resList = coordResultList[id];

      // loop over all results 
      // (= components of the coordinate system x regions) 
      for( UInt iPos = 0; iPos < resList.GetSize(); ++iPos ) {

        // fetch result object
        Result<Double> &  actSol = 
            dynamic_cast<Result<Double>&>(*resList[iPos]);
        EntityIterator it = actSol.GetEntityList()->GetIterator();
        Vector<Double> & actVal = actSol.GetVector();
        actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
        
        // determine coordinate component
        UInt actComp = iPos % dim_;

        // loop over all elements
        Vector<Double> coordDir, locDir(dim_);
        locDir.Init();
        locDir[actComp] = 1.0;
        for ( it.Begin(); !it.IsEnd(); it++ ) {
          Vector<Double> midPoint;
          shared_ptr<ElemShapeMap> esm = GetElemShapeMap(it.GetElem());
          esm->GetGlobMidPoint(midPoint);
          actCosy->Local2GlobalVector(coordDir, locDir, midPoint);
          for( UInt i = 0; i < dim_; ++i ) {
            actVal[it.GetPos()*dim_ + i] = coordDir[i];
          } // components of the local coordinate direction
        } // loop over elements
        ptRes->UpdateResult(resList[iPos]);
      } // loop over components of coordinate system
    } // loop over all coordinate systems
    
    Vector<Double> normal;
    // loop over all surface region results
    for( UInt i = 0; i < surfResultList.GetSize(); ++i ) {
      // fetch result object
      Result<Double> &  actSol = 
          dynamic_cast<Result<Double>&>(*surfResultList[i]);
      EntityIterator it = actSol.GetEntityList()->GetIterator();
      Vector<Double> & actVal = actSol.GetVector();
      
      switch (actSol.GetResultInfo()->resultType) {
      case ELEM_LOC_DIR:
        actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );

        // loop over surface elements
        for ( it.Begin(); !it.IsEnd(); it++ ) {
          const SurfElem * ptElem = it.GetSurfElem();
          shared_ptr<ElemShapeMap> esm = GetElemShapeMap(ptElem);
          Vector<Double> midPoint = 
              Elem::shapes[it.GetElem()->type].midPointCoord;
          esm->CalcNormal( normal, midPoint );


          for(UInt iDim = 0; iDim < dim_; iDim++ ) {
            actVal[it.GetPos()*dim_ + iDim] = normal[iDim];
          }
        } // loop over elements
        break;
        
      case VOLUME:
        actVal.Resize( actSol.GetEntityList()->GetSize() );

        // loop over surface elements
        for ( it.Begin(); !it.IsEnd(); it++ ) {
          const SurfElem * ptElem = it.GetSurfElem();
          shared_ptr<ElemShapeMap> esm = GetElemShapeMap(ptElem);
          actVal[it.GetPos()] = esm->CalcVolume();
        } // loop over elements
        break;
        
      default:
        WARN("Grid info '" << actSol.GetResultInfo()->resultName
             << "is enabled, but calculation is not implemented.");
      }
      
      ptRes->UpdateResult(surfResultList[i]);
    } // loop over surface results

    ptRes->FinishStep();
    ptRes->FinishMultiSequenceStep();
  }

  bool GridCFS::CheckForRegularRegion(RegionIdType reg)
  {
    LOG_DBG(gridcfs) << "CFRR::CheckForRegularRegion(" << reg << ")...";
    
    // determine grid parameter
    RegionData& rd = regionData[reg];
    assert(rd.id == reg);
    assert(rd.type != NOT_SET);

    StdVector<Elem*>& elems = rd.type == VOLUME_REGION ? 
        volElems_[rd.type_idx] : surfElems_[rd.type_idx];
    assert(rd.type == VOLUME_REGION ? 
        volRegionIds_[rd.type_idx] == rd.id : 
        surfRegionIds_[rd.type_idx] == rd.id);

    LOG_DBG2(gridcfs) << "CFRR volume=" << (rd.type == VOLUME_REGION) 
        << " #elems=" << elems.GetSize();
    if(elems.GetSize() == 0) return false;

    // we do not know of the order of the nodes within the element,
    // hence we have to find the diameter vector and compare it.

    shared_ptr<ElemShapeMap> esm = (this->GetElemShapeMap(elems[0]));
    Elem::ShapeType refShape = Elem::GetShapeType(elems[0]->type);

    // only line, quad and hexa can be regular. If you have a
    // regular tetra mesh add the code!
    if(dim_ == 2 && rd.type == SURFACE_REGION && (refShape != Elem::ST_LINE))
      return false;

    if(dim_ == 2 && rd.type == VOLUME_REGION && (refShape != Elem::ST_QUAD))
      return false;

    if(dim_ == 3 && rd.type == SURFACE_REGION && (refShape != Elem::ST_QUAD))
      return false;

    if(dim_ == 3 && rd.type == VOLUME_REGION && (refShape != Elem::ST_HEXA))
      return false;

    Matrix<Double>  coords;
    Vector<Double> diameter; // diameter of reference (first) element

    esm->CalcDiameter(diameter);
    Double dist_tol = diameter.NormL2() * 1e-6;
    LOG_DBG2(gridcfs) << "CFRR diameter=" << diameter.ToString() << " tol=" << dist_tol;

    Vector<Double> test_diameter; // diameter of current element;
    Vector<Double> diff; // difference vector
    for(unsigned int i = 1, in = elems.GetSize(); i < in; i++)
    {
      shared_ptr<ElemShapeMap> testEsm = (this->GetElemShapeMap(elems[i]));
      Elem::ShapeType testShape = Elem::GetShapeType(elems[i]->type);
      if(testShape != refShape) return false;

      testEsm->CalcDiameter(test_diameter);

      diff = diameter-test_diameter;
      LOG_DBG3(gridcfs) << "CFRR test=" << elems[i]->elemNum
                        << " dist=" << diff.NormL2()
                        << " test_diameter=" << test_diameter.ToString();
      if( diff.NormL2() > dist_tol) return false;
    }

    LOG_DBG(gridcfs) << "CFRR: is regular";
    return true; // seems to be regular
  }

  // ======================================================
  // GENERAL GRID INFORMATION
  // ======================================================


  UInt GridCFS::GetNumElemOfType( Elem::FEType type ) {
    return numElemTypes_[type];
  }
  
  UInt GridCFS::GetNumElemOfDim( UInt dim ) {

    UInt numElems = 0;
    switch(dim) {
      case 1:
        numElems = numElemTypes_[Elem::ET_LINE2]
                 + numElemTypes_[Elem::ET_LINE3];
        break;
      case 2:
        numElems = numElemTypes_[Elem::ET_TRIA3]
                 + numElemTypes_[Elem::ET_TRIA6]
                 + numElemTypes_[Elem::ET_QUAD4]
                 + numElemTypes_[Elem::ET_QUAD8]
                 + numElemTypes_[Elem::ET_QUAD9];
        break;
      case 3:
        numElems = numElemTypes_[Elem::ET_TET4]
                 + numElemTypes_[Elem::ET_TET10]
                 + numElemTypes_[Elem::ET_HEXA8]
                 + numElemTypes_[Elem::ET_HEXA20]
                 + numElemTypes_[Elem::ET_HEXA27]
                 + numElemTypes_[Elem::ET_PYRA5]
                 + numElemTypes_[Elem::ET_PYRA13]
                 + numElemTypes_[Elem::ET_WEDGE6]
                 + numElemTypes_[Elem::ET_WEDGE15];
                                 
                                                      
        break;
      default:
        EXCEPTION("Grid can only have dimension up to 3!");
        break;
    }
    return numElems;
  }

  void GridCFS::AddNodes(const UInt numNodes)
  {
#pragma omp critical (GridCFS)
{
    unsigned int old_size = this->numNodes_;
    coords_.Resize(this->numNodes_ + numNodes);
    for(unsigned int  i = old_size; i < coords_.GetSize(); ++i )
      coords_[i].Resize(dim_);

    numNodes_ += numNodes;
}
  }


  void GridCFS::SetNodeCoordinate(const UInt inode, const Vector<Double> & rfPoint)
  {
    if ( inode > numNodes_ ) {
      EXCEPTION( "GridCFS: There are only " << numNodes_
                 << " nodes in the grid. You wanted to set coordinates for "
                 << "node number " << inode );
    }

    if ( (dim_ == 2) && rfPoint.GetSize() > 2 && (abs(rfPoint[2]) > 1.0E-15) ) {
      EXCEPTION( "GridCFS: Dimension of grid is 2D. "
                  << "But you wanted to set the 3D coordinate " << "("
                  << rfPoint[0] << ", " << rfPoint[1] << ", " << rfPoint[2]
                  << ") for node number " << inode);
    }

    UInt idx = inode-1;
    for( UInt i = 0; i < dim_; ++i ) {
      coords_[idx][i] = rfPoint[i];  
    }
    
  }

  UInt GridCFS::GetNumNodes(RegionIdType reg_id) const
  {
    if(reg_id == ALL_REGIONS)
      return numNodes_;

    Integer index = 0;

    // look in volume regions
    index = volRegionIds_.Find(reg_id);
    if(index != -1)
      return numVolElemNodes_[index];

    // look in surface regions
    index = surfRegionIds_.Find(reg_id);
    if(index != -1)
      return numSurfElemNodes_[index];

    EXCEPTION("The region with id '" << reg_id << "' is unknown");
    return 0;
  }


  UInt GridCFS::GetNumNodes( const std::string & nodesName ) const {

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

  UInt GridCFS::GetNumSurfElems(RegionIdType reg_id) const
  {
    if(reg_id != ALL_REGIONS)
    {
      assert(regionData[reg_id].id == reg_id);
      assert(regionData[reg_id].type == SURFACE_REGION);
      return surfElems_[regionData[reg_id].type_idx].GetSize();
    }

    UInt numSurfElems = 0;

    for (UInt i=0; i < surfElems_.GetSize(); i++)
      numSurfElems += surfElems_[i].GetSize();

    return numSurfElems;
  }


  UInt GridCFS::GetNumVolElems(RegionIdType reg_id) const
  {
    if(reg_id != ALL_REGIONS)
    {
      assert(regionData[reg_id].id == reg_id);
      assert(regionData[reg_id].type == VOLUME_REGION);
      return volElems_[regionData[reg_id].type_idx].GetSize();
    }

    UInt numVolElems = 0;

    for (UInt i=0; i < volElems_.GetSize(); i++)
      numVolElems += volElems_[i].GetSize();

    return numVolElems;
  }

  UInt GridCFS::GetNumElems(RegionIdType reg_id) const
  {
    if(reg_id == ALL_REGIONS)
      return numElems_;

    if(regionData[reg_id].type == VOLUME_REGION)
      return GetNumVolElems(reg_id);
    else
      return GetNumSurfElems(reg_id);
  }

  UInt GridCFS::GetNumElems(const StdVector<RegionIdType> & regions) const
  {
    UInt numElems = 0;

    for (UInt i=0; i < regions.GetSize(); i++)
    {
      assert(regionData[regions[i]].id == regions[i]);
      if(regionData[regions[i]].type == SURFACE_REGION)
        numElems += GetNumSurfElems(regions[i]);
      else
        numElems += GetNumVolElems(regions[i]);
    }
    return numElems;
  }

  void GridCFS::AddNamedNodes(const std::string& name, StdVector<unsigned int> & nodeNums)
  {
    // entity names need to be unique, even across different entity types
	  if(nameTypeMap_.find( name) != nameTypeMap_.end())
	  {
      // either the name is used for another entity, or we extend the current type
	    int pos = namedNodeNames_.Find(name);
	    if(pos == -1)
	      EXCEPTION("'" + name + "' is not valid for named nodes, another entity type already uses the name");
	    StdVector<unsigned int> &nN = namedNodes_[pos];
      // now add the new nodeNums (if they are not already in the vector)
      for(unsigned int& n : nodeNums)
        if(!(std::find(nN.begin(), nN.end(), n) != nN.end()))
          nN.Push_back(n);
	  }
	  else
	  {
      namedNodeNames_.Push_back(name);
      namedNodes_.Push_back(nodeNums);
      nameTypeMap_[name] = EntityList::NAMED_NODES;
      entityDim_[name] = 0;
	  }
  }

  void GridCFS::AddNamedElems(const std::string& name, StdVector<UInt> & elemNums)
  {
    // Check if entities with given name exist already - might be of other type.
    // currently no extension of existing named nodes implemented
    if( nameTypeMap_.find( name) != nameTypeMap_.end())
      EXCEPTION( "Entities with name " << name << " are already defined" );

    namedElemNames_.Push_back(name);
    
    // Perform check, that all elements in this list have the same dimension
    UInt elemDim = Elem::shapes[orderedElems_[elemNums[0]-1]->type].dim;
    UInt numElems = elemNums.GetSize();
    for( UInt iElem = 1; iElem < numElems; ++iElem ) {
      if( Elem::shapes[orderedElems_[elemNums[iElem]-1]->type].dim != elemDim ) {
        EXCEPTION( "Element list '" << name << 
                   "' contains elements of different dimensions!")
      }
    }
    
    namedElems_.Push_back(elemNums);
    nameTypeMap_[name] = EntityList::NAMED_ELEMS;
    entityDim_[name] = elemDim;

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


  void GridCFS::GetAdjacentSurfElem( const UInt volElemNum, StdVector<Elem *> & surfEl, const RegionIdType reg_id) {
    Integer index = 0;
    surfEl.Clear();
    index = surfRegionIds_.Find(reg_id);
    if ( index != -1 ) {
      UInt numElems = surfElems_[index].GetSize();
      SurfElem * ptSurfElem;
      for( UInt iElem = 0; iElem <  numElems; ++iElem ) {
        ptSurfElem = dynamic_cast<SurfElem*>(surfElems_[index][iElem]);
        if (ptSurfElem->ptVolElems[0] != nullptr) {
	  if (ptSurfElem->ptVolElems[0]->elemNum == volElemNum)
	    surfEl.Push_back( (Elem *) surfElems_[index][iElem]);
        }
        if (ptSurfElem->ptVolElems[1] != nullptr) {
	  if (ptSurfElem->ptVolElems[1]->elemNum == volElemNum)
	    surfEl.Push_back( (Elem *) surfElems_[index][iElem]);
        }
      }
    } else {
        EXCEPTION( "GridCFS: The surface region with id '" << reg_id
                   << "' was not found in the grid!" ); 

    }
  }

  void GridCFS::GetListOfVolumeRegions( const RegionIdType reg_id, StdVector<RegionIdType> &volRegIds ) {

    // check if region id is a volume anyways, then just return that
    Integer index = 0;
    volRegIds.Clear();

    // look in volume regions
    index = volRegionIds_.Find(reg_id);
    if ( index != -1 ) {
     volRegIds.Resize(1);
     volRegIds[0] = reg_id;
    
    } else {
      // look in surface regions
      index = surfRegionIds_.Find(reg_id);
      if ( index != -1 ) {
        UInt numElems = surfElems_[index].GetSize();
	SurfElem * ptSurfElem;
	Integer iFound;
        for( UInt iElem = 0; iElem <  numElems; ++iElem ) {
          ptSurfElem = dynamic_cast<SurfElem*>(surfElems_[index][iElem]);
	  if (ptSurfElem->ptVolElems[0] != nullptr) {
	    iFound = volRegIds.Find(ptSurfElem->ptVolElems[0]->regionId);
	    if (iFound == -1) // not found
	      volRegIds.Push_back(ptSurfElem->ptVolElems[0]->regionId);

	  } else if(ptSurfElem->ptVolElems[1] != nullptr) {
	    WARN("not implemented");
	  } else {
            EXCEPTION( "GridCFS: The surface region with id '" << reg_id
                   << "' doesn't have a volume attached!" ); 
	  }
	}
      // loop over all elements and save their volRegionIds 
      } else {
        EXCEPTION( "GridCFS: The region with id '" << reg_id
                   << "' was not found in the grid!" ); 

      }
    }


  }

  // ======================================================
  // NODE ACCESS FUNCTIONS
  // ======================================================

  void GridCFS::GetNodesByName( StdVector<UInt> & nodeList,
                                const std::string & name ) {

    // Check if entities with given name exists at all
    if( nameTypeMap_.find( name) == nameTypeMap_.end() ) {
      EXCEPTION( "Entities with name " << name
                 << " not found in nameTypeMap_" );
    }

    // check, which entity type the name belongs to
    EntityList::DefineType defType = nameTypeMap_[name];
    Integer index = -1;

    switch( defType ) {

    case EntityList::REGION:
      GetNodesByRegion( nodeList, region_.Parse(name) );
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
      break;
    }
  }

  void GridCFS::GetNodesByRegion( StdVector<UInt> & nodeList,
                                  const RegionIdType regionId ) {
    std::vector<bool> usedNode(numNodes_ + 1, false);
    nodeList.Clear();
    
    Integer index = 0;
  
    // look in volume regions
    index = volRegionIds_.Find(regionId);
    if ( index != -1 ) {
      UInt numElems = volElems_[index].GetSize();
      for( UInt iElem = 0; iElem <  numElems; ++iElem ) {
        const Elem * el = volElems_[index][iElem];
        const StdVector<UInt> & connect = el->connect;
        UInt numNodes = connect.GetSize();
        for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
          usedNode[connect[iNode]] = true;
        }
      }
      
    } else {
      // look in surface regions
      index = surfRegionIds_.Find(regionId);
      if ( index != -1 ) {
        UInt numElems = surfElems_[index].GetSize();
        for( UInt iElem = 0; iElem <  numElems; ++iElem ) {
          const Elem * el = surfElems_[index][iElem];
          const StdVector<UInt> & connect = el->connect;
          UInt numNodes = connect.GetSize();
          for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
            usedNode[connect[iNode]] = true;
          }
        }
      } else {
        EXCEPTION( "GridCFS: The region with id '" << regionId
                   << "' was not found in the grid!" ); 
      }
    }
    
    UInt count = 0;
#pragma omp parallel for reduction(+:count)
    for (Integer i = 1; i <= (Integer) numNodes_;++i) {
      if (usedNode[i]) {
        count++;
      }
    }
    nodeList.Reserve(count);
    for (UInt i = 1; i <= numNodes_;++i) {
      if (usedNode[i]) {
        nodeList.Push_back(i);
      }
    }
    
  }



  void GridCFS::GetNodeCoordinate( Vector<Double> & rfPoint,
                                   const UInt inode,
                                   bool updated ) const {

    if ( inode > numNodes_ ) {
      EXCEPTION( "GridCFS: There are only " << numNodes_
                 << " nodes in the grid. You requested coordinates for "
                 << "node number " << inode <<". Go check your mesh file!" );
    }

    UInt idx = inode-1;
    rfPoint = coords_[idx];

    if (updated && deltCoords_.GetSize() > 0) {
      rfPoint += deltCoords_[idx];
    }
  }
  

  void GridCFS::GetNodeCoordinates(StdVector<Vector<Double>>& nodeCoords,
                                   const StdVector<UInt>& nodeList,
                                   bool updated ) const {

    nodeCoords.Resize(nodeList.GetSize());
    // check if nodes are available
    for (UInt i = 0; i < nodeList.GetSize(); ++i) {
      if (nodeList[i] > numNodes_ ) {
        EXCEPTION( "GridCFS: There are only " << numNodes_
                   << " nodes in the grid. You requested coordinates for "
                   << "node number " << nodeList[i] <<". Go check your mesh file!" );
      }
      nodeCoords[i] = coords_[nodeList[i]-1]; //node indices in nodeList are 1 based
      if (updated && deltCoords_.GetSize() > 0) {
        nodeCoords[i] += deltCoords_[nodeList[i]-1];
      }
    }
  }

  void GridCFS::GetNodeCoordinate3D( Vector<Double> & rfPoint,
                                   const UInt inode,
                                   bool updated ) const {

    if ( inode > numNodes_ ) {
      EXCEPTION( "GridCFS: There are only " << numNodes_
                 << " nodes in the grid. You requested coordinates for "
                 << "node number " << inode <<". Go check your mesh file!" );
    }

    rfPoint.Resize(3);
    rfPoint.Init();
    UInt idx = inode-1;
    for( UInt i = 0; i < dim_; ++i )
      rfPoint[i] = coords_[idx][i];

    if (updated && deltCoords_.GetSize() > 0) {
      for( UInt i = 0; i < dim_; ++i )
        rfPoint[i] += deltCoords_[idx][i];
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

      if(buildExtendedElemInfo_){
        orderedElems_[idx]->extended = new ExtendedElementInfo();
      }

      orderedElems_[idx]->elemNum = idx+1;
    }

    numElems_ = idx;
  }

  Elem* GridCFS::SearchFistRegionElement(RegionIdType reg) const
  {
    for(unsigned int i = 0, n = orderedElems_.GetSize(); i < n; i++)
      if(orderedElems_[i]->regionId == reg)
        return orderedElems_[i];

    EXCEPTION("no elements found for region id " << reg);
  }

  // Reserve memory for a number of elements without adding them
  void GridCFS::ReserveElems(UInt nElems) {
    orderedElems_.Reserve(orderedElems_.GetCapacity() + nElems);
  }


  void GridCFS::SetElemData(UInt ielem,
                            Elem::FEType type,
                            RegionIdType region,
                            const UInt* connect)
  {
    assert(type != Elem::ET_UNDEF);
    
    if ( ielem > orderedElems_.GetSize() ) {
      EXCEPTION("Element numbers are not contiguous: Element number " << ielem
                << " is out of allowed range. Either you forgot to write out "
                << "some elements or you need to compress the element "
                << "numbering. Please check your mesh!");
    }

    UInt idx=ielem-1;
    Elem* el = orderedElems_[idx];
    el->type = type;
    el->regionId = region;
    UInt d = 2;
    UInt numNodes = Elem::shapes[type].numNodes;

    // set isQuadratic information
    isQuadratic_ |= (Elem::shapes[type].order == 2); 
    
    numElemTypes_[type]++;

    switch(type)
    {
    case Elem::ET_LINE2:
      break;
    case Elem::ET_QUAD4:
      break;
    case Elem::ET_HEXA8:
      d=3;
      break;
    default:
      break;
    }

    if((dim_ == 2) && (d == 3))
    {
      EXCEPTION( "GridCFS: Cannot add 3D element type "
                 << Elem::feType.ToString(type)
                 << " to 2D grid." );
    }

    el->connect.Resize(numNodes);
    bool hasError = false;
    for( UInt i = 0; i < numNodes; i++ ) {
      if( connect[i] > numNodes_ )
        hasError = true;
      el->connect[i] = connect[i];
    }
    if( hasError) {
      EXCEPTION( "Element #" << el->elemNum << " with connectivity ("
                 << el->connect.Serialize() << ") has node number(s) larger "
                 << "than number of nodes in grid (" << numNodes_ << ")" );
    }

    // add correct dimension of element to entityDim_
    if( region != NO_REGION_ID) {
      std::string regionName = region_.ToString(region);
      std::map<std::string, UInt>::iterator it = entityDim_.find(regionName);
      if( it != entityDim_.end() ) {
        if( it->second != Elem::shapes[type].dim ) {
          EXCEPTION( "Region '" << regionName 
                     << "' contains elements of different dimensions!");
        }
      } else {
        entityDim_[regionName] = Elem::shapes[type].dim;
      }
    } 
  }

  void GridCFS::GetElemData(const UInt ielem,
                            Elem::FEType & type,
                            RegionIdType & region,
                            UInt* connect) const
  {
 #ifndef NDEBUG
    if ( ielem > numElems_ ) {
      EXCEPTION( "GridCFS: There are only " << numElems_
                 << " elements in the grid! You requested element number "
                 << ielem << ". Go check your mesh file!" );
    }
    if ( orderedElems_[ielem-1] == nullptr ) {
      EXCEPTION( "Element with Nr. " << ielem << " is not contained in mesh!" );
    }
 #endif

    UInt numNodes;

    type = orderedElems_[ielem-1]->type;
    region = orderedElems_[ielem-1]->regionId;
    numNodes = Elem::shapes[type].numNodes;
    memcpy(connect, &orderedElems_[ielem-1]->connect[0], numNodes*sizeof(UInt));

  }

  void GridCFS::SetElemRegion(UInt ielem, RegionIdType region)
  {
 #ifndef NDEBUG
    if ( ielem > numElems_ ) {
      EXCEPTION( "GridCFS: There are only " << numElems_
                 << " elements in the grid! You requested element number "
                 << ielem << ". Go check your mesh file!" );
    }
    if ( orderedElems_[ielem-1] == nullptr ) {
      EXCEPTION( "Element with Nr. " << ielem << " is not contained in mesh!" );

    }
 #endif

    orderedElems_[ielem-1]->regionId = region;
  }

  void GridCFS::GetElemRegion(UInt ielem, RegionIdType & region)
   {
  #ifndef NDEBUG
     if ( ielem > numElems_ ) {
       EXCEPTION( "GridCFS: There are only " << numElems_
                  << " elements in the grid! You requested element number "
                  << ielem << ". Go check your mesh file!" );
     }
     if ( orderedElems_[ielem-1] == nullptr ) {
       EXCEPTION( "Element with Nr. " << ielem << " is not contained in mesh!" );

     }
  #endif

     region = orderedElems_[ielem-1]->regionId;
   }

  RegionIdType GridCFS::GetElemRegion(UInt ielem)
   {
  #ifndef NDEBUG
     if ( ielem > numElems_ ) {
       EXCEPTION( "GridCFS: There are only " << numElems_
                  << " elements in the grid! You requested element number "
                  << ielem << ". Go check your mesh file!" );
     }
     if ( orderedElems_[ielem-1] == nullptr ) {
       EXCEPTION( "Element with Nr. " << ielem << " is not contained in mesh!" );

     }
  #endif

     return orderedElems_[ielem-1]->regionId;
   }

  void GridCFS::FindElementNeighorhood()
  {
    // TODO: This is from the legacy code -> replace with ShapeMap concept!

    // check if already filled.
    if(orderedElems_[0]->extended->neighborhood != NULL) return;

    // this is expensive but still O(n)
    // for all nodes we store the elements where they participate
    StdVector<StdVector<Elem*> > nodes;
    // we assume that the all nodes are consecutive and start with 0 or 1!
    nodes.Resize(numNodes_ + 1);
    for(unsigned int e = 0; e < numElems_; e++)
    {
      Elem* elem = orderedElems_[e];
      // add this elem to every node
      for(unsigned int n = 0, nn = elem->connect.GetSize(); n < nn; n++)
        nodes[elem->connect[n]].Push_back(elem);
    }

    LOG_DBG3(gridcfs) << "FEN: nodes=" << nodes.ToString();

    // this is a temporary neighborhood, copied to the elements when the size is known.
    StdVector<std::pair<Elem*, int> > neighborhood;
    // this is the temporary pair, as we cannot construct it within StdVector::Push_back()
    std::pair<Elem*, int> pair(nullptr, 0);

    // now add to all elements the neighborhood
    for(unsigned int e = 0; e < numElems_; e++)
    {
      Elem* elem = orderedElems_[e];
      neighborhood.Resize(0); // TODO -> keep capacity

      // process the elements which are connected to our nodes
      for(unsigned int n = 0, nn = elem->connect.GetSize(); n < nn; n++)
      {
        const unsigned int curr_n = elem->connect[n];
        StdVector<Elem*>& list = nodes[curr_n];
        // check for duplicity of the elements
        for(unsigned int c = 0, nc = list.GetSize(); c < nc; c++)
        {
          Elem* cand = list[c];
          if(cand == elem) continue; // we are not a neighbour of ourself
          assert(cand->elemNum != elem->elemNum);

          // if cand is new to neighborhood we add it with counter = 1,
          // otherwise increment the counter
          bool found = false;
          for(unsigned int o = 0; o < neighborhood.GetSize() && ! found; o++)
          {
            if(neighborhood[o].first == cand)
            {
              neighborhood[o].second++;
              found = true;
            }
          }
          // add this element if new
          if(!found)
          {
            pair.first = cand;
            pair.second = 1;
            neighborhood.Push_back(pair);
          }
        }
      }

      // neighborhood is filled now, copy it to elem
      assert(neighborhood.GetSize() > 0); // assume there is no single element case
      elem->extended->neighborhood = new StdVector<std::pair<Elem*, int> >(neighborhood);

      // LOG_DBG2(gridcfs) << "FEN: elem=" << elem->elemNum << " neighbourhood=" << elem->neighborhood->ToString();
    }
  }


  void GridCFS::GetElems( StdVector<Elem*> & elems,
                          const RegionIdType regionId ) {
    LOG_DBG(gridcfs) << "GetElems for region " << region_.ToString(regionId);
    elems.Clear();

    // check if region Id is ALL_REGIONS
    if ( regionId == ALL_REGIONS ) {
      elems.Resize( numElems_ );

      std::copy(orderedElems_.Begin(), orderedElems_.End(),
                elems.Begin());
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
    LOG_DBG(gridcfs) << "GetElems returning '" << elems.GetSize() <<"' elements: " << Elem::ToString(elems);
  }


  StdVector<unsigned int> GridCFS::CalcRegulardGridDiscretization()
  {
    StdVector<unsigned int> grid;

    if(IsGridRegular())
    {
      grid.Resize(3, 0);

      StdVector<double> edges;

      if(volElems_.IsEmpty() || volElems_.First().IsEmpty())
        throw Exception("no volume elements in grid found, maybe dimension mismatch with analysis");

      // take the first vol element of the first vol region. The first ordered element might be a surface element
      GetElemShapeMap(volElems_.First().First(), false)->GetEdgeLength(edges);
      assert(edges.GetSize() == GetDim());

      Matrix<double> m = CalcGridBoundingBox();

      grid[0] = (m[0][1]-m[0][0]) / (edges[0]*.99999); // to avoid rounding errors: 1 / (1/60) gives 59.9999 -> which gives int 59
      grid[1] = (m[1][1]-m[1][0]) / (edges[1]*.99999);
      grid[2] = GetDim() == 3 ? (m[2][1]-m[2][0]) / (edges[2]*.99999) : 1;

      LOG_DBG(gridcfs) << "GRGD: e=" << edges.ToString() << " bb=" << m.ToString() << " -> " << grid.ToString();
    }

    return grid;
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
    elems.Clear();

    Integer index = surfRegionIds_.Find(regionId);
    if ( index != -1 ) {
      UInt numElems = surfElems_[index].GetSize();

      for(UInt i=0; i<numElems; i++)
      {
        elems.Push_back(dynamic_cast<SurfElem*>(surfElems_[index][i]));
      }

    } else {
      EXCEPTION( "GridCFS: The surface region with id '" << regionId
                 << "' (" << GetRegion().ToString(regionId)
                 << ") was not found in the grid!" );
    }
  }


  void GridCFS::GetElemsByName( StdVector<Elem*> & elems,
                                const std::string & elemsName ) {

    LOG_DBG(gridcfs) << "GetElemsByName for name " << elemsName;
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

  void GridCFS::GetElemNumsByName( StdVector<UInt> & elemNums,
                                       const std::string & elemName )
  {
    LOG_DBG(gridcfs) << "GetElemNumsByName for name " << elemName;
    if ( nameTypeMap_.find(elemName) == nameTypeMap_.end() )
    {
      std::cerr << "Available Regions: " << std::endl;
      std::map<std::string, EntityList::DefineType>::iterator eIter = nameTypeMap_.begin();
      for(;eIter != nameTypeMap_.end(); ++eIter){
        std::cerr << eIter->first << std::endl;
      }
      throw Exception("There are no entities with name '" + elemName + "' in the grid");
    }
    
    Integer idx = -1;
    EntityList::DefineType defType = nameTypeMap_[elemName];

    switch ( defType )
    {
      case EntityList::REGION:
        
        idx = volRegionIds_.Find( region_.Parse(elemName) );
        if ( idx != -1 )
        {
          UInt numElems = volElems_[idx].GetSize();
          elemNums.Resize( numElems );
          for ( UInt i=0; i<numElems; ++i )
          {
            elemNums[i] = volElems_[idx][i]->elemNum;
          }
        }
        else
        {
          idx = surfRegionIds_.Find( region_.Parse(elemName) );
          if ( idx != -1 )
          {
            UInt numElems = surfElems_[idx].GetSize();
            elemNums.Resize( numElems );
            for ( UInt i=0; i<numElems; ++i )
            {
              elemNums[i] = surfElems_[idx][i]->elemNum;
            }
          }
          else
          {
            EXCEPTION( "The region with name '" << elemName
                       << "' was not found in the grid" );
          }
          
        }
        break;
        
      case EntityList::NAMED_ELEMS:
        
        idx = namedElemNames_.Find( elemName );
        if ( idx != -1 )
        {
          elemNums.Resize( namedElems_[idx].GetSize() );
          std::copy( namedElems_[idx].Begin(),
                     namedElems_[idx].End(),
                     elemNums.Begin() );
        }
        else
        {
          EXCEPTION( "There are no named elements called '" << elemName
                     << "' in the grid" );
        }
        break;
        
      default:
        EXCEPTION( "GetElemNumsByName(" + elemName + ") cannot be called for named nodes" );
        break;
    }
  }
  

  void GridCFS::GetElemNodes( StdVector<UInt> & connect,
                              const UInt iElem ) {

    if ( iElem > numElems_ ) {
      EXCEPTION( "GridCFS: There are only " << numElems_
                 << " elements in the grid! You requested element number "
                 << iElem << ". Go check your mesh file!" );
    }

 #ifndef NDEBUG
    if ( orderedElems_[iElem-1] == nullptr ) {
      EXCEPTION( "Element with Nr. " << iElem << " is not contained in mesh!" );
    }
 #endif

    connect = orderedElems_[iElem-1]->connect;

  }


  void GridCFS::GetElemNodesCoord( Matrix<Double> & coordMat,
                                   const StdVector<UInt> & connect,
                                   bool updated ) {

    LOG_DBG2(gridcfs) << "GetElemNodeCoord() for connect list: " << connect.ToString();
    coordMat.Resize(dim_, connect.GetSize());

    if( updated == true && deltCoords_.GetSize() != 0 ) {
      for (UInt k = 0, cs = connect.GetSize() ; k < cs; k++) {
        for (UInt actDim=0; actDim < dim_; actDim++) {
          coordMat[actDim][k] = coords_[connect[k]-1][actDim];
          coordMat[actDim][k] += deltCoords_[connect[k]-1][actDim];
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


  void GridCFS::GetElemsNextToNodes( StdVector<const Elem*> & elemList,
                                     const StdVector<UInt> & nodeList,
                                     const StdVector<RegionIdType> & regionIds) {
    elemList.Clear();
    for (UInt n = 0; n < nodeList.GetSize(); n++) {
      StdVector<const Elem*> nElemList;
      GetElemsNextToNode(nElemList, nodeList[n], regionIds);
      for (UInt iE = 0; iE < nElemList.GetSize(); iE++) {
        bool found = false;
        for (UInt i = 0; i < elemList.GetSize() && !found; i++) {
          found = elemList[i] == nElemList[iE];
        }
        if (!found) {
          elemList.Push_back(nElemList[iE]);
        }
      }
    }
  }

  void GridCFS::GetElemsNextToNode( StdVector<const Elem*> & elemList,
                                     const UInt & node) {
    SetNodesToElemsMap();
    
    const UInt maxIdx = nodeElemMapIndices_[node + 1];
    elemList.Clear();
    
    for (UInt idx = nodeElemMapIndices_[node]; idx < maxIdx; idx++) {
      const Elem* elem = GetElem(nodeElemMap_[idx]);
      elemList.Push_back(elem);
    }
  }

  void GridCFS::GetElemsNextToNode( StdVector<const Elem*> & elemList,
                                     const UInt & node,
                                     const StdVector<RegionIdType>& regionIds) {
    SetNodesToElemsMap();
    
    const UInt maxIdx = nodeElemMapIndices_[node + 1];
    const UInt nRegions = regionIds.GetSize();
    elemList.Clear();
    
    for (UInt idx = nodeElemMapIndices_[node]; idx < maxIdx; idx++) {
      const Elem* elem = GetElem(nodeElemMap_[idx]);
      for (UInt iR = 0; iR < nRegions; iR++) {
        if (regionIds[iR] == elem->regionId) {
          elemList.Push_back(elem);
        }
      }
    }
  }

  void GridCFS::GetNumOfElemsNextToNodes( UInt & num,
                                     const UInt & node,
                                     const StdVector<RegionIdType>& regionIds) {
    StdVector<const Elem*> elemList;
    GetElemsNextToNode(elemList, node, regionIds);
    num = elemList.GetSize();
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
                                    StdVector<const Elem*> & elemList,
                                    bool onlyLinNodes) {

    std::set<UInt> elemNodes;
    std::set<UInt>::iterator it;
    UInt iElem, iNode, numElemCorners;

    // First, create a set with node numbers of elements
    for ( iElem = 0; iElem < elemList.GetSize(); iElem++ ) {
      StdVector<UInt> const & connecth = elemList[iElem]->connect;
      ElemShape & actShape = Elem::shapes[elemList[iElem]->type];
      if (onlyLinNodes == true)
        numElemCorners = actShape.numNodes;
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

  void GridCFS::GetNodesOfElemList( StdVector<UInt> & nodeList,
                                    StdVector<Elem*> & elemList,
                                    bool onlyLinNodes) {

    std::set<UInt> elemNodes;
    std::set<UInt>::iterator it;
    UInt iElem, iNode, numElemCorners;

    // First, create a set with node numbers of elements
    for ( iElem = 0; iElem < elemList.GetSize(); iElem++ ) {
      StdVector<UInt> const & connecth = elemList[iElem]->connect;
      ElemShape & actShape = Elem::shapes[elemList[iElem]->type];
      if (onlyLinNodes == true)
        numElemCorners = actShape.numNodes;
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


  void GridCFS::SetNodeOffset( const StdVector<UInt>& nodes,
                               const Vector<Double>& offsets ) {

    // Check if node offsets were already set
    if( deltCoords_.GetSize() == 0 ) {
      deltCoords_.Resize( coords_.GetSize() );
      for( UInt i = 0; i < coords_.GetSize(); ++i ) {
        deltCoords_[i].Resize(dim_);
        deltCoords_[i].Init();
      }
    }

    // Set delta coordinates
    for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
      Vector<Double> actOffset(dim_);
      for( UInt iDim = 0; iDim < dim_; iDim++ ) {
        actOffset[iDim] = offsets[iNode*dim_ + iDim];
      }
      deltCoords_[nodes[iNode]-1]= actOffset;
    }
  }


  bool GridCFS::HasNodalOffset() {

    if( deltCoords_.GetSize() != 0 ) {
      return true;
    } else {
      return false;
    }

  }
  void GridCFS::ClearNodeToElemConnectivity() {
    mappedNodeToElems_ = false;
    nodeElemMapIndices_.Clear(false);
    nodeElemMap_.Clear(false);
  }
  
  void GridCFS::SetNodesToElemsMap()
  {
    #pragma omp critical (CoefFunctionAccumulator)
    {
      if(!mappedNodeToElems_) {
	nodeElemMapIndices_.Resize(GetNumNodes()+2);
        nodeElemMapIndices_.Init(0);
        for (UInt e = 0; e < numElems_; e++) {
          Elem* elem = orderedElems_[e];
          for (UInt n = 0; n < elem->connect.GetSize(); n++) {
            nodeElemMapIndices_[elem->connect[n]]++;
          }
        }
        UInt idx = 0;
        for (UInt n = 1; n < GetNumNodes() + 2; n++) {
          UInt add = nodeElemMapIndices_[n];
          nodeElemMapIndices_[n] = idx;
          idx += add;
        }
        UInt dummy = GetNumElems() + 2;
        nodeElemMap_.Resize(idx);
        nodeElemMap_.Init(dummy);
        for (UInt e = 0; e < numElems_; e++) {
          Elem* elem = orderedElems_[e];
          for (UInt n = 0; n < elem->connect.GetSize(); n++) {
            UInt node = elem->connect[n];
            UInt sIdx;
            for (sIdx = nodeElemMapIndices_[node]; nodeElemMap_[sIdx] != dummy; sIdx++) {}
            nodeElemMap_[sIdx] = e + 1;
          }
        }
        mappedNodeToElems_ = true;
      }
    }
  }

  void GridCFS::TriggerAutoLayerGeneration() {
    // if no param object is present, just leave
    if (!param_)
      return;

    PtrParamNode layerGenNode = param_->Get("domain")->Get("layerGenerationList", ParamNode::PASS);
    // if autoLayerGeneration is not specified, leave
    if (!layerGenNode) 
      return;
    // number of layer generations specified
    UInt numGenLayers = layerGenNode->GetChildren().size();
    // leave if there is no layer generations specified
    if (numGenLayers < 1)
      return;

    // initialize variables...
    // specified region to act on
    std::string surfRegionName;
    // ID of surface region
    RegionIdType surfRegionId;
    // the layer gen node
    PtrParamNode actLayerGenNode;
    // region on grid
    std::string regionNameOnGrid;
    // nodeIds of a region
    StdVector<UInt> surfRegionNodeIds;
    // bools to check if the specified surface region is on the grid
    bool surfRegionFound = false;

    // check if specified surface region actually exist..
    // loop over specified layer generations
    for (UInt iLayers = 0; iLayers < numGenLayers; ++iLayers) {
      actLayerGenNode = layerGenNode->GetChildren()[iLayers];
      // loop over surface regions on grid and compare
      for (UInt iSurfRegion = 0; iSurfRegion < this->GetNumSurfRegions(); ++iSurfRegion) {
        actLayerGenNode->Get("sourceSurfRegion")->GetValue("name", surfRegionName);
        regionNameOnGrid = this->GetRegion().ToString(surfRegionIds_[iSurfRegion]);
        // if correct surface region is there, take the next step...
        if (surfRegionName == regionNameOnGrid) {
          surfRegionFound = true;
          // check if the surfaceRegion contains nodes at all
          this->GetNodesByRegion(surfRegionNodeIds, surfRegionIds_[iSurfRegion]);
          if (surfRegionNodeIds.IsEmpty()) {
            EXCEPTION("Surface region '" << surfRegionName << "' does not contain nodes!");
          }
          // now we know everything is fine, so start layer generation
          surfRegionId = GetRegionId(surfRegionName);
          this->GenerateExternalLayer(surfRegionId, actLayerGenNode);
        }
      }
      // throw exception if a specified surface is not found
      if (!surfRegionFound)
        EXCEPTION("SurfaceRegion " << surfRegionName << " not found on the specified Grid!");
      // reset the bool to start with next layer
      surfRegionFound = false;
    }
  }


  void GridCFS::GenerateExternalLayer(const RegionIdType surfRegionId, const PtrParamNode layerGenNode) {
    // declare variables...
    // layer generation parameters
    Double elemHeight = 0;
    Double numLayers = 0;
    // Id of region nodes and coordinates of nodes
    StdVector<UInt> surfRegionNodeIds;
    StdVector<Vector<Double>> surfNodeCoords;
    UInt numSurfNodes;
    // elements on regions
    StdVector<SurfElem*> surfRegionElems;
    UInt numSurfElems;
    StdVector<StdVector<UInt>> surfElemConnectivity;
    // name, id, and elems of new region
    std::string layerName;
    RegionIdType layerRegionId;
    // filepaths for reading / writing geometry from / to file
    std::string readFilePath, writeFilePath;
    // string to store the simple geometry type, if specified
    std::string analyticGeomType;

    // parameters of elements in new layer...
    // connectivity the current element that is created within the layer
    StdVector<UInt> layerConnectivity;
    // type of created element 
    Elem::FEType currLayerElemType;
    // type of element on the surface region
    Elem::FEType currSurfElemType;
    // pointers to created elements in the layer
    StdVector<StdVector<Elem*>> addedElems;
    // corresponding IDs
    StdVector<StdVector<UInt>> addedElemIds;
    // number of nodes of a element on the surface region
    UInt numNodesInSurfElement;
    // and of the newly created element
    UInt numNodesInLayerElement;
    // store the node connectivity of the surface region and every subsequent isosurface
    StdVector<StdVector<UInt>> connectNodeIdx;
    // helper variable to pass layer connectivity to function
    UInt* ptrLayerConnectivity;
    // coordinates of the current node (where the new node is created upon)
    StdVector<Vector<Double>> currNodeCoords;
    // IDs within the new layer
    StdVector<StdVector<UInt>> allLayerNodeIds;
    // name and ID of newly created surface regions (each new iso surface within the layer)
    // we create these new iso-surface regions for being able to compute surface geometry
    // later on
    std::string newSurfRegionName;
    RegionIdType newSurfRegionId;
    StdVector<StdVector<SurfElem*>> addedSurfElems;
    StdVector<StdVector<UInt>> addedSurfElemIds;

    // get children of layer generation node to retrieve info
    PtrParamNode layerParamsNode, geomParamsNode, writeNode;
    layerParamsNode = layerGenNode->Get("extrusionParameters");
    geomParamsNode = layerGenNode->Get("surfGeometry");

    // extract parameters from layerGenNode
    layerParamsNode->GetValue("numLayers", numLayers);
    layerParamsNode->GetValue("elemHeight", elemHeight);

    // get, node Ids, coordinates and elements
    this->GetNodesByRegion(surfRegionNodeIds, surfRegionId);
    this->GetNodeCoordinates(surfNodeCoords, surfRegionNodeIds, false);
    this->GetSurfElems(surfRegionElems, surfRegionId);
    numSurfElems = this->GetNumElems(surfRegionId);
    numSurfNodes = surfRegionNodeIds.GetSize();

    // extract specified name of new region and add to grid
    layerGenNode->GetValue("name", layerName);
    layerRegionId = this->AddRegion(layerName, VOLUME_REGION);

    // if an analytic geometry type is specified, we use a simplified computation
    if (geomParamsNode->Get("analyticApproximation", ParamNode::PASS)) {
      ComputeGeometryOfSimpleType(surfRegionId, layerGenNode);
    }
    else if (geomParamsNode->Get("fromInputFile", ParamNode::PASS)) {
      // if a read path is specified, read from file
      ReadGeometryFromFile(layerGenNode);
    }
    else if (geomParamsNode->Get("cgalApproximation", ParamNode::PASS)) {
      // if cgalApproximation is specified, use cgal to approximate a generic convex boundary
      ComputeGeometryOnSurfaceRegionNodes(layerGenNode);
    }
    // check for specified write files...
    if (layerGenNode->Get("outputFile", ParamNode::PASS)) {
      layerGenNode->Get("outputFile")->GetValue("name", writeFilePath);
      // if a write file path is specified, write the geometry to the file
      WriteGeometryToFile(layerGenNode);
    }

    // prepare computation of new nodes
    currNodeCoords = surfNodeCoords;
    allLayerNodeIds = StdVector<StdVector<UInt>>(numLayers+1);
    allLayerNodeIds[0] = surfRegionNodeIds;
    UInt newNodeId;
    // compute new nodes iteratively...
    for (UInt iLayers = 1; iLayers <= numLayers; ++iLayers) {
      allLayerNodeIds[iLayers].Resize(numSurfNodes);
      for (UInt iNodes = 0; iNodes < numSurfNodes; ++iNodes) {
        // compute new node coordinates
        currNodeCoords[iNodes] += geometryRegionMap_[surfRegionId]->normalVectors_[iNodes] * elemHeight;
        // add new node to grid
        this->AddNode( currNodeCoords[iNodes], newNodeId );
        // WARN("Added Node with ID " << newNodeId << " on layer " << iLayers);
        // store iD
        allLayerNodeIds[iLayers][iNodes] = newNodeId;
      }
    }

    // add new elems to grid...
    // as the new elems are prismatic, there will be one layer per linear element
    if (this->IsQuadratic() == true) {
      EXCEPTION("Layer Generation for quadratic elements not implemented yet. Use linear elements instead.");
    } else {
      // prepare computation
      surfElemConnectivity.Resize(numSurfElems);
      connectNodeIdx.Resize(numSurfElems);
      addedElems.Resize(numLayers);
      addedElemIds.Resize(numLayers);
      addedSurfElems.Resize(numLayers);
      addedSurfElemIds.Resize(numLayers);
      // iteratively create layer (volume) elements, iso-surface regions and surface elements
      for (UInt iLayers = 0; iLayers < numLayers; ++iLayers) {
        addedElems[iLayers].Resize(numSurfElems);
        addedElemIds[iLayers].Resize(numSurfElems);
        addedSurfElems[iLayers].Resize(numSurfElems);
        addedSurfElemIds[iLayers].Resize(numSurfElems);
        
        // add a new surface region for every new iso surface
        newSurfRegionName = layerName + "_IsoSurface_" + std::to_string(iLayers+1);
        newSurfRegionId = this->AddRegion(newSurfRegionName, SURFACE_REGION);

        // add new elements...
        for (UInt iSurfElems = 0; iSurfElems < numSurfElems; ++iSurfElems) {
          // create new elements. They will be assigned to the orderedElems_ AddVolumeElems()
          // and AddSurfaceElems() functions. So they will be deleted in the GridCFS destructor
          addedElems[iLayers][iSurfElems] = new Elem;
          addedSurfElems[iLayers][iSurfElems] = new SurfElem;
        }
        // add new elements of current layer to the grid and obtain element IDs
        AddVolumeElems( layerRegionId, addedElems[iLayers], addedElemIds[iLayers]);
        AddSurfaceElems( newSurfRegionId, addedSurfElems[iLayers], addedSurfElemIds[iLayers]);

        // next, we need to set all the necessary information for the new elements...
        for (UInt iSurfElems = 0; iSurfElems < numSurfElems; ++iSurfElems) {
          // get type of current surface element
          currSurfElemType = surfRegionElems[iSurfElems]->type;
          // check for type of the surface element and assign type of layer element accordingly
          switch (currSurfElemType) {
            case Elem::FEType::ET_TRIA3:
                currLayerElemType = Elem::FEType::ET_WEDGE6;
                numNodesInSurfElement = 3;
                break;
            case Elem::FEType::ET_QUAD4:
                currLayerElemType = Elem::FEType::ET_HEXA8;
                numNodesInSurfElement = 4;
                break;
            default:
                EXCEPTION("Layer Generation for surface element type "<< currSurfElemType <<" not implemented yet!" <<
                          "use linear triangles or quadrangles instead!");
          }
          // set the number of nodes in the current layer element
          numNodesInLayerElement = numNodesInSurfElement*2;

          // get the connectivity of the corresponding surface element and store for later use
          // we only need to do this once for each surface element on the interface
          if (iLayers == 0) {
            surfElemConnectivity[iSurfElems].Resize(numNodesInSurfElement);
            surfElemConnectivity[iSurfElems] = surfRegionElems[iSurfElems]->connect;
            connectNodeIdx[iSurfElems].Resize(numNodesInSurfElement);
            // find the indices of the connection list of the surface elements as we need them 
            // to assign the connections in the new layer elements
            for (UInt iNode = 0; iNode < numNodesInSurfElement; ++iNode) {
              connectNodeIdx[iSurfElems][iNode] = allLayerNodeIds[0].Find(surfElemConnectivity[iSurfElems][iNode]);
            }
            // check if the new elements have more nodes than elements in the original grid and set 
            // the maximum number of nodes occurring in any global element
            maxNumElemNodes_ = (numNodesInLayerElement > maxNumElemNodes_) ? numNodesInLayerElement : maxNumElemNodes_;
          } // if (iLayers == 0)

          // assign connectivity to current layer element...
          layerConnectivity.Resize(numNodesInLayerElement);
          for (UInt iNode = 0; iNode < numNodesInSurfElement; ++iNode) {
            layerConnectivity[iNode] = allLayerNodeIds[iLayers][connectNodeIdx[iSurfElems][iNode]];
            layerConnectivity[iNode+numNodesInSurfElement] = allLayerNodeIds[iLayers+1][connectNodeIdx[iSurfElems][iNode]];
          }
          // convert from StdVector to UInt* as SetElemData() only takes UInt*
          // point on the volume-element connectivity
          ptrLayerConnectivity  = &layerConnectivity[0];
          // assign type, region, and connectivity to element
          this->SetElemData(addedElemIds[iLayers][iSurfElems], currLayerElemType, layerRegionId, ptrLayerConnectivity);
          // point on the surface-element connectivity
          ptrLayerConnectivity  = &layerConnectivity[numNodesInSurfElement];
          this->SetElemData(addedSurfElemIds[iLayers][iSurfElems], currSurfElemType, newSurfRegionId, ptrLayerConnectivity);
        }
        // assign surface region to the volume region in order to find the connected regions lateron
        volumeSurfaceRegionMap_[layerRegionId].Push_back(newSurfRegionId);
      }
      // finally, also add the interface surface region to the map
      volumeSurfaceRegionMap_[layerRegionId].Insert(0, surfRegionId);
    }
    // Depending on the direction of layer generation and orientation of the interface surface elements,
    // the volume elements might not be oriented correctly. Hence, we check for negative Jacobi determinants 
    // and try correcing the orientation if not.
    CorrectElementConnectivities(layerRegionId);
    // finally, as the grid now contains new nodes and elements, reset mappedNodeToElems_ so that 
    // SetNodesToElemsMap() will be called once more when needed in future.
    mappedNodeToElems_ = false;
  };


  void GridCFS::ReadGeometryFromFile(const PtrParamNode layerGenNode) {
    // get parameters from the layer generation node...
    // filepath, delimiter character and comment character
    std::string readFilePath, delimiter, commentCharacter;
    // which info is stored in which column
    UInt nodeNumberCol,
         xComp_nVecCol, yComp_nVecCol, zComp_nVecCol, 
         xComp_tMinVecCol, yComp_tMinVecCol, zComp_tMinVecCol,
         xComp_tMaxVecCol, yComp_tMaxVecCol, zComp_tMaxVecCol, 
         minCurvCol, maxCurvCol;
    PtrParamNode layerParamsNode = layerGenNode->Get("extrusionParameters");
    PtrParamNode readNode = layerGenNode->Get("surfGeometry")->Get("fromInputFile");
    readNode->GetValue("name", readFilePath);
    readNode->GetValue("delimiter", delimiter);
    readNode->GetValue("commentCharacter", commentCharacter);
    readNode->GetValue("nodeNumberCol", nodeNumberCol);
    readNode->GetValue("xComp_nVecCol", xComp_nVecCol);
    readNode->GetValue("yComp_nVecCol", yComp_nVecCol);
    readNode->GetValue("zComp_nVecCol", zComp_nVecCol);
    readNode->GetValue("xComp_tMinVecCol", xComp_tMinVecCol);
    readNode->GetValue("yComp_tMinVecCol", yComp_tMinVecCol);
    readNode->GetValue("zComp_tMinVecCol", zComp_tMinVecCol);
    readNode->GetValue("xComp_tMaxVecCol", xComp_tMaxVecCol);
    readNode->GetValue("yComp_tMaxVecCol", yComp_tMaxVecCol);
    readNode->GetValue("zComp_tMaxVecCol", zComp_tMaxVecCol);
    readNode->GetValue("minCurvCol", minCurvCol);
    readNode->GetValue("maxCurvCol", maxCurvCol);

    // check for valid columns
    if(xComp_nVecCol==0 || yComp_nVecCol==0 || zComp_nVecCol==0 || 
       xComp_tMinVecCol==0 || yComp_tMinVecCol==0 || zComp_tMinVecCol==0 ||
       xComp_tMaxVecCol==0 || yComp_tMaxVecCol==0 || zComp_tMaxVecCol==0 || 
       minCurvCol==0 || maxCurvCol==0){
      EXCEPTION("Wrong parameters for reading geometry file: Column indices need to be one based.");
    }
    // check for valid delimiter / comment characters
    if(commentCharacter.size() > 1 || delimiter.size() > 1){
      EXCEPTION("Wrong parameters for reading geometry file: Comment and delimiter strings need to be single characters!");
    }
    // check for correct filepath
    std::ifstream readFile(readFilePath.c_str());
    if(!readFile.good()){
      EXCEPTION("Wrong parameters for reading geometry file: Could not find file \"" + readFilePath + "\"!");
    }
    else if (!readFile) {
      EXCEPTION("Error reading geometry file from the auto layer generation node.")
    }

    // now we know the parameters fit and the file exists, so initialize the geometryRegionMap_
    std::string surfRegionName;
    RegionIdType surfRegionId;
    StdVector<UInt> nodeIds;
    UInt numNodes;
    // get the surface node Ids
    layerGenNode->Get("sourceSurfRegion")->GetValue("name", surfRegionName);
    surfRegionId = GetRegionId(surfRegionName);
    this->GetNodesByRegion(nodeIds, surfRegionId);
    numNodes = nodeIds.GetSize();
    // add new entry in the geometryRegionMap_ to store geometry
    geometryRegionMap_[surfRegionId] = shared_ptr<NodeGeometry>(new NodeGeometry(numNodes));
    // set the node Ids to the struct
    geometryRegionMap_[surfRegionId]->nodeIds_ = nodeIds;

    // variables for file reading
    std::string currentLine;
    string::size_type startPosition;   // position to start reading for each line
    StdVector<std::string> readColumn; // a single column
    std::istringstream iss;            // stream for reading
    std::string readToken;             // a single token
    readFile >> std::ws;               // remove whitespace

    // variables for storing the data
    UInt readNodeId;
    UInt storedNodeIdx = 0;
    UInt guess = 0;

    // now read the file...
    while (std::getline(readFile, currentLine)) {
      // skip empty or commented lines
      if (currentLine.empty() || currentLine[0] == commentCharacter[0]) {
        continue;
      }
      // ignore leading whitespace
      startPosition = 0;
      while (startPosition < currentLine.size() && std::isspace(currentLine[startPosition], std::locale()))
        startPosition++;
      currentLine.erase(0,startPosition);

      std::istringstream iss(currentLine);
      readColumn.Clear();
      // store content line by line into vector
      while (std::getline(iss, readToken, delimiter[0])) {
        readColumn.push_back(readToken);
      }
      readNodeId = std::stoul(readColumn[0]);
      // try to find the correct index to store wrong ordered data, with an initial guess of correctly ordered data
      storedNodeIdx = geometryRegionMap_[surfRegionId]->nodeIds_.Find(readNodeId, guess, false);
      ++guess;
      // store the information accordingly
      geometryRegionMap_[surfRegionId]->normalVectors_[storedNodeIdx].Resize(3);
      geometryRegionMap_[surfRegionId]->normalVectors_[storedNodeIdx][0] = std::stod(readColumn[xComp_nVecCol-1]);
      geometryRegionMap_[surfRegionId]->normalVectors_[storedNodeIdx][1] = std::stod(readColumn[yComp_nVecCol-1]);
      geometryRegionMap_[surfRegionId]->normalVectors_[storedNodeIdx][2] = std::stod(readColumn[zComp_nVecCol-1]);
      geometryRegionMap_[surfRegionId]->minPrincipalVectors_[storedNodeIdx].Resize(3);
      geometryRegionMap_[surfRegionId]->minPrincipalVectors_[storedNodeIdx][0] = std::stod(readColumn[xComp_tMinVecCol-1]);
      geometryRegionMap_[surfRegionId]->minPrincipalVectors_[storedNodeIdx][1] = std::stod(readColumn[yComp_tMinVecCol-1]);
      geometryRegionMap_[surfRegionId]->minPrincipalVectors_[storedNodeIdx][2] = std::stod(readColumn[zComp_tMinVecCol-1]);
      geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[storedNodeIdx].Resize(3);
      geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[storedNodeIdx][0] = std::stod(readColumn[xComp_tMaxVecCol-1]);
      geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[storedNodeIdx][1] = std::stod(readColumn[yComp_tMaxVecCol-1]);
      geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[storedNodeIdx][2] = std::stod(readColumn[zComp_tMaxVecCol-1]);
      geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[storedNodeIdx] = std::stod(readColumn[minCurvCol-1]);
      geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[storedNodeIdx] = std::stod(readColumn[maxCurvCol-1]);
    }
    readFile.close();
  };

  void GridCFS::WriteGeometryToFile(const PtrParamNode layerGenNode) {
    // get parameters from the layer generation node...
    // filepath, delimiter character and comment character
    std::string writeFilePath, delimiter, commentCharacter;
    // get children of the layer gen node
    PtrParamNode layerParamsNode = layerGenNode->Get("extrusionParameters");
    PtrParamNode writeNode = layerGenNode->Get("outputFile");
    
    writeNode->GetValue("name", writeFilePath);
    writeNode->GetValue("delimiter", delimiter);
    writeNode->GetValue("commentCharacter", commentCharacter);

    // check for valid delimiter / comment characters
    if(commentCharacter.size() > 1 || delimiter.size() > 1){
      EXCEPTION("Wrong parameters for reading geometry file: Comment and delimiter strings need to be single characters!");
    }

    // check if directory-path for the file exists...
    // this is done the same way as in SinglePDE::ReadSensorArrayResults()
    // search for last Slash in fileName
    int idx_lastSlash = writeFilePath.find_last_of("/");
    // if idx_lastSlash = -1 -> "/" not found, else position of the last slash
    // if there is a "/" in the filename -> save directory is not "." -> check if it exists
    if ( idx_lastSlash != -1){
      // get directory name
      std::string directoryName;
      directoryName = writeFilePath.substr(0,idx_lastSlash);
      // ensure errno is cleared and call mkdir with the directory name
      errno = 0;
      int mkdir_call;
      #ifndef WIN32
        mkdir_call = mkdir( directoryName.c_str(), S_IRWXU | S_IRGRP | S_IWGRP | S_IROTH );
      #else
        mkdir_call = _mkdir( directoryName.c_str() );
      #endif
      if ( mkdir_call == -1 && errno == EEXIST ){
        // directory exists, do nothing
        errno = 0;
      } else if ( mkdir_call == 0 ){
        // directory didn't exist but was created, do nothing
      } else{
        // directory didn't exist, and couldn't be created -> raise exception
        EXCEPTION("The directory: '" << directoryName << "' to save the geometry file doesn't exist and couldn't be created! Please create it by hand!" );
      }
    } else {
      // no slash in filename -> do nothing
    }

    // everything seems fine so far, so gather the info for storing...
    // header names
    std::vector<std::string> headerStrings = 
      {"nodeNumber", 
       "nVecX", "nVecY", "nVecZ", 
       "tMinVecX", "tMinVecY", "tMinVecZ", 
       "tMaxVecX", "tMaxVecY", "tMaxVecZ", 
       "minCurv", "maxCurv", 
       "xCoord", "yCoord", "zCoord"};
    // the surface id
    std::string surfRegionName;
    RegionIdType surfRegionId;
    UInt numNodes;
    layerGenNode->Get("sourceSurfRegion")->GetValue("name", surfRegionName);
    surfRegionId = GetRegionId(surfRegionName);
    // node coordinates
    StdVector<Vector<Double>> nodeCoords;
    this->GetNodeCoordinates(nodeCoords, geometryRegionMap_[surfRegionId]->nodeIds_, false);
    // and the number of nodes
    numNodes = geometryRegionMap_[surfRegionId]->numNodes_;

    // now write the file...
    // create output stream
    std::ofstream  writeFile(writeFilePath, std::ios::out);
    if (!writeFile.is_open()) {
                EXCEPTION("Failed to open the geometry file for writing. Try to change the specified name.")
    }
    // Ensure that no precision is lost
    writeFile.precision(15);

    // write header
    writeFile << commentCharacter << " --- Nodal Geometry ---" << std::endl;
    writeFile << commentCharacter;
    for (UInt iColumn = 0; iColumn < headerStrings.size()-1; ++iColumn) {
      writeFile << headerStrings[iColumn] << delimiter;
    }
    writeFile << headerStrings[headerStrings.size()-1] << std::endl;
    writeFile << "# ---------------------------------------------------------------------------";
    writeFile << std::endl;

    // write data
    for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
      writeFile << geometryRegionMap_[surfRegionId]->nodeIds_[iNodes] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][0] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][1] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][2] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][0] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][1] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][2] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][0] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][1] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][2] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iNodes] << delimiter;
      writeFile << geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iNodes] << delimiter;
      writeFile << nodeCoords[iNodes][0] << delimiter;
      writeFile << nodeCoords[iNodes][1] << delimiter;
      writeFile << nodeCoords[iNodes][2] << std::endl;
    }
    writeFile.close();
  };

  void GridCFS::ComputeGeometryOfSimpleType(const RegionIdType& surfRegionId, const PtrParamNode layerGenNode) {
    // declare variables...
    // Node IDs and coordinates of the surface region
    StdVector<UInt> nodeIds;
    StdVector<Vector<Double>> nodeCoords;
    this->GetNodesByRegion(nodeIds, surfRegionId);
    this->GetNodeCoordinates(nodeCoords, nodeIds, false);
    // number of nodes on surface
    UInt numNodes = nodeIds.GetSize();
    // add new entry in the geometryRegionMap_ to store geometry
    geometryRegionMap_[surfRegionId] = shared_ptr<NodeGeometry>(new NodeGeometry(numNodes));
    // set the node Ids to the struct
    geometryRegionMap_[surfRegionId]->nodeIds_ = nodeIds;

    // get type of geometry
    std::string geomType = layerGenNode->Get("surfGeometry")->Get("analyticApproximation")->GetChildren()[0]->GetName();
    ParamNodeList children = layerGenNode->Get("surfGeometry")->Get("analyticApproximation")->Get(geomType)->GetChildren();

    // set origin offset by MP expression. Default is defined in the xsd scheme at (0,0,0).
    Vector<Double> origin(3);
    MathParser* mp = domain->GetMathParser();
    UInt mpHandle = mp->GetNewHandle(true);
    std::string mpExpression;
    layerGenNode->Get("surfGeometry")->Get("analyticApproximation")->Get(geomType)->GetValue("origin_x", mpExpression);
    mp->SetExpr(mpHandle, mpExpression);
    origin[0] = mp->Eval(mpHandle);
    mpHandle = mp->GetNewHandle(true);
    layerGenNode->Get("surfGeometry")->Get("analyticApproximation")->Get(geomType)->GetValue("origin_y", mpExpression);
    mp->SetExpr(mpHandle, mpExpression);
    origin[1] = mp->Eval(mpHandle);
    mpHandle = mp->GetNewHandle(true);
    layerGenNode->Get("surfGeometry")->Get("analyticApproximation")->Get(geomType)->GetValue("origin_z", mpExpression);
    mp->SetExpr(mpHandle, mpExpression);
    origin[2] = mp->Eval(mpHandle);

    // compute geometry based on the passed type
    if (geomType == "sphere") {

      // curvatures can be assigned via the sphere radius
      Double r = sqrt(pow(nodeCoords[0][0]-origin[0],2) + pow(nodeCoords[0][1]-origin[1],2) + pow(nodeCoords[0][2]-origin[2],2));
      Double curv = 1 / r;
      // to compute the tangential vectors, we compute the numerical differential where we vary theta and phi by a small portion
      Double phi, theta, eps, nor, x, y, z, xVar, yVar, zVar;
      eps = 1e-6; // variation of the angle
      for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
        x = nodeCoords[iNodes][0] - origin[0];
        y = nodeCoords[iNodes][1] - origin[1];
        z = nodeCoords[iNodes][2] - origin[2];
        // assign normal vectors
        geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].Resize(this->dim_);
        geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][0] = x;
        geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][1] = y;
        geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][2] = z;
        // compute angles for getting tangential vectors
        if (z/r >= 1.0)
          theta = 0;
        else if (z/r <= -1.0)
          theta = M_PI;
        else
          theta = acos(z/r);
        if (x/sqrt(pow(x,2) + pow(y,2)) >= 1.0 || sqrt(pow(x,2) + pow(y,2)) < 1e-15)
          phi = 0;
        else if (x/sqrt(pow(x,2) + pow(y,2)) <= -1.0)
          phi = M_PI;
        else
          phi = acos(x/sqrt(pow(x,2) + pow(y,2)));
        if (y < 0)
          phi *= -1;
        // first tangential vector
        xVar = r * cos(phi) * sin(theta+eps);
        yVar = r * sin(phi) * sin(theta+eps);
        zVar = r * cos(theta+eps);
        nor = sqrt(pow(x-xVar,2) + pow(y-yVar,2) + pow(z-zVar,2));
        // assign
        geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes].Resize(this->dim_);
        geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][0] = (x - xVar)/nor;
        geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][1] = (y - yVar)/nor;
        geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][2] = (z - zVar)/nor;
        // compute second tangential vector via cross product
        geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].CrossProduct(
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes], 
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes]);
        // finally, assign curvatures
        geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iNodes] = curv;
        geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iNodes] = curv;
      }
    }
    else if (geomType == "cylinder") {
      // get the specified h axis
      std::string axisDirection;
      layerGenNode->Get("surfGeometry")->Get("analyticApproximation")->GetChildren()[0]->GetValue("axisDirection", axisDirection);
      if (axisDirection == "x") {
        // the maximum curvature can be assigned via the cylinder radius
        Double r = sqrt(pow(nodeCoords[0][1]-origin[1],2) + pow(nodeCoords[0][2]-origin[2],2));
        Double curv = 1 / r;
        Double x, y, z, nor;
        for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
          x = 0;
          y = nodeCoords[iNodes][1] - origin[1];
          z = nodeCoords[iNodes][2] - origin[2];
          nor = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
          // assign normal vectors
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][0] = x/nor;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][1] = y/nor;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][2] = z/nor;
          // minimum principal vector points into h axis
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][0] = 1.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][1] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][2] = 0.0;
          // compute maximum principal vector via cross product
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].CrossProduct(
            geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes], 
            geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes]);
          // finally, assign curvatures
          geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iNodes] = 0;
          geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iNodes] = curv;
        }
      }
      else if (axisDirection == "y") {
        // the maximum curvature can be assigned via the cylinder radius
        Double r = sqrt(pow(nodeCoords[0][0]-origin[0],2) + pow(nodeCoords[0][2]-origin[2],2));
        Double curv = 1 / r;
        Double x, y, z, nor;
        for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
          x = nodeCoords[iNodes][0] - origin[0];
          y = 0;
          z = nodeCoords[iNodes][2] - origin[2];
          nor = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
          // assign normal vectors
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][0] = x/nor;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][1] = y/nor;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][2] = z/nor;
          // minimum principal vector points into h axis
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][0] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][1] = 1.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][2] = 0.0;
          // compute maximum principal vector via cross product
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].CrossProduct(
            geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes], 
            geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes]);
          // finally, assign curvatures
          geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iNodes] = 0;
          geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iNodes] = curv;
        }
      }
      else if (axisDirection == "z") {
        // the maximum curvature can be assigned via the cylinder radius
        Double r = sqrt(pow(nodeCoords[0][0]-origin[0],2) + pow(nodeCoords[0][1]-origin[1],2));
        Double curv = 1 / r;
        Double x, y, z, nor;
        for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
          x = nodeCoords[iNodes][0] - origin[0];
          y = nodeCoords[iNodes][1] - origin[1];
          z = 0;
          nor = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
          // assign normal vectors
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][0] = x/nor;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][1] = y/nor;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][2] = z/nor;
          // minimum principal vector points into h axis
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][0] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][1] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][2] = 1.0;
          // compute maximum principal vector via cross product
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].CrossProduct(
            geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes], 
            geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes]);
          // finally, assign curvatures
          geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iNodes] = 0;
          geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iNodes] = curv;
        }
      }
    }
    else if (geomType == "plane") {
      // get the specified normal direction
      std::string normalDirection;
      layerGenNode->Get("surfGeometry")->Get("analyticApproximation")->GetChildren()[0]->GetValue("normalDirection", normalDirection);
      if (normalDirection == "x") {
        // use distance vector to origin for obtaining the sign of the normal vector
        int nVecSign = (nodeCoords[0][0] > origin[0]) ? 1 : -1;
        for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][0] = 1.0 * nVecSign;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][1] = 0.0;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][2] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][0] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][1] = 1.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][2] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][0] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][1] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][2] = 1.0;
          geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iNodes] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iNodes] = 0.0;
        }
      }
      else if (normalDirection == "y") {
        // use distance vector to origin for obtaining the sign of the normal vector
        int nVecSign = (nodeCoords[0][1] > origin[1]) ? 1 : -1;
        for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][0] = 0.0;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][1] = 1.0 * nVecSign;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][2] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][0] = 1.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][1] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][2] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][0] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][1] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][2] = 1.0;
          geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iNodes] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iNodes] = 0.0;
        }
      }
      else if (normalDirection == "z") {
        // use distance vector to origin for obtaining the sign of the normal vector
        int nVecSign = (nodeCoords[0][2] > origin[2]) ? 1 : -1;
        for (UInt iNodes = 0; iNodes < numNodes; ++iNodes) {
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][0] = 0.0;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][1] = 0.0;
          geometryRegionMap_[surfRegionId]->normalVectors_[iNodes][2] = 1.0 * nVecSign;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][0] = 1.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][1] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iNodes][2] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes].Resize(this->dim_);
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][0] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][1] = 1.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iNodes][2] = 0.0;
          geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iNodes] = 0.0;
          geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iNodes] = 0.0;
        }
      }
    } else {
      EXCEPTION("You passed an invalid geometry type.");
    }
  };

  // =======================================================================
  // CGAL Functions
  // =======================================================================
#if defined(USE_CGAL) && defined(USE_EIGEN)
  void GridCFS::ComputeGeometryOnSurfaceRegionNodes(const PtrParamNode layerGenNode) {
    // declare variables...
    PtrParamNode cgalParamsNode =  layerGenNode->Get("surfGeometry")->Get("cgalApproximation");
    // set parameters for Monge fitting
    UInt degreePolyFit, degreeMongeCoeff;
    cgalParamsNode->GetValue("degreePolyFit", degreePolyFit);
    cgalParamsNode->GetValue("degreeMongeCoeff", degreeMongeCoeff);
    if (degreeMongeCoeff > degreePolyFit || degreeMongeCoeff > 4) {
      EXCEPTION("Parameter 'degreeMongeCoeff' must be < 'degreePolyFit' and < 4!");
    }
    // get the surface id
    RegionIdType surfRegionId;
    std::string surfRegionName;
    layerGenNode->Get("sourceSurfRegion")->GetValue("name", surfRegionName);
    surfRegionId = GetRegionId(surfRegionName);
    // a monge form that is computed for every node
    MongeForm mongeForm;
    // minimum required points for fitting the monge base (so interpolation is performed instead of approximation)
    UInt minNumPoints = (degreePolyFit + 1) * (degreePolyFit + 2) / 2;
    // all node IDs of the surface region
    StdVector<UInt> nodeIds;
    this->GetNodesByRegion(nodeIds, surfRegionId);
    // Node IDs and coordinates of the surface region
    StdVector<UInt> surfRegionNodeIds;
    StdVector<Vector<Double>> surfNodeCoords;
    this->GetNodesByRegion(surfRegionNodeIds, surfRegionId);
    this->GetNodeCoordinates(surfNodeCoords, surfRegionNodeIds, false);

    // number of nodes on surface
    UInt numNodes = nodeIds.GetSize();
    // vector that holds the coordinates of the point set for which the form is computed 
    StdVector<Vector<Double>> currNodeCoords;
    // a vector in Point_3 format. This data format is used by CGAL, so I need to convert the currNodeCoords.
    std::vector<DPoint> coordsPoint_3;
    // corresponding node IDs of the current point set
    StdVector<UInt> currNodeIds;
    // surface elements surrounding the vertex node on which the monge base is fitted
    StdVector<const Elem*> currElemsNextToNode;
    // StdVector needed to call GetElemsNextToNodes()
    StdVector<RegionIdType> searchRegionIds;
    searchRegionIds.Push_back(surfRegionId);

    // variables needed to ensure consistent orientation of the monge base 
    // so that normal vectors point outwards the (presumably convex) surface
    Vector<Double> innerPoint = Vector<Double>(3,0);
    // compute the average of all points in the current point cloud to obtain a 
    // point inside the convex layer (needed for orienting the normal vectors)
    for (UInt iSurfNodes = 0; iSurfNodes < numNodes; ++iSurfNodes) {
      innerPoint[0] += surfNodeCoords[iSurfNodes][0];
      innerPoint[1] += surfNodeCoords[iSurfNodes][1];
      innerPoint[2] += surfNodeCoords[iSurfNodes][2];
    }
    innerPoint /= Double(numNodes);
    DPoint innerPoint_3(innerPoint[0], innerPoint[1], innerPoint[2]);
    // variables needed to put the vertex onto the [0] position for each point set
    UInt zeroId;
    UInt vertexIdx;
    // vector to extract current geometry vector from vector_3 format before storing
    Vector<Double> tempVec = Vector<Double>(3,0);
    // temporarily store double value
    Double tempVar;

    // add new entry in the geometryRegionMap_ to store geometry
    geometryRegionMap_[surfRegionId] = shared_ptr<NodeGeometry>(new NodeGeometry(numNodes));
    // set the node Ids to the struct
    geometryRegionMap_[surfRegionId]->nodeIds_ = nodeIds;
    
    // compute surface normals. We need to iterate ofer every surface node and 
    // gather at least the 1-ring neighbourhood for approximation
    for (UInt iSurfNodes = 0; iSurfNodes < numNodes; ++iSurfNodes) {
      currNodeIds.Clear();
      currNodeIds.Push_back(nodeIds[iSurfNodes]);
      // gather enough points for monge fitting
      while (currNodeIds.GetSize() < minNumPoints) {
        // get elements 
        this->GetElemsNextToNodes(currElemsNextToNode, currNodeIds, searchRegionIds);
        // get node ids of elements
        this->GetNodesOfElemList(currNodeIds, currElemsNextToNode);
      }
      // assure that the vertex ID is the first entry in the vector so that CGAL uses it as vertex
      zeroId = currNodeIds[0];
      vertexIdx = currNodeIds.Find(nodeIds[iSurfNodes]);
      currNodeIds[0] = currNodeIds[vertexIdx];
      currNodeIds[vertexIdx] = zeroId;
      // get coordinates of the current nodeset
      this->GetNodeCoordinates(currNodeCoords, currNodeIds, false);
      // first, we need to switch the format of out point-coordinate representation
      this->ConvertVectorToPoint_3Format(coordsPoint_3, currNodeCoords);
      // then, set up the monge form
      this->SetUpMongeForm(mongeForm, degreePolyFit, degreeMongeCoeff, coordsPoint_3);

      // vector from inner point to current vertex in the DVector format that is used by CGAL.
      DVector innerVec(innerPoint_3, coordsPoint_3[0]);
      // orient the monge base accordingly
      mongeForm.comply_wrt_given_normal(innerVec);
      // extract geometry from MongeForm and store...
      // extract normal vector and store
      tempVec[0] = mongeForm.normal_direction().x();
      tempVec[1] = mongeForm.normal_direction().y();
      tempVec[2] = mongeForm.normal_direction().z();
      geometryRegionMap_[surfRegionId]->normalVectors_[iSurfNodes] = tempVec;
      // we imply to have only convex surfaces... 
      // hence, we invert the sign of the CGAL curvatures and exchange the min/max 
      // curvatures and direction vectors. This way we should get only curvatures >= 0...
      // store maximum principal directions
      tempVec[0] = -mongeForm.minimal_principal_direction().x();
      tempVec[1] = -mongeForm.minimal_principal_direction().y();
      tempVec[2] = -mongeForm.minimal_principal_direction().z();
      geometryRegionMap_[surfRegionId]->maxPrincipalVectors_[iSurfNodes] = tempVec;
      // minimum principal directions
      tempVec[0] = -mongeForm.maximal_principal_direction().x();
      tempVec[1] = -mongeForm.maximal_principal_direction().y();
      tempVec[2] = -mongeForm.maximal_principal_direction().z();
      geometryRegionMap_[surfRegionId]->minPrincipalVectors_[iSurfNodes] = tempVec;
      // maximum principal curvatures
      tempVar = -mongeForm.principal_curvatures(1);
      geometryRegionMap_[surfRegionId]->maxPrincipalCurvatures_[iSurfNodes] = tempVar;
      // minimum principal curvatures
      tempVar = -mongeForm.principal_curvatures(0);
      if (tempVar >= 0)
        geometryRegionMap_[surfRegionId]->minPrincipalCurvatures_[iSurfNodes] = tempVar;
      else
        EXCEPTION("Only convex surfaces are supported for the layerGeneration!");
    }
  };
#else
  void GridCFS::ComputeGeometryOnSurfaceRegionNodes(const PtrParamNode layerGenNode) {
    // this function is a dummy for a future implementation that does not require CGAL
    EXCEPTION("Missing dependencies for GridCFS::ComputeGeometryOnSurfaceRegionNodes!")
  };
#endif



  // =======================================================================
  // Helper Methods
  // =======================================================================

  void GridCFS::CreateSurfaceElements( StdVector<Elem*>& elems,
                                       StdVector<SurfElem*>& surfElems ) {

    LOG_DBG2(gridcfs) << "Starting to map surface elements";

    // 1.) Create vector of vector of elems
    StdVector<StdVector<UInt> > elemNrPerNode;
    UInt nrNodes, iRegion, iElem;
    Elem * ptVolElem = nullptr;
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
    LOG_DBG(gridcfs) << "There are " << elems.GetSize() 
                     << " surface elements to be mapped";
    UInt numElems = elems.GetSize();
    surfElems.Resize(numElems);
    for( UInt iEl = 0; iEl < numElems; ++iEl ) {
      oldElem = elems[iEl];
     // create new surface element
      myElem = new SurfElem();
      myElem->elemNum = oldElem->elemNum;
      myElem->type = oldElem->type;
      myElem->regionId = oldElem->regionId;
      myElem->connect = oldElem->connect;
      if(buildExtendedElemInfo_){
        myElem->extended = new ExtendedElementInfo;
      }
      surfElems[iEl] = myElem;

      // delete old volume element
      delete oldElem;
    }

    // 3.) Iterate over all surface elements and look for each
    //     element, if all of its nodes can be assigned to one or
    //     two neighbours

    UInt surfNodeNr = 0;
    UInt elemsFound = 0;
    UInt elemsAssigned = 0;

    
    for( UInt iEl = 0; iEl < numElems; ++iEl ) {
      elemsAssigned = 0;

      myElem = surfElems[iEl];

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
            myElem->ptVolElems[0] = ptVolElem;
          }
          else {
            myElem->ptVolElems[1] = ptVolElem;
          }

          elemsAssigned++;
        }
      } // loop over element numbers of first node

      // sanity check (avoid the impossible ;-)
      if ( elemsAssigned > 2 ) {
        WARN( "Found " << elemsAssigned
                   << " volume elements for surface element no. "
                   << myElem->elemNum );
      }

    } // loop over surface elements


    // The following code is not needed anymore, as the ElemShapeMap::CalcNormal
    // ALWAYS calculates a normal pointing OUT of the first volume neighbor, which
    // is excatly the functionality previoulsy implemented with the normalSign and
    // the old, unoriented version of CalcNormal().
//    // 4.) Iterate over all surface elements and calculate surface
//    //     flag by comparing the directed and the undirected surface
//    //     normal. If both differ, the surfaceNormalSign = -1, otherwise 1.
//    Vector<Double> normalUndefSign, normalDefSign;
//    Double sign;
//
//    for( surfElIt = surfElems.begin();
//         surfElIt != surfElems.end();
//         surfElIt++ ) {
//
//      myElem = surfElIt->second;
//
//      // check, if each surface element has at least one volume neighbour
//      if ( myElem->ptVolElems[0] == NULL ) {
//        //  EXCEPTION( "Pointer to first volume element is NULL for surface"
//        //                    << " element no. "
//        //                    << surfElems_[iRegion][iSurfElem]->elemNum << ".\n"
//        //                    << "Please check your mesh-file!" );
//        //         }
//        myElem->normalSign = 0;
//      } else {
//        
//        shared_ptr<ElemShapeMap> esm(new LagrangeElemShapeMap(this));
//        esm->SetElem(myElem );
//        LocPoint lp = Elem::shapes[myElem->type].midPointCoord;
//        esm->CalcNormal( normalUndefSign, lp );
//        esm->CalcNormalOutOfVol( normalDefSign, lp, *myElem->ptVolElems[0] );
//
//        // Check if all entries have the same sign by calculating
//        // a scalar product between both vectors.
//        // If it is positive, they point in the same direction,
//        // otherwise an angle of 180 lies in between.
//        
//        sign = normalUndefSign * normalDefSign;
//
//        if ( sign > 0.0 ) {
//          myElem->normalSign = 1;
//        } else {
//          myElem->normalSign = -1;
//        }
//        
//        if(myElem->normalSign != 1 ) {
//          EXCEPTION("SIGN IS DIFFERENT FROM 1!");
//        }
//
//      }
//    }

  }

  void GridCFS::ToInfo(PtrParamNode in) 
  { 
    PtrParamNode gridNode = in->Get("grids")->Get("grid",ParamNode::APPEND);
    gridNode->Get("gridId")->SetValue(gridId_); 
    gridNode->Get("dimensions")->SetValue(GetDim()); 
    gridNode->Get("elements")->SetValue(GetNumElems()); 
    gridNode->Get("nodes")->SetValue(GetNumNodes()); 
    
    // we only have this info when doing homogenization
    if (progOpts->DoDetailedInfo()) {
      in->Get("hull_volume")->SetValue(CalcHullVolume());
      in->Get("structure_volume")->SetValue(CalcVolumeOfAllRegions());
    }

    StdVector<unsigned int> reg = CalcRegulardGridDiscretization();
    if(!reg.IsEmpty()) {
      in->Get("nx")->SetValue(reg[0]);
      in->Get("ny")->SetValue(reg[1]);
      in->Get("nz")->SetValue(reg[2]);
    }

    PtrParamNode list = in->Get("regions"); 

    for(unsigned int i = 0; i < regionData.GetSize(); i++ )
    { 
      PtrParamNode in_ = list->Get("region", ParamNode::APPEND);
      RegionData& rd = regionData[i];
      in_->Get("name")->SetValue(rd.name);
      in_->Get("id")->SetValue(rd.id);
      in_->Get("type")->SetValue(rd.type == VOLUME_REGION ? "volume" : "surface");
      in_->Get("regular")->SetValue(rd.regular);
      in_->Get("hom")->SetValue(rd.homogeneous);
      in_->Get("nodes")->SetValue(GetNumNodes(rd.id));
      in_->Get("elems")->SetValue(GetNumElems(rd.id));
      in_->Get("isQuadratic")->SetValue(IsQuadratic());
      if(progOpts->DoDetailedInfo()) {
        in_->Get("vol")->SetValue(this->CalcVolumeOfRegion(rd.id,true));
      }
    }

    list = in->Get("namedNodes");
    for(unsigned int i = 0; i < namedNodes_.GetSize(); i++)
    {
      PtrParamNode pn = list->GetByVal("nodes", "name", namedNodeNames_[i],ParamNode::APPEND);
      pn->Get("count")->SetValue(namedNodes_[i].GetSize());
      if(namedNodes_[i].GetSize() == 1)
       pn->Get("coord")->SetValue("(" + coords_[namedNodes_[i][0]-1].ToString(TS_PLAIN,",") + ")");
    }

    list = in->Get("namedElements");
    for(unsigned int i = 0; i < namedElems_.GetSize(); i++)
      list->GetByVal("elements", "name", namedElemNames_[i], ParamNode::APPEND)->Get("count")->SetValue(namedElems_[i].GetSize());
    
    // coordinate systems
    if(domain) // does not exist for cfstool
      domain->ToInfo(in);

    // in the cfstool case progOpts is not set
    if(progOpts != NULL && progOpts->DoExportGrid())
    {
      PtrParamNode nl = in->Get("grid/nodeList");
      // Setup large array to resolve zero to many node names
      StdVector<StdVector<unsigned int> > node_names(coords_.GetSize() + 1); // 1-based!
      for(unsigned int s = 0, sn = namedNodeNames_.GetSize(); s < sn; s++)
      {
        const StdVector<unsigned int>& nn = namedNodes_[s];
        for(unsigned int n = 0; n < nn.GetSize(); n++)
          node_names[nn[n]-1].Push_back(s); // fuck mixed 1-based and 0-based :(
      }

      for(unsigned int n = 0, nn = coords_.GetSize(); n < nn; n++)
      {
        PtrParamNode node = nl->Get("node", ParamNode::APPEND);
        node->Get("id", ParamNode::APPEND)->SetValue(n+1); // 1-based!
        node->Get("x", ParamNode::APPEND)->SetValue(coords_[n][0]);
        node->Get("y", ParamNode::APPEND)->SetValue(coords_[n][1]);
        node->Get("z", ParamNode::APPEND)->SetValue(dim_ > 2 ? coords_[n][2] : 0.0); // unfortunately there is garbage set :(

        const StdVector<unsigned int>& ni = node_names[n];
        node->Get("names", ParamNode::APPEND)->SetValue(ni.GetSize());
        for(unsigned int c = 0; c < ni.GetSize(); c++)
          node->Get("name_" + lexical_cast<std::string>(c), ParamNode::APPEND)->SetValue(namedNodeNames_[ni[c]]);
      }
      node_names.Clear();
      // same game for the element names as for the nodes
      StdVector<StdVector<unsigned int> > elem_names(GetNumElems() + 1); // 1-based!
      for(unsigned int s = 0, sn = namedElemNames_.GetSize(); s < sn; s++)
      {
        const StdVector<unsigned int>& nn = namedElems_[s];
        for(unsigned int n = 0; n < nn.GetSize(); n++)
          elem_names[nn[n]-1].Push_back(s);
      }

      PtrParamNode rl = in->Get("grid/regionList");
      for(unsigned int r = 0; r < regionData.GetSize(); r++)
      {
        RegionData& rd = regionData[r];
        SetElementBarycenters(rd.id, false);
        PtrParamNode reg = rl->Get("region", ParamNode::APPEND);
        reg->Get("name")->SetValue(rd.name);
        const StdVector<Elem*>& elems = rd.type == VOLUME_REGION ? volElems_[rd.type_idx] : surfElems_[rd.type_idx];
        for(unsigned int e = 0, n = elems.GetSize(); e < n; e++)
        {
          Elem* elem = elems[e];
          PtrParamNode el = reg->Get("element", ParamNode::APPEND);
          el->Get("id", ParamNode::APPEND)->SetValue(elem->elemNum);
          el->Get("type", ParamNode::APPEND)->SetValue(Elem::feType.ToString(elem->type));
//          Point& bc = elem->barycenter;
//          el->Get("x")->SetValue(bc.data[0], ParamNode::APPEND);
//          el->Get("y")->SetValue(bc.data[1], ParamNode::APPEND);
//          el->Get("z")->SetValue(bc.data[2], ParamNode::APPEND);

          const StdVector<unsigned int>& con = elem->connect;
          el->Get("nodes", ParamNode::APPEND)->SetValue(con.GetSize());
          for(unsigned int c = 0; c < con.GetSize(); c++)
            el->Get("node_" + lexical_cast<std::string>(c), ParamNode::APPEND)->SetValue(con[c]);

          const StdVector<unsigned int>& ni = elem_names[n];
          el->Get("names", ParamNode::APPEND)->SetValue(ni.GetSize());
          for(unsigned int c = 0; c < ni.GetSize(); c++)
            el->Get("name_" + lexical_cast<std::string>(c), ParamNode::APPEND)->SetValue(namedElemNames_[ni[c]]);

        }

      }
    }
  }

  Double GridCFS::CalcHullVolume(bool updated)
  {
    double s = CalcVolumeOfAllRegions(updated);
    // Volume of the bounding box of the grid
    Double cube_vol = 1.0;
    Matrix<Double> m = CalcGridBoundingBox();
    for( UInt d = 0; d < dim_; d++ )
    {
      cube_vol *= m[d][1] - m[d][0];
    }
    LOG_DBG(gridcfs) << "Volume of rectangular dense mesh: " << s;

    if( std::abs(s - cube_vol) / s < 1e-5 ) {
      assert(s >= 0);
      return s;
    }

    // From here on we have either a sparse or non rectangular mesh.
    // We calculate the volume of the mesh using its (possibly virtual)
    // vertices.
    Vector<Double> p(dim_), q(dim_), r(dim_);
    StdVector<Vector<Double> > points, verts;
    points.Reserve(namedNodeNames_.GetSize() * dim_);
    verts.Reserve(pow(2,dim_) * dim_);
    // This stores the indices of named nodes, which we take into account in the following
    // calculations. E.g., if there are too little nodes in namedNodes[i], we skip these nodes.
    assert(namedNodeNames_.GetSize() > 3);
    StdVector<UInt> pointsToNamedNodesMap;
    pointsToNamedNodesMap.Reserve(namedNodeNames_.GetSize());
    for( UInt i=0; i < namedNodeNames_.GetSize(); i++ )
    {
      if( namedNodeNames_[i] == "center" ) continue;
      if( namedNodes_[i].GetSize() < dim_ ) continue;

      pointsToNamedNodesMap.Push_back(i);
      const StdVector<UInt>& nodes = namedNodes_[i];

      if( dim_ == 2 )
      {
        // In 2D the vertices are calculated as intersections of the
        // lines defined by the first and the last node in a group of
        // named nodes with the same name.
        GetNodeCoordinate(p, nodes[0], false);
        points.Push_back(p);
        GetNodeCoordinate(p, nodes[nodes.GetSize()-1], false);
        points.Push_back(p);
      } else {
        // In 3D the vertices are calculated as intersections of the
        // planes defined by the first, the last node and a third one
        // in a group of named nodes with the same name.
        GetNodeCoordinate(p, nodes[0], false);
        points.Push_back(p);
        GetNodeCoordinate(p, nodes[round(sqrt(nodes.GetSize()))-1], false);
        points.Push_back(p);
        GetNodeCoordinate(p, nodes[nodes.GetSize()-1], false);
        points.Push_back(p);
      }
    }

    // Make sure we have enough points
    if( (dim_ == 2 && points.GetSize() < 4) || (dim_ == 3 && points.GetSize() < 9) ) {
      return -1;
    }

    // Calculate vertices
    Double det, det1, det2;
    if( dim_ == 2 ) {
      // We apply a line-line intersection algorithm using determinants.
      for( UInt i=0; i < points.GetSize()/2.0; i++ )
      {
        for( UInt j=i+1; j < points.GetSize()/2.0; j++ )
        {
          det = (points[2*i][0] - points[2*i+1][0]) * (points[2*j][1] - points[2*j+1][1]) - (points[2*j][0] - points[2*j+1][0]) * (points[2*i][1] - points[2*i+1][1]);
          if( abs(det) > 1e-10 ) {
            det1 = points[2*i][0] * points[2*i+1][1] - points[2*i+1][0] * points[2*i][1];
            det2 = points[2*j][0] * points[2*j+1][1] - points[2*j+1][0] * points[2*j][1];
            r[0] = ( det1 * (points[2*j][0] - points[2*j+1][0]) - det2 * (points[2*i][0] - points[2*i+1][0]) ) / det;
            r[1] = ( det1 * (points[2*j][1] - points[2*j+1][1]) - det2 * (points[2*i][1] - points[2*i+1][1]) ) / det;
            verts.Push_back(r);
          }
        }
      }

      // Calculate the volume as absolute value of the determinant
      // of the matrix, which is defined by two vectors connecting
      // vertices and thus spanning a parallelogram.
      // Note: Due to Cavalieri's principle it does not matter if on of
      // this vectors is a diagonal of our original parallelogram
      p = verts[1] - verts[0];
      q = verts[2] - verts[0];
      det = p[0] * q[1] - p[1] * q[0];

      LOG_DBG(gridcfs) << "CGV: 2D Volume of sparse and/or non rectangular mesh: " << det;
      return std::abs(det);
    } else {
      Vector<Double> n(dim_), n1(dim_), n2(dim_), n3(dim_), p1(dim_), p2(dim_), p3(dim_);
      StdVector<Vector<Double> > normals;
      for( UInt i=0; i < points.GetSize()/3.0; i++ )
      {
        // Calculate the normal of each plane. If the normal is zero the
        // three points lie on a straight line. We then replace the
        // second point.
        Double radicand = 0;
        UInt j = 2;
        UInt k = pointsToNamedNodesMap[i];
        while( radicand < 1e-8 && j < namedNodes_[k].GetSize()-1 ) {
          p = points[3*i+1] - points[3*i];
          q = points[3*i+2] - points[3*i];
          // normal vector
          n[0] = p[1]*q[2] - p[2]*q[1];
          n[1] = p[2]*q[0] - p[0]*q[2];
          n[2] = p[0]*q[1] - p[1]*q[0];
	  
          LOG_DBG3(gridcfs) << "CGV: p=" << p.ToString() << " q=" << q.ToString();
	  
          radicand = pow(n[0],2) + pow(n[1],2) + pow(n[2],2);

          GetNodeCoordinate(points[3*i+1], namedNodes_[k][j], false);
          j++;
        }
        if (radicand < 1e-8) {
          EXCEPTION("Could not calculate normal of bounding plane. The nodes of "
                    << namedNodeNames_[k] << " seem to lie on a straight line "
                    << "and thus do not define a plane.");
        }
        // Normalize normal
        n = n / sqrt( radicand );
        normals.Push_back(n);
      }

      // Calculate the intersection point of three planes
      // Let the planes be specified by a point $p_i$ and a unit normal vector $n_i$.
      // Then the unique point of intersection is given by
      // x = |n1 n2 n3|^(-1) [ (x1 * n1) (n2 x n3) + (x2 * n2) (n3 x n1) + (x3 * n3) (n1 x n2) ]
      Double fac1, fac2, fac3;
      Vector<Double> cross1(dim_), cross2(dim_), cross3(dim_);
      for( UInt i=0; i < normals.GetSize(); i++ )
      {
        n1 = normals[i];
        p1 = points[3*i];
        for( UInt j=i+1; j < normals.GetSize(); j++ )
        {
          n2 = normals[j];
          p2 = points[3*j];
          for( UInt k=j+1; k < normals.GetSize(); k++ )
          {
            n3 = normals[k];
            p3 = points[3*k];
            det = n1[0]*n2[1]*n3[2] + n2[0]*n3[1]*n1[2] + n3[0]*n1[1]*n2[2] - n1[0]*n3[1]*n2[2] - n2[0]*n1[1]*n3[2] - n3[0]*n2[1]*n1[2];
            if( det != 0 ) {
              fac1 = p1[0]*n1[0] + p1[1]*n1[1] + p1[2]*n1[2];
              cross1[0] = n2[1]*n3[2] - n2[2]*n3[1];
              cross1[1] = n2[2]*n3[0] - n2[0]*n3[2];
              cross1[2] = n2[0]*n3[1] - n2[1]*n3[0];
              fac2 = p2[0]*n2[0] + p2[1]*n2[1] + p2[2]*n2[2];
              cross2[0] = n3[1]*n1[2] - n3[2]*n1[1];
              cross2[1] = n3[2]*n1[0] - n3[0]*n1[2];
              cross2[2] = n3[0]*n1[1] - n3[1]*n1[0];
              fac3 = p3[0]*n3[0] + p3[1]*n3[1] + p3[2]*n3[2];
              cross3[0] = n1[1]*n2[2] - n1[2]*n2[1];
              cross3[1] = n1[2]*n2[0] - n1[0]*n2[2];
              cross3[2] = n1[0]*n2[1] - n1[1]*n2[0];

              r[0] = (fac1 * cross1[0] + fac2 * cross2[0] + fac3 * cross3[0]) / det;
              r[1] = (fac1 * cross1[1] + fac2 * cross2[1] + fac3 * cross3[1]) / det;
              r[2] = (fac1 * cross1[2] + fac2 * cross2[2] + fac3 * cross3[2]) / det;
              verts.Push_back(r);
            }
          }
        }
      }
      assert(!verts.IsEmpty());
      // We have to find three linearly independent vectors
      det = 0.0;
      for( UInt i=verts.GetSize()-1; i > 2 && abs(det) < 1e-10 ; i-- )
      {
        p = verts[1] - verts[0];
        q = verts[2] - verts[0];
        r = verts[i] - verts[0];
        det = p[0]*q[1]*r[2] + q[0]*r[1]*p[2] + r[0]*p[1]*q[2] - p[0]*r[1]*q[2] - q[0]*p[1]*r[2] - r[0]*q[1]*p[2];
        LOG_DBG2(gridcfs) << "CGV: " << i << " det=" << det << " p=" << p << " q=" << q << " r=" << r;
      }

      LOG_DBG(gridcfs) << "3D Volume of sparse and/or non rectangular mesh: " << det;
      assert(det != 0);
      return std::abs(det);
    }
    assert(false);
    return -1;
  }

  double GridCFS::CalcVolumeOfRegion(const RegionIdType regionId, bool updated)
  {
    StdVector<Elem*> elems;
    GetElems(elems,regionId);

    double volume = 0.0;

    // according to Jens, NACS does not use the scaling with depth_ in this function
    // only CalcVolumeOfEntityList uses CalcVolume(true) in GridNACS
    for(unsigned int i = 0, n = elems.GetSize(); i < n; i++ )
      volume += GetElemShapeMap(elems[i], updated)->CalcVolume();

    return volume;
  }
  
  Double GridCFS::CalcVolumeOfEntityList( shared_ptr<EntityList> ent,
                                          bool updated ) {
    Double volume = 0.0;
    // get elements of entity list
    if( ent->GetType() == EntityList::ELEM_LIST ||
        ent->GetType() == EntityList::SURF_ELEM_LIST ) {
      EntityIterator it = ent->GetIterator();
      
      // loop over all elements
      for( ; !it.IsEnd(); it++ ) {
        
        const Elem * ptEl = it.GetElem();
        
        
        shared_ptr<ElemShapeMap> esm = GetElemShapeMap( ptEl, updated );
        // sum up element contribution
        // enable scaling with depth_ for 2d plane case as it is done in NACS
        volume += esm->CalcVolume(true);
      }
    } else {
      EXCEPTION( "CalcVolumeOfEntityList only possible for element "
          << "and surface element list" );
    }
    return volume;
  }
  
  void GridCFS::CalcBoundingBoxOfRegion (const RegionIdType regId,
                                         Matrix<Double> & minMax,
                                         CoordSystem* cSys){
    if(!cSys)
      cSys = domain->GetCoordSystem();

    assert(cSys != nullptr);
    StdVector<Elem*> elemssd;

    minMax.Resize(dim_,2);
    Double largeVal = 1e33;
     for(UInt i = 0; i < dim_; i++ ) {
       minMax[i][0] =   largeVal;
       minMax[i][1] = - largeVal;
     }

    GetElems(elemssd, regId );

    for (UInt actEl=0; actEl< elemssd.GetSize(); ++actEl) {
      StdVector<UInt> & connecth = elemssd[actEl]->connect;

      Matrix<Double> ptCoord;
      GetElemNodesCoord(ptCoord, connecth,  false );

      Vector<Double> globCoord(dim_), locCoord(dim_);
      // loop over nodes
      for (UInt i=0; i< ptCoord.GetNumCols(); i++) {

        // convert global coordinate to local coordinate
        for( UInt iDim = 0; iDim < dim_ ; ++iDim )  {
          globCoord[iDim] = ptCoord[iDim][i];
        }
        cSys->Global2LocalCoord(locCoord, globCoord);

        // determine min / max of propagation region
        for( UInt iDim = 0; iDim < dim_ ; ++iDim )  {
          if ( locCoord[iDim] < minMax[iDim][0] )
            minMax[iDim][0] = locCoord[iDim];
          if ( locCoord[iDim] > minMax[iDim][1] )
            minMax[iDim][1] = locCoord[iDim];

        }
      } // loop over nodes
    }
  }


  void GridCFS::AddNode( const Vector<Double> & coord, UInt & inode )
  {
    if(!isInitialized_)
      EXCEPTION("Cannot add node to uninitialized grid!");

    if(coord.GetSize() != 3)
      EXCEPTION("Node to be added has wrong dimension!");
#pragma omp critical (GridCFS)
{
    coords_.Push_back(coord);
    inode = ++numNodes_;

    if (deltCoords_.GetSize() > 0) {
      Vector<Double> zero(dim_);
      zero.Init();
      deltCoords_.Push_back(zero);
    }
}
  }

  void GridCFS::AddNodes( const StdVector< Vector<Double> > & coords,
                          StdVector< UInt > & inodes)
  {
    if(!isInitialized_)
      EXCEPTION("Cannot add nodes to uninitialized grid!");
#pragma omp critical (GridCFS)
{
    UInt i, n;

    n=coords.GetSize();
    inodes.Resize(n);

    for(i=0; i<n; i++)
    {
      coords_.Push_back(coords[i]);
      numNodes_++;
      inodes[i] = numNodes_;
    }

    if (deltCoords_.GetSize() > 0) {
      Vector<Double> zero(dim_);
      zero.Init();
      for (i = 0; i < n; ++i)
        deltCoords_.Push_back(zero);
    }
}
  }

  void GridCFS::AddSurfaceElems( const RegionIdType regionid,
                                 const StdVector< SurfElem* > & surfelems,
                                 StdVector< UInt > & elemids)
  {
    UInt i, n;
    UInt numNodes;

    if(!isInitialized_)
      EXCEPTION("Cannot add surface elements to uninitialized grid!");

    Integer regionIdx = surfRegionIds_.Find(regionid);

    if(regionIdx == -1)
      EXCEPTION("Surface regionid not found!");

    n=surfelems.GetSize();
    elemids.Resize(n);

    StdVector<UInt> surfRegionNodes;
    GetNodesByRegion( surfRegionNodes, regionid);
    for(i=0; i<n; i++)
    {
      // TODO: a check should be added to avoid insertions
      // of already existing elements
      surfelems[i]->regionId = regionid;
      numElems_++;
      surfelems[i]->elemNum = numElems_;

      orderedElems_.Push_back(surfelems[i]);
      surfElems_[regionIdx].Push_back(surfelems[i]);
      elemids[i] = numElems_;

      numNodes = surfelems[i]->connect.GetSize();

      // Loop over all nodes an check, if they are already contained in
      // the list of nodes
      const StdVector<UInt> & connect = surfelems[i]->connect;
      std::set<UInt> s;
      for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
        s.insert(connect[iNode]);
       // if( surfRegionNodes.Find(connect[iNode]) == -1 ) {
        //surfRegionNodes.Push_back(connect[iNode]);
       // }
      }
      std::set<UInt>::iterator iter = s.begin();
      for(;iter!=s.end();iter++){
        surfRegionNodes.Push_back(*iter);
      }
    }
    
    // in the end store back the number of surface element nodes
    numSurfElemNodes_[regionIdx] = surfRegionNodes.GetSize();
  }

  void GridCFS::AddVolumeElems( const RegionIdType regionid,
                                const StdVector< Elem* > & volelems,
                                StdVector< UInt > & elemids)
  {
    UInt i, n;
    UInt numNodes;

    if(!isInitialized_)
      EXCEPTION("Cannot add volume elements to uninitialized grid!");

    Integer regionIdx = volRegionIds_.Find(regionid);

    if(regionIdx == -1)
      EXCEPTION("Volume regionid not found!");


    n=volelems.GetSize();
    elemids.Resize(n);

    StdVector<UInt> volRegionNodes;
    GetNodesByRegion( volRegionNodes, regionid);

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

      numNodes = volelems[i]->connect.GetSize();

      // Loop over all nodes an check, if they are already contained in
      // the list of nodes
      const StdVector<UInt> & connect = volelems[i]->connect;
      for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
        if( volRegionNodes.Find(connect[iNode]) == -1 ) {
          volRegionNodes.Push_back(connect[iNode]);
        }
      }
    }
    // in the end store back the number of surface element nodes
    numVolElemNodes_[regionIdx] = volRegionNodes.GetSize();
  }


  void GridCFS::ClearRegion( const RegionIdType regionid )
  {
    StdVector<Elem*> newOrderedElems;
    StdVector<Elem*> elems;
    UInt numElems;
    UInt i, n;

    // look in volume regions
    Integer index = volRegionIds_.Find(regionid);
    if ( index != -1 ) {
      n = volElems_[index].GetSize();

      for(i=0; i<n; i++)
      {
        orderedElems_[volElems_[index][i]->elemNum-1] = nullptr;
        delete volElems_[index][i];
      }

      volElems_[index].Clear();
      numVolElemNodes_[index] = 0;
    } else {
      // look in surface regions
      index = surfRegionIds_.Find(regionid);
      if ( index != -1 ) {
        n = surfElems_[index].GetSize();

        for(i=0; i<n; i++)
        {
          orderedElems_[surfElems_[index][i]->elemNum-1] = nullptr;
          delete surfElems_[index][i];
        }

        surfElems_[index].Clear();
        numSurfElemNodes_[index] = 0;
      } else {
        EXCEPTION("GridCFS: The region with id '" << regionid
                  << "' was not found in the grid!");
      }
    }


    numElems = 0;
    for(i=0; i<numElems_; i++)
    {
      if(orderedElems_[i] != nullptr)
      {
        newOrderedElems.Push_back(orderedElems_[i]);
        //        std::cout << "Clear Region: " << orderedElems_[i]->elemNum << " -> " << numElems << std::endl;
        orderedElems_[i]->elemNum = ++numElems;
      }
    }

    orderedElems_ = newOrderedElems;
    numElems_ = numElems;
  }

  void GridCFS::DeleteNamedNodes(const std::string &name) {
    // check if the entry in the nameTypeMap is set
    if (nameTypeMap_.find(name) == nameTypeMap_.end()) {
      EXCEPTION("Node list '" << name << "' does not exist.");
    }
    // check if the specified node list exists
    Integer idx = namedNodeNames_.Find(name);
    if (idx != -1) {
      std::set<UInt> sortedNodes;
      // put the nodes into a std::set to get an ordered list
      sortedNodes.insert(namedNodes_[idx].Begin(), namedNodes_[idx].End());

      if (*sortedNodes.begin() == coords_.GetSize()-sortedNodes.size()+1) {
        // is a block at the end of the vector, so just resize
        numNodes_ = coords_.GetSize() - sortedNodes.size();
        coords_.Resize(numNodes_);
        deltCoords_.Resize(numNodes_);
      } else {
        // Please implement the general case if you need it. But beware!
        // It's complicated, because you need to renumber all nodes and
        // therefore change the connectivity as well!
        EXCEPTION("You try to delete named nodes which are not a continuous block at the end of the node array." << std::endl <<
                  "This is not allowed as it would change the mesh connectivity!" << std::endl <<
                  "Check the declaration of your ncInterfaces! ()");
      }
      namedNodeNames_.Erase( (UInt) idx );
      namedNodes_.Erase( (UInt) idx );
      nameTypeMap_.erase(name);
    }
    else {
      EXCEPTION("Cannot delete '" << name << "': not a node list");
    }
  }
  
  void GridCFS::CorrectElementConnectivities(RegionIdType regionId) {
    // define variables
    Matrix<Double> jacobian;
    std::set<const Elem*> corrElems, failedElems;
    Double jacDet = 0;
    UInt numElems;
    StdVector<Elem*> elems;
    // if no region id is passed, check all elements on the mesh.
    // otherwise, check only the specified region elements.
    if (regionId == -1) {
      numElems = numElems_;
      elems = orderedElems_;
    } else {
      GetElems(elems, regionId);
      numElems = elems.GetSize();
    }
    for (UInt iElems = 0; iElems < numElems; ++iElems) {
      Elem* el = elems[iElems];
      shared_ptr<ElemShapeMap> esm = GetElemShapeMap( el, false );

      jacDet = esm->CalcJDet( jacobian, Elem::shapes[el->type].midPointCoord);
      if( jacDet < 0 ) {
        try {
        el->CorrectConnectivity(*this);
        // at this point, we can be sure that the element connectivity
        // was adjusted correctly
        corrElems.insert(el);
        } catch (Exception& ex) {
          // at this point, the correction failed e.g. due to a missing
          // implementation or a totally weird element connectivity.
          //before we give up, lets try a brute force attack
          WARN("Trying to correct connectivity by permutating array. This can be costly! Recheck the mesh!");
          bool success = false;
          do{
            shared_ptr<ElemShapeMap> esmTMP = GetElemShapeMap( el, false );
            jacDet = esmTMP->CalcJDet( jacobian, Elem::shapes[el->type].midPointCoord);
            if(jacDet > 0){
              success = true;
              break;
            }
          }while( std::next_permutation(el->connect.Begin(),el->connect.End()) );
          if(!success){
            failedElems.insert(el);
          }else{
            corrElems.insert(el);
          }
        }
        //we need to reset the element in case of precached maps
        esm->SetElem(el,false);
      }
    }    
    // if some elements were successfully reoriented, issue warning
    if(corrElems.size() > 0 ) {
      std::stringstream out;
      if(corrElems.size() > 10){
        out << "A total number of " << corrElems.size() << " elements"
            << " had a wrong orientation and were reoriented."
            << "\n Usually this does not lead to further errors.";
      }else{
        out << "The following elements have a wrong orientation and "
            << "were re-oriented:\n";
        std::set<const Elem*>::iterator it = corrElems.begin();
        for( ; it != corrElems.end(); it++  ) {
          out << (*it)->elemNum << ", ";
        }
        out << "\n\nPlease check your mesh!\n";
      }
        WARN( out.str().c_str() );
    }

    // if some elements could not be reoriented, generate exception
    if(failedElems.size() > 0 ) {
      std::stringstream out;
      out << "The following elements have a wrong orientation and "
          << "could NOT be re-oriented:\n";
      std::set<const Elem*>::iterator it = failedElems.begin();
      for( ; it != failedElems.end(); it++  ) {
        out << (*it)->elemNum << ", ";
      }
      out << "\n\nPlease check your mesh!\n";
      EXCEPTION( out.str() );
    }
  }

  void GridCFS::makeNameNodesFromLines()
  {
    if(!param_ || !param_->Has("domain/surfRegionList")) return;

    ParamNodeList list = param_->Get("domain/surfRegionList")->GetList("surfRegion");
    std::map<std::string, std::string> excludeSurf;
    StdVector<UInt> nodeList;
    StdVector<UInt> nodeListTmp1, nodeListTmp2;
    for(UInt i = 0; i < list.GetSize(); i++)
    {
      if (list[i]->Get("makeNamedNodes")->As<bool>())
      {
        const std::string nameTmp = list[i]->Get("name")->As<std::string>();
        excludeSurf[nameTmp] = list[i]->Get("excludeSurface")->As<std::string>();
      }
    }
    std::map<std::string, std::string>::const_iterator namedIter = excludeSurf.begin();
    for ( ; namedIter != excludeSurf.end(); ++namedIter)
    {
      const std::string& surfName = namedIter->first;
      GetNodesByRegion(nodeList, region_.Parse(namedIter->first));
      std::string nodeRegName = "nodes_" + surfName;
      if (excludeSurf[surfName] != "")
      {
        UInt iTmp = region_.Parse(excludeSurf[surfName]);
        GetNodesByRegion(nodeListTmp2, surfRegionIds_[iTmp]);
        for (UInt iTmp1 = 0; iTmp1 < nodeList.GetSize(); ++iTmp1)
        {
          const UInt& currNode = nodeList[iTmp1];
          if (nodeListTmp2.Find(currNode) == -1)
            nodeListTmp1.Push_back(currNode);
        }
        nodeList = nodeListTmp1;
        nodeListTmp1.Clear();
      }
      AddNamedNodes(nodeRegName, nodeList);
    }
  }
  

  void GridCFS::MapMidSideNodes() {


    // Leave if grid consists just of linear order elements
    if( !isQuadratic_)
      return;

    // Loop over all elements
    UInt numElems = orderedElems_.GetSize();
    for( UInt iEl = 0; iEl < numElems; ++iEl ) { 
      const Elem * elem = orderedElems_[iEl];
      const ElemShape & sh = Elem::shapes[elem->type];
      
      // Loop over all mid-side nodes of element
      for( UInt iNode = sh.numVertices; iNode < sh.numNodes; ++iNode ) {
        const UInt nodeNum = elem->connect[iNode];
        // Obtain local coordinate from ElemShape and store it
        const LocPoint & lp = LocPoint(sh.nodeCoords[iNode]);
        std::pair<const Elem*, LocPoint> entry(elem, lp);
        midNodeProjections_[nodeNum].Push_back(entry);
      } // loop: nodes
    } // loop: elements
    
    // In the end free any non-needed memory
    boost::unordered_map<UInt, NodeElemMatch>::iterator it;
    it = midNodeProjections_.begin();
    for( ; it != midNodeProjections_.end(); ++it ) {
      it->second.Trim();
    }
  }

  void GridCFS::ExportGrid(PtrParamNode out)
  {
    PtrParamNode nl = out->Get("nodeList");
    // Setup large array to resolve zero to many node names
    StdVector<StdVector<unsigned int> > node_names(coords_.GetSize() + 1); // 1-based!
    for(unsigned int s = 0, sn = namedNodeNames_.GetSize(); s < sn; s++)
    {
      const StdVector<unsigned int>& nn = namedNodes_[s];
      for(unsigned int n = 0; n < nn.GetSize(); n++)
        node_names[nn[n]-1].Push_back(s); // fuck mixed 1-based and 0-based :(
    }

    StdVector<std::string>& block = nl->GetFastBulkBlock();
    block.Resize(coords_.GetSize());
    for(unsigned int n = 0, nn = coords_.GetSize(); n < nn; n++)
    {
      std::stringstream ss;

      ss << "<node id=\"" << (n + 1)
         << "\" x=\"" << coords_[n][0]
         << "\" y=\"" << coords_[n][1]
         << "\" z=\"" << (dim_ > 2 ? coords_[n][2] : 0.0) << "\"";

      const StdVector<unsigned int>& ni = node_names[n];
      ss << " names=\"" << ni.GetSize() << "\"";
      for(unsigned int c = 0; c < ni.GetSize(); c++)
        ss << " name_" << c << "=\"" << namedNodeNames_[ni[c]] << "\"";
      ss << "/>";
      block[n] = ss.str();
    }
    node_names.Clear();
    // same game for the element names as for the nodes
    StdVector<StdVector<unsigned int> > elem_names(GetNumElems() + 1); // 1-based!
    for(unsigned int s = 0, sn = namedElemNames_.GetSize(); s < sn; s++)
    {
      const StdVector<unsigned int>& nn = namedElems_[s];
      for(unsigned int n = 0; n < nn.GetSize(); n++)
        elem_names[nn[n]-1].Push_back(s);
    }

    PtrParamNode rl = out->Get("regionList");
    for(unsigned int r = 0; r < regionData.GetSize(); r++)
    {
      RegionData& rd = regionData[r];
      SetElementBarycenters(rd.id, false);
      PtrParamNode reg = rl->Get("region", ParamNode::APPEND);
      reg->Get("name")->SetValue(rd.name);
      const StdVector<Elem*>& elems = rd.type == VOLUME_REGION ? volElems_[rd.type_idx] : surfElems_[rd.type_idx];

      StdVector<std::string>& block = reg->GetFastBulkBlock();
      block.Resize(elems.GetSize());
      for(unsigned int e = 0, n = elems.GetSize(); e < n; e++)
      {
        const Elem* elem = elems[e];
        const Point& bc = elem->extended->barycenter;
        const StdVector<unsigned int>& con = elem->connect;
        const StdVector<unsigned int>& ni = elem_names[n];

        std::stringstream ss;
        ss << "<element id=\"" << (e + 1)
           << "\" type=\"" << Elem::feType.ToString(elem->type)
           << "\" x=\"" << bc.data[0]
           << "\" y=\"" << bc.data[1]
           << "\" z=\"" << bc.data[2]
           << "\" nodes=\"" << con.GetSize() << "\"";

        for(unsigned int c = 0; c < con.GetSize(); c++)
          ss << " node_" << c << "=\"" << con[c] << "\"";

        ss << " name=\"" << ni.GetSize() << "\"";
        for(unsigned int c = 0; c < ni.GetSize(); c++)
          ss << " name_" << c << "=\"" << namedElemNames_[ni[c]] << "\"";
        ss << "/>";
        block[e] = ss.str();
      }
    }
  }

} // end namespace
