// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// This reader reads meshes created with Gmsh (http://geuz.org/gmsh/) in
// the formats msh1 (legacy), msh2 (ascii), msh2 (binary, with little
// and big endian formats).
//
// Example for XML file:
//
// <input>
//   <gmsh fileName="box_msh1.msh" coordSysId="myCSys" scaleFac="0.5">
//     <region name="air" physicalEntity="12456" linearize="true"/>
//   </gmsh>
// </input>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#include <boost/algorithm/string/trim.hpp>
#include <filesystem>

#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Domain/Domain.hh"

#include "GmshHelper.hh"
#include "SimInputGmsh.hh"

namespace CoupledField {
    
  SimInputGmsh::SimInputGmsh(std::string fileName, PtrParamNode inputNode,
                             PtrParamNode infoNode) :
    SimInput(fileName, inputNode, infoNode),
    numNodesPerElem_(27),
    dim_(0),
    coordSysId_("default"),
    scaleFac_(1.0)
  {
    GmshHelper::InitElemNodeMap();

    PtrParamNode pNode;
    ParamNodeList pNodeList;

    pNode = myParam_->Get("coordSysId", ParamNode::PASS);
    if(pNode) {
      coordSysId_ = myParam_->Get("coordSysId", ParamNode::PASS)->As<std::string>();
    }
    
    pNode = myParam_->Get("scaleFac", ParamNode::PASS);
    if(pNode) {
      scaleFac_ = myParam_->Get("scaleFac", ParamNode::PASS)->As<Double>();
    }

    readOnlySomeRegions_ = false;
    pNodeList = myParam_->GetList("region");
    for(UInt i=0, n= pNodeList.GetSize(); i<n; i++)
    {
      std::string regionName =
        pNodeList[i]->Get("name", ParamNode::PASS)->As<std::string>();

      bool linearize = false;
      pNodeList[i]->GetValue("linearize", linearize, ParamNode::PASS);
      
      Integer physEntity = 0;
      pNodeList[i]->GetValue("physicalEntity", physEntity, ParamNode::PASS);

      physEntities2RegionNames_[physEntity] = regionName;
      linearizeRegions_[regionName] = linearize;
      
      //      std::cout << "region " << regionName << " linearize " << linearize << " physicalEntity " << physEntity << std::endl;
      readOnlySomeRegions_ = true;
    }



  }

  SimInputGmsh::~SimInputGmsh()
  {
  }

  void SimInputGmsh::GetPhysEnt2RegionMapFromXML() 
  {
  }

  void SimInputGmsh::InitModule()
  {
    // Tetrahedron:
    //                   edge 1: nodes 1 -> 2                
    //    v                   2:       1 -> 3                                            
    //    |                   3:       1 -> 4                
    //    |                   4:       2 -> 3                
    //    |                   5:       2 -> 4                
    //    3                   6:       3 -> 4                
    //    |\
    //    | \            face 1: edges  1 -3  5  nodes 1 2 4
    //    |__\2_____u         2:       -1  2 -4        1 3 2
    //    1\ /                3:       -2  3 -6        1 4 3
    //      \4                4:        4 -5  6        2 3 4

    // Indices
    UInt i,j;
    char line[256];

    // Converters for endianess
    EndianSwapper<UInt> uiswap;
    EndianSwapper<Integer> iswap;
    EndianSwapper<Double> dswap;

    Double versionNumber = 1;
    UInt fileType = 0;

    std::ifstream in(fileName_.c_str(), std::ios::binary);
    if ( !std::filesystem::exists( fileName_ ))
    {
      EXCEPTION("Input file does not exist:\n" << fileName_ << std::endl);
    } else if(!in.good()) {
      EXCEPTION("Could not open '" << fileName_ << "'!");
    }

    while( in.good() ) {
      in.getline(line, sizeof(line), '\n');
      std::string sline(line);
      boost::trim(sline);

      // std::cout << line << std::endl;

      if( sline == "$MeshFormat" ) 
      {
        // Variable to determine endianess
        UInt oneBinary = 1;
        UInt dataSize = 0;

        // std::cout << "Starting $MeshFormat Block" << std::endl;
        in >> versionNumber;
        unsigned majorVer = static_cast<unsigned>(versionNumber);
        if(majorVer > 2) 
        {
          EXCEPTION("Gmsh input files with version " << versionNumber
                    << " are not supported at this time!\n" <<
                    "The Gmsh reader currently supports up to version 2.")
        }
        in >> fileType;
        in >> dataSize;

        
        if(fileType == 1)
        {
          // first read newline
          in.read((char*) &oneBinary, 1);
          in.read((char*) &oneBinary, sizeof(oneBinary));
        }
        
        //        std::cout << versionNumber << " " << fileType << " "  << dataSize << std::endl;

        if(oneBinary != 1)
        {
          uiswap.swap = true;
          iswap.swap = true;
          dswap.swap = true;
        }        

        continue;
      }


      if( sline == "$EndMeshFormat" ) 
      {
        //        std::cout << "Ending $MeshFormat Block" << std::endl;
      }

      if( sline == "$NOD" ||
          sline == "$NOE" ||
          sline == "$Nodes") {
        UInt numNodes = 0;
        int id;
        in >> numNodes;

        //        std::cout << "reading " << numNodes << " nodes:" << std::endl;
        int idx = nodeCoords_.GetSize() / 3;
        nodeCoords_.Resize( nodeCoords_.GetSize() + numNodes * 3 );

        if(fileType == 1) 
        {
          // Read newlines until real data
          in.read((char*)(&id), 1);

          // Read node data as one block
          std::vector<NodeStruct> dummy(numNodes);

          for(UInt i=0;i<numNodes;i++) {
            in.read((char*)(&dummy[i].nodeId), sizeof(dummy[i].nodeId));
            in.read((char*)(&dummy[i].x), sizeof(dummy[i].x));
            in.read((char*)(&dummy[i].y), sizeof(dummy[i].y));
            in.read((char*)(&dummy[i].z), sizeof(dummy[i].z));

            nodeCoords_[idx*3+0] = dswap.EndianSwapBytes(dummy[i].x);
            nodeCoords_[idx*3+1] = dswap.EndianSwapBytes(dummy[i].y);
            nodeCoords_[idx*3+2] = dswap.EndianSwapBytes(dummy[i].z);
            
            nodeNumMap_[uiswap.EndianSwapBytes(dummy[i].nodeId)] = idx+1;

            idx++;
          }
        }
        else 
        {  
          double x,y,z;

          for(UInt i=0;i<numNodes;i++) {
            in >> id >>  x >> y >> z;
            
            //            std::cout << x << " " << y << " "  << z << std::endl;
            
            nodeCoords_[idx*3+0] = x;
            nodeCoords_[idx*3+1] = y;
            nodeCoords_[idx*3+2] = z;
            
            nodeNumMap_[id] = idx+1;
            
            idx++;
          }
        }
      }        


      if( sline == "$EndNodes" ) 
      {
        //        std::cout << "Ending $Nodes Block" << std::endl;
      }

      if( sline == "$ELM" || sline == "$Elements") {
        // read elements
        UInt numElements = 0;
        in >> numElements;
        //        std::cout << ".msh reading elements" << std::endl << numElements << std::endl;
        
        UInt idx = elementTypes_.size();
        
        // If no physical entity to region name map has been given in XML
        // generate default region names.
        bool genDefaultRegions = physEntities2RegionNames_.empty();        

        elementTypes_.resize(idx + numElements);
        elementPhysicsTypes_.resize(idx + numElements);
        connectivity_.resize(idx*numNodesPerElem_ + numElements*numNodesPerElem_);
        std::fill(&connectivity_[idx*numNodesPerElem_],
                  &connectivity_[idx*numNodesPerElem_ + numElements*numNodesPerElem_-1],
                  0);
          
        UInt id, type, regPhys, regElem, numElementNodes, numTags;
        Elem::FEType feType;
        std::vector<UInt> tags;

        // Read newlines until real data
        if(fileType == 1) 
        {
          in.read((char*)(&id), 1);
        }
        
        // Read newlines until real data
        if(fileType == 1) 
        {
          UInt elemCounter = 0;
          UInt numFollowingElements = 0;
          std::vector<UInt> dummy;
          UInt elemRecordLength = 0;
          
          while(elemCounter < numElements) 
          {
            // Read header
            in.read((char*)(&type), sizeof(type));
            in.read((char*)(&numFollowingElements), sizeof(numFollowingElements));
            in.read((char*)(&numTags), sizeof(numTags));

            type = uiswap.EndianSwapBytes(type);
            numFollowingElements = uiswap.EndianSwapBytes(numFollowingElements);
            numTags = uiswap.EndianSwapBytes(numTags);

            GmshHelper::ElemType2FEType(type, feType, numElementNodes);

            UInt feDim = Elem::shapes[feType].dim;
            dim_ = dim_ > feDim ? dim_ : feDim;

            elemRecordLength = 1 + numTags + numElementNodes;
            dummy.resize(numFollowingElements * elemRecordLength);
            
            in.read((char*)&dummy[0], dummy.size()*sizeof(dummy[0]));

//             std::cout << "idx " << idx  << " elemRecordLength " << elemRecordLength << std::endl;
//             std::cout << "ELEMHEAD numFollowingElements " << numFollowingElements << " numTags " << numTags << " numElementNodes " << numElementNodes  << std::endl;

            for( i=0; i<numFollowingElements; i++ ) {
              elementTypes_[idx] = feType;
              elementPhysicsTypes_[idx] = uiswap.EndianSwapBytes( dummy[i*elemRecordLength + 1] );

              if( genDefaultRegions &&
                  physEntities2RegionNames_.find(elementPhysicsTypes_[idx]) ==
                  physEntities2RegionNames_.end()) 
              {
                // Add a new default region name
                std::stringstream sstr;
                
                sstr << "physical_entity_" << elementPhysicsTypes_[idx];
                physEntities2RegionNames_[elementPhysicsTypes_[idx]] = sstr.str();
              }

//               std::copy(&dummy[(i+1)*elemRecordLength-numElementNodes],
//                         &dummy[(i+1)*elemRecordLength],
//                         &connectivity_[idx*numNodesPerElem_]);

              typedef GmshHelper::ElemNodeMapType::left_map::const_iterator left_const_iterator;

              // Iterate over nodes element and reorder connectivity
              for( left_const_iterator left_iter = GmshHelper::elemNodeMap_[feType].left.begin(), iend = GmshHelper::elemNodeMap_[feType].left.end();
                   left_iter != iend;
                   ++left_iter)
              {
                UInt idx2 = idx*numNodesPerElem_ + left_iter->second;
                connectivity_[idx2] = uiswap.EndianSwapBytes(dummy[(i+1)*elemRecordLength-numElementNodes + left_iter->first]);
                //                connectivity_[idx2] = nodeNumMap_[connectivity_[idx2]];

                // std::cout << left_iter->first << " --> " << left_iter->second << std::endl;
                //                std::cout << connectivity_[idx2] << " --> " << (left_iter->first+1) << std::endl;
              }

              idx++;              
            }

            elemCounter += numFollowingElements;

          }
        }
        else
        {
          for( i=0; i<numElements; i++ ) {
              
              if(versionNumber == 1)
              {
                in >> id >> type >> regPhys >> regElem >> numElementNodes;
                //                std::cout << "id " << id << " type " << type << " regPhys " << regPhys << std::endl;
              }
              else
              {
                in >> id >> type >> numTags;
                tags.resize(numTags);
                
                //                std::cout << "id " << id << " type " << type << " numTags " << numTags << std::endl;
                
                for( j=0; j<numTags; j++ ) {
                  in >> tags[j];                
                }
                
                //                std::cout << "tags[0] " << tags[0] << " tags[1] " << tags[1] << " tags[2] " << tags[2] << std::endl;
                
                // region!
                regPhys = tags[0];
              }
              
              GmshHelper::ElemType2FEType(type, feType, numElementNodes);
              UInt feDim = Elem::shapes[feType].dim;
              dim_ = dim_ > feDim ? dim_ : feDim;
              
              elementTypes_[idx] = feType;
              elementPhysicsTypes_[idx] = regPhys;

              if( genDefaultRegions &&
                  physEntities2RegionNames_.find(regPhys) ==
                  physEntities2RegionNames_.end()) 
              {
                // Add a new default region name
                std::stringstream sstr;
                
                sstr << "physical_entity_" << regPhys;
                physEntities2RegionNames_[regPhys] = sstr.str();
              }
              
              UInt buf[128];
              UInt elemNode;
              for( elemNode = 0; elemNode < numElementNodes; elemNode++)
              {
                in >> buf[elemNode];
              }
              
              typedef GmshHelper::ElemNodeMapType::left_map::const_iterator left_const_iterator;
              // Iterate over nodes element and reorder connectivity
              elemNode=0;
              for( left_const_iterator left_iter = GmshHelper::elemNodeMap_[feType].left.begin(), iend = GmshHelper::elemNodeMap_[feType].left.end();
                   left_iter != iend;
                   ++left_iter, elemNode++)
              {
                UInt idx2 = idx*numNodesPerElem_ + left_iter->second;

                connectivity_[idx2] = buf[elemNode];
                //                connectivity_[idx2] = nodeNumMap_[dummy];

                //                std::cout << dummy << std::endl;
                // std::cout << left_iter->first << " --> " << left_iter->second << std::endl;
                //                std::cout << connectivity_[idx2] << " --> " << (left_iter->first+1) << std::endl;

              }

              idx++;
            }          
          }
      }

      if( sline == "$EndElements" ) 
      {
        //        std::cout << "Ending $Elements Block" << std::endl;
        //        EXCEPTION("Done reading elements");
      }

      if( sline == "$PhysicalNames" ) 
      {
        UInt numPhysicalNames = 0;
        std::stringstream sstr;
        std::map<std::string, UInt> physNamesToDim;
        std::map<std::string, UInt> physNamesToId;
        std::map<std::string, UInt>::const_iterator it, end;
        
        in >> numPhysicalNames;

        for(UInt i=0; i<=numPhysicalNames; i++)
        {
          UInt physDim = 3;
          UInt physId = 0;
          std::string physName;
       
          in.getline(line,sizeof(line));
          // std::cout << line << std::endl;
          sline = line;
          boost::trim(sline);

          if(!sline.length())
            continue;
          
          sstr.clear();
          sstr.str(sline);

          
          if( versionNumber >= 2)
          {
            sstr >> physDim;

            dim_ = dim_ > physDim ? dim_ : physDim;
            
//             if(!sstr.good())
//               continue;
          }
          
          sstr >> physId >> physName;
          //          if(!sstr.good())
          //            continue;
          
          // remove " from physName
          boost::trim_if(physName, boost::is_any_of("\" \t"));

          // std::cout << physDim << " " << physId << " #" << physName  << "#" << std::endl;

          std::map<UInt, std::string>::iterator regIt, regEnd;
          regIt = physEntities2RegionNames_.begin();
          regEnd = physEntities2RegionNames_.end();
          bool readRegion = false;
          
          // Erase wrong physical id- physical name pairs given in XML file and 
          // replace them with correct ones from Gmsh file.
          for( ; regIt != regEnd; ) {
            if(regIt->second == physName) {
              physEntities2RegionNames_.erase(regIt++);
              readRegion = true;
              continue;
            } else {
              ++regIt;
            }
          }
          
          if(readOnlySomeRegions_ && !readRegion)
            continue;

          physNamesToDim[physName] = physDim;
          physNamesToId[physName] = physId;
        }

        it = physNamesToDim.begin();
        end = physNamesToDim.end();
        
        if(dim_ < 2) dim_ = 2;
        
        for( ; it != end; it++) {
          std::string physName = it->first;
          UInt physId = physNamesToId[it->first];
          UInt physDim = it->second;
          
          switch(physDim) {
            case 0:
              physEntities2NamedNodes_[physId] = physName;
              break;
            
            case 1:
              if(dim_ == 2) {
                physEntities2RegionNames_[physId] = physName;
              } else {
                if(physName != "GmshNamedElems") 
                {
                  physEntities2RegionNames_[physId] = "GmshNamedElems";
                  physEntities2NamedElems_[physId] = physName;
                }
              }
              break;
            
            case 2:
            case 3:
              physEntities2RegionNames_[physId] = physName;
              break;
          }
        }
        
        //        std::cout << "Ending $Elements Block" << std::endl;
        //        EXCEPTION("Done reading elements");
      }

      if( sline == "$EndPhysicalNames" )
      {
        //        std::cout << "Ending $Elements Block" << std::endl;
        //        EXCEPTION("Done reading elements");
      }
    }
  }
  

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputGmsh::GetDim() {
    return dim_;
  }
  
  void SimInputGmsh::GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                              std::map<UInt, UInt>& numSteps,
                                              bool isHistory) {
    // At the moment only the grid may be read.
    analysis.clear();
    numSteps.clear();
  }

  void SimInputGmsh::LinearizeElem(Elem::FEType* elemType) {
    static std::map<Elem::FEType, Elem::FEType> elemTypeMap;

    if(!elemTypeMap.size()) {
      elemTypeMap[Elem::ET_POINT]   = Elem::ET_POINT;
      elemTypeMap[Elem::ET_LINE2]   = Elem::ET_LINE2;
      elemTypeMap[Elem::ET_LINE3]   = Elem::ET_LINE2;
      elemTypeMap[Elem::ET_TRIA3]   = Elem::ET_TRIA3;
      elemTypeMap[Elem::ET_TRIA6]   = Elem::ET_TRIA3;
      elemTypeMap[Elem::ET_QUAD4]   = Elem::ET_QUAD4;
      elemTypeMap[Elem::ET_QUAD8]   = Elem::ET_QUAD4;
      elemTypeMap[Elem::ET_QUAD9]   = Elem::ET_QUAD4;
      elemTypeMap[Elem::ET_TET4]    = Elem::ET_TET4;
      elemTypeMap[Elem::ET_TET10]   = Elem::ET_TET4;
      elemTypeMap[Elem::ET_HEXA8]   = Elem::ET_HEXA8;
      elemTypeMap[Elem::ET_HEXA20]  = Elem::ET_HEXA8;
      elemTypeMap[Elem::ET_HEXA27]  = Elem::ET_HEXA8;
      elemTypeMap[Elem::ET_PYRA5]   = Elem::ET_PYRA5;
      elemTypeMap[Elem::ET_PYRA13]  = Elem::ET_PYRA5;
      elemTypeMap[Elem::ET_PYRA14]  = Elem::ET_PYRA5;
      elemTypeMap[Elem::ET_WEDGE6]  = Elem::ET_WEDGE6;
      elemTypeMap[Elem::ET_WEDGE15] = Elem::ET_WEDGE6;
      elemTypeMap[Elem::ET_WEDGE18] = Elem::ET_WEDGE6;
    }

    *elemType = elemTypeMap[*elemType];
  }

  void SimInputGmsh::ReadMesh(Grid *mi)
  {
    UInt idx;
    std::map<UInt, UInt> nodeNumMap2;
    std::vector<UInt> readElemIndices;

    std::map<UInt, RegionIdType> regionIdMap;
    std::map<UInt, std::string>::iterator regIt, regEnd;
    regIt = physEntities2RegionNames_.begin();
    regEnd = physEntities2RegionNames_.end();

    UInt numElements = 0;

    // Go over all elements and check if any complete quadratic elements
    // are present. Gmsh at least until version 2.8.3 has the bug that
    // complete pyramids will be generated, even if incomplete elements
    // have been switched on.
    bool incompleteElems = true;
    for( UInt i=0, ne = elementTypes_.size(); i < ne; i++ ) {
      Elem::FEType feType = elementTypes_[i];
      
      if(feType == Elem::ET_QUAD9 ||
         feType == Elem::ET_WEDGE18 ||
         feType == Elem::ET_HEXA27) 
      {
        incompleteElems = false;
        break;
      }
    }
    
    for( ; regIt != regEnd; regIt++) {
      RegionIdType actRegionId;
      const Enum<RegionIdType> regEnum = mi->GetRegion();
      RegionIdType def = NO_REGION_ID;
      actRegionId = regEnum.Parse(regIt->second, def);
      if(actRegionId == NO_REGION_ID) {
        actRegionId = mi->AddRegion(regIt->second);
      }
      
      regionIdMap[regIt->first] = actRegionId;

      // Go over all elements
      for( UInt i=0, ne = elementTypes_.size(); i < ne; i++ ) {
        // If current element does not belong to any region which
        // needs to read, go on...
        if(elementPhysicsTypes_[i] != regIt->first)
          continue;
        
        numElements++;
        readElemIndices.push_back(i);

        // If elements in this region need to be linearized do so
        if(linearizeRegions_[regIt->second])
          LinearizeElem(&elementTypes_[i]);

        // Make incomplete pyramids out of complete ones, if otherwise
        // only incomplete elements are used.
        if(incompleteElems && (elementTypes_[i] == Elem::ET_PYRA14))
        {
          elementTypes_[i] = Elem::ET_PYRA13;
        }        

        // Iterate over all nodes of element and build up map of nodes
        // which must be transferred to the Grid.

        UInt numElemNodes = Elem::shapes[elementTypes_[i]].numNodes;
        //        std::cout << "numElemNodes " << numElemNodes << std::endl;
        for(UInt j=0; j<numElemNodes; j++) {
          UInt nodeNum = connectivity_[numNodesPerElem_*i + j];
          if(nodeNumMap2.find(nodeNum) != nodeNumMap2.end())
            continue;
          
          nodeNumMap2[nodeNum] = nodeNumMap_[nodeNum];

          //          std::cout << "nodeNum " << nodeNum << " nodeNumMap2[nodeNum] " << nodeNumMap2[nodeNum] << std::endl;

        }
      }
    }

    // If a different coordinate system than the default one was specified
    // we map the nodal coordinates into this coordinate system.
    if(coordSysId_ != "default" || scaleFac_ != 1.0) 
    {  
      CoordSystem* coordSys = domain->GetCoordSystem(coordSysId_);
      TransformNodes(*coordSys, scaleFac_);
    }

    // Add nodes to grid and remap node numbers
    UInt numNodesInGrid = mi->GetNumNodes();
    mi->AddNodes( nodeNumMap2.size() );
    Vector<Double> p(3);
    std::map<UInt, UInt>::iterator it, end;
    it = nodeNumMap2.begin();
    end = nodeNumMap2.end();
    
    for( UInt i=numNodesInGrid + 1; it != end; it++, i++ ) {
      //      std::cout << "node it->second " << it->second << " it->first " << it->first << std::endl;
      idx = (it->second-1)*3;
      p[0] = nodeCoords_[idx + 0];
      p[1] = nodeCoords_[idx + 1];
      p[2] = nodeCoords_[idx + 2];
      it->second = i;
      mi->SetNodeCoordinate(it->second, p );
    }


    UInt numKnownElems = 0;
    for( UInt i=0; i < numElements; i++ ) {
      Elem::FEType feType = elementTypes_[readElemIndices[i]];
      if(feType == Elem::ET_UNDEF || feType == Elem::ET_POINT) {
          continue;
      }
      numKnownElems++;
    }

    // Add elements to grid and check which ones are named elements
    std::map<UInt, std::string >::const_iterator peIt, peEnd;
    std::map<UInt, StdVector<UInt> > namedElems;

    UInt numElemsInGrid = mi->GetNumElems() + 1;
    mi->AddElems( numKnownElems );
    for( UInt i=0; i < numElements; i++ ) {
      Elem::FEType feType = elementTypes_[readElemIndices[i]];
      if(feType == Elem::ET_UNDEF || feType == Elem::ET_POINT) {
        continue;
      }

      UInt* conn = &connectivity_[numNodesPerElem_*readElemIndices[i]];

      UInt numElemNodes = Elem::shapes[feType].numNodes;
      for(UInt j=0; j<numElemNodes; j++) {
        conn[j] = nodeNumMap2[conn[j]];
      }

      UInt ePhysType = elementPhysicsTypes_[readElemIndices[i]];
      
      mi->SetElemData( numElemsInGrid, feType,
                       regionIdMap[ePhysType],
                       conn );

      // Check for named elements and add them to a map
      peIt = physEntities2NamedElems_.begin();
      peEnd = physEntities2NamedElems_.end();
      for( ; peIt != peEnd; peIt++ ) {
        if(ePhysType != peIt->first) continue;
        // std::cout << "Named Elems: " << peIt->first << " " << peIt->second << std::endl;
   
        namedElems[peIt->first].Push_back( numElemsInGrid );
      }
                       
      numElemsInGrid++;
    }

    // Add named nodes to grid
    std::map<UInt, StdVector<UInt> > namedNodes;
    std::map<UInt, StdVector<UInt> >::iterator nnIt, nnEnd;

    peIt = physEntities2NamedNodes_.begin();
    peEnd = physEntities2NamedNodes_.end();

    for( ; peIt != peEnd; peIt++ ) {
      // std::cout << "Named Nodes: " << peIt->first << " " << peIt->second << std::endl;
      
//      elementTypes_[idx] = feType;
//      elementPhysicsTypes_[idx] = regPhys;
      
      for(UInt i=0, n=elementPhysicsTypes_.size(); i<n; i++) {
        if(elementPhysicsTypes_[i] != peIt->first) continue;
        
        UInt* conn = &connectivity_[numNodesPerElem_*i];
        // std::cout << "We have found element physics type: " << elementPhysicsTypes_[i] << std::endl;
        // std::cout << "Node connectivity: " << conn[0] << std::endl;

        if(nodeNumMap2.find(conn[0]) != nodeNumMap2.end()) {
          // std::cout << "Node is " << conn[0] << " maps to " << nodeNumMap2[conn[0]] << std::endl;
          namedNodes[peIt->first].Push_back(nodeNumMap2[conn[0]]);
        }        
      }
    }
    
    nnIt = namedNodes.begin();
    nnEnd = namedNodes.end();
    for( ; nnIt != nnEnd; nnIt++ ) {
      mi->AddNamedNodes(physEntities2NamedNodes_[nnIt->first], nnIt->second);
    }

    // Add named elements to grid
    std::map<UInt, StdVector<UInt> >::iterator neIt, neEnd;

    neIt = namedElems.begin();
    neEnd = namedElems.end();
    for( ; neIt != neEnd; neIt++ ) {
      mi->AddNamedElems(physEntities2NamedElems_[neIt->first], neIt->second);
    }

  }

  void SimInputGmsh::TransformNodes(CoordSystem& coordSys, double scaleFac)
  {
    UInt numNodes = nodeCoords_.GetSize() / 3;
    Vector<Double> p, globPoint;
    p.Resize(3);
    globPoint.Resize(3);

    for(UInt i=0; i<numNodes; i++) 
    {
      UInt idx = i*3;
      p[0] = nodeCoords_[idx + 0];
      p[1] = nodeCoords_[idx + 1];
      p[2] = nodeCoords_[idx + 2];
      coordSys.Global2LocalCoord(globPoint, p);
      
      nodeCoords_[idx + 0] = globPoint[0] * scaleFac;
      nodeCoords_[idx + 1] = globPoint[1] * scaleFac;
      nodeCoords_[idx + 2] = globPoint[2] * scaleFac;
    }
  }
}
