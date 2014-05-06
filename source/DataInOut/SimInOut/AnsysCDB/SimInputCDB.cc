// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdarg>
#include <cctype>
#include <cstdlib>

#include <boost/algorithm/string/trim.hpp>

#include "SimInputCDB.hh"

// TODO: - Workbench FE Modeler generiert Gitter ohne Oberflächenelemente nur node Components. In diesem Fall müssen die Volumenelemente in eine extra Dummyregion gesteckt werden.
// TODO: - Herausfinden, wie man Oberflächenelemente erzeugen kann.

namespace CoupledField {

  SimInputCDB::SimInputCDB(std::string fileName, PtrParamNode inputNode,
                           PtrParamNode infoNode) :
    SimInput(fileName, inputNode, infoNode ),
    strict_(true),
    degen_(false),
    fSize_(0),
    numNodeGroups_(0),
    numElemGroups_(0),
    numRegions_(0),
    ansysSubTypes_(0),
    ansysElTypes_(0),
    ans2NacsET_(NULL),
    maxNumElemNodes_(0)
  {
    capabilities_.insert(SimInput::MESH);

    if(myParam_->Has("strict")) 
    {
      strict_ = myParam_->Get("strict")->As<bool>();
    }
    
    if(myParam_->Has("degenerate")) 
    {
      degen_ = myParam_->Get("degenerate")->As<bool>();
    }
  }


  SimInputCDB::~SimInputCDB()
  {
    CloseCDBFile();
  }


  void SimInputCDB::InitModule()
  {
    std::string line;
    std::ostringstream errMsg;
    numNodeGroups_ = 0;
    numElemGroups_ = 0;
    
    numNBlocks_ = numEBlocks_ = numCMBlocks_ = 0;
    numNCmnds_ = numENCmnds_ = numEselCmnds_ = numNselCmnds_ = numCMCmnds_ = 0;
    numSFECmnds_ = 0;
    
    if (!strict_) {
      EXCEPTION("cdb based interface currently working only in strict mode!\n"
                << "You will have to renumber your model (numcmp) eventually!");
    }
    
    // first we establish cdb file and test for blocked or unblocked cdb file format
    OpenCDBFile(fileName_);
    
    PreScanCDBFile();
    
    //--------------------------------------------------------------------
    // in case that a different file instead of cdb is used
    // for the ElementTypeMap (how to write that?)
    // ( 2-file case )
    // cdb file must be reread completely, see below
    // otherwise, cdb file must be repositioned to the very beginning
    //--------------------------------------------------------------------
    //
    // as ET commands are located before node and element definitions,
    // lets handle these first
    InitAnsys2NacsTypes();
    SetAnsys2NacsElementTypeMap();
    
    // nodes and elements
    if ( linePtsNBlocks_.size() > 0 ) {
      ReadCoordinatesBlocked();
    } else if ( linePtsNCmnds_.size() > 0 ) {
      ReadCoordinatesUnBlocked();
    }
    if ( linePtsEBlocks_.size() > 0 ) {
      ReadElementsBlocked();
    } else if ( linePtsENCmnds_.size() > 0 ) {
      ReadElementsUnBlocked();
    }

    // region and group definitions
    ReadRegionsAndGroups();
  }
    
  void SimInputCDB::ReadMesh(Grid *mi)
  {

    mi_ = mi;
    UInt dim = mi_->GetDim();

    // Get nodes and add them to grid
    UInt numNodes = nodalCoords_.size()/3;
    mi_->AddNodes(numNodes);

    std::vector< double > coords;

    for(UInt i=0; i<numNodes; i++) {
      Vector<Double> loc(3);
      loc[0] = nodalCoords_[i*3+0];
      loc[1] = nodalCoords_[i*3+1];
      loc[2] = nodalCoords_[i*3+2];
      mi_->SetNodeCoordinate(i+1, loc);
    }

    UInt numElems = topology_.size();
    bool hasVolElems = false;
    bool hasSurfElems = false;
    if(regionNames_.empty())
    {
      for(UInt i=0; i<numElems; i++)
      {
        Elem::FEType feType = Elem::FEType(elemTypes_[i+1]);
        UInt elemDim = Elem::shapes[feType].dim;

        if(hasVolElems == false && elemDim == dim) 
        {
          hasVolElems = true;
        }
        
        if(hasSurfElems == false && elemDim != dim) 
        {
          hasSurfElems = true;
        }
      }

      if(hasVolElems) 
      {
        regionNames_.push_back("volElems");
      }
      if(hasSurfElems) 
      {
        regionNames_.push_back("surfElems");
      }
    }

    StdVector<Integer> ids;
    mi_->AddRegions(regionNames_, ids);

    mi_->AddElems(numElems);

    for(UInt i=0; i<numElems; i++)
    {
      Integer regionId = NO_REGION_ID;
      
      if(elemRegionMap_.find(i+1) != elemRegionMap_.end()) 
      {
        regionId = ids[elemRegionMap_[i+1]];
      }
      
      Elem::FEType feType = Elem::FEType(elemTypes_[i+1]);
      UInt elemDim = Elem::shapes[feType].dim;
      if(regionId == NO_REGION_ID) 
      {
        if(elemDim == dim) 
        {
          regionId = 0;
        }
        else
        {
          regionId = 1;
        }
      }

      mi_->SetElemData(i+1, feType, regionId, &topology_[i+1][0]);
    }

    // Add named nodes to grid.
    for(UInt i=0, n=nodeGroupNames_.size(); i<n; i++) 
    {
      std::string name = nodeGroupNames_[i];
      StdVector<UInt>& nodes = nodeGroupData_[i];
      
      mi_->AddNamedNodes(name, nodes);
    }
    
    // Add named elems to grid.
    for(UInt i=0, n=elemGroupNames_.size(); i<n; i++) 
    {
      std::string name = elemGroupNames_[i];
      StdVector<UInt>& nodes = elemGroupData_[i];
      
      mi_->AddNamedElems(name, nodes);
    }
  }

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputCDB::GetDim()
  {
    return dim_;
  }
  
  UInt SimInputCDB::GetNumNodes()
  {
    return nodalCoords_.size()/3;
  }
    
  UInt SimInputCDB::GetNumElems(const Integer dim)
  {
    UInt numElems = topology_.size();
    return numElems;
  }
  
  UInt SimInputCDB::GetNumRegions()
  {
    EXCEPTION("Not implemented!");
    return 0;
  }

  UInt SimInputCDB::GetNumNamedNodes()
  {
    EXCEPTION("Not implemented!");
    return 0;
  }

  UInt SimInputCDB::GetNumNamedElems()
  {
    EXCEPTION("Not implemented!");
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputCDB::GetAllRegionNames( StdVector<std::string> & regionNames )
  {
    EXCEPTION("Not implemented!");
  }
    
  void SimInputCDB::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                       const UInt dim )
  {
    EXCEPTION("Not implemented!");
  }

  void SimInputCDB::GetNodeNames( StdVector<std::string> &nodeNames )
  {
    EXCEPTION("Not implemented!");
  }

  void SimInputCDB::GetElemNames( StdVector<std::string> & elemNames )
  {
    EXCEPTION("Not implemented!");
  }

  void SimInputCDB::DegenerateElement(const Elem::FEType elemTypeIn,
                                      Elem::FEType& elemTypeOut,
                                      std::vector<UInt>& elemNodes)
  {
    static std::vector<UInt> newElemNodes;
    static std::map<Elem::FEType, std::vector<UInt> > idxMap;
    elemTypeOut = elemTypeIn;

    switch(elemTypeIn)
    {
    case Elem::ET_TRIA3:
      elemTypeOut = Elem::ET_QUAD4;
      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2};
        std::copy(elemIdxMap, elemIdxMap+4, std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_TRIA6:
      elemTypeOut = Elem::ET_QUAD8;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2, 3, 4, 2, 5};
        std::copy(elemIdxMap, elemIdxMap+8, std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_TET4:
      elemTypeOut = Elem::ET_HEXA8;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2, 3, 3, 3, 3};
        std::copy(elemIdxMap, elemIdxMap+8,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_TET10:
      elemTypeOut = Elem::ET_HEXA20;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2,
                             3, 3, 3, 3,
                             4, 5, 2, 6,
                             3, 3, 3, 3,
                             7, 8, 9, 9};
        std::copy(elemIdxMap, elemIdxMap+20,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_PYRA5:
      elemTypeOut = Elem::ET_HEXA8;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 3,
                             4, 4, 4, 4};
        std::copy(elemIdxMap, elemIdxMap+8,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_PYRA13:
      elemTypeOut = Elem::ET_HEXA20;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 3,
                             4, 4, 4, 4,
                             5, 6, 7, 8,
                             4, 4, 4, 4,
                             9, 10, 11, 12};
        std::copy(elemIdxMap, elemIdxMap+20,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_WEDGE6:
      elemTypeOut = Elem::ET_HEXA8;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2,
                             3, 4, 5, 5};
        std::copy(elemIdxMap, elemIdxMap+8,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    case Elem::ET_WEDGE15:
      elemTypeOut = Elem::ET_HEXA20;

      if(idxMap[elemTypeIn].empty())
      {
        UInt elemIdxMap[] = {0, 1, 2, 2,
                             3, 4, 5, 5,
                             6, 7, 2, 8,
                             9, 10, 5, 11,
                             12, 13, 14, 14};
        std::copy(elemIdxMap, elemIdxMap+20,  std::back_inserter(idxMap[elemTypeIn]));
      }
      break;

    default:
      return;
    }

    newElemNodes.resize(Elem::shapes[elemTypeOut].numNodes);
    std::vector<UInt>::const_iterator it, end;
    it = idxMap[elemTypeIn].begin();
    end = idxMap[elemTypeIn].end();

    for(UInt i=0 ; it != end; it++, i++)
    {
      newElemNodes[i] = elemNodes[*it];
    }

    std::copy(newElemNodes.begin(), newElemNodes.end(), elemNodes.begin());
  }

  void SimInputCDB::ResortNodes(std::vector<UInt>& elemNodes)
  {
     std::set<UInt> nodeSet;
     UInt pos;

     nodeSet.insert(elemNodes[0]);
     pos = 1;

     for(UInt i=1, n=elemNodes.size(); i<n; i++) {
       if(nodeSet.find(elemNodes[i]) != nodeSet.end()) {
         continue;
       }

       elemNodes[pos++] = elemNodes[i];
       nodeSet.insert(elemNodes[i]);
     }
  }
  void SimInputCDB::PreScanCDBFile() {

    std::string line;

    // for large file sizes, report file size rounded to 100 MB
#if(WIN32 || __MINGW32__)
    Double tmpVal = (Double) (inSize_/(1024.0*1024.0*1024.0));
#else
    Double tmpVal = (Double) (fSize_/(1024.0*1024.0*1024.0));
#endif
    int tmp = (int) (tmpVal*10.0);
    tmpVal = tmp*0.1;
    std::cout << "Initial scan of input file";
    if (tmpVal > 0.1)
      std::cout << " ("    << tmpVal << " GB)";

    std::cout << std::endl;

#if(WIN32 || __MINGW32__)
    __int64 pos = _ftelli64( inStream_ );
#else
    unsigned long pos = inFile_.tellg();
#endif

    UInt count = 0;
    if (pos != 0) 
      EXCEPTION("File scan: initial position not at beginning of file");

    while(GetNextLine(line)) {
      count++;
      if (line.substr(0,6) == "NBLOCK" || line.substr(0,6) == "nblock") {
        numNBlocks_++;
        linePtsNBlocks_.push_back(pos);
      }
      if (line.substr(0,6) == "EBLOCK" || line.substr(0,6) == "eblock") {
        numEBlocks_++;
        linePtsEBlocks_.push_back(pos);
      }
      if (line.substr(0,2) == "N," || line.substr(0,2) == "n,") {
        numNCmnds_++;
        linePtsNCmnds_.push_back(pos);
      }
      if (line.substr(0,3) == "EN," || line.substr(0,3) == "en,") {
        if (line.find("ATTR") != line.npos || line.find("attr") != line.npos) {
          numENCmnds_++;
          linePtsENCmnds_.push_back(pos);
        }
      }
      if (line.substr(0,7) == "CMBLOCK" || line.substr(0,7) == "cmblock") {
        numCMBlocks_++;
        linePtsCMBlocks_.push_back(pos);
      }
      if (line.substr(0,4) == "ESEL" || line.substr(0,4) == "esel") {
        numEselCmnds_++;
        linePtsEselCmnds_.push_back(pos);
      }
      if (line.substr(0,4) == "NSEL" || line.substr(0,4) == "nsel") {
        numNselCmnds_++;
        linePtsNselCmnds_.push_back(pos);
      }
      if (line.substr(0,3) == "CM," || line.substr(0,3) == "cm,") {
        numCMCmnds_++;
        linePtsCMCmnds_.push_back(pos);
      }
      if (line.substr(0,3) == "ET," || line.substr(0,3) == "et,") {
        lineETCmnds_.push_back(line);
      }
      if (line.substr(0,5) == "KEYOP" || line.substr(0,5) == "keyop") {
        lineKEYOPTCmnds_.push_back(line);
      }
      if (line.find("Pressure Using Surface Effect Elements") != line.npos) {
        linePtsPSECmnds_.push_back(pos);
      }
      if (line.find("Construct Weak Springs") != line.npos) {
        linePtsPSEStop_ = pos;
      }
      if (line.substr(0,4) == "SFE," || line.substr(0,4) == "sfe,") {
        numSFECmnds_++;
        linePtsSFECmnds_.push_back(pos);
      }
      pos = GetFilePosition();
      if (count%10000000 ==0) {
        std::cout << "Lines read " << count << "\r"; 
        std::cout.flush();
      }
      
    }
    std::cout << std::endl;

    if (numNBlocks_>0)
      std::cout << "Number of node blocks found " << numNBlocks_ << std::endl;
    if (numEBlocks_>0) 
      std::cout << "Number of element blocks found " << numEBlocks_ << std::endl;
    if (numCMBlocks_>0)
      std::cout << "Number of cm blocks found " << numCMBlocks_ << std::endl;
    if (numNCmnds_>0)
      std::cout << "Number of N commands found " << numNCmnds_ << std::endl;
    if (numENCmnds_>0)
      std::cout << "Number of EN commands found " << numENCmnds_ << std::endl;
    if (numEselCmnds_>0)
      std::cout << "Number of ESEL commands found " << numEselCmnds_ << std::endl;
    if (numCMCmnds_>0)
      std::cout << "Number of CM commands found " << numCMCmnds_ << std::endl;
    if (lineETCmnds_.size()>0)
      std::cout << "Number of ET commands found " << lineETCmnds_.size() << std::endl;
    if (linePtsPSECmnds_.size()>0)
      std::cout << "Number of Pressure Surface groups found " << linePtsPSECmnds_.size() << std::endl;
    if (numSFECmnds_>0)
      std::cout << "Number of SFE commands found " << numSFECmnds_ << std::endl;

    std::cout << "Finished scanning input file" << std::endl;

  }

  void SimInputCDB::InitAnsys2NacsTypes() {

    ansysSubTypes_ = 12;
    ansysElTypes_  = 281;

    ans2NacsET_ = new UInt[ansysSubTypes_*(ansysElTypes_+1)];
    for (UInt i=0; i<=ansysSubTypes_*ansysElTypes_; i++)
      ans2NacsET_[i] = Elem::ET_UNDEF;

    std::vector<UInt> etList;

    // NACS type 2=Elem::ET_LINE2
    etList.push_back(1); etList.push_back(3); etList.push_back(4); etList.push_back(8);
    etList.push_back(10); etList.push_back(23); etList.push_back(24); etList.push_back(32);
    etList.push_back(33); etList.push_back(44); etList.push_back(54); etList.push_back(68);
    etList.push_back(125); etList.push_back(129); etList.push_back(138); etList.push_back(180);
    etList.push_back(188);
    etList.push_back(189); etList.push_back(200); etList.push_back(208); etList.push_back(251);
    for (UInt i=0; i<etList.size(); i++) {
      ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_LINE2;
    }
    // ansys element 200,2 must also be accounted for
    ans2NacsET_[(200-1)*ansysSubTypes_+3] = Elem::ET_LINE2;
    // ansys element 153,0 must also be accounted for
    ans2NacsET_[(153-1)*ansysSubTypes_+1] = Elem::ET_LINE2;

    // NACS type 3=Elem::ET_LINE3
    etList.clear();
    etList.push_back(161); etList.push_back(209); 
    for (UInt i=0; i<etList.size(); i++) {
      ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_LINE3;
    }
    // ansys element 200,1 and 200,3 must also be accounted for
    ans2NacsET_[(200-1)*ansysSubTypes_+2] = Elem::ET_LINE3;
    ans2NacsET_[(200-1)*ansysSubTypes_+4] = Elem::ET_LINE3;
    // ansys element 153,1
    ans2NacsET_[(153-1)*ansysSubTypes_+2] = Elem::ET_LINE3;

    // NACS type 4 = TRIA3 (currently only 200,5)
    //etList.clear();
    //for (UInt i=0; i<etList.size(); i++) {
    //  ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = 5;
    //}
    ans2NacsET_[(200-1)*ansysSubTypes_+5] = Elem::ET_TRIA3;

    // NACS type 5 = TRIA6
    etList.clear();
    etList.push_back(35); etList.push_back(146); 
    for (UInt i=0; i<etList.size(); i++) {
      ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_TRIA6;
    }
    // ansys element 200,5 must also be accounted for
    ans2NacsET_[(200-1)*ansysSubTypes_+6] = Elem::ET_TRIA6;

    // NACS type 6 = QUAD4
    etList.clear();
    etList.push_back(13); etList.push_back(25); etList.push_back(29); etList.push_back(42);
    etList.push_back(55); etList.push_back(63); etList.push_back(67); etList.push_back(75);
    etList.push_back(79); etList.push_back(106); etList.push_back(115); etList.push_back(130);
    etList.push_back(141); etList.push_back(162); etList.push_back(181); etList.push_back(182);
    etList.push_back(192); etList.push_back(202); etList.push_back(252);
    for (UInt i=0; i<etList.size(); i++) {
      ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_QUAD4;
    }
    // ansys element 200,6 must also be accounted for
    ans2NacsET_[(200-1)*ansysSubTypes_+7] = Elem::ET_QUAD4;
    // ansys element 154,0 must also be accounted for
    ans2NacsET_[(154-1)*ansysSubTypes_+1] = Elem::ET_QUAD4;

    // NACS type 7 = QUAD8
    etList.clear();
    etList.push_back(53); etList.push_back(77); etList.push_back(78); etList.push_back(82);
    etList.push_back(83); etList.push_back(88); etList.push_back(108); etList.push_back(118);
    etList.push_back(121); etList.push_back(145); etList.push_back(183); etList.push_back(190);
    etList.push_back(223); etList.push_back(230); etList.push_back(281);
    for (UInt i=0; i<etList.size(); i++) {
      ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_QUAD8;
    }
    // ansys element 200,7 must also be accounted for
    ans2NacsET_[(200-1)*ansysSubTypes_+8] = Elem::ET_QUAD8;
    // ansys element 154,1 must also be accounted for
    ans2NacsET_[(154-1)*ansysSubTypes_+2] = Elem::ET_QUAD8;

    // NACS type 9 - TET4 
    //etList.clear();
    //for (UInt i=0; i<etList.size(); i++) {
    //  ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_TET4;
    //}
    // ansys element 200,8
    ans2NacsET_[(200-1)*ansysSubTypes_+9] = Elem::ET_TET4;

    // NACS type 10 - TET10
    etList.clear();
    etList.push_back(87); etList.push_back(92); etList.push_back(98); etList.push_back(119);
    etList.push_back(123); etList.push_back(127); etList.push_back(148); etList.push_back(168);
    etList.push_back(187); etList.push_back(227); etList.push_back(232);
    for (UInt i=0; i<etList.size(); i++) {
      ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_TET10;
    }
    // ansys element 200,9
    ans2NacsET_[(200-1)*ansysSubTypes_+10] = Elem::ET_TET10;

    // NACS type 11 - HEX8
    etList.clear();
    etList.push_back(5); etList.push_back(30); etList.push_back(45); etList.push_back(46);
    etList.push_back(62); etList.push_back(65); etList.push_back(69); etList.push_back(70);
    etList.push_back(80); etList.push_back(96); etList.push_back(97); etList.push_back(107);
    etList.push_back(142); etList.push_back(164); etList.push_back(185); etList.push_back(195);
    for (UInt i=0; i<etList.size(); i++) {
      ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_HEXA8;
    }
    // ansys element 200,10
    ans2NacsET_[(200-1)*ansysSubTypes_+11] = Elem::ET_HEXA8;

    // NACS type 12 - HEX20
    etList.clear();
    etList.push_back(89); etList.push_back(90); etList.push_back(95); etList.push_back(117);
    etList.push_back(120); etList.push_back(122); etList.push_back(128); etList.push_back(147);
    etList.push_back(186); etList.push_back(191); etList.push_back(226); etList.push_back(231);
    for (UInt i=0; i<etList.size(); i++) {
      ans2NacsET_[(etList[i]-1)*ansysSubTypes_+1] = Elem::ET_HEXA20;
    }
    // ansys element 200,11
    ans2NacsET_[(200-1)*ansysSubTypes_+12] = Elem::ET_HEXA20;

    // NACS type -1 --> should be ignored.....
    // ans2NacsET_[(14-1)*ansysSubTypes_+1] = Elem::ET_UNDEF;

//----------------------------------------------------------------------------------------------
//    dimensions currently not needed
//    ansETDim_ = new UInt[ansysElTypes_+1];
//    for (UInt i=0; i<=ansysElTypes_; i++)
//      ansETDim_[i] = 0;
//    // dimensions: 2d
//    etList.clear();
//    etList.push_back(1); etList.push_back(3); etList.push_back(9); etList.push_back(12);
//    etList.push_back(13); etList.push_back(23); etList.push_back(25); etList.push_back(29);
//    etList.push_back(32); etList.push_back(35); etList.push_back(42); etList.push_back(53);
//    etList.push_back(54); etList.push_back(55); etList.push_back(63); etList.push_back(67);
//    etList.push_back(75); etList.push_back(77); etList.push_back(78); etList.push_back(79);
//    etList.push_back(81); etList.push_back(82); etList.push_back(83); etList.push_back(109);
//    etList.push_back(110); etList.push_back(118); etList.push_back(121); etList.push_back(129);
//    etList.push_back(131); etList.push_back(132); etList.push_back(141); etList.push_back(145);
//    etList.push_back(146); etList.push_back(150); etList.push_back(151); etList.push_back(153);
//    etList.push_back(162); etList.push_back(169); etList.push_back(181); etList.push_back(182);
//    etList.push_back(183); etList.push_back(190); etList.push_back(192); etList.push_back(193);
//    etList.push_back(223); etList.push_back(230);
//    for (UInt i=0; i<etList.size(); i++) {
//      ansETDim_[etList[i]] = 2;
//    }

//    // dimensions: 3d
//    etList.clear();
//    etList.push_back(4); etList.push_back(5); etList.push_back(8); etList.push_back(24);
//    etList.push_back(30); etList.push_back(33); etList.push_back(43); etList.push_back(44);
//    etList.push_back(45); etList.push_back(47); etList.push_back(52); etList.push_back(62);
//    etList.push_back(65); etList.push_back(69); etList.push_back(70); etList.push_back(80);
//    etList.push_back(87); etList.push_back(90); etList.push_back(92); etList.push_back(95);
//    etList.push_back(96); etList.push_back(97); etList.push_back(98); etList.push_back(111);
//    etList.push_back(115); etList.push_back(117); etList.push_back(119); etList.push_back(120);
//    etList.push_back(122); etList.push_back(123); etList.push_back(127); etList.push_back(128);
//    etList.push_back(130); etList.push_back(136); etList.push_back(138); etList.push_back(139);
//    etList.push_back(142); etList.push_back(147); etList.push_back(148); etList.push_back(152);
//    etList.push_back(154); etList.push_back(156); etList.push_back(160); etList.push_back(161);
//    etList.push_back(164); etList.push_back(166); etList.push_back(168); etList.push_back(170);
//    etList.push_back(185); etList.push_back(186); etList.push_back(187); etList.push_back(188);
//    etList.push_back(189); etList.push_back(194); etList.push_back(195); etList.push_back(226);
//    etList.push_back(227); etList.push_back(231); etList.push_back(232); etList.push_back(236);
//    etList.push_back(237);
//    for (UInt i=0; i<etList.size(); i++) {
//      ansETDim_[etList[i]] = 3;
//    }
//----------------------------------------------------------------------------------------------

  }

  void SimInputCDB::SetAnsys2NacsElementTypeMap()
  {
    std::string line;
    std::ostringstream errMsg;
    UInt ansysType = 0, ansElemType = 0;
    Integer ansysSubType = -1;
    UInt nacsElemType = 0;
    std::vector< std::string > tokens(32);

    ansNacsMap_.clear();

    UInt numETCmnds = lineETCmnds_.size();

    for (UInt i=0; i<numETCmnds; i++) {

      // ET, ITYPE, Ename, KOP1, KOP2, KOP3, KOP4, KOP5, KOP6, INOPR
      // Defines a local element type from the element library.
      line = lineETCmnds_[i];
      UInt numTok = SplitLine(line, tokens);

      ansElemType = ansysType = nacsElemType = 0;
      ansysSubType = -1;

      // Ignore all  ASCII characters  of category name  and just  account for
      // digits on ET lines. E.g. et,2,MESH200
      size_t pos = 0;
      while(!std::isdigit(tokens[2][pos]))
      {
        pos++;
      }

      ansElemType = std::strtoul(tokens[1].c_str(), NULL, 0);
      ansysType = std::strtoul(tokens[2].substr(pos).c_str(), NULL, 0);

      if (numTok > 3) {
        if (ansysType==200 || ansysType==154 ) {
          ansysSubType = std::strtoul(tokens[3].c_str(), NULL, 0);
        }
      }

      // If the subtype of the element could not be determined yet, we have to
      // look it up in the KEYOPT lines.
      if(ansysSubType < 0) 
      {        
        if(ansysType==200 || ansysType==154) 
        {
          // KEYOPT, ITYPE, KNUM, VALUE
          // Sets element key options.
          // ITYPE: Element type number as defined on the ET command.
          // KNUM: Number of the KEYOPT to be defined (KEYOPT(KNUM)).
          // VALUE: Value of this KEYOPT.
          
          UInt itype = 0;
          Integer knum = -1;
          Integer value = -1;
          
          for (UInt j=0, n=lineKEYOPTCmnds_.size(); j<n; j++) {
            line = lineKEYOPTCmnds_[j];
            numTok = SplitLine(line, tokens);
            itype = std::strtoul(tokens[1].c_str(), NULL, 0);

            if(ansElemType != itype) continue;
            
            knum = std::strtoul(tokens[2].c_str(), NULL, 0);
            value = std::strtoul(tokens[3].c_str(), NULL, 0);
            break;
          }
          
          if(knum < 0 || value < 0) 
          {
            errMsg << "ANSYS element type is 200 or 154 and no subtype could " 
                   << "be determined!";
            if(strict_) {
              EXCEPTION(errMsg.str());
            } else {
              std::cerr << errMsg.str() << std::endl;
            }
          }
        
          if(knum == 1) 
          {
            ansysSubType = value+1;
          }
        }
        else
        {
          ansysSubType = 1;
        }        
      }

      // check for valid ANSYS element type
      if(ansElemType==0 || ansysType==0 ) {
        errMsg << "ANSYS element type pointer or type is 0. A serious problem "
            << "occurred during processing the interface output "
            << "from ANSYS. Please contact SIMetris support!";
        if(strict_) {
          EXCEPTION(errMsg.str());
        } else {
          std::cerr << errMsg.str() << std::endl;
        }
      }

      if (ansysType > 0 && ansysType <= ansysElTypes_) {
        nacsElemType = ans2NacsET_[(ansysType-1)*ansysSubTypes_+ansysSubType];
      } else {
        errMsg << "Error: invalid ansys element type encountered in cdb interface: "
               << ansysType;
        EXCEPTION(errMsg.str());
      }

      // check for valid NACS element type
      if( nacsElemType == Elem::ET_UNDEF || nacsElemType == Elem::ET_POINT ) {
        errMsg << "Warning:\nAn ANSYS element type which is mapped to NACS element type UNDEFINED or POINT"
            << " has been found.\nThe ANSYS type, which has no NACS counterpart, is " << ansElemType
	    << " and the selected ANSYS element is " << ansysType << ".\n";
        if(strict_) {
          std::cout << errMsg.str() << std::endl;
        } else {
          std::cerr << errMsg.str() << std::endl;
        }
      }

      ansNacsMap_[ansElemType]=nacsElemType;

    }
  }

  void SimInputCDB::ReadCoordinatesBlocked()
  {
    // read node coordinates defined in nblock commands
    // cf. Programmer's Manual for Mechanical APDL - NBLOCK Command

    std::string line;
    std::ostringstream errMsg;
    double x, y, z;
    UInt fileNodeNum = 0;
    // UInt numFields = 0;
    UInt nodeNum = 1;
    UInt maxNodeNum = 0;
    UInt numNodes = 0;
    std::vector< std::string > tokens(32);

    // loop on all found nblock commands
    for (UInt ib=0; ib<linePtsNBlocks_.size(); ib++) {

      if (!GetLine(line,linePtsNBlocks_[ib]))
        EXCEPTION("Read error in GetLine");

      if ((line.find("NBLOCK") == line.npos) && (line.find("nblock") == line.npos)) {
        errMsg << "NBLOCK " << ib << " expected to start at line no. " << linePtsNBlocks_[ib]
               << " not found!\nLine read was " << line;
        EXCEPTION(errMsg);
      }
      std::cout << "Processing node block " << ib+1 << std::endl;

      // read number of fields in blocked node lines from NBLOCK line
      UInt numTok = SplitLine(line, tokens);
      // numFields = std::strtoul(tokens[1].c_str(), NULL, 0);

      // retrieve format by splitting a line of the form (3i8,6e20.13)
      std::vector<int> offsets;
      GetNextLine(line);
      UInt numInts = DecodeBlockFormatLine(line, offsets);

      // read-in procedure for blocked case
      UInt numNodesInBlock = 0;

      GetNextLine(line);
      while(line.substr(0,2) != "N," &&
            line.substr(0,2) != "n," &&
            line.substr(0,2) != "-1" &&
            line.substr(0,1) != "!") {

        // End of a typical node block in a .cdb file:
        //     5706       0       0 9.3694972500000E+02-3.4946417960000E+02
        // N,R5.3,LOC,       -1,

        // End of a typical node block in an .inp file:
        //    5706    9.369497250E+002   -3.494641796E+002    0.000000000E+000
        // ! end of nblock command
        // !

        numTok = SplitLine(line, tokens, "", &offsets);

        fileNodeNum = x = y = z = 0;

        fileNodeNum = std::strtoul(tokens[0].c_str(), NULL, 0);
        UInt numActualFloats = numTok - numInts;
        if(numActualFloats > 0) 
        {
          x = std::strtod(tokens[numInts + 0].c_str(), NULL);
          if(numActualFloats > 1) 
          {
            y = std::strtod(tokens[numInts + 1].c_str(), NULL);
            if(numActualFloats > 2) 
            {
              z = std::strtod(tokens[numInts + 2].c_str(), NULL);
            }
          }
        }

        StoreSingleNode(fileNodeNum,x,y,z, nodeNum, numNodes, maxNodeNum);
	numNodesInBlock++;
        if (numNodes%1000000 == 0) {
          std::cout << "Nodes read " << numNodes << "\r";
          std::cout.flush();
        }

        GetNextLine(line);
      }
      std::cout << std::endl;

      errMsg << "Node array written by ANSYS seems to have gaps!";
      if(maxNodeNum != numNodes && strict_) {
        EXCEPTION(errMsg.str());
      } else if(maxNodeNum != numNodes) {
        std::cerr << errMsg.str() << std::endl;

        std::cerr << "maxNodeNum " << maxNodeNum << std::endl;
        std::cerr << "numNodes " << numNodes << std::endl;
      }
      std::cout << "nblock " << ib+1 << " (" << numNodesInBlock << " nodes found in block)" << std::endl;
    }
    std::cout << "Finished reading nodes (" << numNodes << " read)" << std::endl;
  }

  void SimInputCDB::ReadCoordinatesUnBlocked()
  {
    std::string line;
    std::ostringstream errMsg;
    double x, y, z;
    UInt fileNodeNum = 0;
    UInt nodeNum = 1;
    UInt maxNodeNum = 0;
    UInt numNodes = 0;
    std::vector< std::string > tokens(32);

    // loop on all found n commands
    for (UInt ib=0; ib<linePtsNCmnds_.size(); ib++) {

      GetLine(line,linePtsNCmnds_[ib]);

      // process line for current node
      fileNodeNum = x = y = z = 0;

      SplitLine(line, tokens);

      fileNodeNum = std::strtoul(tokens[3].c_str(), NULL, 0);
      x = std::strtod(tokens[6].c_str(), NULL);
      y = std::strtod(tokens[7].c_str(), NULL);
      z = std::strtod(tokens[8].c_str(), NULL);

      StoreSingleNode(fileNodeNum,x,y,z, nodeNum, numNodes, maxNodeNum);
      if (numNodes%1000000 == 0) {
        std::cout << "Nodes read " << numNodes << "\r";
        std::cout.flush();
      }
    }
    std::cout << std::endl;

    errMsg << "Node array written by ANSYS seems to have gaps!";
    if(maxNodeNum != numNodes && strict_)
    {
      EXCEPTION(errMsg.str());
    }
    else if(maxNodeNum != numNodes)
    {
      std::cerr << errMsg.str() << std::endl;

      std::cerr << "maxNodeNum " << maxNodeNum << std::endl;
      std::cerr << "numNodes " << numNodes << std::endl;
    }
    std::cout << "Finished reading nodes (" << numNodes << " read)" << std::endl;
  }

  void SimInputCDB::StoreSingleNode(UInt fileNodeNum,
                                    double x, double y, double z,
                                    UInt &nodeNum,
                                    UInt &numNodes,
                                    UInt &maxNodeNum) {

    std::ostringstream errMsg;

    // Sollten die Knotennummern in der Knotendatei nicht mit der
    // internen Knotennummerierung uebereinstimmen eine Meldung ausgeben

    if(fileNodeNum != nodeNum) {
      errMsg << "Nodes seem to be not consecutive. Node nr. in file "
          << fileNodeNum << " does not match internal node nr. "
          << nodeNum << ". Please check your model!";
      if(strict_) {
        EXCEPTION(errMsg.str());
      } else {
        std::cerr << errMsg.str() << std::endl;
      }
      errMsg.clear(); errMsg.str("");
    }

    // Hier abfragen, ob die Knotennummern konsekutiv sind
    if(!nodeNumsMap_.empty() &&
        (nodeNumsMap_.rbegin()->first + 1) != nodeNum)
    {
      errMsg << "Nodes are not consecutive: " << nodeNumsMap_.rbegin()->first
          << " -> " << nodeNum;
      if(strict_) {
        EXCEPTION(errMsg.str());
      } else {
        std::cerr << errMsg.str() << std::endl;
      }
      errMsg.clear(); errMsg.str("");
    }

    // Dieser Fehler ist fatal! Er hat nichts mit strict oder relaxed zu tun.
    if(nodeNumsMap_[nodeNum] != 0)
      EXCEPTION("Node " << nodeNum << " has been referenced before!\n");

    nodalCoords_.push_back(x);
    nodalCoords_.push_back(y);
    nodalCoords_.push_back(z);

    maxNodeNum = maxNodeNum < nodeNum ? nodeNum : maxNodeNum;
    numNodes++;
    nodeNumsMap_[nodeNum] = numNodes;
    nodeNum++;
  }

  void SimInputCDB::ReadElementsBlocked() {
    std::string line;
    std::vector< std::string > tokens(32);
    std::ostringstream errMsg;

    // Read element types
    UInt elemType, elemMat;
    // UInt ansElemNum = 0;
    UInt elemNum = 1;
    UInt elemDim = 0;
    UInt dim = 0;
    std::vector<UInt> elemNodes(40);
    std::vector<UInt> lineContent(40);
    std::set<UInt> elemNodeSet;
    UInt numElemNodes;
    // maxNonSolidNodes is the maximum number of element nodes per non-solid eblock
    // if maxNonSolidNodes > 10, 2 records have to be read per element
    UInt maxNonSolidNodes = 0;
    bool isSolidEBlock = false;
    UInt numSkipElems = 0, firstSkippedElem = 0;

    // loop on all found eblock commands
    for (UInt ib=0; ib<linePtsEBlocks_.size(); ib++) {

      if (!GetLine(line,linePtsEBlocks_[ib])) {
        errMsg << "EBLOCK " << ib << " expected to start at line no. "
               << linePtsEBlocks_[ib] << " not found!\nLine read was " << line;
        EXCEPTION(errMsg);
      }

      // do we have a solid eblock ?
      if (line.find("solid") != line.npos || line.find("SOLID") != line.npos )
        isSolidEBlock = true;

      UInt numTok = SplitLine(line, tokens);

      // in case of a  nonsolid eblock, we need to know  the maximum number of
      // element nodes per element
      if (!isSolidEBlock) {
        maxNonSolidNodes = std::strtoul(tokens[1].c_str(), NULL, 0);
      }

      // check format line
      std::vector<int> offsets;
      GetNextLine(line);
      DecodeBlockFormatLine(line, offsets);

      UInt numElemsInBlock = 0;

      GetNextLine(line);
      while(line.find("-1") == line.npos) {
        std::fill(elemNodes.begin(), elemNodes.end(), 0);
        elemNodeSet.clear();

        numTok = SplitLine(line, tokens, "", &offsets);

        if (isSolidEBlock) {
          elemMat = std::strtoul(tokens[0].c_str(), NULL, 0);
          elemType = std::strtoul(tokens[1].c_str(), NULL, 0);
          // ansElemNum = std::strtoul(tokens[10].c_str(), NULL, 0);
          numElemNodes = std::strtoul(tokens[8].c_str(), NULL, 0);

          for (UInt i=0, n=offsets.size()-11; i<n; i++) {
            elemNodes[i] = std::strtoul(tokens[i + 11].c_str(), NULL, 0);
          }

          if (numElemNodes > 8) {
            GetNextLine(line);
            numTok = SplitLine(line, tokens, "", &offsets);

            for (UInt i=0; i<numTok; i++) {
              elemNodes[i+8] = std::strtoul(tokens[i].c_str(), NULL, 0);
            }
          }
	} else {
          elemMat = std::strtoul(tokens[3].c_str(), NULL, 0);
          elemType = std::strtoul(tokens[1].c_str(), NULL, 0);
          // ansElemNum = std::strtoul(tokens[0].c_str(), NULL, 0);

          numElemNodes = 0;
          for (UInt i=0, n=10; i<n; i++) {
            elemNodes[i] = std::strtoul(tokens[i + 5].c_str(), NULL, 0);
	    numElemNodes++;
          }

          if (maxNonSolidNodes > 10) {
            GetNextLine(line);
            numTok = SplitLine(line, tokens, "", &offsets);

            for (UInt i=0, n=numTok; i<n; i++) {
              elemNodes[i+10] = std::strtoul(tokens[i].c_str(), NULL, 0);
              numElemNodes++;
            }
          }
	}

        if (ansNacsMap_[elemType] != Elem::ET_UNDEF) {
          elemTypes_[elemNum] = ansNacsMap_[elemType];
          // Determine number of element nodes by inserting nodes into a set.
          elemNodeSet.insert(elemNodes.begin(), elemNodes.end());
          elemNodeSet.erase((UInt) 0);
          // reset numElemNodes as set elemNodeSet might contain fewer node than elemNodes
          numElemNodes = elemNodeSet.size();

          maxNumElemNodes_ = numElemNodes > maxNumElemNodes_ ? numElemNodes : maxNumElemNodes_;

          Elem::FEType newFEType = DegenTypeToNativeType(elemTypes_[elemNum], numElemNodes);

          if(degen_) {
            DegenerateElement((Elem::FEType)elemTypes_[elemNum], newFEType, elemNodes);
          } else {
            ResortNodes(elemNodes);
          }

          elemMaterials_.push_back(elemMat);
          elemTypes_[elemNum] = newFEType;
          elemNumsMap_[elemNum] = elemNum;
          elemDim = Elem::shapes[newFEType].dim;
          dim = dim < elemDim ? elemDim : dim;

          std::copy(elemNodes.begin(), elemNodes.begin() + Elem::shapes[newFEType].numNodes,
              std::back_inserter(topology_[elemNum]));

          elemNum++;

	} else {
          if (numSkipElems == 0)
            firstSkippedElem = elemNum;

          numSkipElems++;
        }

        numElemsInBlock++;
        GetNextLine(line);
        if (elemNum%1000000 == 0) {
          std::cout << "Elements read " << elemNum << "\r";
          std::cout.flush();
        }

      }
      std::cout << std::endl;
      std::cout << "Element block " << ib << "(" << numElemsInBlock << " elements read)" << std::endl;
      if (numSkipElems > 0)
        std::cout << numSkipElems << " ANSYS elements have been skipped, as no NACS counterpart exists" 
                  << std::endl;
    }

    // store dimension to settings
    dim_ = dim;
    std::cout << "Finished reading elements (" << elemNum-1 << " read)" << std::endl;

    if (numSkipElems>0 && elemNum > firstSkippedElem)
      std::cout << "WARNING: non skipped elements found after first skipped element (= " 
                << firstSkippedElem << ").\nThis might lead to incorrect mesh transfer.\n" 
		<< "Please check carefully!" << std::endl;

  }

  void SimInputCDB::ReadElementsUnBlocked() {
    std::string line;
    std::ostringstream errMsg;
    std::vector< std::string > tokens(32);

    // Read element types
    UInt elemType;
    UInt elemNum = 1;
    // UInt ansElemNum = 0;
    UInt elemDim = 0;
    UInt dim = 0;
    std::vector<UInt> elemNodes(40);
    std::vector<UInt> lineContent(40);
    std::set<UInt> elemNodeSet;
    UInt numElemNodes;
    UInt numSkipElems = 0, firstSkippedElem = numENCmnds_+1;

    // loop on all found en commands

    UInt lCount = 0;
    for (UInt ib=0; ib<linePtsENCmnds_.size(); ib++) {

      GetLine(line,linePtsENCmnds_[ib]);
      lCount++;
      UInt numTok = SplitLine(line, tokens);

      numElemNodes = std::strtoul(tokens[3].c_str(), NULL, 0);
      elemType = std::strtoul(tokens[5].c_str(), NULL, 0);
      // ansElemNum = std::strtoul(tokens[9].c_str(), NULL, 0);

      GetNextLine(line);
      lCount++;
      numTok = SplitLine(line, tokens);

      std::fill(elemNodes.begin(), elemNodes.end(), 0);
      elemNodeSet.clear();

      for (UInt i=0, n=numTok-3; i<n; i++) {
        elemNodes[i] = std::strtoul(tokens[i+3].c_str(), NULL, 0);
      }

      if (numElemNodes > 8) {
        GetNextLine(line);
        lCount++;
        numTok = SplitLine(line, tokens);

        for (UInt i=0, n=numTok-3; i<n; i++) {
          elemNodes[i+8] = std::strtoul(tokens[i+3].c_str(), NULL, 0);
        }
      }

      if (numElemNodes > 16) {
        GetNextLine(line);
        lCount++;
        numTok = SplitLine(line, tokens);

        for (UInt i=0, n=numTok-3; i<n; i++) {
          elemNodes[i+16] = std::strtoul(tokens[i+3].c_str(), NULL, 0);
        }
      }

      if (ansNacsMap_[elemType] != Elem::ET_UNDEF) {
        elemTypes_[elemNum] = ansNacsMap_[elemType];

        // Determine number of element nodes by inserting nodes into a set.
        elemNodeSet.insert(elemNodes.begin(), elemNodes.end());
        elemNodeSet.erase((UInt) 0);
        numElemNodes = elemNodeSet.size();

        maxNumElemNodes_ = numElemNodes > maxNumElemNodes_ ? numElemNodes : maxNumElemNodes_;

        Elem::FEType newFEType = DegenTypeToNativeType(elemTypes_[elemNum], numElemNodes);

        if(degen_) {
          DegenerateElement((Elem::FEType)elemTypes_[elemNum], newFEType, elemNodes);
        } else {
          ResortNodes(elemNodes);
        }

        elemTypes_[elemNum] = newFEType;
        elemNumsMap_[elemNum] = elemNum;
        elemDim = Elem::shapes[newFEType].dim;
        dim = dim < elemDim ? elemDim : dim;

        std::copy(elemNodes.begin(), elemNodes.begin() + Elem::shapes[newFEType].numNodes,
            std::back_inserter(topology_[elemNum]));

        elemNum++;

      } else {
        if (numSkipElems ==0)
          firstSkippedElem = elemNum;

        numSkipElems++;
      }
      if (elemNum%1000000 == 0) {
        std::cout << "Elements read " << elemNum << "\r";
        std::cout.flush();
      }

    }
    std::cout << std::endl;

    // store dimension to settings
    dim_ = dim;
    std::cout << "Finished reading elements (" << elemNum-1 << " read)" << std::endl;
    if (numSkipElems > 0)
      std::cout << numSkipElems << " ANSYS elements have been skipped, as no NACS counterpart exists" 
                << std::endl;

    // sanity check
    if (elemNum > firstSkippedElem)
      std::cout << "WARNING: non skipped elements found after first skipped element (= " 
                << firstSkippedElem << ").\nThis might lead to incorrect mesh transfer.\n" 
		<< "Please check carefully!" << std::endl;
  }

  void SimInputCDB::ReadRegionsAndGroups() {

    std::string line;
    std::string regnam,regtype;

    // std::cout << "Reading blocked " << linePtsCMBlocks_.size() << std::endl;
    // process all cmblocks
    if (linePtsCMBlocks_.size() > 0) 
      ReadBlockedRegionsAndGroups();

    // std::cout << "Reading unblocked " << linePtsCMCmnds_.size() << std::endl;
    if (linePtsCMCmnds_.size() > 0)
      ReadUnBlockedRegionsAndGroups();

    // in case of WB files, look for surface elements used for pressure loading
    // this is due to that fact that these are not contained in any named selection
    // std::cout << "Reading surface elems " << linePtsPSECmnds_.size() << std::endl;
    if (linePtsPSECmnds_.size() > 0)
      GenElGroupFromPrsSrfElems();

    if (numSFECmnds_ > 0)
      GenElGroupFromSurfForceElems();

  }

  void SimInputCDB::ReadBlockedRegionsAndGroups() {

    std::string line;
    std::ostringstream errMsg;
    std::string regnam,regtype;
    std::vector< std::string > tokens(32);

    std::vector<std::string>::const_iterator it, end;

    UInt numdata;
    int *dataVal;

    // loop on all found cmblock commands
    for (UInt ib=0; ib<linePtsCMBlocks_.size(); ib++) {

      if (!GetLine(line,linePtsCMBlocks_[ib])) {
        errMsg << "CMBLOCK " << ib << " expected to start at line no. " << linePtsCMBlocks_[ib]
               << " not found!\nLine read was " << line;
        EXCEPTION(errMsg);
      }

      // read region name, component type number of elems from CMBLOCK line
      UInt numTok = SplitLine(line, tokens, "!", NULL, true);

      regnam = tokens[1];
      regtype = tokens[2];      
      numdata = std::strtoul(tokens[3].c_str(), NULL, 0);

      dataVal = new int[numdata];

      // Don't skip format specs
      std::vector<int> offsets;
      GetNextLine(line);
      DecodeBlockFormatLine(line, offsets);

      // data spec
      UInt count = 0;
      while (count < numdata) {
        GetNextLine(line);
        numTok = SplitLine(line, tokens, "", &offsets);

        for (UInt i=0; i<numTok; i++) {
          dataVal[count] = std::strtoul(tokens[i].c_str(), NULL, 0);
          count++;
        }
      }

      if (regtype == "NODE" || regtype == "node") {
        StoreNodeGroup(regnam,numdata,dataVal);
      } else if (regtype == "ELEM" || regtype == "elem") {
        // test dimension of element type of first element in component
        // if this equals model dimension -> region
        // otherwise -> element group
        UInt dim = dim_;

        if (Elem::shapes[Elem::FEType(elemTypes_[dataVal[0]])].dim == dim) {
          StoreRegion(regnam,numdata,dataVal);
        } else {
          StoreElemGroup(regnam,numdata,dataVal);
        }
      }
      delete dataVal;
    }

  }

  void SimInputCDB::ReadUnBlockedRegionsAndGroups() {

    std::stringstream sstr;
    std::ostringstream errMsg;
    std::string line;
    std::string regnam,regtype;

    std::vector<std::string>::const_iterator it, end;

    for (UInt ib=0; ib<linePtsCMCmnds_.size(); ib++) {

      if (!GetLine(line,linePtsCMCmnds_[ib])) {
        errMsg << "CM command " << ib << " expected at line no. " << linePtsCMCmnds_[ib]
               << " not found!\nLine read was " << line;
        EXCEPTION(errMsg);
      }

      // possible region/group name
      sstr.clear();sstr.str("");
      UInt pos1 = (UInt) line.find(",")+1, pos2 = (UInt) line.rfind(",");
      sstr << line.substr(pos1,pos2-pos1);
      sstr >> regnam;

      // elem or node component
      if (line.find("ELEM",line.rfind(",")+1) != line.npos || 
                    line.find("elem",line.rfind(",")+1) != line.npos) {

        if (linePtsEselCmnds_.size() > 0)
          ReadUnBlockedElemRegionOrGroup(ib,regnam);

      } else if (line.find("NODE",line.rfind(",")+1) != line.npos || 
                    line.find("node",line.rfind(",")+1) != line.npos) {

        if (linePtsNselCmnds_.size() > 0)
          ReadUnBlockedNodeGroup(ib,regnam);

      }

    }
  }

  void SimInputCDB::GenElGroupFromPrsSrfElems() {

    std::stringstream sstr;
    std::ostringstream errMsg;
    std::string line;
    std::string regnam,regtype;

    std::vector<std::string>::const_iterator it, end;

    UInt elemNum;
    StdVector<UInt> elemNumbers;

    for (UInt ips=0; ips<linePtsPSECmnds_.size(); ips++) {

      if (!GetLine(line,linePtsPSECmnds_[ips])) {
        errMsg << "pressure comment line" << ips << " expected at line no. " << linePtsCMCmnds_[ips]
               << " not found!\nLine read was " << line;
        EXCEPTION(errMsg);
      }

      // find related eblock command
      for (UInt ib=1; ib<linePtsEBlocks_.size(); ib++) {
        if ( linePtsEBlocks_[ib] >= linePtsPSECmnds_[ips] &&
               linePtsEBlocks_[ib] < linePtsPSEStop_ ) {

          // process/get element numbers of eblock command
          if (!GetLine(line,linePtsEBlocks_[ib])) {
            errMsg << "eblock line" << ib << " expected at line no. " << linePtsEBlocks_[ib]
                   << " not found!\nLine read was " << line;
            EXCEPTION(errMsg);
          }
          // skipp format specs
          GetNextLine(line);

	  elemNumbers.Clear();
          GetNextLine(line);
          while(line.find("-1") == line.npos) {
            sstr.clear(); sstr.str("");
            sstr << line.substr(0,8);
            sstr >> elemNum;
	    elemNumbers.Push_back(elemNum);
            GetNextLine(line);
          }
          std::cout << "Found " << elemNumbers.GetSize() << " elements for eblock" << std::endl;
          
          GenerateSurfElGroup(ib, elemNumbers);
        }
      }
    }
  }

  void SimInputCDB::GenElGroupFromSurfForceElems() {

    std::stringstream sstr;
    std::ostringstream errMsg;
    std::string line;
    std::string regnam,regtype;

    std::vector<std::string>::const_iterator it, end;

    UInt sfeElemNum;
    StdVector<UInt> elemNumbers;
    std::string actSfeID="", lastSfeID="", strElemNum="";
    UInt ib=0;

    // in case that there is an format error (ANSYS pre 14.5), we need to know about
    // the last "volume" element to get the correct surface element numbers
    UInt dim = dim_;
    UInt numTotalElems = elemNumsMap_.size();
    UInt lastVolElem = numTotalElems;
    while ( (Elem::shapes[Elem::FEType(elemTypes_[lastVolElem])].dim < dim) && (lastVolElem > 0) ) {
      lastVolElem--;
    }
    std::cout << "Last volume element found is " << lastVolElem << std::endl;

    for (UInt ips=0; ips<numSFECmnds_; ips++) {

      GetLine(line,linePtsSFECmnds_[ips]);

      size_t pos1 = line.find(",");
      size_t pos2 = line.find(",",pos1+1);
      size_t pos3 = line.find(",",pos2+1);
      size_t pos4 = line.find(",",pos3+1);
      size_t pos5 = line.find(",",pos4+1);
      sstr.clear(); sstr.str("");

      // the following is to handle possible integer format overflows
      strElemNum = line.substr(pos1+1,pos2-pos1-1);
      if (strElemNum.substr(0,1) == "*") {
        sfeElemNum = lastVolElem + ips + 1;
      } else {
        sstr << strElemNum;
        sstr >> sfeElemNum;
      }

      sstr.clear(); sstr.str("");
      sstr << line.substr(pos5+1);
      sstr >> actSfeID;

      if (actSfeID != lastSfeID) {
        // generate SurfElGroup for lastSfeID (in case that one is already defined - ib>0)
        if (ib > 0) {
          GenerateSurfElGroup(ib, elemNumbers);
        }
        // reset data
	elemNumbers.Clear();
        // trigger new group
        lastSfeID = actSfeID;
	ib++;
      }

      // add element to new SurfElGroup
      elemNumbers.Push_back(sfeElemNum);
    }
    // do not forget the last defined SurfElGroup
    GenerateSurfElGroup(ib, elemNumbers);
  }

  void SimInputCDB::GenerateSurfElGroup(UInt elblock, StdVector<UInt> elemNumbers) {

    std::vector<bool> foundNodeGroup;
    std::set<UInt> elGrpNodeSet,nodeGrpNodeSet;

    for (UInt el=0; el < elemNumbers.GetSize(); el++) {
      UInt elemNum = elemNumbers[el];
      for (UInt elNode=0; elNode < topology_[elemNum].size() ; elNode++) {
        elGrpNodeSet.insert(topology_[elemNum][elNode]);
      }
    }

    foundNodeGroup.resize(numNodeGroups_);
    for (UInt nblock=0; nblock < numNodeGroups_; nblock++) {
      foundNodeGroup[nblock] = false;
      nodeGrpNodeSet.clear();
      for (UInt grpNode=0; grpNode<nodeGroupData_[nblock].GetSize() ; grpNode++) {
        nodeGrpNodeSet.insert(nodeGroupData_[nblock][grpNode]);
      }
      if (elGrpNodeSet == nodeGrpNodeSet)
        foundNodeGroup[nblock] = true;

    }

    UInt numFound = 0, grpFound = 0;
    for (UInt nblock=0; nblock < foundNodeGroup.size(); nblock++) {
      if (foundNodeGroup[nblock]) {
        std::cout << "Node group " << nodeGroupNames_[nblock] << " (grp. no. " << nblock 
		  << ") has been found" << std::endl;
        numFound++;
	grpFound = nblock;
      }
    }
    if (numFound != 1) {
        std::cout << "Case encountered which can not be handled: not excactly one node group has been found" << std::endl;
    } else {
      StoreElemGroup(nodeGroupNames_[grpFound], elemNumbers.GetSize(), elemNumbers);
      nodeGroupNames_[grpFound] = nodeGroupNames_[grpFound] + "_NODES";
    }
  }

//  void SimInputCDB::GenerateSurfElGroupV1(UInt elblock, std::vector<int> elemNumbers) {
//
//    std::vector<bool> foundNodeGroup;
//    bool foundMissingNode;
//    std::set<UInt> elGrpNodeSet,nodeGrpNodeSet;
//
//    for (UInt el=0; el < elemNumbers.size(); el++) {
//      UInt elemNum = elemNumbers[el];
//      for (UInt elNode=0; elNode < topology_[elemNum].size() ; elNode++) {
//        elGrpNodeSet.insert(topology_[elemNum][elNode]);
//      }
//    }
//
//    foundNodeGroup.resize(numNodeGroups_);
//    for (UInt nblock=0; nblock < numNodeGroups_; nblock++) {
//      foundNodeGroup[nblock] = true;
//
//      for (UInt el=0; el < elemNumbers.size(); el++) {
//        UInt elemNum = elemNumbers[el];
//        for (UInt elNode=0; elNode < topology_[elemNum].size() ; elNode++) {
//          if (foundNodeGroup[nblock]) {
//            bool nodeFound = false;
//            for (UInt grpNode=0; grpNode<nodeGroupData_[nblock].size() && !nodeFound; grpNode++) {
//              if ( topology_[elemNum][elNode] == nodeGroupData_[nblock][grpNode] ) {
//                nodeFound = true;
//              }
//            }
//            foundNodeGroup[nblock] = nodeFound;
//          }
//        }
//
//      }
//    }
//
//    UInt numFound = 0, grpFound;
//    for (UInt nblock=0; nblock < foundNodeGroup.size(); nblock++) {
//      if (foundNodeGroup[nblock]) {
//        std::cout << "Node group " << nodeGroupNames_[nblock] << " (grp. no. " << nblock 
//		  << ") has been found" << std::endl;
//        numFound++;
//	grpFound = nblock;
//      }
//    }
//    if (numFound != 1) {
//        std::cout << "Case encountered which can not be handled: not excactly one node group has been found" << std::endl;
//    } else {
//      StoreElemGroup(nodeGroupNames_[grpFound], elemNumbers.size(), elemNumbers);
//      nodeGroupNames_[grpFound] = nodeGroupNames_[grpFound] + "_NODES";
//    }
//  }

  void SimInputCDB::ReadUnBlockedElemRegionOrGroup(UInt numCmCmd, std::string regnam) {

    std::stringstream sstr;
    std::ostringstream errMsg;
    std::string line;

    std::vector<std::string>::const_iterator it, end;

    UInt elemNum, matNum;
    UInt numdata;
    StdVector<UInt> dataVal;

    dataVal.Clear();

    // look for all esel command between current cm command and aforegoing (if any)
    UInt cmLinePos = linePtsCMCmnds_[numCmCmd], cmLinePrevPos = 0;
    if (numCmCmd >= 1)
      cmLinePrevPos = linePtsCMCmnds_[numCmCmd-1];

    UInt numEselStart = 0, numEselStop = 0;
    for (UInt eb=0; eb<linePtsEselCmnds_.size(); eb++) {
      if (linePtsEselCmnds_[eb] < cmLinePos)
        numEselStop = linePtsEselCmnds_[eb];
    }
    for (UInt eb=linePtsEselCmnds_.size()-1; eb>0; eb--) {
      if (linePtsEselCmnds_[eb] > cmLinePrevPos)
        numEselStart = linePtsEselCmnds_[eb];
    }
    if (linePtsEselCmnds_[0] > cmLinePrevPos)
      numEselStart = linePtsEselCmnds_[0];

    if (numEselStop == 0) {
      return;
    }

    bool addGroupOrRegion = false;
    for (UInt eb=numEselStart; eb<=numEselStop; eb++) {
      GetLine(line,eb);

      // support for either esel,s|a,mat,, (wb) or esel,s|a,elem,,
      if (line.find("ESEL,S,MAT",0,10) != line.npos ||
            line.find("esel,s,mat",0,10) != line.npos ||
              line.find("ESEL,A,MAT",0,10) != line.npos ||
                line.find("esel,a,mat",0,10) != line.npos) {
        addGroupOrRegion = true;
        sstr.clear();sstr.str("");
        sstr << line.substr(line.rfind(",")+1);
        sstr >> matNum;
        for (elemNum=0; elemNum< elemMaterials_.size(); elemNum++) {
          if (elemMaterials_[elemNum] == matNum) {
            dataVal.Push_back(elemNum+1);
          }
        }
      }

      if (line.find("ESEL,S,ELEM",0,11) != line.npos ||
            line.find("esel,s,elem",0,11) != line.npos ||
              line.find("ESEL,A,ELEM",0,11) != line.npos ||
                line.find("esel,a,elem",0,11) != line.npos) {
        addGroupOrRegion = true;
        size_t pos1,pos2,pos3;
        UInt npos1=0,npos2=0,npos3=0;
	pos1 = line.find(",",11); npos1 = (UInt) pos1;
	pos1 = line.find(",",pos1+1); npos1 = (UInt) pos1;

	pos2 = line.find(",",npos1+1);
        if (pos2 != line.npos)
          npos2 = (UInt) pos2;

	pos3 = line.find(",",npos2+1);
        if (pos3 != line.npos)
          npos3 = (UInt) pos3;

	UInt eStart=0,eStop=0,eIncr=1;
        sstr.clear();sstr.str("");
	if (npos2>0)
          sstr << line.substr(npos1+1,npos2-npos1-1);
	else
          sstr << line.substr(npos1+1);
        sstr >> eStart;

	eStop = eStart;
	if (npos2>0) {
          sstr.clear();sstr.str("");
	  if (npos3>0)
            sstr << line.substr(npos2+1,npos3-npos2-1);
	  else
            sstr << line.substr(npos2+1);

          sstr >> eStop;
	  if (npos3>0) {
            sstr.clear();sstr.str("");
            sstr << line.substr(npos3+1);
            sstr >> eIncr;
          }
        }
	for (UInt elem=eStart; elem<=eStop; elem+=eIncr)
          dataVal.Push_back(elem);
      }

    }
    if (addGroupOrRegion) {
      // test dimension of element type of first element in component
      // if this equals model dimension -> region
      // otherwise -> element group
      UInt dim = dim_;
      
      numdata = dataVal.GetSize();
      if (Elem::shapes[Elem::FEType(elemTypes_[dataVal[0]])].dim == dim) {
        StoreRegion(regnam,numdata,dataVal);
      } else {
        StoreElemGroup(regnam,numdata,dataVal);
      }
    }
  }

  void SimInputCDB::ReadUnBlockedNodeGroup(UInt numCmCmd, std::string regnam) {

    std::stringstream sstr;
    std::ostringstream errMsg;
    std::string line;

    std::vector<std::string>::const_iterator it, end;

    UInt numdata;
    StdVector<UInt> dataVal;

    dataVal.Clear();

    // look for all esel command between current cm command and aforegoing (if any)
    UInt cmLinePos = linePtsCMCmnds_[numCmCmd], cmLinePrevPos = 0;
    if (numCmCmd >= 1)
      cmLinePrevPos = linePtsCMCmnds_[numCmCmd-1];

    UInt numNselStart = 0, numNselStop = 0;
    for (UInt eb=0; eb<linePtsNselCmnds_.size(); eb++) {
      if (linePtsNselCmnds_[eb] < cmLinePos)
        numNselStop = linePtsNselCmnds_[eb];
    }
    for (UInt eb=linePtsNselCmnds_.size()-1; eb>0; eb--) {
      if (linePtsNselCmnds_[eb] > cmLinePrevPos)
        numNselStart = linePtsNselCmnds_[eb];
    }
    if (linePtsNselCmnds_[0] > cmLinePrevPos)
      numNselStart = linePtsNselCmnds_[0];
    if (numNselStop == 0) {
      return;
    }

    bool addGroup = false;
    for (UInt nb=numNselStart; nb<=numNselStop; nb++) {
      GetLine(line,nb);

      if (line.find("NSEL,S,NODE",0,11) != line.npos ||
            line.find("nsel,s,node",0,11) != line.npos ||
              line.find("NSEL,A,NODE",0,11) != line.npos ||
                line.find("nsel,a,node",0,11) != line.npos) {
        addGroup = true;
        size_t pos1,pos2,pos3;
        UInt npos1=0,npos2=0,npos3=0;
	pos1 = line.find(",",11); npos1 = (UInt) pos1;
	pos1 = line.find(",",pos1+1); npos1 = (UInt) pos1;

	pos2 = line.find(",",npos1+1);
        if (pos2 != line.npos)
          npos2 = (UInt) pos2;

	pos3 = line.find(",",npos2+1);
        if (pos3 != line.npos)
          npos3 = (UInt) pos3;

	UInt nStart=0,nStop=0,nIncr=1;
        sstr.clear();sstr.str("");
	if (npos2>0)
          sstr << line.substr(npos1+1,npos2-npos1-1);
	else
          sstr << line.substr(npos1+1);
        sstr >> nStart;

	nStop = nStart;
	if (npos2>0) {
          sstr.clear();sstr.str("");
	  if (npos3>0)
            sstr << line.substr(npos2+1,npos3-npos2-1);
	  else
            sstr << line.substr(npos2+1);

          sstr >> nStop;
	  if (npos3>0) {
            sstr.clear();sstr.str("");
            sstr << line.substr(npos3+1);
            sstr >> nIncr;
          }
        }
	for (UInt node=nStart; node<=nStop; node+=nIncr)
          dataVal.Push_back(node);
      }
    }

    if (addGroup) {
      numdata = dataVal.GetSize();
      StoreNodeGroup(regnam,numdata,dataVal);
    }
  }

  void SimInputCDB::StoreNodeGroup(std::string grpname,
                                   UInt numdata,
                                   int* dataVal) {

    nodeGroupNames_.push_back(grpname);

    StdVector<UInt> tempVec;

    if (numdata==1) {
      tempVec.Push_back(dataVal[0]);
    } else {
      for (UInt i=0; i<numdata-1; i++) {
        int n1 = dataVal[i];
        int n2 = dataVal[i+1];
        // following cases:
        //   n1, n2 positive -> add dataVal[n1] only
        //   n1 positive, n2 negative -> add n1,......,abs(n2) only
        //   n1 negative -> end of previous n1,n2 case
        // --> n1 must always be positive, in order to add values

        if (n1 > 0) {
          if (n2 < 0) {
            for(int n=n1; n<=abs(n2); n++) {
              tempVec.Push_back(n);
            } 
          } else {
            tempVec.Push_back(n1);
	    if (i==numdata-2)
              tempVec.Push_back(n2);
          }
        }
      }
    }
    nodeGroupData_.push_back(tempVec);
    std::cout << "NodeGroup " << numNodeGroups_+1 << " contains "
              << nodeGroupData_[numNodeGroups_].GetSize() << " entries"
              << std::endl;
    numNodeGroups_++;
  }

  void SimInputCDB::StoreNodeGroup(std::string grpname, 
                                   UInt numdata,
                                   StdVector<UInt> dataVal) {

    nodeGroupNames_.push_back(grpname);

    nodeGroupData_.push_back(dataVal);
    std::cout << "NodeGroup " << numNodeGroups_+1 << " contains "
              << nodeGroupData_[numNodeGroups_].GetSize() << " entries"
              << std::endl;
    numNodeGroups_++;
  }

  void SimInputCDB::StoreElemGroup(std::string grpname,
                                   UInt numdata,
                                   int* dataVal) {

    elemGroupNames_.push_back(grpname);

    StdVector<UInt> tempVec;

    if (numdata==1) {
      tempVec.Push_back(dataVal[0]);
    } else {
      for (UInt i=0; i<numdata-1; i++) {
        int n1 = dataVal[i];
        int n2 = dataVal[i+1];
        // following cases:
        //   n1, n2 positive -> add dataVal[n1] only
        //   n1 positive, n2 negative -> add n1,......,abs(n2) only
        //   n1 negative -> end of previous n1,n2 case
        // --> n1 must always be positive, in order to add values

        if (n1 > 0) {
          if (n2 < 0) {
            for(int n=n1; n<=abs(n2); n++) {
              tempVec.Push_back(n);
            } 
          } else {
            tempVec.Push_back(n1);
	    if (i==numdata-2)
              tempVec.Push_back(n2);
          }
        }
      }
    }
    elemGroupData_.push_back(tempVec);
    std::cout << "ElemGroup " << numElemGroups_+1 << " contains "
              << elemGroupData_[numElemGroups_].GetSize() << " entries"
              << std::endl;
    numElemGroups_++;
  }

  void SimInputCDB::StoreElemGroup(std::string grpname,
                                   UInt numdata,
                                   StdVector<UInt> dataVal) {

    elemGroupNames_.push_back(grpname);

    elemGroupData_.push_back(dataVal);
    std::cout << "ElemGroup " << numElemGroups_+1 << " contains "
              << elemGroupData_[numElemGroups_].GetSize() << " entries"
              << std::endl;
    numElemGroups_++;
  }

  void SimInputCDB::StoreRegion(std::string grpname,
                                UInt numdata,
                                int* dataVal) {
    std::cout << "Region name: " << grpname << std::endl;
    regionNames_.push_back(grpname);

    StdVector<UInt> tempVec;

    if (numdata==1) {
      tempVec.Push_back(dataVal[0]);
    } else {
      for (UInt i=0; i<numdata-1; i++) {
        int n1 = dataVal[i];
        int n2 = dataVal[i+1];
        // following cases:
        //   n1, n2 positive -> add dataVal[n1] only
        //   n1 positive, n2 negative -> add n1,......,abs(n2) only
        //   n1 negative -> end of previous n1,n2 case
        // --> n1 must always be positive, in order to add values

        if (n1 > 0) {
          if (n2 < 0) {
            for(int n=n1; n<=abs(n2); n++) {
              tempVec.Push_back(n);
            } 
          } else {
            tempVec.Push_back(n1);
	    if (i==numdata-2)
              tempVec.Push_back(n2);
          }
        }
      }
    }
    regionData_.push_back(tempVec);
    std::cout << "Region " << grpname << " contains "
              << regionData_[numRegions_].GetSize() << " entries"
              << std::endl;
    for (UInt i=0; i<tempVec.GetSize(); i++) {
      // If an element is referenced in more than one domain region
      // throw an error, since this is not allowed.
      if (elemRegionMap_.find(tempVec[i]) != elemRegionMap_.end()) {
        EXCEPTION("Element nr. " << tempVec[i]
            << " is contained in more than one region. This affects "
            << "region " << regionNames_[numRegions_] << " and region "
            << regionNames_[elemRegionMap_[tempVec[i]]] );
      }
      elemRegionMap_[tempVec[i]] = numRegions_;
//      elemNumsOrig_[numRegions_].push_back(tempVec[i]);
    }
    numRegions_++;
  }

  void SimInputCDB::StoreRegion(std::string grpname,
                                UInt numdata,
                                StdVector<UInt> dataVal) {
    regionNames_.push_back(grpname);

    regionData_.push_back(dataVal);
    std::cout << "Region " << grpname << " contains "
              << regionData_[numRegions_].GetSize() << " entries"
              << std::endl;
    for (UInt i=0; i<dataVal.GetSize(); i++) {
      // If an element is referenced in more than one domain region
      // throw an error, since this is not allowed.
      if (elemRegionMap_.find(dataVal[i]) != elemRegionMap_.end()) {
        EXCEPTION("Element nr. " << dataVal[i]
            << " is contained in more than one region. This affects "
            << "region " << regionNames_[numRegions_] << " and region "
            << regionNames_[elemRegionMap_[dataVal[i]]] );
      }
      elemRegionMap_[dataVal[i]] = numRegions_;
//      elemNumsOrig_[numRegions_].push_back(dataVal[i]);
    }
    numRegions_++;
  }

  Elem::FEType SimInputCDB::DegenTypeToNativeType(UInt type, UInt numNodes)
  {
    Elem::FEType ret = (Elem::FEType) type;

    switch((Elem::FEType) type)
    {
    case Elem::ET_QUAD4: // rectangle
      switch(numNodes)
      {
      case 3:
        ret = Elem::ET_TRIA3;
        break;

      case 4:
        ret = Elem::ET_QUAD4;
        break;
      }
      break;

    case Elem::ET_QUAD8: // quad. rectangle
      switch(numNodes)
      {
      case 6:
        ret = Elem::ET_TRIA6;
        break;

      case 8:
        ret = Elem::ET_QUAD8;
        break;
      }
      break;

    case Elem::ET_HEXA8: // hexa
      switch(numNodes)
      {
      case 4:
        ret = Elem::ET_TET4;
        break;

      case 5:
        ret = Elem::ET_PYRA5;
        break;

      case 6:
        ret = Elem::ET_WEDGE6;
        break;

      case 8:
        ret = Elem::ET_HEXA8;
        break;
       }
       break;

     case Elem::ET_HEXA20: // quad. hexa
      switch(numNodes)
      {
        case 10:
          ret = Elem::ET_TET10;
          break;

        case 13:
          ret = Elem::ET_PYRA13;
          break;

        case 15:
          ret = Elem::ET_WEDGE15;
          break;

        case 20:
          ret = Elem::ET_HEXA20;
          break;
       }
       break;

     default:
       break;
    }

    return ret;
  }

#if(WIN32 || __MINGW32__)
  void SimInputCDB::OpenCDBFile(std::string fn)
  {
    std::string filename=fn.c_str();

    if (fopen_s( &inStream_, fn.c_str(), "rb" ) != 0 ) {
      EXCEPTION("Can't open " << filename);
    }

    int result;
    // Determine size of file
    result = _fseeki64( inStream_, 0, SEEK_SET);
    if ( result )
      EXCEPTION(" fseek failed positioning to beginning of file " << filename);

    __int64 pos1, pos2;

    pos1 = _ftelli64( inStream_ );

    result = _fseeki64( inStream_, 0, SEEK_END);
    if ( result )
      EXCEPTION(" fseek failed positioning to end of file " << filename);

    pos2 = _ftelli64( inStream_ );
    inSize_ = pos2 - pos1;

    // start from the beginning
    result = _fseeki64( inStream_, 0, SEEK_SET);
    if ( result )
      EXCEPTION(" fseek failed positioning to beginning of file " << filename);

  }

  void SimInputCDB::CloseCDBFile() { inFile_.close(); }

  bool SimInputCDB::GetNextLine(std::string& line)
  {
    __int64 pos = _ftelli64( inStream_ );

    if ( pos >= inSize_ )
      return false;

    UInt const buflen = 512;
    char buffer[buflen];

    if ( fgets ( buffer, buflen, inStream_ ) == NULL )
      EXCEPTION("fgets read error");

    line = buffer;
    return true;
  }

  bool SimInputCDB::GetLine(std::string& line, __int64 pos)
  {
    int const buflen = 512;
    char buffer[buflen];

    int result;

    result = _fseeki64( inStream_, pos, SEEK_SET);

    if ( result )
      EXCEPTION("fseek failed to find pos " << pos << " of cdb file" );

    if ( fgets ( buffer, buflen, inStream_ ) == NULL )
      EXCEPTION("fgets read error");

    line = buffer;

    return true;
  }

  __int64 SimInputCDB::GetFilePosition() {
      return _ftelli64( inStream_ );
  }

#else

  void SimInputCDB::OpenCDBFile(std::string fn)
  {
    std::string filename=fn.c_str();

    inFile_.clear();
    inFile_.open(fn.c_str(), std::ios::binary);

    if (!inFile_) {
      EXCEPTION("Can't open " << filename);
    }

    // Let's remember which files we opened, so that we can delete them in the end.
    // fileNames_.push_back(fn);

    // Determine size of file
    inFile_.seekg(0,std::ios::end);
    fSize_ = inFile_.tellg();

    // start from the beginning
    inFile_.seekg(0,std::ios::beg);
  }

  void SimInputCDB::CloseCDBFile() { inFile_.close(); }

  bool SimInputCDB::GetNextLine(std::string& line)
  {
    if (inFile_.tellg() >= fSize_)
      return false;

    UInt const buflen = 512;
    char buffer[buflen];

    inFile_.getline(buffer,buflen);
    line = buffer;
    //     std::cout << "CDB file completely read in" << std::endl;
    return true;
  }

  bool SimInputCDB::GetLine(std::string& line, std::streampos pos)
  {
    UInt const buflen = 512;
    char buffer[buflen];

    inFile_.seekg(pos,std::ios::beg);
    bool success = true;
    if (inFile_.tellg() >= fSize_)
      success = false;

    inFile_.getline(buffer,buflen);
    if (inFile_.bad())
      EXCEPTION("Read error in GetLine");

    line = buffer;

    return success;
  }

  unsigned long SimInputCDB::GetFilePosition() {
      return inFile_.tellg();
  }
#endif

  UInt SimInputCDB::SplitLine(const std::string& line,
                              std::vector< std::string >& tokens,
                              const std::string& addSplitChars,
                              std::vector<int>* chunkSizes,
                              bool trim,
                              const std::string& trimChars) const
  {   
    UInt i=0;

    // If  chunkSize   is  greater   zero  lets   split  by   constant  length
    // bits. Otherwise, we split along commas for APDL.
    if(chunkSizes) 
    {
      int* offsets = &(*chunkSizes)[0];
      boost::offset_separator f(offsets, offsets+(*chunkSizes).size());
      boost::tokenizer<boost::offset_separator> tok(line, f);

      for(boost::tokenizer<boost::offset_separator>::iterator it=tok.begin();
          it != tok.end();
          ++it, i++)
      {
        tokens[i] = *it;
      }
    }    
    else
    {
      typedef boost::tokenizer<boost::char_separator<char> > Tok;
      std::stringstream sstr;
      sstr << "," << addSplitChars;      
      boost::char_separator<char> sep(sstr.str().c_str());

      Tok t(line, sep);
      Tok::iterator it, end;

      it = t.begin();
      end = t.end();
      
      for( ; it != end; it++, i++) {
        tokens[i] = *it;        
      }      
    }
  
    if(trim) 
    {
      if(trimChars == "") 
      {
        for(UInt n=0; n<i; n++) {
          boost::trim(tokens[n]);
        }   
      }
      else
      {
        for(UInt n=0; n<i; n++) {
          boost::trim_if(tokens[n], boost::is_any_of(trimChars));
        }   
      }
    }
    
    return i;
  }
  
  UInt SimInputCDB::DecodeBlockFormatLine(const std::string& line,
                                          std::vector<int>& chunkSizes) const
  {
    std::vector< std::string > tokens(32);
    
    chunkSizes.clear();
    
    // retrieve format by splitting a line of the form (3i8,6e20.13) or (19i8)
    UInt numTok = SplitLine(line, tokens);

    std::string intLenStr = tokens[0];
    std::string floatLenStr = "";
    if(numTok > 1) 
    {
      floatLenStr = tokens[1];
    }
    
    numTok = SplitLine(intLenStr, tokens, "(iI)");
    UInt numInts = std::strtoul(tokens[0].c_str(), NULL, 0);
    UInt intWidth = std::strtoul(tokens[1].c_str(), NULL, 0);

    // Generate an offset vector for the different fields.
    for(UInt i=0; i<numInts; i++) 
    {
      chunkSizes.push_back(intWidth);
    }
    
    if(floatLenStr != "") 
    {
      numTok = SplitLine(floatLenStr, tokens, "eE.");
      UInt numFloats = std::strtoul(tokens[0].c_str(), NULL, 0);
      UInt floatWidth = std::strtoul(tokens[1].c_str(), NULL, 0);
    
      for(UInt i=0; i<numFloats; i++) 
      {
        chunkSizes.push_back(floatWidth);
      }      
    }
    
    return numInts;
  }
  
}
