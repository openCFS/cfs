// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <string>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "SimInputMESH.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

// declare class specific logging stream
DEFINE_LOG(mesh, "mesh")

namespace CoupledField {


  SimInputMESH::SimInputMESH(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode) :
    SimInput(fileName, inputNode, infoNode )
  {
    capabilities_.insert( SimInput::MESH);
  }

  SimInputMESH::~SimInputMESH()
  {
    inFile_.close();
  }

  void SimInputMESH::InitModule()
  {
    inFile_.open( fileName_.c_str(), std::ios::binary );
    if (!inFile_.good()) 
      throw Exception("Cannot open mesh file '" + fileName_ + "'");

    // parse header
    // [Info]
    // Version 1
    // ...
    // Num 2d-line      : 16
    std::string line;
    while(std::getline(inFile_, line)) 
     if(line.find("[Info]") != std::string::npos)
       break;    
    if(!inFile_.good() && line.find("[Info]"))   
      throw Exception("Input file is no valid openCFS mesh file: " + fileName_);
    // first 9 lines: strict "Key Value" format
    for(int i = 0; i < 9; i++) 
    {
      std::getline(inFile_, line);
      std::istringstream ss(line);
      std::string key; unsigned int val;
      if (!(ss >> key >> val))
        throw Exception("Unexpected format in [Info] header line " + std::to_string(i+1) + ": " + line);
      header_[key] = val;
      LOG_DBG(mesh) << "IM: header '" << line << "' -> " << key << " = " << val;
    }

    // skip lines in form of Num 2d-line,quad : 0 until empty line
    while(line.size() > 1) // empty or Windows empty line
      std::getline(inFile_, line);

    dim_ = header_.at("Dimension");
    maxNumNodes_ = header_.at("NumNodes");
    maxNumElems_ = header_.at("Num3DElements") + header_.at("Num2DElements") + header_.at("Num1DElements");

    LOG_DBG(mesh) << "IM: header of file '" << fileName_ << "' nodes=" << maxNumNodes_ << " elems=" << maxNumElems_ << " pos=" << inFile_.tellg(); 
  }

  void SimInputMESH::ReadMesh(Grid *grid)
  {
    LOG_DBG(mesh) << "RM: pos=" << inFile_.tellg();
    mi_ = grid;

    if(grid->GetDim() != dim_)
      throw Exception("Mesh '" + fileName_ + "' has dimension " + std::to_string(dim_) + " but grid has " + std::to_string(grid->GetDim()));

    // the header is read and we continue
    // [Nodes]
    // #NodeNr x-coord y-coord z-coord
    // 1  0.0  0.0  0.0
    // 2  0.25  0.0  0.0    
    std::string line;
    std::getline(inFile_, line); // consume "[Nodes]"
    assert(line == "[Nodes]");
    std::getline(inFile_, line); // consume "#NodeNr x-coord y-coord z-coord"
    assert(line[0] == '#');

    grid->AddNodes(maxNumNodes_); // create empty nodes, by reading we overwrite

    Vector<double> loc(3);

    for (unsigned int i = 0; i < maxNumNodes_; i++)
    {
      std::getline(inFile_, line);
      const char* p = line.data();
      char* end;

      // always expect one character spacing
      // parse the node number and coords the classic way as std::strtod() is not that stable
      [[maybe_unused]] unsigned int nodeNr = std::strtoul(p, &end, 10); p = end + 1; 
      loc[0] = std::strtod(p, &end); p = end + 1;
      loc[1] = std::strtod(p, &end); p = end + 1;
      loc[2] = std::strtod(p, &end);
      
      assert(nodeNr == i+1); // shall be as expected, otherwise the elements will fail
      grid->SetNodeCoordinate(i+1, loc); // for GridCFS simply sets coords_
    }
    
    // [1D Elements]
    // #ElemNr  ElemType  NrOfNodes  Level
    // #Node1 Node2 ... NodeNrOfNodes
    //
    // [2D Elements]
    // #ElemNr  ElemType  NrOfNodes  Level
    // #Node1 Node2 ... NodeNrOfNodes
    // 1 6 4 mech
    // 1 2 7 6
    
    // create empty elements, by reading we overwrite
    grid->AddElems(maxNumElems_); 

    // we parse the Blocks [1D Elements], [2D Elements], [3D Elements] - old test with reverse order :(
    for (int i = 0; i < 3; i++)
      ParseElementBlock();

    auto entities = ParseNamedEntities("[Node BC]",header_.at("NumNodeBC"));
    for(auto& pair : entities)
      grid->AddNamedNodes(pair.first, pair.second);

    entities = ParseNamedEntities("[Save Nodes]",header_.at("NumSaveNodes"));
    for(auto& pair : entities)
      grid->AddNamedNodes(pair.first, pair.second);

    entities = ParseNamedEntities("[Save Elements]",header_.at("NumSaveElements"));
    for(auto& pair : entities)
      grid->AddNamedElems(pair.first, pair.second);
  }

  void SimInputMESH::ParseElementBlock()
  {
    Grid* grid = mi_;
    // [2D Elements]
    // #ElemNr  ElemType  NrOfNodes  Level
    // #Node1 Node2 ... NodeNrOfNodes
    // 1 6 4 mech
    // 1 2 7 6

    std::string line;
    // skip empty lines - handle also Windows
    while(std::getline(inFile_, line) && line.size() < 2)
      if(!inFile_.good())
        throw Exception("Unexpected end of file searching [*D Elements] block in " + fileName_);
    boost::algorithm::trim_right(line);
    // new files are [1D Elements], ..., [3D Elements] - older files are reversed
    if(line != "[1D Elements]" && line != "[2D Elements]" && line != "[3D Elements]")
      throw Exception("Unexpected mesh file format. Got '" + line + "' expected [1/2/3D Elements]");
    std::string key = std::string("Num") + line[1] + "DElements";
    assert(header_.find(key) != header_.end());
    unsigned int num = header_.at(key);
    LOG_DBG(mesh) << "PEB: block=" << line << " num=" << num << " at pos=" << inFile_.tellg();

    std::getline(inFile_, line); // consume "#ElemNr  ElemType  NrOfNodes  Level"
    assert(line[0] == '#');
    std::getline(inFile_, line); // consume "Node1 Node2 ... NodeNrOfNodes"
    assert(line[0] == '#');
   
    // with all regions at hand we can now parse the block again
    Vector<unsigned int> connectivity; 
    for(unsigned int i = 0; i < num; i++)
    {
      std::getline(inFile_, line); // 1 6 4 mech
      boost::algorithm::trim_right(line); // remove possible Windows endings, fast if nothing to do
      const char* p = line.data();

      unsigned int elemNr = 0, elemType = 0, numNodes = 0;
      p = std::from_chars(p, p+20, elemNr).ptr + 1; // skip single space
      p = std::from_chars(p, p+20, elemType).ptr + 1;
      p = std::from_chars(p, p+20, numNodes).ptr + 1;
      std::string region(p, line.data() + line.size() - p);

      // read the connectivity line
      std::getline(inFile_, line); // 1 2 7 6
      p = line.data();
      connectivity.Resize(numNodes);
      for(unsigned int n = 0; n < numNodes; n++) 
        p = std::from_chars(p, p+20, connectivity[n]).ptr + 1; // skip single space

      RegionIdType regid = grid->GetRegionId(region, true); // silent for new region  
      if(regid == NO_REGION_ID) 
        regid = grid->AddRegion(region); // populate just the name
      
      // we could slightly speed up SetElemData() by handling the entityDim_ ourselves 
      Elem::FEType cfs_elem = AnsysType2ElemType(elemType);
      grid->SetElemData(elemNr, cfs_elem, regid, connectivity.GetPointer());

    }
  } 

  StdVector<std::pair<std::string, StdVector<unsigned int>>> SimInputMESH::ParseNamedEntities(const std::string& label, unsigned int num)
  {
    std::string line;
    // skip empty lines - handle also Windows
    while(std::getline(inFile_, line) && line.size() < 2)
      if(!inFile_.good())
        throw Exception("Unexpected end of file searching [*D Elements] block in " + fileName_);
    assert(line == label); // e.g. "[Save Elements]"
    std::getline(inFile_, line); // consume "#NodeNr Name"
    assert(line[0] == '#');

    // we parse stuff like 
    // 1 bottom
    // 2 bottom
    // 78 top

    // we cannot use a map as we want to preserve the order of the names for tests
    StdVector<std::pair<std::string, StdVector<unsigned int>>> entities;
    std::string lastName;
    unsigned last_idx = 1234;
    
    for (unsigned int i = 0; i < num; i++)
    {
      std::getline(inFile_, line);
      boost::algorithm::trim(line); // old tests have left padding
      const char* p = line.data();
      unsigned int nr = 0; // 1-based
      p = std::from_chars(p, p+20, nr).ptr + 1; // skip single space
      std::string name(p, line.data() + line.size() - p);
      boost::algorithm::trim_left(name); // handle legacy test case
      if (name != lastName) 
      {
        // a change of name does not guarantee that the name is new
        int idx = -1;
        for(unsigned i = 0; i < entities.GetSize(); i++)
          if(entities[i].first == name) { 
            idx = i; 
            last_idx = i;
            break; 
          }        
        if(idx == -1)
        {
          entities.Push_back(std::make_pair(name, StdVector<unsigned int>()));
          entities.Last().second.Reserve(num-i); // too big if we have many names
          last_idx = entities.GetSize() - 1;
        }
        lastName = name;
      }
      entities[last_idx].second.Push_back(nr);
    }
    return entities;
  }

  Elem::FEType SimInputMESH::AnsysType2ElemType(unsigned int itype)
  {
    switch( itype ) 
    {
    case 101:
      return Elem::ET_LINE3;
    case 100:
      return Elem::ET_LINE2;
    case 4:
      return Elem::ET_TRIA3;
    case 5:
      return Elem::ET_TRIA6;
    case 6:
      return Elem::ET_QUAD4;
    case 7:
      return Elem::ET_QUAD8;
    case 107:
      return Elem::ET_QUAD9;
    case 8:
      return Elem::ET_TET4;
    case 9:
      return Elem::ET_TET10;
    case 10:
      return Elem::ET_HEXA8;
    case 11:
      return Elem::ET_HEXA20;
    case 111:
      return Elem::ET_HEXA27;
    case 12:
      return Elem::ET_PYRA5;
    case 13:
      //    (*warning) << "Pyram. with quadratic shape functions" << 
      //    "do not work well for some cases "<< 
      //    "(i.e. electric field intensity). Please verify your results."; 
      return Elem::ET_PYRA13;
    case 14:
      return Elem::ET_WEDGE6;
    case 15:
      return Elem::ET_WEDGE15;
    default:
      break;
    }
    
    // This place should never be reached!
    return Elem::ET_UNDEF;
  }

}
