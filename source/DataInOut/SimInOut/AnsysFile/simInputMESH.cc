// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iterator>
#include <string>

#include "DataInOut/simInput.hh"
#include "Domain/grid.hh"
#include "General/exception.hh"
#include "Utils/Point.hh"
#include "Utils/StdVector.hh"
#include "simInputMESH.hh"
#include "stdarg.h"

using namespace std;

namespace CoupledField {


    SimInputMESH::SimInputMESH(std::string fileName, PtrParamNode inputNode) :
        SimInput(fileName, inputNode), numElemsDim_(3)
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
    if ( !inFile_.good() ) {
      EXCEPTION("I am unable to open mesh file " << fileName_);
    }
    
    ReadHeader();
  }
  
  void SimInputMESH::ReadHeader() {
    string line;
    getline(inFile_, line, '\n');
    if(line != "[Info]"){
      EXCEPTION("The mesh file does not start with the correct header!");
    }
    boost::trim(line);
    while(!line.empty()){
      getline(inFile_, line, '\n');
      boost::trim(line);
      if(boost::starts_with(line, "Dimension ")){
        dim_ = boost::lexical_cast<UInt>(line.substr(10));
      }else if(boost::starts_with(line, "Num")){
        line = line.substr(3);
        if(line == " "){
          continue;
        }
        if(boost::starts_with(line, "Nodes ")){
          numNodes_ = boost::lexical_cast<UInt>(line.substr(6));
        }else if(boost::starts_with(line, "3DElements ")){
          numElemsDim_[2] = boost::lexical_cast<UInt>(line.substr(11));          
        }else if(boost::starts_with(line, "2DElements ")){
          numElemsDim_[1] = boost::lexical_cast<UInt>(line.substr(11));          
        }else if(boost::starts_with(line, "1DElements ")){
          numElemsDim_[0] = boost::lexical_cast<UInt>(line.substr(11));          
        }else if(boost::starts_with(line, "NodeBC ")){
          numNodeBC_ = boost::lexical_cast<UInt>(line.substr(7));
        }else if(boost::starts_with(line, "SaveNodes ")){
          numSaveNodes_ = boost::lexical_cast<UInt>(line.substr(10));
        }else if(boost::starts_with(line, "SaveElements ")){
          numNamedElems_ = boost::lexical_cast<UInt>(line.substr(13));
        }
      }
    }
    // we have read the header and are now positioned at [Nodes]
  }
    
  void SimInputMESH::ReadMesh(Grid *mi) {
    string section;
    
    UInt numElems = GetNumElems();
    UInt numNodes = GetNumNodes();

    if(mi){ // we can be called with NULL (from XMLSkeleton), we then read the mesh but do not write anything to the grid
      mi_ = mi;
      mi_->AddElems(numElems);
      mi_->AddNodes(numNodes);
    }
    
    StdVector<StdVector<UInt> > nodeIndices;
    StdVector<StdVector<UInt> > elemIndices;
    while(true){
      section = ForwardToSection(); // Forward to the next section
      if(section.empty()){
        break;
      }else if(section == "[Nodes]"){ // read Nodes
        int ibuf;
        double x,y,z;
        for(UInt i = 1; i <= numNodes; i++){
          inFile_ >> ibuf >> x >> y >> z; // nodes are number(ignored) x y z coordinate (numNodes times)
          if(mi){
            mi_->SetNodeCoordinate(i, Point(x, y, z));
          }
        }        
      }else if(section.length() > 2 && section.substr(2) == "D Elements]"){ // read elements
        int d = lexical_cast<int>(section[1]); // get dimension from section heading 
        UInt numElemsD = GetNumElems(d);
        UInt eNum, eType, eNodes, node;
        string reg, lastreg;
        RegionIdType regionId(NO_REGION_ID);
        vector<UInt> nodes;
        for(UInt i = 0; i < numElemsD; i++){
          inFile_ >> eNum >> eType >> eNodes >> reg; // first line for an element is elemNum elemType noNodes region
          if(reg != lastreg){ // if we have a different region compared to the last element, fetch a new regionId. 
            // Note the ObtainRegionId returns the correct regionId if we already had this region, and also adds the region to the grid
            regionId = ObtainRegionId(reg, d, mi_);
            lastreg = reg;
          }
          for(UInt j = 0; j < eNodes; j++){ // second line is noNodes nodeidxs
            inFile_ >> node;
            nodes.push_back(node);
          }
          if(mi){
            mi_->SetElemData(eNum, AnsysType2ElemType(eType), regionId, &nodes[0]);
          }
          nodes.resize(0);
        }          
      }else if(section == "[Node BC]" || section == "[Save Nodes]" || section == "[Save Elements]"){
        // these sections always are node/element index and a name
        // we have to collect all names for nodes and elements and add them to the grid in one call
        bool elems = section == "[Save Elements]";
        StdVector<string>& names = elems ? elemNames_ : nodeNames_; // this is a vector of the node/elem names
        StdVector<StdVector<UInt> >& indices = elems ? elemIndices : nodeIndices; // this is a vector of vectors of UInt which correspond to the node/elem name above
        UInt num = elems ? numNamedElems_ : (section == "[Save Nodes]" ? numSaveNodes_ : numNodeBC_);
        string name, lastname;
        int nameidx(-1);
        UInt idx;
        for(UInt i = 0; i < num; i++){
          inFile_ >> idx >> name;
          if(name != lastname){ // if name is the same as last time, simply keep nameidx
            nameidx = names.Find(name); // can return -1
            if(nameidx < 0){
              names.Push_back(name);
              nameidx = names.GetSize()-1;
            }
            lastname = name;
            if(nameidx >= (int)indices.GetSize()){ // if new name added, also add new index vector
              indices.Push_back(StdVector<UInt>());
            }
          }
          indices[nameidx].Push_back(idx);
        }
      }
      if(!inFile_.good()){
        EXCEPTION("Error while reading section " << section << "! Malformed data!")
      }
      string tmp;
      getline(inFile_, tmp, '\n');
      boost::trim(tmp);
      if(!tmp.empty()){
        EXCEPTION("Error while reading section " << section << "! Did not find an empty line after the data. Check whether the Number of (elements/nodes) given in the header is correct!");
      }
    } // while(true)
    if(mi){
      unsigned int n = nodeNames_.GetSize();
      for(unsigned int i = 0; i < n; i++){
        mi_->AddNamedNodes(nodeNames_[i], nodeIndices[i]);
      }
      n = elemNames_.GetSize();
      for(unsigned int i = 0; i < n; i++){
        mi_->AddNamedElems(elemNames_[i], elemIndices[i]);
      }
    }
  }

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputMESH::GetDim() {
    return(dim_);
  }
  
  UInt SimInputMESH::GetNumNodes(){
    return(numNodes_);
  }
    
  UInt SimInputMESH::GetNumElems(const Integer dim){
    UInt numElems = 0;
    // 1.) return number of all elements
    if(dim == 0){
      if(GetDim() == 3){
        numElems += GetNumElems(3);
      }
      numElems += GetNumElems(2);
      numElems += GetNumElems(1);
    }else if ( dim >=1 && dim <= 3 ){
      numElems = numElemsDim_[dim-1];
    }else{
      EXCEPTION("Dimension " << dim << " is out of range!");
    }
    return(numElems);
  }
  
  UInt SimInputMESH::GetNumRegions(){
    return(regions_.GetSize());
  }

  UInt SimInputMESH::GetNumNamedNodes(){
    return(numNodeBC_ + numSaveNodes_);
  }

  UInt SimInputMESH::GetNumNamedElems(){
    return(numNamedElems_);
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================
    
  void SimInputMESH::GetAllRegionNames( StdVector<std::string> & regionNames ){
    unsigned int n = regions_.GetSize();
    for(unsigned int i = 0; i < n; i++){
      regionNames.Push_back(regions_[i].name);
    }
  }
  
  void SimInputMESH::GetRegionNamesOfDim( StdVector<std::string> & regionNames, const UInt dim ) {
    regionNames.Clear();
    
    if(regions_.GetSize() == 0){ // read mesh if necessary, as this is also called from XMLSkeletonConf, everything else is used after InitModule and ReadMesh
      ReadMesh(NULL);
      // seeking is not really necessary, but if someone once gets the idea to create a grid after having called this from XMLSkeleton, it would still work
      inFile_.seekg(0, ios_base::beg);
    }
    
    unsigned int n = regions_.GetSize();
    for(unsigned int i = 0; i < n; i++){
      RegionInfo& ri = regions_[i];
      if(ri.dim == dim){
        regionNames.Push_back(ri.name);
      }
    }
  }

  void SimInputMESH::GetNodeNames( StdVector<std::string> &nodeNames ) {
    nodeNames = nodeNames_;
  }

  void SimInputMESH::GetElemNames( StdVector<std::string> & elemNames ) {
    elemNames = elemNames_;
  }

  RegionIdType SimInputMESH::ObtainRegionId(const std::string & regionName, const UInt dim, Grid* grd) {
    size_t n = regions_.GetSize();
    for(size_t i = 0; i < n; i++){
      if(regions_[i].name == regionName){
        return(grd ? regions_[i].id : i);
      }
    }
    RegionInfo ri;
    ri.name = regionName;
    ri.dim = dim;
    ri.id = grd ? grd->AddRegion(regionName) : NO_REGION_ID;
    regions_.Push_back(ri);
    return(grd ? regions_.GetSize()-1 : ri.id);
  }

  Elem::FEType SimInputMESH::AnsysType2ElemType( const UInt itype ) {
    switch( itype ) {
    case 101:
      return Elem::LINE3;
    case 100:
      return Elem::LINE2;
    case 4:
      return Elem::TRIA3;
    case 5:
      return Elem::TRIA6;
    case 6:
      return Elem::QUAD4;
    case 7:
      return Elem::QUAD8;
    case 107:
      return Elem::QUAD9;
    case 8:
      return Elem::TET4;
    case 9:
      return Elem::TET10;
    case 10:
      return Elem::HEXA8;
    case 11:
      return Elem::HEXA20;
    case 111:
      return Elem::HEXA27;
    case 12:
      return Elem::PYRA5;
    case 13:
      //    (*warning) << "Pyram. with quadratic shape functions" << 
      //    "do not work well for some cases "<< 
      //    "(i.e. electric field intensity). Please verify your results."; 
      return Elem::PYRA13;
    case 14:
      return Elem::WEDGE6;
    case 15:
      return Elem::WEDGE15;
    default:
      EXCEPTION("Undefined element in mesh: " << itype);
      break;
    }
    
    // This place should never be reached!
    return Elem::UNDEF;
  }
  
  string SimInputMESH::ForwardToSection(){
    string line, tmp;
    do{
      if(!getline(inFile_, line, '\n')){
        return("");
      }
      boost::trim(line);
    }while(line.empty() || (line.length() > 2 && line[0] != '[' && line.substr(line.length()-1) != "]"));
    // skip all comment lines
    while(inFile_.peek() == '#'){
      getline(inFile_, tmp, '\n');
    }
    return(line);
  }
  
}
