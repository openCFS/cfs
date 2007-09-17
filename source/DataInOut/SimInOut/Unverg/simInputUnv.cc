
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "simInputUnv.hh"
#include "unv_if.hh"

extern const char *nodeDataTypesStr[30];
extern const char *elemDataTypesStr[10];

namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simInputUNV, "SimInputUNV")

  SimInputUnv::SimInputUnv( std::string fileName, ParamNode * inputNode ) 
    : SimInput(fileName, inputNode) {
    capabilities_.insert( SimInput::MESH);
    capabilities_.insert( SimInput::MESH_RESULTS);
  }
  
  void SimInputUnv::InitModule()
  {
    // read in mesh data, in order to provide information for
    // dim, number nodes, number of elements, number of regions etc.
  }

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputUnv::GetDim() {
   return 2;
  }
  
  UInt SimInputUnv::GetNumNodes(){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }
    
  UInt SimInputUnv::GetNumElems(const Integer dim){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }
  
  UInt SimInputUnv::GetNumRegions(){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }

  UInt SimInputUnv::GetNumNamedNodes(){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }

  UInt SimInputUnv::GetNumNamedElems(){
    CoupledField::Warning("SimInputUnv::ReadMesh() not implemented");
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputUnv::GetAllRegionNames( StdVector<std::string> & regionNames ){
    std::cerr << "SimInputUnv::ReadMesh() not implemented\n";
  }

  void SimInputUnv::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                   const UInt dim )
  {
    EXCEPTION("SimInputUnv::GetRegionNamesOfDim() not implemented");
  }
    
  void SimInputUnv::GetNodeNames( StdVector<std::string> & nodeNames )
  {
    EXCEPTION("SimInputUnv::GetNodeNames() not implemented");
  }
  
  void SimInputUnv::GetElemNames( StdVector<std::string> & elemNames )
  {
    EXCEPTION("SimInputUnv::GetElemNames() not implemented");
  }

  void
  SimInputUnv::ReadMesh( Grid* mi )
  {
    //!! TODO also extract surface meshes
    mi_ = mi;
    
    LOG_TRACE(simInputUNV) << "reading base mesh from file " << fileName_ << ":";
    
    long n,set;
    CapaInterfaceC capaInterface;
    GDataInfo      datainfo;
    
    if(capaInterface.ReadUniversalfile(fileName_.c_str()) < 0)
    {
      EXCEPTION("Could not read universal file " << fileName_);
      return;
    }
    
    LOG_TRACE(simInputUNV) << "Finished reading";
    
    int nnumtimesteps=0;
    int nnodes=capaInterface.GetNumOfNodes();
    int nelems=capaInterface.GetNumOfElements();
    int dim=capaInterface.GetDimension();
    int nset55=capaInterface.GetNumOfNodeDataSets();
    int nset56=capaInterface.GetNumOfElemDataSets();

    capaInterface.GetDataInfo(datainfo);
 #if 0
    for (set=0; set<nset55; set++) {
      std::cout << "*** node data set: " << set << std::endl;
      std::cout << " datainfo.types55[set] " << datainfo.types55[set] << " => " << nodeDataTypesStr[datainfo.types55[set]] << std::endl;
      int setid = datainfo.Nsetinfo[0][set].idx;
      int edata = datainfo.Nsetinfo[0][set].ndata;
      std::cout << " idx   = " << setid <<std::endl;
      std::cout << " ndata = " << edata <<std::endl;

    }
    for (set=0; set<nset56; set++) {
      std::cout << "*** element data set: " << set << std::endl;
      std::cout << " datainfo.types56[set] " << datainfo.types56[set] << " => " << elemDataTypesStr[datainfo.types56[set]] << std::endl;
      int setid = datainfo.Esetinfo[0][set].idx;
      int edata = datainfo.Esetinfo[0][set].ndata;
      std::cout << " idx   = " << setid <<std::endl;
      std::cout << " ndata = " << edata <<std::endl;
    }
 #endif

    LOG_TRACE(simInputUNV) << "Number of nodes            : "<<nnodes;
    LOG_TRACE(simInputUNV) << "Number of elements         : "<<nelems;
    LOG_TRACE(simInputUNV) << "Number of node data sets   : "<<nset55;
    LOG_TRACE(simInputUNV) << "Number of element data sets: "<<nset56;
    

    // get coordinates
    LOG_TRACE(simInputUNV) << "reading vertex coordinates";

    Point  p;
    long   node,elem;
    bool   isDimPlausible = false;
    //    ocs::Vec3d oldNormal;

    // Check if dim is really sane!
    for (node=1; node<=nnodes; node++) {
      capaInterface.GetPos(node, &p[0]);

      if(p[0] != 0)
      {
        dim = 3;
        break;
      }
    }

    mi_->AddNodes(nnodes);

    for (node=1; node<=nnodes; node++) {
      capaInterface.GetPos(node, &p[0]);

      // In the 2D case CFS++ writes the 3D coordinates in the
      // following form: 0 x y
      if(dim == 2)
      {
        p[0] = p[1];
        p[1] = p[2];
        p[2] = 0;
      }

      mi_->SetNodeCoordinate(node, p);
    }

    

    // get elements
    long elemColor, unvElemType, numElemNodes, elemNodes[32], dummy;
    UInt elemNodesUInt[32];
    std::vector< FEType >  element_types(nelems);
    std::vector< UInt > num_vertices_of_elements(nelems);
    std::vector< UInt > element_partition(nelems);
    std::set<UInt> partitions;
    FEType eType;
    std::stringstream strBuffer;


    mi_->AddElems(nelems);

    long totalidxs = 0;
    for (n=0; n<nelems; ++n) {
      capaInterface.GetElemColor(n+1, elemColor);
      capaInterface.GetElemType(n+1, unvElemType);
      capaInterface.GetElemNodes(n+1, numElemNodes, elemNodes);
      std::copy(elemNodes, elemNodes+32, elemNodesUInt);

      eType = UnvType2ElemType(unvElemType);

      if(dim == 2)
      {
        switch(eType)
        {
        case ET_TRIA6:
          dummy = elemNodes[1];
          elemNodes[1] = elemNodes[2];
          elemNodes[2] = elemNodes[4];
          elemNodes[4] = elemNodes[3];
          elemNodes[3] = dummy;
          break;
        case ET_QUAD8:
          dummy = elemNodes[1];
          elemNodes[1] = elemNodes[2];
          elemNodes[2] = elemNodes[4];
          elemNodes[4] = elemNodes[3];
          elemNodes[3] = elemNodes[6];
          elemNodes[6] = elemNodes[5];
          elemNodes[5] = elemNodes[4];
          elemNodes[4] = dummy;
          break;
        }
      }        

      mi_->SetElemData(n+1, eType, elemColor-1,
                       elemNodesUInt);

      partitions.insert(elemColor);
    }

    std::vector<std::string> regionNames;
    std::vector<UInt> regionDims;
    std::vector<RegionIdType> regionIds;
    std::set<UInt>::iterator pit, pend;

    pit = partitions.begin();
    pend = partitions.end();

    for(; pit != pend ; pit++)
    {
      strBuffer.str("");
      strBuffer.clear();
      strBuffer << "Partition" << (*pit);
      regionNames.push_back(strBuffer.str());
      regionDims.push_back(3);
    }
    mi_->AddRegions(regionNames, regionIds);

    std::string attrname;
    double data1[32], data2[32];

    for(UInt k=0; k<datainfo.numtype55; k++)
    {
      attrname = nodeDataTypesStr[datainfo.types55[k]];
      std::cout << attrname;

      for (UInt set=0; set<datainfo.n55[k]; set++)
      {
        std::cout << " idx " << datainfo.Nsetinfo[k][set].idx
                  << " time " << datainfo.Nsetinfo[k][set].time
                  << " ndata " << datainfo.Nsetinfo[k][set].ndata << std::endl;

        long ndata = datainfo.Nsetinfo[k][set].ndata;
        for (UInt node=1; node<=nnodes; node++) {
          capaInterface.GetNodeData(node,
                                    datainfo.Nsetinfo[k][set].idx,
                                    data1,
                                    data2);

          for (UInt n=0; n<ndata; n++) {
            std::cout << data1[n] << " ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
    }

    for(UInt k=0; k<datainfo.numtype56; k++)
    {
      attrname = elemDataTypesStr[datainfo.types56[k]];
      std::cout << attrname;

      for (UInt set=0; set<datainfo.e56[k]; set++)
      {
        std::cout << " idx " << datainfo.Esetinfo[k][set].idx
                  << " time " << datainfo.Esetinfo[k][set].time
                  << " edata " << datainfo.Esetinfo[k][set].ndata << std::endl;

        long edata=datainfo.Esetinfo[k][set].ndata;
        for (UInt elem=1; elem<=nelems; elem++) {
          capaInterface.GetElemData(elem, datainfo.Esetinfo[k][set].idx, data1, data2);

          for (UInt n=0; n<edata; n++) {
            std::cout << data1[n] << " ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
    }

    return;
  }

  FEType SimInputUnv::UnvType2ElemType( const UInt elemType ) {

    switch (elemType) {
    case 91:  return ET_TRIA3;
    case 94:  return ET_QUAD4;
    case 92:  return ET_TRIA6;
    case 95:  return ET_QUAD8;
    case 111: return ET_TET4;  // tetraeder 1.ord
    case 112: return ET_WEDGE6;  // prism     1.ord
    case 115: return ET_HEXA8;  // hexaeder  1.ord
    case 113: return ET_WEDGE15;  // prism     2.ord
    case 116: return ET_HEXA20;  // hexaeder  2.ord
    }

    // This place should never be reached!
    return ET_UNDEF;
  }

  bool SimInputUnv::ReadData()
  {
  }
  
}
