// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
  DEFINE_LOG(simInputUNV, "SimInputUnv")

  SimInputUnv::SimInputUnv( std::string fileName, ParamNode * inputNode ) 
    : SimInput(fileName, inputNode) {
    capabilities_.insert( SimInput::MESH);
    capabilities_.insert( SimInput::MESH_RESULTS);

    std::vector<std::string> axisMapTemp;
    axisMapTemp.push_back("x");
    axisMapTemp.push_back("y");
    axisMapTemp.push_back("z");

    myParam_->Get("mapXTo", axisMapTemp[0], false);
    myParam_->Get("mapYTo", axisMapTemp[1], false);
    myParam_->Get("mapZTo", axisMapTemp[2], false);

    // Prepare mapping of coordinates
    axisMap_.Resize(3);
    for(UInt i=0; i<3; i++)
    {
      const std::string axis = axisMapTemp[i];
      
      if(axis == "x")
      {
        axisMap_[i] = 0;
      } else if(axis == "y")
      {
        axisMap_[i] = 1;
      } else if(axis == "z")
      {
        axisMap_[i] = 2;
      } else if(axis == "zero")
      {
        axisMap_[i] = 3;
      }
    }
    
    // std::cout << "UNV axis map: " << axisMap_[0] << " " << axisMap_[1] << " "  << axisMap_[2] << std::endl;
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
    // Converting of Unv files does not work properly
    return 3;
  }

  void SimInputUnv::GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                              std::map<UInt, UInt>& numSteps,
                                              bool isHistory) {
    // At the moment only the grid may be read.
    analysis.clear();
    numSteps.clear();
  }

  UInt SimInputUnv::GetNumNodes(){
    CoupledField::Warning("SimInputUnv::GetNumNodes() not implemented");
    return 0;
  }
    
  UInt SimInputUnv::GetNumElems(const Integer dim){
    CoupledField::Warning("SimInputUnv::GetNumElems() not implemented");
    return 0;
  }
  
  UInt SimInputUnv::GetNumRegions(){
    CoupledField::Warning("SimInputUnv::GetNumRegions() not implemented");
    return 0;
  }

  UInt SimInputUnv::GetNumNamedNodes(){
    CoupledField::Warning("SimInputUnv::GetNumNamedNodes() not implemented");
    return 0;
  }

  UInt SimInputUnv::GetNumNamedElems(){
    CoupledField::Warning("SimInputUnv::GetNumNamedElems() not implemented");
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputUnv::GetAllRegionNames( StdVector<std::string> & regionNames ){
    std::cerr << "SimInputUnv::GetAllRegionNames() not implemented\n";
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
    
		// COMPWARNING: unused variable long set;
    long n;
    CapaInterfaceC capaInterface;
    GDataInfo      datainfo;
    
    if(capaInterface.ReadUniversalfile(fileName_.c_str()) < 0)
    {
      EXCEPTION("Could not read universal file " << fileName_);
      return;
    }
    
    LOG_TRACE(simInputUNV) << "Finished reading";
    
    // COMPWARNING: unused variable int nnumtimesteps=0;
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

    Point  p, pTemp;
    long   node; // COMPWARNING: unused variable long elem;
    // COMPWARNING: unused variable bool   isDimPlausible = false;
    //    ocs::Vec3d oldNormal;

    mi_->AddNodes(nnodes);

    for (node=1; node<=nnodes; node++) {
      capaInterface.GetPos(node, &pTemp[0]);
      // unused variable Double d;
      // In the 2D case CFS++ writes the 3D coordinates in the
      // following form: 0 x y
      // if(dim == 2)
      // {
      //  p[0] = p[1];
      //  p[1] = p[2];
      //  p[2] = 0;
      // }

      p[0] = axisMap_[0] == 3 ? 0.0 : pTemp[axisMap_[0]];
      p[1] = axisMap_[1] == 3 ? 0.0 : pTemp[axisMap_[1]];
      p[2] = axisMap_[2] == 3 ? 0.0 : pTemp[axisMap_[2]];
      
      // Check if dim is really sane!
      //if(p[2] != 0)
      //{
      //  dim = 3;
      //  break;
      //}
               
      mi_->SetNodeCoordinate(node, p);
    }

    

    // get elements
    long elemColor, unvElemType, numElemNodes, elemNodes[32], dummy;
    UInt elemNodesUInt[32];
    std::vector< Elem::FEType >  element_types(nelems);
    std::vector< UInt > num_vertices_of_elements(nelems);
    std::vector< UInt > element_partition(nelems);
    std::set<UInt> partitions;
    Elem::FEType eType;
    std::stringstream strBuffer;


    mi_->AddElems(nelems);

    // COMPWARNING: unused variable long totalidxs = 0;
    for (n=0; n<nelems; ++n) {
      capaInterface.GetElemColor(n+1, elemColor);
      capaInterface.GetElemType(n+1, unvElemType);
      capaInterface.GetElemNodes(n+1, numElemNodes, elemNodes);

      eType = UnvType2ElemType(unvElemType);

      if(dim == 2)
      {
        switch(eType)
        {
        case Elem::ET_TRIA6:
          dummy = elemNodes[1];
          elemNodes[1] = elemNodes[2];
          elemNodes[2] = elemNodes[4];
          elemNodes[4] = elemNodes[3];
          elemNodes[3] = dummy;
          break;
        case Elem::ET_QUAD8:
          dummy = elemNodes[1];
          elemNodes[1] = elemNodes[2];
          elemNodes[2] = elemNodes[4];
          elemNodes[4] = elemNodes[3];
          elemNodes[3] = elemNodes[6];
          elemNodes[6] = elemNodes[5];
          elemNodes[5] = elemNodes[4];
          elemNodes[4] = dummy;
          break;
				default:
					break;
        }
      }        

      for (UInt i=0; i<32; i++) {
        elemNodesUInt[i] = (UInt) elemNodes[i];
      }

      mi_->SetElemData(n+1, eType, elemColor-1,
                       elemNodesUInt);

      partitions.insert(elemColor);

    }

    StdVector<std::string> regionNames;
    StdVector<RegionIdType> regionIds;
    std::set<UInt>::iterator pit, pend;

    pit = partitions.begin();
    pend = partitions.end();

    for(; pit != pend ; pit++)
    {
      strBuffer.str("");
      strBuffer.clear();
      strBuffer << "Partition" << (*pit);
      regionNames.Push_back(strBuffer.str());
    }
    mi_->AddRegions(regionNames, regionIds);

    std::string attrname;
    double data1[32], data2[32];

    for(UInt k=0; (int) k<datainfo.numtype55; k++)
    {
      attrname = nodeDataTypesStr[datainfo.types55[k]];
      std::cout << attrname;

      for (UInt set=0; (int) set<datainfo.n55[k]; set++)
      {
        std::cout << " idx " << datainfo.Nsetinfo[k][set].idx
                  << " time " << datainfo.Nsetinfo[k][set].time
                  << " ndata " << datainfo.Nsetinfo[k][set].ndata << std::endl;

        long ndata = datainfo.Nsetinfo[k][set].ndata;
        for (UInt node=1; (int) node<=nnodes; node++) {
          capaInterface.GetNodeData(node,
                                    datainfo.Nsetinfo[k][set].idx,
                                    data1,
                                    data2);

          for (Integer n=0; n<ndata; n++) {
            std::cout << data1[n] << " ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
    }

    for(UInt k=0; (int) k<datainfo.numtype56; k++)
    {
      attrname = elemDataTypesStr[datainfo.types56[k]];
      std::cout << attrname;

      for (UInt set=0; (int) set<datainfo.e56[k]; set++)
      {
        std::cout << " idx " << datainfo.Esetinfo[k][set].idx
                  << " time " << datainfo.Esetinfo[k][set].time
                  << " edata " << datainfo.Esetinfo[k][set].ndata << std::endl;

        long edata=datainfo.Esetinfo[k][set].ndata;
        for (UInt elem=1; (int) elem<=nelems; elem++) {
          capaInterface.GetElemData(elem, datainfo.Esetinfo[k][set].idx, data1, data2);

          for (Integer n=0; n<edata; n++) {
            std::cout << data1[n] << " ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
    }

    return;
  }

  Elem::FEType SimInputUnv::UnvType2ElemType( const UInt elemType ) {

    switch (elemType) {
    case 91:  return Elem::ET_TRIA3;
    case 94:  return Elem::ET_QUAD4;
    case 92:  return Elem::ET_TRIA6;
    case 95:  return Elem::ET_QUAD8;
    case 111: return Elem::ET_TET4;  // tetraeder 1.ord
    case 112: return Elem::ET_WEDGE6;  // prism     1.ord
    case 115: return Elem::ET_HEXA8;  // hexaeder  1.ord
    case 113: return Elem::ET_WEDGE15;  // prism     2.ord
    case 116: return Elem::ET_HEXA20;  // hexaeder  2.ord
    }

    // This place should never be reached!
    return Elem::ET_UNDEF;
  }

  bool SimInputUnv::ReadData()
  {
		return false;
  }
  
}
