// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <string>
#include <algorithm>

#include "SimInputRefElems.hh"

namespace CoupledField {

  SimInputRefElems::SimInputRefElems(std::string fileName, PtrParamNode inputNode,
                                     PtrParamNode infoNode) :
    SimInput(fileName, inputNode, infoNode)
  {
    capabilities_.insert( SimInput::MESH);
  }
  
  
  SimInputRefElems::~SimInputRefElems()
  {
  }
  
  
  void SimInputRefElems::InitModule()
  {
    
    std::cout << "Got file name '" << fileName_ << "'" << std::endl;

    UInt idx = fileName_.find( ".refelem");
    
    //    std::cout << "POS: " << idx << "" << std::endl;

    std::string eTypeStr = fileName_.substr(0, idx);
    
    // std::cout << "SUBSTR: " << eTypeStr << "" << std::endl;

    feType_ = Elem::feType.Parse(eTypeStr);

    // Elem::feType.ToString(elemType);
    std::cout << "FETYPE: " << Elem::feType.ToString(feType_) << "" << std::endl;

    dim_ = Elem::shapes[feType_].dim;
    maxNumNodes_ = Elem::shapes[feType_].numNodes;
    maxNumElems_ = 1;

    // EXCEPTION("HALT");
  }
  
  void SimInputRefElems::ReadMesh(Grid *mi)
  {
    mi_ = mi;

    std::vector<double> coords;
    std::vector<UInt> connect;

    mi_->AddElems(maxNumElems_);
    mi_->AddNodes(maxNumNodes_);

    UInt i, idx;

    coords.resize(maxNumNodes_ * 3);
    connect.resize(maxNumNodes_);
    
    switch(feType_) 
    {
    case Elem::ET_LINE3:
      idx = 2;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
    case Elem::ET_LINE2:
      idx = 0;
      coords[idx*3+0] = -1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 1;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      break;

    case Elem::ET_TRIA6:
      idx = 3;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 4;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 5;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
    case Elem::ET_TRIA3:
      idx = 0;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 1;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 2;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      break;

    case Elem::ET_QUAD9:
      idx = 8;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
    case Elem::ET_QUAD8:
      idx = 4;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = -1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 5;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 6;
      coords[idx*3+0] = 0.0;
      coords[idx*3+1] = 1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 7;
      coords[idx*3+0] = -1.0;
      coords[idx*3+1] = 0.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
    case Elem::ET_QUAD4:
      idx = 0;
      coords[idx*3+0] = -1.0;
      coords[idx*3+1] = -1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 1;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = -1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 2;
      coords[idx*3+0] = 1.0;
      coords[idx*3+1] = 1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      idx = 3;
      coords[idx*3+0] = -1.0;
      coords[idx*3+1] = 1.0;
      coords[idx*3+2] = 0.0;
      connect[idx] = idx+1;
      break;

    case Elem::ET_TET10:
      idx = 4;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 5;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 6;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 7;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;

      idx = 8;
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;

      idx = 9;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;

    case Elem::ET_TET4:
      idx = 0;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 1;
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 2;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

      idx = 3;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
      break;

    case Elem::ET_HEXA27:
      idx = 20;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 21;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 22;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 23;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 24;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 25;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 26;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

    case Elem::ET_HEXA20:
      idx = 8;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 9;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 10;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 11;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                          
      idx = 12;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 13;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 14;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 15;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                          
      idx = 16;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 17;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 18;          
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                          
      idx = 19;          
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

    case Elem::ET_HEXA8:
      idx = 0;
      coords[idx*3+0] = -1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 1;            
      coords[idx*3+0] = 1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 2;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 3;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 4;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 5;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 6;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 7;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
      break;

    case Elem::ET_PYRA14:
      idx = 13;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

    case Elem::ET_PYRA13:
      idx = 5;          
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 6;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 7;           
      coords[idx*3+0] = 0;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 8;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 9;           
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;
                            
      idx = 10;            
      coords[idx*3+0] = -0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;
                             
      idx = 11;            
      coords[idx*3+0] = -0.5;
      coords[idx*3+1] = -0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;
                             
      idx = 12;            
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = -0.5;
      coords[idx*3+2] = 0.5;
      connect[idx] = idx+1;

    case Elem::ET_PYRA5:
      idx = 0;
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 1;            
      coords[idx*3+0] = -1;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 2;           
      coords[idx*3+0] = -1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 3;           
      coords[idx*3+0] = 1;
      coords[idx*3+1] = -1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 4;           
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
      break;

    case Elem::ET_WEDGE18:
      idx = 15;           
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                            
      idx = 16;             
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                           
      idx = 17;             
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

    case Elem::ET_WEDGE15:
      idx = 6;           
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                            
      idx = 7;             
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                           
      idx = 8;             
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;
                           
      idx = 9;             
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                           
      idx = 10;            
      coords[idx*3+0] = 0.5;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                           
      idx = 11;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0.5;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;
                           
      idx = 12;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                           
      idx = 13;            
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;
                           
      idx = 14;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 0;
      connect[idx] = idx+1;

    case Elem::ET_WEDGE6:
      idx = 0;
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 1;            
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 2;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = -1;
      connect[idx] = idx+1;    
                            
      idx = 3;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 4;            
      coords[idx*3+0] = 1;
      coords[idx*3+1] = 0;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
                            
      idx = 5;            
      coords[idx*3+0] = 0;
      coords[idx*3+1] = 1;
      coords[idx*3+2] = 1;
      connect[idx] = idx+1;    
      break;

    default:
      break;
    }

    for(i=1;i<=maxNumNodes_;i++) {
      Vector<Double> coord(3);
      coord[2] = 0;
      for(UInt j=0; j<dim_; j++) 
      {
        coord[j] = coords[(i-1)*3+j];
      }
      
      mi_->SetNodeCoordinate( i, coord );
    }

    StdVector<std::string> names;
    StdVector<Integer> ids;
    names.Resize(1);
    names[0] = Elem::feType.ToString(feType_);
    mi_->AddRegions(names, ids);

    mi_->SetElemData(1, feType_, ids[0], &connect[0]);

  }
  
  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputRefElems::GetDim() {
    return dim_;
  }
  
  UInt SimInputRefElems::GetNumNodes(){
    return maxNumNodes_;
  }
    
  UInt SimInputRefElems::GetNumElems(const Integer dim){    
    return maxNumElems_;
  }
  
  UInt SimInputRefElems::GetNumRegions(){
    return 1;
  }

  UInt SimInputRefElems::GetNumNamedNodes(){
    return 0;
  }

  UInt SimInputRefElems::GetNumNamedElems(){
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputRefElems::GetAllRegionNames( StdVector<std::string> & regionNames ){
    regionNames.Resize(1);
    regionNames[0] = Elem::feType.ToString(feType_);
  }
    
  void SimInputRefElems::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                              const UInt dim ) {
    if(Elem::shapes[feType_].dim == dim)
    {
      regionNames.Resize(1);
      regionNames[0] = Elem::feType.ToString(feType_);
    } else 
    {
      regionNames.Resize(0);
    }
    
  }

  void SimInputRefElems::GetNodeNames( StdVector<std::string> &nodeNames ) {
  }

  void SimInputRefElems::GetElemNames( StdVector<std::string> & elemNames ) {
  }

  // ======================================================
  // ENTITY ACCESS
  // ======================================================
  
  void SimInputRefElems::GetCoordinates( std::vector< double > & nodeCoords ) {
  }

   
  void SimInputRefElems::GetElements( std::vector< std::vector<UInt> > & elems,
                                  std::vector< std::vector<Elem::FEType> > & elemTypes,
                                  std::vector< std::vector<UInt> > & elemNums,                                
                                  std::vector<RegionIdType> & regionIds,
                                  const UInt dim ) {
  }

  void SimInputRefElems::GetNamedNodes(StdVector<StdVector<UInt> > &nodes,
                                   StdVector<std::string> &nodeNames )
  {
  }

  void SimInputRefElems::GetNamedElems(StdVector<StdVector<UInt> > & elems,
                                   StdVector<std::string> & elemNames)
  {  
  }
  
}
