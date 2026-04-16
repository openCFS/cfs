// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <algorithm>
#include <complex>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "SimInputUnv.hh"
#include "Domain/Results/BaseResults.hh"
#include "Domain/Results/ResultInfo.hh"
#include "unv_if.hh"

extern const char *nodeDataTypesStr[30];
extern const char *elemDataTypesStr[10];

namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simInputUNV, "SimInputUnv")

  SimInputUnv::SimInputUnv( std::string fileName, PtrParamNode inputNode,
                            PtrParamNode infoNode) 
    : SimInput(fileName, inputNode, infoNode ) {
    capabilities_.insert( SimInput::MESH);
    capabilities_.insert( SimInput::MESH_RESULTS);

    std::vector<std::string> axisMapTemp;
    axisMapTemp.push_back("x");
    axisMapTemp.push_back("y");
    axisMapTemp.push_back("z");

    myParam_->GetValue("mapXTo", axisMapTemp[0], ParamNode::PASS);
    myParam_->GetValue("mapYTo", axisMapTemp[1], ParamNode::PASS);
    myParam_->GetValue("mapZTo", axisMapTemp[2], ParamNode::PASS);

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
    dim_ = 3;
  }

  // ======================================================
  // GENERAL MESH INFORMATION
  // ======================================================
  UInt SimInputUnv::GetDim() {
   return dim_;
  }

  UInt SimInputUnv::GetNumNodes(){
    return numNodes_;
  }
    
  UInt SimInputUnv::GetNumElems(const Integer dim){
    return numElems_;
  }
  
  UInt SimInputUnv::GetNumRegions(){
    return regionNames_.GetSize();
  }

  UInt SimInputUnv::GetNumNamedNodes(){
    return 0;
  }

  UInt SimInputUnv::GetNumNamedElems(){
    return 0;
  }
  
  // ======================================================
  // ENTITY NAME ACCESS
  // ======================================================

  void SimInputUnv::GetAllRegionNames( StdVector<std::string> & regionNames ){
    regionNames.Clear();
    std::copy(regionNames_.Begin(), regionNames_.End(), regionNames.Begin());
  }

  void SimInputUnv::GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                   const UInt dim )
  {
    regionNames.Clear();
    if (dim == dim_)
      std::copy(regionNames_.Begin(), regionNames_.End(), regionNames.Begin());
    // TODO: determine surface regions
  }
    
  void SimInputUnv::GetNodeNames( StdVector<std::string> & nodeNames )
  {
    nodeNames.Clear();
  }
  
  void SimInputUnv::GetElemNames( StdVector<std::string> & elemNames )
  {
    elemNames.Clear();
  }

  void
  SimInputUnv::ReadMesh( Grid* mi )
  {
    //!! TODO also extract surface meshes
    mi_ = mi;
    
    LOG_TRACE(simInputUNV) << "reading base mesh from file " << fileName_ << ":";
    
    unsigned long n;
    
    if(capaIf_.ReadUniversalfile(fileName_.c_str()) < 0)
    {
      EXCEPTION("Could not read universal file " << fileName_);
      return;
    }
    
    LOG_TRACE(simInputUNV) << "Finished reading";
    
    analysis_ = capaIf_.GetAnalysisType();
    numNodes_ = capaIf_.GetNumOfNodes();
    numElems_ = capaIf_.GetNumOfElements();
    int nset55 = capaIf_.GetNumOfNodeDataSets();
    int nset56 = capaIf_.GetNumOfElemDataSets();
    dim_ = capaIf_.GetDimension();
    // if dim was not read from unv file, initialize it to 2
    // (we will correct it later, if dim is 3, actually)
    if (dim_ == 0) dim_ = 2;

    capaIf_.GetDataInfo(datainfo_);
#ifndef NDEBUG
    long set;
    for (int i=0; i<datainfo_.numtype55; ++i) {
      std::cout << "*** node data type: "
                << nodeDataTypesStr[datainfo_.types55[i]]
                << "\tdatasets: " << datainfo_.n55[i] << std::endl;
      for (set=0; set<datainfo_.n55[i]; ++set) {
        std::cout << "****** dataset " << (set+1) 
                  << ":\tindex: " << datainfo_.Nsetinfo[i][set].idx
                  << "\ttime/freq: " << datainfo_.Nsetinfo[i][set].time
                  << "\tsize: " << datainfo_.Nsetinfo[i][set].ndata
                  << std::endl; 
      }
    }
    for (int i=0; i<datainfo_.numtype56; ++i) {
      std::cout << "*** element data type: "
                << elemDataTypesStr[datainfo_.types56[i]]
                << "\tdatasets: " << datainfo_.e56[i] << std::endl;
      for (set=0; set<datainfo_.e56[i]; ++set) {
        std::cout << "****** dataset " << (set+1)
                  << ":\tindex: " << datainfo_.Esetinfo[i][set].idx
                  << "\ttime/freq: " << datainfo_.Esetinfo[i][set].time
                  << "\tsize: " << datainfo_.Esetinfo[i][set].ndata
                  << std::endl;
      }
    }
#endif

    LOG_TRACE(simInputUNV) << "Number of nodes            : " << numNodes_;
    LOG_TRACE(simInputUNV) << "Number of elements         : " << numElems_;
    LOG_TRACE(simInputUNV) << "Number of node data sets   : " << nset55;
    LOG_TRACE(simInputUNV) << "Number of element data sets: " << nset56;
    

    // get coordinates
    LOG_TRACE(simInputUNV) << "reading vertex coordinates";

    Vector<Double>  p(3), pTemp(3);
    unsigned long   node;

    mi_->AddNodes(numNodes_);

    for (node=1; node<=numNodes_; ++node) {
      capaIf_.GetPos(node, &pTemp[0]);
      // unused variable Double d;
      // In the 2D case openCFS writes the 3D coordinates in the
      // following form: 0 x y
      // if(dim_ == 2)
      // {
      //  p[0] = p[1];
      //  p[1] = p[2];
      //  p[2] = 0;
      // }

      p[0] = axisMap_[0] == 3 ? 0.0 : pTemp[axisMap_[0]];
      p[1] = axisMap_[1] == 3 ? 0.0 : pTemp[axisMap_[1]];
      p[2] = axisMap_[2] == 3 ? 0.0 : pTemp[axisMap_[2]];
      
      // Check if dim is really sane!
      if(p[2] != 0.0)
      {
        dim_ = 3;
        //break; // why break here?
      }
               
      mi_->SetNodeCoordinate(node, p);
    }

    

    // get elements
    RegionIdType regionId;
    long elemColor, unvElemType, numElemNodes, elemNodes[32], dummy;
    UInt elemNodesUInt[32];
    //std::vector< Elem::FEType >  element_types(nelems);
    //std::vector< UInt > num_vertices_of_elements(nelems);
    //std::vector< UInt > element_partition(nelems);
    std::map<long, RegionIdType> color2Region;
    std::map<long, RegionIdType>::iterator findIt;
    Elem::FEType eType;
    std::stringstream strBuffer;

    mi_->AddElems(numElems_);

    for (n=0; n<numElems_; ++n) {
      capaIf_.GetElemColor(n+1, elemColor);
      capaIf_.GetElemType(n+1, unvElemType);
      capaIf_.GetElemNodes(n+1, numElemNodes, elemNodes);

      eType = UnvType2ElemType(unvElemType);

      if(dim_ == 2)
      {
        switch(eType)
        {
          case Elem::ET_QUAD4:
            dummy = elemNodes[1];
            elemNodes[1] = elemNodes[3];
            elemNodes[3] = dummy;
            break;
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
      
      // make sure that regionIds are used in consecutive order
      findIt = color2Region.find(elemColor);
      if (findIt == color2Region.end()) {
        strBuffer.str("");
        strBuffer.clear();
        strBuffer << "Region" << (elemColor+1);
        regionNames_.Push_back(strBuffer.str());
        regionId = mi_->AddRegion(strBuffer.str(), false);
        color2Region[elemColor] = regionId;
      }
      else {
        regionId = findIt->second;
      }

      mi_->SetElemData(n+1, eType, regionId, elemNodesUInt);

    }

  }

  // =========================================================================
  // GENERAL SOLUTION INFORMATION
  // =========================================================================

  void SimInputUnv::GetNumMultiSequenceSteps(
      std::map<UInt, BasePDE::AnalysisType>& analysis,
      std::map<UInt, UInt>& numSteps,
      bool isHistory )
  {
    analysis.clear();
    numSteps.clear();
    
    if (analysis_ == CapaInterfaceC::NoData)
      return;
    
    try {
      analysis[1] = AnalysisCAPA2CFS(analysis_);
    } catch (Exception &ex) {
      WARN(ex.GetMsg() << " Only mesh will be read from UNV file.");
      return;
    }

    UInt maxNumSteps = 0;
    
    for (int i=0; i<datainfo_.numtype55; ++i) {
      if ( (UInt)datainfo_.n55[i] > maxNumSteps )
        maxNumSteps = (UInt)datainfo_.n55[i];
    }
    for (int i=0; i<datainfo_.numtype56; ++i) {
      if ( (UInt)datainfo_.e56[i] > maxNumSteps )
        maxNumSteps = (UInt)datainfo_.e56[i];
    }
    int nHistSteps = 0;
    capaIf_.GetHistXData(nHistSteps, NULL);
    if (nHistSteps > (int) maxNumSteps)
      maxNumSteps = nHistSteps;
    
    numSteps[1] = maxNumSteps;
  }
  
  void SimInputUnv::GetResultTypes(
      UInt sequenceStep, 
      StdVector<shared_ptr<ResultInfo> >& infos,
      bool isHistory)
  {
    if (sequenceStep != 1) {
      EXCEPTION("Invalid sequence step: " << sequenceStep);
    }
    
    int i, numDofs;
    
    infos.Clear();
    
    for ( i=0; i<datainfo_.numtype55; ++i ) {
      shared_ptr<ResultInfo> ptInfo( new ResultInfo() );
      try {
        ptInfo->resultType = NodeResultCAPA2CFS(datainfo_.types55[i]);
        ptInfo->resultName = SolutionTypeEnum.ToString(ptInfo->resultType);
        ptInfo->dofNames.Clear();
        if (datainfo_.Nsetinfo[i][0].ndata == 1) {
          ptInfo->dofNames.Push_back("-");
        } else {
          char dof='x';
          for (int j=0; j<datainfo_.Nsetinfo[i][0].ndata; ++j)
            ptInfo->dofNames.Push_back(std::string(1, dof++));
        }
        ptInfo->unit = MapSolTypeToUnit(ptInfo->resultType);
        ptInfo->complexFormat = REAL_IMAG;
        ptInfo->definedOn = ResultInfo::NODE;
        ptInfo->entryType = (datainfo_.Nsetinfo[i][0].ndata == 1 ?
                             ResultInfo::SCALAR : ResultInfo::VECTOR);
        infos.Push_back(ptInfo);
      } catch (Exception &ex) {
        continue;
      }
    }
    
    for ( i=0; i<datainfo_.numtype56; ++i ) {
      shared_ptr<ResultInfo> ptInfo( new ResultInfo() );
      try {
        ptInfo->resultType = ElemResultCAPA2CFS(datainfo_.types56[i]);
        ptInfo->resultName = SolutionTypeEnum.ToString(ptInfo->resultType);
        ptInfo->dofNames.Clear();
        if (datainfo_.Esetinfo[i][0].ndata == 1) {
          ptInfo->dofNames.Push_back("-");
        } else {
          char dof='x';
          for (int j=0; j<datainfo_.Esetinfo[i][0].ndata; ++j)
            ptInfo->dofNames.Push_back(std::string(1, dof++));
        }
        ptInfo->unit = MapSolTypeToUnit(ptInfo->resultType);
        ptInfo->complexFormat = REAL_IMAG;
        ptInfo->definedOn = ResultInfo::ELEMENT;
        ptInfo->entryType = (datainfo_.Esetinfo[i][0].ndata == 1 ?
                             ResultInfo::SCALAR : ResultInfo::VECTOR);
        infos.Push_back(ptInfo);
      } catch (Exception &ex) {
        continue;
      }
    }
    
    for ( i=0; i<datainfo_.numtype58; ++i ) {
      shared_ptr<ResultInfo> ptInfo( new ResultInfo() );
      try {
        ptInfo->resultType = NodeResultCAPA2CFS(datainfo_.types58[i]);
        ptInfo->resultName = SolutionTypeEnum.ToString(ptInfo->resultType);
        ptInfo->dofNames.Clear();
        numDofs = capaIf_.GetHistNumData(i);
        if (numDofs == 1) {
          ptInfo->dofNames.Push_back("-");
        } else {
          char dof='x';
          for (int j=0; j<numDofs; ++j)
            ptInfo->dofNames.Push_back(std::string(1, dof++));          
        }
        ptInfo->unit = MapSolTypeToUnit(ptInfo->resultType);
        ptInfo->complexFormat = REAL_IMAG;
        ptInfo->definedOn = ResultInfo::NODE;
        ptInfo->entryType = (numDofs == 1 ? ResultInfo::SCALAR
                             : ResultInfo::VECTOR);
        infos.Push_back(ptInfo);
      } catch (Exception &ex) {
        continue;
      }
    }
  }

  void SimInputUnv::GetStepValues( UInt sequenceStep,
                                   shared_ptr<ResultInfo> info,
                                   std::map<UInt, Double>& steps,
                                   bool isHistory)
  {
    if (sequenceStep != 1) {
      EXCEPTION("Invalid sequence step: " << sequenceStep);
    }
    
    steps.clear();
    
    int capaDataType, i, j;
    if (info->definedOn == ResultInfo::NODE) {
      capaDataType = NodeResultCFS2CAPA(info->resultType);
      for ( i=0; i<datainfo_.numtype55; ++i) {
        if (datainfo_.types55[i] == capaDataType) {
          for ( j=0; j<datainfo_.n55[i]; ++j) {
            steps[j+1] = datainfo_.Nsetinfo[i][j].time;
          }
          break;
        }
      }
      if (steps.size() == 0) { // try history data
        for ( i=0; i<datainfo_.numtype58; ++i ) {
          if (datainfo_.types58[i] == capaDataType) {
            int numSteps = 0;
            capaIf_.GetHistXData(numSteps, NULL);
            Double *stepValues = new Double[numSteps];
            capaIf_.GetHistXData(numSteps, stepValues);
            for ( j=0; j<numSteps; ++j ) {
              steps[j+1] = stepValues[j];
            }
            delete[] stepValues;
            break;
          }
        }
      }
    }
    else if (info->definedOn == ResultInfo::ELEMENT) {
      capaDataType = ElemResultCFS2CAPA(info->resultType);
      for ( i=0; i<datainfo_.numtype56; ++i) {
        if (datainfo_.types56[i] == capaDataType) {
          for ( j=0; j<datainfo_.e56[i]; ++j) {
            steps[j+1] = datainfo_.Esetinfo[i][j].time;
          }
          break;
        }
      }
    }
    else { return; }
  }

  void SimInputUnv::GetResultEntities( UInt sequenceStep,
                                       shared_ptr<ResultInfo> info,
                                       StdVector<shared_ptr<EntityList> >& list,
                                       bool isHistory )
  {
    if (sequenceStep != 1) {
      EXCEPTION("Invalid sequence step: " << sequenceStep);
    }

    list.Clear();
    
    // TODO: determine the regions that a result lives on
    // There is no easy way to tell, if a result lives only on some of the
    // regions. So we just return all regions in any case.
    UInt numRegions = regionNames_.GetSize();
    EntityList::ListType listType;
    
    switch (info->definedOn) {
      case ResultInfo::NODE:
        listType = EntityList::NODE_LIST;
        break;
      case ResultInfo::ELEMENT:
        listType = EntityList::ELEM_LIST;
        break;
      default:
        EXCEPTION("UNV files support only results per node/element.");
    }

    for ( UInt i=0; i<numRegions; ++i) {
      list.Push_back( mi_->GetEntityList(listType, regionNames_[i]));
    }
  }

  void SimInputUnv::GetResult( UInt sequenceStep,
                               UInt stepNum,
                               shared_ptr<BaseResult> result,
                               bool isHistory )
  {
    if (sequenceStep != 1) {
      EXCEPTION("Invalid sequence step: " << sequenceStep);
    }

    shared_ptr<ResultInfo> resInfo = result->GetResultInfo();
    UInt numResEntities = result->GetEntityList()->GetSize();

    int capaDataType, i, setIdx;
    UInt j=0, iDof, resVecSize;
    UInt numDofs = resInfo->dofNames.GetSize();
    Double data1[16], data2[16];
    
    if (resInfo->definedOn == ResultInfo::NODE) {
      capaDataType = NodeResultCFS2CAPA(resInfo->resultType);
      
      for ( i=0; i<datainfo_.numtype55; ++i) {
        if (datainfo_.types55[i] == capaDataType) {
          if ((stepNum == 0) || ((int)stepNum > datainfo_.n55[i])) {
            EXCEPTION("Step " << stepNum << " of " << resInfo->resultName
                      << " does not exist.");
          }

#if 0
          if (result->GetEntityList()->GetSize() != numNodes_) {
            EXCEPTION("Wrong number of nodes for result '"
                      << resInfo->resultName << "': should be " << numNodes_);
          }
#endif
          
          if ((int)numDofs != datainfo_.Nsetinfo[i][stepNum-1].ndata) {
            EXCEPTION("Wrong number of DOFs for result '"
                      << resInfo->resultName << "': should be "
                      << datainfo_.Nsetinfo[i][stepNum-1].ndata);
          }
          
          setIdx = datainfo_.Nsetinfo[i][stepNum-1].idx;
          resVecSize = numDofs * numResEntities;
          
          if (analysis_ == CapaInterfaceC::CwData && 
              result->GetEntryType() == BaseMatrix::COMPLEX ) {

            Vector<Complex> &resVec
                = dynamic_cast< Result<Complex>& >(*result).GetVector();
            resVec.Resize(resVecSize);
              std::cout << "i " << i << " stepNum " << stepNum << std::endl;
            
            EntityIterator eit = result->GetEntityList()->GetIterator();
            for ( j=0; j<numResEntities; ++j ) {
              UInt nodeNumber = eit.GetNode();
//              setIdx = datainfo_.Nsetinfo[i][stepNum-1].idx;
//              int setIdx2 = datainfo_.Nsetinfo[i][(stepNum-1)/2+1].idx;
              capaIf_.GetNodeData(nodeNumber, setIdx, data1, data2);
              for ( iDof=0; iDof<numDofs; ++iDof ) {
                resVec[j*numDofs+iDof] = Complex(data2[iDof], 0);
              }
              eit++;
            }

                
          } else {
          if (analysis_ == CapaInterfaceC::ComplexCwData) {
            if ( result->GetEntryType() != BaseMatrix::COMPLEX ) {
              EXCEPTION("Must provide storage of type Complex for result '"
                        << resInfo->resultName << "'.");
            }
            
            Vector<Complex> &resVec
                = dynamic_cast< Result<Complex>& >(*result).GetVector();
            resVec.Resize(resVecSize);
            
            EntityIterator eit = result->GetEntityList()->GetIterator();
            for ( j=0; j<numResEntities; ++j ) {
              UInt nodeNumber = eit.GetNode();
              capaIf_.GetNodeData(nodeNumber, setIdx, data1, data2);
              for ( iDof=0; iDof<numDofs; ++iDof ) {
                resVec[j*numDofs+iDof] = Complex( data1[iDof], data2[iDof] );
              }
              eit++;
            }
            
          } else {
            if ( result->GetEntryType() != BaseMatrix::DOUBLE ) {
              EXCEPTION("Must provide storage of type Double for result '"
                        << resInfo->resultName << "'.");
            }
            
            Vector<Double> &resVec
                = dynamic_cast< Result<Double>& >(*result).GetVector();
            resVec.Resize(resVecSize);
            
            EntityIterator eit = result->GetEntityList()->GetIterator();
            for ( j=0; j<numResEntities; ++j ) {
              UInt nodeNumber = eit.GetNode();
              capaIf_.GetNodeData(nodeNumber, setIdx, data1, data2);
              for ( iDof=0; iDof<numDofs; ++iDof ) {
                resVec[j*numDofs+iDof] = data1[iDof];
              }
              eit++;
            }
          }
          }
          break;
        }
        
      }
      
      if (datainfo_.numtype55 == 0) {
        for ( i=0; i<datainfo_.numtype58; ++i ) {
          if (datainfo_.types58[i] == capaDataType) {
	    UInt histNumData = capaIf_.GetHistNumData(i);
            if (histNumData != numDofs) {
              EXCEPTION("Wrong number of DOFs for result '"
                        << resInfo->resultName << "': should be "
                        << capaIf_.GetHistNumData(i));
            }
            break;
          }
        }
        if (i == datainfo_.numtype58) {
          EXCEPTION("result '" << resInfo->resultName
                    << "' not found in UNV file.");
        }
        
        shared_ptr<EntityList> nodeList = result->GetEntityList();
        EntityIterator it = nodeList->GetIterator();
        resVecSize = numDofs * nodeList->GetSize();
        
        if (analysis_ == CapaInterfaceC::ComplexCwData) {
          if (result->GetEntryType() != BaseMatrix::COMPLEX) {
            EXCEPTION("Must provide storage of type Complex for result '"
                      << resInfo->resultName << "'.");
          }
          
          Vector<Complex> &resVec
              = dynamic_cast< Result<Complex>& >(*result).GetVector();
          resVec.Resize(resVecSize);
          
          for ( it.Begin(); !it.IsEnd(); it++ ) {
            if (capaIf_.GetNodeHistData(it.GetNode(), capaDataType, stepNum,
                                        data1, data2) == 0) {
              for ( iDof=0; iDof<numDofs; ++iDof ) {
                resVec[it.GetPos()*numDofs+iDof]
                    = Complex( data1[iDof], data2[iDof] );
              }
            } else {
              EXCEPTION("Could not read result '" << resInfo->resultName
                  << "' in step " << stepNum << " at node " << j+1);
            }
          }
        } else {
          if (result->GetEntryType() != BaseMatrix::DOUBLE) {
            EXCEPTION("Must provide storage of type Double for result '"
                      << resInfo->resultName << "'.");
          }
          
          Vector<Double> &resVec
              = dynamic_cast< Result<Double>& >(*result).GetVector();
          resVec.Resize(resVecSize);
          
          for ( it.Begin(); !it.IsEnd(); it++ ) {
            if (capaIf_.GetNodeHistData(it.GetNode(), capaDataType, stepNum,
                                        data1, data2) == 0) {
              for ( iDof=0; iDof<numDofs; ++iDof ) {
                resVec[it.GetPos()*numDofs+iDof] = data1[iDof];
              }
            } else {
              EXCEPTION("Could not read result '" << resInfo->resultName
                  << "' in step " << stepNum << " at node " << j+1);
            }
          }
        }

      }
      
    }
    else if (resInfo->definedOn == ResultInfo::ELEMENT) {
      capaDataType = ElemResultCFS2CAPA(resInfo->resultType);
      for ( i=0; i<datainfo_.numtype56; ++i) {
        if (datainfo_.types56[i] == capaDataType) {
          if ((stepNum == 0) || ((int)stepNum > datainfo_.e56[i])) {
            EXCEPTION("Step " << stepNum << " of " << resInfo->resultName
                      << " does not exist.");
          }

#if 0          
          if (result->GetEntityList()->GetSize() != numElems_) {
            EXCEPTION("Wrong number of elements for result '"
                      << resInfo->resultName << "': should be " << numElems_);
          }
#endif
          
          if ((int)numDofs != datainfo_.Esetinfo[i][stepNum-1].ndata) {
            EXCEPTION("Wrong number of DOFs for result '"
                      << resInfo->resultName << "': should be "
                      << datainfo_.Esetinfo[i][stepNum-1].ndata);
          }
          
          setIdx = datainfo_.Esetinfo[i][stepNum-1].idx;
          resVecSize = numDofs * numResEntities;
          
          if (analysis_ == CapaInterfaceC::ComplexCwData) {
            if (result->GetEntryType() != BaseMatrix::COMPLEX) {
              EXCEPTION("Must provide storage of type Complex for result '"
                        << resInfo->resultName << "'.");
            }
            
            Vector<Complex> &resVec
                = dynamic_cast< Result<Complex>& >(*result).GetVector();
            resVec.Resize(resVecSize);
            
            EntityIterator eit = result->GetEntityList()->GetIterator();
            for ( j=0; j<numResEntities; ++j ) {
              UInt elemNum = eit.GetElem()->elemNum;
              capaIf_.GetElemData(elemNum, setIdx, data1, data2);
              for ( iDof=0; iDof<numDofs; ++iDof ) {
                resVec[j*numDofs+iDof] = Complex( data1[iDof], data2[iDof] );
              }
              eit++;
            }
          } else {
            if (result->GetEntryType() != BaseMatrix::DOUBLE) {
              EXCEPTION("Must provide storage of type Double for result '"
                        << resInfo->resultName << "'.");
            }
            
            Vector<Double> &resVec
                = dynamic_cast< Result<Double>& >(*result).GetVector();
            resVec.Resize(resVecSize);
            
            EntityIterator eit = result->GetEntityList()->GetIterator();
            for ( j=0; j<numResEntities; ++j ) {
              UInt elemNum = eit.GetElem()->elemNum;
              capaIf_.GetElemData(elemNum, setIdx, data1, data2);
              for ( iDof=0; iDof<numDofs; ++iDof ) {
                resVec[j*numDofs+iDof] = data1[iDof];
              }
              eit++;
            }
          }
          break;
        }
      }
    }
    else {
      EXCEPTION("UNV files support only results per node/element.");
    }

    /*std::string attrname;
    double data1[32], data2[32];

    std::cout.precision(7);
    std::cout.width(10);
    
    for(UInt k=0; (int) k<datainfo_.numtype55; k++)
    {
      attrname = nodeDataTypesStr[datainfo_.types55[k]];
      std::cout << attrname;

      for (UInt set=0; (int) set<datainfo_.n55[k]; set++)
      {
        std::cout << " idx " << datainfo_.Nsetinfo[k][set].idx
                  << " time " << datainfo_.Nsetinfo[k][set].time
                  << " ndata " << datainfo_.Nsetinfo[k][set].ndata << std::endl;

        long ndata = datainfo_.Nsetinfo[k][set].ndata;
        for (UInt node=1; (int) node<=nnodes; node++) {
          capaIf_.GetNodeData(node,
                                    datainfo_.Nsetinfo[k][set].idx,
                                    data1,
                                    data2);

          for (UInt n=0; n<ndata; n++) {
            std::cout << data1[n] << " " << data2[n] << " ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
    }

    for(UInt k=0; (int) k<datainfo_.numtype56; k++)
    {
      attrname = elemDataTypesStr[datainfo_.types56[k]];
      std::cout << attrname;

      for (UInt set=0; (int) set<datainfo_.e56[k]; set++)
      {
        std::cout << " idx " << datainfo_.Esetinfo[k][set].idx
                  << " time " << datainfo_.Esetinfo[k][set].time
                  << " edata " << datainfo_.Esetinfo[k][set].ndata << std::endl;

        long edata=datainfo_.Esetinfo[k][set].ndata;
        for (UInt elem=1; (int) elem<=nelems; elem++) {
          capaIf_.GetElemData(elem, datainfo_.Esetinfo[k][set].idx, data1, data2);

          for (UInt n=0; n<edata; n++) {
            std::cout << data1[n] << " ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << std::endl;
    }

    return;*/
  }

  // =========================================================================
  // GENERAL HELPER FUNCTIONS
  // =========================================================================
  
  Elem::FEType SimInputUnv::UnvType2ElemType( const UInt elemType ) {

    switch (elemType) {
      case 11:
      case 21:
      case 22:
        return Elem::ET_LINE2; // line first oder
      case 23:
      case 24: 
        return Elem::ET_LINE3; // line second oder
      case 41:
      case 51:
      case 61:
      case 71:
      case 81:
      case 91:
        return Elem::ET_TRIA3; // triangle first oder
      case 42:
      case 52:
      case 62:
      case 72:
      case 82:
      case 92:
        return Elem::ET_TRIA6; //triangle second order
      case 44:
      case 54:
      case 64:
      case 74:
      case 84:
      case 94:
        return Elem::ET_QUAD4; //quadrilateral first order
      case 45:
      case 55:
      case 65:
      case 75:
      case 85:
      case 95:
        return Elem::ET_QUAD8; //quadrilateral second order
      case 111:
        return Elem::ET_TET4;  // tetrahedron first order
      case 118:
        return Elem::ET_TET10;  // tetrahedron second order
      case 101:
      case 112:
        return Elem::ET_WEDGE6;  // wedge first order
      case 102:
      case 113:
        return Elem::ET_WEDGE15;  // wedge second order
      case 104:
      case 115:
        return Elem::ET_HEXA8;  // brick first order
      case 105:
      case 116:
        return Elem::ET_HEXA20;  // brick second order
    }

    // This place should never be reached!
    return Elem::ET_UNDEF;
  }

  BasePDE::AnalysisType SimInputUnv::AnalysisCAPA2CFS(int capaType) {
    switch (capaType) {
      case CapaInterfaceC::StaticData:
        return BasePDE::STATIC;
      case CapaInterfaceC::ComplexCwData:
      case CapaInterfaceC::CwData:
        return BasePDE::HARMONIC;
      case CapaInterfaceC::ModalData:
        return BasePDE::EIGENFREQUENCY;
      case CapaInterfaceC::TimeData:
        return BasePDE::TRANSIENT;
      default:
        EXCEPTION("Unknown CAPA analysis type: " << capaType);
    }
  }

  SolutionType SimInputUnv::NodeResultCAPA2CFS(int capaType) {
    switch (capaType) {
      case CapaInterfaceC::Displacements:
        return MECH_DISPLACEMENT;
      case CapaInterfaceC::Velocities:
        return MECH_VELOCITY;
      case CapaInterfaceC::Accelerations:
        return MECH_ACCELERATION;
      case CapaInterfaceC::ElectricPot:
      case CapaInterfaceC::ElectricSPot:
        return ELEC_POTENTIAL;
      case CapaInterfaceC::ElectricCharge:
        return ELEC_CHARGE;
      case CapaInterfaceC::VelocityPot:
        return ACOU_POTENTIAL;
      case CapaInterfaceC::VelocityPotDeriv1:
        return ACOU_POTENTIAL_DERIV_1;
      case CapaInterfaceC::VelocityPotDeriv2:
        return ACOU_POTENTIAL_DERIV_2;
      case CapaInterfaceC::MagScalarPot:
      case CapaInterfaceC::MagVectorPot:
        return MAG_POTENTIAL;
      case CapaInterfaceC::Temperature:
        return HEAT_TEMPERATURE;
      default:
        EXCEPTION("Unknown CAPA node result type: " << capaType);
    }
  }

  int SimInputUnv::NodeResultCFS2CAPA(SolutionType cfsType) {
    switch (cfsType) {
      case MECH_DISPLACEMENT:
        return CapaInterfaceC::Displacements;
      case MECH_VELOCITY:
        return CapaInterfaceC::Velocities;
      case MECH_ACCELERATION:
        return CapaInterfaceC::Accelerations;
      case ELEC_POTENTIAL:
        return CapaInterfaceC::ElectricSPot;
      case ELEC_CHARGE:
        return CapaInterfaceC::ElectricCharge;
      case ACOU_POTENTIAL:
        return CapaInterfaceC::VelocityPot;
      case ACOU_POTENTIAL_DERIV_1:
        return CapaInterfaceC::VelocityPotDeriv1;
      case ACOU_POTENTIAL_DERIV_2:
        return CapaInterfaceC::VelocityPotDeriv2;
      case MAG_POTENTIAL:
        return CapaInterfaceC::MagVectorPot;
      case HEAT_TEMPERATURE:
        return CapaInterfaceC::Temperature;
      default:
        EXCEPTION("result type '" << SolutionTypeEnum.ToString(cfsType)
                  << "' has no corresponding CAPA result type.");
    }
  }

  SolutionType SimInputUnv::ElemResultCAPA2CFS(int capaType) {
    switch (capaType) {
      case CapaInterfaceC::Strains:
        return MECH_STRAIN;
      case CapaInterfaceC::Stresses:
        return MECH_STRESS;
      case CapaInterfaceC::ElectricField:
        return ELEC_FIELD_INTENSITY;
      case CapaInterfaceC::MagneticField:
        return MAG_FIELD_INTENSITY;
      case CapaInterfaceC::EddyCurrent:
        return MAG_EDDY_CURRENT_DENSITY;
      default:
        EXCEPTION("Unknown CAPA element result type: " << capaType);
    }
  }

  int SimInputUnv::ElemResultCFS2CAPA(SolutionType cfsType) {
    switch (cfsType) {
      case MECH_STRAIN:
        return CapaInterfaceC::Strains;
      case MECH_STRESS:
        return CapaInterfaceC::Stresses;
      case ELEC_FIELD_INTENSITY:
        return CapaInterfaceC::ElectricField;
      case MAG_FIELD_INTENSITY:
        return CapaInterfaceC::MagneticField;
      case MAG_EDDY_CURRENT_DENSITY:
        return CapaInterfaceC::EddyCurrent;
      default:
        EXCEPTION("result type '" << SolutionTypeEnum.ToString(cfsType)
                  << "' has no corresponding CAPA result type.");
    }
  }
  
} // namespace CoupledField
