// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SimInputensight.cc
 *       \brief    <Description>
 *
 *       \date     21.02.2014
 *       \author   ahueppe
 */
//================================================================================================

#include "SimInputEnsight.hh"

#include "Domain/Results/ResultInfo.hh"
#include "Domain/Results/BaseResults.hh"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#pragma clang diagnostic ignored "-Warray-parameter"

#include "vtkEnSight6Reader.h"
#include "vtkEnSight6BinaryReader.h"
#include "vtkEnSightGoldReader.h"
#include "vtkEnSightGoldBinaryReader.h"
#include "vtkEnSightReader.h"

#include "vtkDataArrayCollection.h"
#include "vtkDataArrayCollectionIterator.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkArrayIteratorTemplate.h"
#include "vtkFloatArray.h"

#pragma clang diagnostic pop

#include <boost/numeric/conversion/cast.hpp>

// This is due to the OLAS New operator!!!
#undef New



#include <boost/tokenizer.hpp>

namespace CoupledField {

SimInputEnsight::SimInputEnsight(std::string fileName,
                                     PtrParamNode inputNode,
                                     PtrParamNode infoNode )
 : SimInputVTKBased(fileName, inputNode,  infoNode),
   numAvailSteps_(0){
  capabilities_.insert( MESH );
  capabilities_.insert( MESH_RESULTS );
  neutralTime_ = 0.0;

}

SimInputEnsight::~SimInputEnsight(){

}

void SimInputEnsight::SetTimeValue(Double val){

  vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

  if(val < 0.0)
  {
    reader->SetTimeValue(0);
  } else
  {
    reader->SetTimeValue(val);
  }
}

void SimInputEnsight::GetTimeValues(){
  // Get the time sets from the reader.
  std::set<Double> timeSteps;
  std::set<Double>::const_iterator tsIt, tsEnd;
  vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

  this->SetTimeValue(-1.0);
  reader_->Modified();
  reader_->Update();

  vtkDataArrayCollection* timeSets = reader->GetTimeSets();

  // Iterate through the time sets.
  vtkDataArrayCollectionIterator* daciter = vtkDataArrayCollectionIterator::New();
  daciter->SetCollection(timeSets);
  for(daciter->GoToFirstItem(); !daciter->IsDoneWithTraversal();
      daciter->GoToNextItem())
  {
    // Each time set is stored in one message.
    vtkDataArray* da = daciter->GetDataArray();
    for(int i=0; i < da->GetNumberOfTuples(); ++i){
      Double tmp = da->GetTuple1(i);
      // apparently, the read-in values do not have double precision...
      Double test = this->dround(tmp,8);
      timeSteps.insert(test);
    }
  }
  daciter->Delete();
  timeValues_.assign( timeSteps.begin(),timeSteps.end() );
  numAvailSteps_ = timeValues_.size();
}

void SimInputEnsight::CreateReader()
{
  // Create new EnSight reader
  vtkGenericEnSightReader *reader = vtkEnSightReader::New();
  reader->SetCaseFileName(this->fileName_.c_str());
  if (reader->DetermineEnSightVersion(1) == -1) {
    EXCEPTION("Can not open EnSight case file: " << this->fileName_.c_str());
  }
  reader_ = reader;
  reader->ReadAllVariablesOn();
  this->SetTimeValue(-1.0);
  reader_->Modified();
  reader_->Update();
  this->dim_ = GetDim();

  // now Fill the basic CFS_result definition
  FillResultMap();

  // fill the user supplied result names and figure out how the result is definied
  ParamNodeList variableDef = this->myParam_->Get("variableList")->GetList("variable");
  for (UInt aNode = 0; aNode < variableDef.GetSize(); aNode++) {
    // try to parse the element name to a result name
    // extract the string, fill the names
    PtrParamNode varNode = variableDef[aNode];
    std::string cfsVarName = varNode->Get("CFSVarName")->As<std::string>();
    std::string varName = varNode->Get("EnsightVarName")->As<std::string>();
    bool staticVar = varNode->Get("static")->As<bool>();
    SolutionType cfsSol = SolutionTypeEnum.Parse(cfsVarName);
    this->cfsEnsightResMap_[cfsSol].isComplex = false;
    this->cfsEnsightResMap_[cfsSol].isStatic = staticVar;
    std::string dofName = varNode->Get("dof")->As<std::string>();

    if (dofName == "") {
      reader->SetPointArrayStatus(varName.c_str(), 1);
      this->cfsEnsightResMap_[cfsSol].dofNameMap.Push_back(varName);
      this->cfsEnsightResMap_[cfsSol].isValid_.Push_back(false);
    }
    else {
      if (this->cfsEnsightResMap_[cfsSol].dofNameMap.GetSize() == 0) {
        this->cfsEnsightResMap_[cfsSol].dofNameMap.Resize(this->dim_);
        this->cfsEnsightResMap_[cfsSol].isValid_.Resize(this->dim_);
      }
      if (dofName == "x") {
        reader->SetPointArrayStatus(varName.c_str(), 1);
        this->cfsEnsightResMap_[cfsSol].dofNameMap[0] = varName;
        this->cfsEnsightResMap_[cfsSol].isValid_[0] = false;
      }
      else if (dofName == "y") {
        reader->SetPointArrayStatus(varName.c_str(), 1);
        this->cfsEnsightResMap_[cfsSol].dofNameMap[1] = varName;
        this->cfsEnsightResMap_[cfsSol].isValid_[1] = false;
      }
      else if (dofName == "z") {
        if (this->dim_ == 2) {
          EXCEPTION("You passed a z dof for an Ensight file for a 2D mesh");
        }
        reader->SetPointArrayStatus(varName.c_str(), 1);
        this->cfsEnsightResMap_[cfsSol].dofNameMap[2] = varName;
        this->cfsEnsightResMap_[cfsSol].isValid_[2] = false;
      }
      else {
        EXCEPTION("Unknown dof name: " << dofName << " detected for Ensight variable: " << varName);
      }
    }
  }
  ValidateResultDefinition();
}

void SimInputEnsight::ValidateResultDefinition(){
  //lets check if the supplied information is valid and
  //adapt the resultInfor definition
  vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

  UInt numVTKResults = reader->GetNumberOfVariables ();
  std::map<SolutionType,ResultDef>::iterator resIter = cfsEnsightResMap_.begin();
  while(resIter != cfsEnsightResMap_.end()){

    // in case of vectorial results we need to check if every component is defined
    StdVector<std::string>& aVec = resIter->second.dofNameMap;
    StdVector<bool>& aValid = resIter->second.isValid_;

    unsigned int numDofNames = aVec.GetSize();
    if (numDofNames == 0) {
      resIter++;
      continue;
    }
    else {
      // check if there are the same names occurring in different dofs
      for (UInt iDof1 = 0; iDof1 < numDofNames; ++iDof1) {
        for (UInt iDof2 = iDof1+1; iDof2 < numDofNames; ++iDof2) {
          if (aVec[iDof1] == aVec[iDof2]) {
            EXCEPTION("The same variable name " << aVec[iDof1] << " is used for different dofs of the same vector!");
          }
        }
      }
    }

    for (UInt k = 0; k < numDofNames; k++) {
      bool wasVector = false;
      // for each user supplied string we look through the ensight file
      for (UInt aVar = 0; aVar < numVTKResults; aVar++) {
        int varType = reader->GetVariableType(aVar);
        std::string ensightName = reader->GetDescription(aVar);
        if (aVec[k] == ensightName) {
          aValid[k] = true;
          // determine if element or node result
          switch (varType) {
          case vtkEnSightReader::COMPLEX_VECTOR_PER_NODE:
            resIter->second.isComplex = true;
            continue;
          case vtkEnSightReader::VECTOR_PER_NODE:
            resIter->second.res->definedOn = ResultInfo::NODE;
            wasVector = true;
            continue;
          case vtkEnSightReader::COMPLEX_SCALAR_PER_NODE:
            resIter->second.isComplex = true;
          case vtkEnSightReader::SCALAR_PER_NODE:
            resIter->second.res->definedOn = ResultInfo::NODE;
            if (wasVector)
              EXCEPTION("May not mix SCALAR_PER_ELEMENT and VECTOR_PER_ELEMENT in one result. VarName: " << ensightName);
            continue;
          case vtkEnSightReader::COMPLEX_VECTOR_PER_ELEMENT:
            resIter->second.isComplex = true;
            continue;
          case vtkEnSightReader::VECTOR_PER_ELEMENT:
            resIter->second.res->definedOn = ResultInfo::ELEMENT;
            wasVector = true;
            continue;
          case vtkEnSightReader::COMPLEX_SCALAR_PER_ELEMENT:
            resIter->second.isComplex = true;
            continue;
          case vtkEnSightReader::SCALAR_PER_ELEMENT:
            resIter->second.res->definedOn = ResultInfo::ELEMENT;
            if (wasVector)
              EXCEPTION("May not mix SCALAR_PER_ELEMENT and VECTOR_PER_ELEMENT in one result. VarName: " << ensightName);
            continue;
          default:
            EXCEPTION("VTK variable type not supported. Must be (COMPLEX)_(SCALAR|VECTOR)_PER_(NODE|ELEMENT).");
          }
        }
      }
    }
    if (aValid.Find(false) == -1) {
      availResults_.insert(resIter->first);
    }
    resIter++;
  }
  reader->ReadAllVariablesOff();
  this->SetTimeValue(-1.0);
  reader_->Modified();
  reader_->Update();
}

void SimInputEnsight::GetResultTypes( UInt sequenceStep,
                                           StdVector<shared_ptr<ResultInfo> >& infos,
                                           bool isHistory ){
  std::set<SolutionType>::iterator aRes = availResults_.begin();
  while(aRes != availResults_.end()){
    this->cfsEnsightResMap_[*aRes].res->isStatic = cfsEnsightResMap_[*aRes].isStatic;
    infos.Push_back(this->cfsEnsightResMap_[*aRes].res);
    aRes++;
  }
}

void SimInputEnsight::FillResultMap(SolutionType sType, ResultInfo::EntryType eType) {
  this->cfsEnsightResMap_[sType].res.reset( new ResultInfo);
  shared_ptr<ResultInfo> resInfoPtr = this->cfsEnsightResMap_[sType].res;
  resInfoPtr->resultType = sType;
  resInfoPtr->resultName = SolutionTypeEnum.ToString(sType);
  if (eType == ResultInfo::TENSOR) {
    if(dim_ == 3) {
      resInfoPtr->dofNames = "xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"; 
    } else {
      resInfoPtr->dofNames = "xx", "xy", "yx", "yy";
    }
  } else if (eType == ResultInfo::VECTOR) {
    if(dim_ == 3) {
      resInfoPtr->dofNames = "x", "y", "z"; 
    } else {
      resInfoPtr->dofNames = "x", "y";
    }
  } else {
    resInfoPtr->dofNames = "";
  }
  resInfoPtr->unit = MapSolTypeToUnit(sType);
  resInfoPtr->definedOn = ResultInfo::NODE;
  resInfoPtr->entryType = eType;
}

void SimInputEnsight::FillResultMap(){
  FillResultMap(FLUIDMECH_VELOCITY, ResultInfo::VECTOR);
  FillResultMap(FLUIDMECH_PRESSURE, ResultInfo::SCALAR);
  FillResultMap(FLUIDMECH_PRESSURE_TIME_DERIV_1, ResultInfo::SCALAR);
  FillResultMap(FLUIDMECH_PRESSURE_DERIV_2, ResultInfo::SCALAR);
  FillResultMap(FLUIDMECH_PRESSURE_DERIV_1, ResultInfo::VECTOR);
  FillResultMap(FLUIDMECH_VORTICITY, ResultInfo::VECTOR);
  FillResultMap(MEAN_FLUIDMECH_VELOCITY, ResultInfo::VECTOR);
  FillResultMap(FLUIDMECH_STRESS, ResultInfo::VECTOR);
  FillResultMap(FLUIDMECH_DENSITY, ResultInfo::SCALAR);
  FillResultMap(FLUIDMECH_TEF, ResultInfo::SCALAR);
  FillResultMap(FLUIDMECH_TEMP, ResultInfo::SCALAR);
  FillResultMap(FLUIDMECH_TKE, ResultInfo::SCALAR);
  FillResultMap(FLUIDMECH_TDR, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_1, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_2, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_3, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_4, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_5, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_6, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_7, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_8, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_9, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_10, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_11, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_12, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_13, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_14, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_15, ResultInfo::SCALAR);
  FillResultMap(OPT_RESULT_16, ResultInfo::SCALAR);
}


void SimInputEnsight::GetResultEntities( UInt sequenceStep,
                                      shared_ptr<ResultInfo> info,
                                      StdVector<shared_ptr<EntityList> >& list,
                                      bool isHistory){
  this->SetTimeValue(-1);
  reader_->Modified();
  reader_->Update();

  //here we assume, that all regions contain results in all timesteps
  // this is a little dangerous as the situation can be different if multiple timesets and or filesets
  // are involved. We will come to this later...

  SolutionType solT = info->resultType;

  if(info->definedOn != cfsEnsightResMap_[solT].res->definedOn){
    EXCEPTION("Requested result " << SolutionTypeEnum.ToString(solT) << " on " << ResultInfo::EntityUnknownTypeEnum_.ToString(info->definedOn)
               << " but it is actually given on " << ResultInfo::EntityUnknownTypeEnum_.ToString(cfsEnsightResMap_[solT].res->definedOn) )
  }

  //we only cover volume regions right now.
  //const StdVector<std::string>& aRegions =  this->regionNamesOfDim_[dim_];
  EntityList::ListType listType;
  switch( info->definedOn ) {
  case ResultInfo::NODE:
  //case ResultInfo::PFEM:
    listType = EntityList::NODE_LIST;
    break;
  case ResultInfo::ELEMENT:
    listType = EntityList::ELEM_LIST;
    break;
  default:
    EXCEPTION( "Only results defined on nodes and elements "
               << "can be read in from Ensight file up to now" );
  }

  StdVector<RegionIdType> regions;
  mi_->GetVolRegionIds(regions);

  list.Clear();
  for( UInt i = 0; i < regions.GetSize(); i++ ) {
    list.Push_back( mi_->GetEntityList( listType, mi_->GetRegion().ToString(regions[i]) ) );
  }
}

template<typename T>
void SimInputEnsight::GetNodeResult( UInt sequenceStep,
    UInt stepValue,
    shared_ptr<BaseResult> result){

  vtkPointData* pointData = NULL;
  vtkDataSet* ds = NULL;
  bool scalarForVector;
  StdVector< vtkDataArray * > CurArrays;

  //obtain ensight reader
  //vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

  SolutionType solT = result->GetResultInfo()->resultType;
  UInt numDofs = result->GetResultInfo()->dofNames.GetSize();
  UInt numEntities = result->GetEntityList()->GetSize();
  UInt resVecSize =  numEntities * numDofs;

  Vector<T> & resVec = dynamic_cast<Result<T>& >(*result).GetVector();
  resVec.Resize( resVecSize );
  resVec.Init();

  shared_ptr<EntityList> resList = result->GetEntityList();
  EntityIterator iter = resList->GetIterator();

  //try to obtain region ID from entity list name.
  RegionIdType id = mi_->GetRegion().Parse(resList->GetName());
  UInt blockId = this->regionAssoc_[id];

  ds = vtkDataSet::SafeDownCast(reader_->GetOutput()->GetBlock(blockId));
  pointData = ds->GetPointData();

  UInt numNodeArrays = pointData->GetNumberOfArrays();
  if(numDofs>1 && pointData->IsArrayAnAttribute(vtkDataSetAttributes::VECTORS) == -1){
    scalarForVector = false;
  }else{
    scalarForVector = true;
  }

  //get grip of VTK arrays
  CurArrays.Resize(numNodeArrays);
  StdVector<std::string>::iterator nameIter = cfsEnsightResMap_[solT].dofNameMap.Begin();
  for (UInt array=0; array < numNodeArrays; array++,nameIter++){
    CurArrays[array] = pointData->GetArray(nameIter->c_str());
  }

  StdVector<T> myValues(numDofs);
  while(!iter.IsEnd()){
    //we assume region nodes are sorted continuously
    if(!scalarForVector){
      myValues.Init();
      CurArrays[0]->GetTuple(iter.GetPos(),myValues.GetPointer());
      for(UInt j=0; j<numDofs; j++) {
        resVec[iter.GetPos()*numDofs+j] = myValues[j];
      }
    }else{
      T myVal=0;
      for(UInt j=0; j<numDofs; j++) {
        CurArrays[j]->GetTuple(iter.GetPos(),&myVal);
        resVec[iter.GetPos()*numDofs+j] = myVal;
      }
    }
    iter++;
  }
}

template<typename T>
void SimInputEnsight::GetElemResult( UInt sequenceStep,
    UInt stepValue,
    shared_ptr<BaseResult> result){

  vtkCellData* cellData = NULL;
  vtkDataSet* ds = NULL;
  bool scalarForVector = false;
  StdVector< vtkDataArray * > CurArrays;

  //obtain ensight reader
  //vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

  SolutionType solT = result->GetResultInfo()->resultType;
  UInt numDofs = result->GetResultInfo()->dofNames.GetSize();
  UInt numEntities = result->GetEntityList()->GetSize();
  UInt resVecSize =  numEntities * numDofs;

  Result<T>* myResult = dynamic_cast<Result<T>* >(result.get());
  Vector<T> & resVec = myResult->GetVector();
  resVec.Resize( resVecSize );
  resVec.Init();

  shared_ptr<EntityList> resList = result->GetEntityList();
  ElemList* aEList = dynamic_cast<ElemList*>(resList.get());

  //try to obtain region ID from entity list name.
  RegionIdType id = mi_->GetRegion().Parse(resList->GetName());
  UInt blockId = this->regionAssoc_[id];

  ds = vtkDataSet::SafeDownCast(reader_->GetOutput()->GetBlock(blockId));
  cellData = ds->GetCellData();

  UInt numNodeArrays = cellData->GetNumberOfArrays();
  if(numNodeArrays == 0){
    std::cerr << "No cell arrays are defined for step #" << stepValue << std::endl;
    return;
  }

  // get grip of VTK arrays
  CurArrays.Resize(numNodeArrays);
  StdVector<std::string>::iterator nameIter = cfsEnsightResMap_[solT].dofNameMap.Begin();
  for (UInt array = 0; array < numNodeArrays; array++, nameIter++) {
    CurArrays[array] = cellData->GetArray(nameIter->c_str());
  }

  // check if we have a vector result that is passed by multiple Ensight scalar results
  if (numDofs > 1 && cfsEnsightResMap_[solT].dofNameMap.GetSize() == numDofs) {
    // Ensight scalar results and CFS vector result
    scalarForVector = true;
  }

  //fill index of cell arrays
  StdVector<T> myValues(numDofs);
  for(UInt aElem=0;aElem<numEntities;++aElem){
    const Elem* curE = aEList->GetElem(aElem);
    UInt idx =  aElem;
    UInt vtkIdx = this->globElemLocElem_[curE->elemNum];

    //we assume region nodes are sorted continuously
    if(!scalarForVector){
      myValues.Init();
      CurArrays[0]->GetTuple(vtkIdx,myValues.GetPointer());
      for(UInt j=0; j<numDofs; j++) {
        resVec[idx*numDofs+j] = myValues[j];
      }
    }else{
      T myVal=0;
      for(UInt j=0; j<numDofs; j++) {
        CurArrays[j]->GetTuple(vtkIdx,&myVal);
        resVec[idx*numDofs+j] = myVal;
      }
    }
  }
}

void SimInputEnsight::GetResult( UInt sequenceStep,
                             UInt stepValue,
                             shared_ptr<BaseResult> result,
                             bool isHistory){
  bool elemResult = false;


  vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

  SolutionType solT = result->GetResultInfo()->resultType;
  StdVector<std::string>& actVTKVar = cfsEnsightResMap_[solT].dofNameMap;

  //activate only the desired quantities

  if(this->timeValues_.size() < stepValue){
   WARN("Requested step not available in time value map setting to first step");
   this->SetTimeValue(-1.0);
  }else{
    this->SetTimeValue(this->timeValues_[stepValue]);
  }

  reader->ReadAllVariablesOff();
  //disable all arrays and enable only the requested
  DisableAllArrays();
  if(cfsEnsightResMap_[solT].res->definedOn == ResultInfo::NODE){
   for(UInt i=0;i<actVTKVar.GetSize();i++)
     reader->SetPointArrayStatus(actVTKVar[i].c_str(), 1);
  }else if (cfsEnsightResMap_[solT].res->definedOn == ResultInfo::ELEMENT){
    elemResult = true;
    for(UInt i=0;i<actVTKVar.GetSize();i++)
      reader->SetCellArrayStatus(actVTKVar[i].c_str(), 1);
  }
  reader_->Modified();
  reader_->Update();

  if( result->GetEntryType() == BaseMatrix::DOUBLE ) {
    if(elemResult){
      GetElemResult<Double>(sequenceStep,stepValue,result);
    }else{
      GetNodeResult<Double>(sequenceStep,stepValue,result);
    }
  }else{
    EXCEPTION("Complex valued results are not supported");
  }
}

void SimInputEnsight::DisableAllArrays(){
  vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);
  UInt numArray = 0;
  //first cell arrays
  numArray = reader->GetNumberOfCellArrays();
  for(UInt i=0;i<numArray;i++){
    const char* name = reader->GetCellArrayName(i);
    reader->SetCellArrayStatus(name,0);
  }
  numArray = reader->GetNumberOfPointArrays();
  for(UInt i=0;i<numArray;i++){
    const char* name = reader->GetPointArrayName(i);
    reader->SetPointArrayStatus(name,0);
  }
}

Double SimInputEnsight::dround(Double a, int ndigits) {
  if(a == 0)
    return 0.0;

  int    exp_base10 = round(log10(a));
  double man_base10 = a*std::pow(10.0,-exp_base10);
  double factor     = std::pow(10.0,-ndigits+1);
  double truncated_man_base10 = man_base10 - fmod(man_base10,factor);
  double rounded_remainder    = fmod(man_base10,factor)/factor;

  rounded_remainder = rounded_remainder > 0.5 ? 1.0*factor : 0.0;

  return (truncated_man_base10 + rounded_remainder)*std::pow(10.0,exp_base10) ;
}

// Explicit template instantiation
template void SimInputEnsight::GetNodeResult<Double>(UInt,UInt,shared_ptr<BaseResult>);
template void SimInputEnsight::GetElemResult<Double>(UInt,UInt,shared_ptr<BaseResult>);

}


