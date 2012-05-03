// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGrid.cc 
 *       \brief    Implementation file for the Grid interpolation CoefFunction
 *
 *       \date     Jan. 23, 2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#include "CoefFunctionGrid.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "FeBasis/FeSpace.hh"
#include "DataInOut/ResultHandler.hh"

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tr1/type_traits.hpp>


namespace CoupledField{

template<class DATA_TYPE>
CoefFunctionGridBase<DATA_TYPE>::CoefFunctionGridBase(PtrParamNode configNode){
  this->ReadXmlNode(configNode);

  dependType_ = GENERAL;


  UInt dofDim = 0;
  if(resultInfo_->entryType == ResultInfo::SCALAR){
    dofDim = 1;
    dimType_ = SCALAR;
  }else if(resultInfo_->entryType == ResultInfo::VECTOR){
    dofDim = resultInfo_->dofNames.GetSize();
    dimType_ = VECTOR;
  }else if(resultInfo_->entryType == ResultInfo::TENSOR){
    dofDim = resultInfo_->dofNames.GetSize();
    dimType_ = TENSOR;
  }
  Grid* curGrid = domain->GetGrid(gridId_);
  UInt dim = curGrid->GetDim();
  //right now we hardcode this here as we have always nodal results
  //one can think about setting this operator externally but we come to this later
  if(dofDim == 1){
    if(dim == 2)
      this->myOperator_.reset(new IdentityOperator<FeH1,2,1,DATA_TYPE>());
    else if(dim == 3)
      this->myOperator_.reset(new IdentityOperator<FeH1,3,1,DATA_TYPE>());
  }else if(dofDim == 2){
    if(dim == 2)
      this->myOperator_.reset(new IdentityOperator<FeH1,2,2,DATA_TYPE>());
    else if(dim == 3)
      this->myOperator_.reset(new IdentityOperator<FeH1,3,2,DATA_TYPE>());
  }else if(dofDim == 3){
    if(dim == 2)
      this->myOperator_.reset(new IdentityOperator<FeH1,2,3,DATA_TYPE>());
    else if(dim == 3)
      this->myOperator_.reset(new IdentityOperator<FeH1,3,3,DATA_TYPE>());
  }else{
    EXCEPTION("Dimension of result not supported. Is it tensorial???")
  }

  // THE FOLLOWING IS NOT USED RIGHT NOW! WILL BE IN CASE OF TRANSIENTS
  // obtain handle from internal variable coefficient function
 //mp_ = domain->GetMathParser();
 //mHandleStep_ = mp_->GetNewHandle(true);

 // // register callback mechanism if expression changes
 // mp_->AddExpChangeCallBack(
 //     boost::bind(&CoefFunctionGridBase<DATA_TYPE>::Recalculate, this ),
 //     mHandleStep_ );

  //bool complex = std::tr1::is_same<DATA_TYPE,Complex>::value;
  //std::string value;
  //if(complex)
  //  value = "f";
  //else
  //  value = "t";
  //
  //mp_->SetExpr(mHandleStep_, value);



}

template<class DATA_TYPE>
CoefFunctionGridBase<DATA_TYPE>::~CoefFunctionGridBase(){

}

template<class DATA_TYPE>
void CoefFunctionGridBase<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                                            const LocPointMapped& lpm ){
  EXCEPTION("TENSORIAL RESULTS NOT SUPPORTED FOR INTERPOLATION FROM EXTERNAL GRIDS");

}

template<class DATA_TYPE>
void CoefFunctionGridBase<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                                            const LocPointMapped& lpm ){
  //1. call to grid to find the element
  //2. compute solution of this element
  //3. apply operator
  assert(dimType_ != TENSOR);
  if(dimType_ == SCALAR)
    CoefMat.Resize(1);
  else if (dimType_ == VECTOR)
    CoefMat.Resize(domain->GetGrid(gridId_)->GetDim());

  CoefMat.Init();
#ifdef USE_INTERPOLATION
  //const Elem* targetElem = lpm.ptEl;
  Vector<Double> globCoord;
  Vector<Double> localCoord;
  Vector<DATA_TYPE> elemSol;

  lpm.shapeMap->Local2Global(globCoord,lpm.lp);
  StdVector< Vector<Double> > locCoords;
  StdVector<const Elem*> elems = domain->GetGrid(gridId_)->GetElemsAtGlobalCoord(globCoord,locCoords);
  const Elem* sourceElem=NULL;
  for(UInt i =0;i<elems.GetSize();i++){
    if(srcRegions_.find(domain->GetGrid(gridId_)->GetRegion().ToString(elems[i]->regionId)) != srcRegions_.end()){
      sourceElem = elems[i];
      localCoord = locCoords[i];
      break;
    }
  }
  if(!sourceElem){
    WARN("Could not find element for coord " << globCoord);
    return;
  }

  shared_ptr<ElemShapeMap> esm = domain->GetGrid(gridId_)->GetElemShapeMap( sourceElem, true );
  LocPoint myLp = localCoord;
  LocPointMapped lpmS;
  lpmS.Set(myLp,esm);

  feFunct_->GetElemSolution( elemSol, sourceElem);
  BaseFE * ptFe = feFunct_->GetFeSpace()->GetFe(sourceElem->elemNum);
  myOperator_->ApplyOp(CoefMat,lpmS,ptFe,elemSol);
#endif
}

template<class DATA_TYPE>
void CoefFunctionGridBase<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                                           const LocPointMapped& lpm ){

  assert(dimType_ == SCALAR );

  //const Elem* targetElem = lpm.ptEl;
  Vector<Double> globCoord;
  Vector<Double> localCoord;
  Vector<DATA_TYPE> elemSol;
  Vector<DATA_TYPE> ptSol(1);
  ptSol.Init();

#ifdef USE_INTERPOLATION
  lpm.shapeMap->Local2Global(globCoord,lpm.lp);
  StdVector< Vector<Double> > locCoords;
  StdVector<const Elem*> elems = domain->GetGrid(gridId_)->GetElemsAtGlobalCoord(globCoord,locCoords);
  //const Elem* sourceElem = domain->GetGrid(gridId_)->GetElemAtGlobalCoord(targetElem,globCoord,localCoord);
  const Elem* sourceElem = NULL;
  for(UInt i =0;i<elems.GetSize();i++){
    if(srcRegions_.find(domain->GetGrid(gridId_)->GetRegion().ToString(elems[i]->regionId)) != srcRegions_.end()){
      sourceElem = elems[i];
      localCoord = locCoords[i];
      break;
    }

  }
  if(!sourceElem){
    return;
  }

  shared_ptr<ElemShapeMap> esm = domain->GetGrid(gridId_)->GetElemShapeMap( sourceElem, true );
  LocPoint myLp = localCoord;
  LocPointMapped lpmS;
  lpmS.Set(myLp,esm);
  feFunct_->GetElemSolution( elemSol, sourceElem);
  BaseFE * ptFe = feFunct_->GetFeSpace()->GetFe(sourceElem->elemNum);
  myOperator_->ApplyOp(ptSol,lpmS,ptFe,elemSol);
#endif
  CoefMat = ptSol[0];
}



template<class DATA_TYPE>
std::string CoefFunctionGridBase<DATA_TYPE>::ToString() const {
  return "ToSting is not implemented";
}

template<class DATA_TYPE>
UInt CoefFunctionGridBase<DATA_TYPE>::GetVecSize() const {
  return domain->GetGrid(gridId_)->GetDim();
}

template<class DATA_TYPE>
void CoefFunctionGridBase<DATA_TYPE>::GetTensorSize( UInt& numRows, UInt& numCols ) const {
  EXCEPTION("Not tensor valued data available yet");
}


template<class DATA_TYPE>
void CoefFunctionGridBase<DATA_TYPE>::Recalculate() {
  //no need to recalculate right now
 // UInt aSeqStep = domain->GetDriver()->GetActSequenceStep();
 // UInt stepnum = GetStepNum();
 // feFunct_ = domain->GetResultHandler()->GetFeFunction<DATA_TYPE>(inputId_,aSeqStep,stepnum,solType_,srcRegions_);
}

template<class T>
struct map_data_compare : public std::binary_function<typename T::value_type, 
                                                      typename T::mapped_type, 
                                                      bool>
{
public:
    bool operator() (typename T::value_type &pair, 
                     typename T::mapped_type i) const
    {
        return pair.second == i;
    }
};


template<class DATA_TYPE>
UInt CoefFunctionGridBase<DATA_TYPE>::GetStepNum(){
  Double aTimeFreq = this->mp_->Eval(mHandleStep_);
  UInt step = 0;
  typedef std::map<UInt,Double> MapType;



  MapType::iterator stepIter = std::find_if( stepValueMap_.begin(), stepValueMap_.end(), std::bind2nd(map_data_compare<MapType>(), aTimeFreq) );

  if(stepIter==stepValueMap_.end()){
    //ok there we have to trigger temporal interpolation
    //for now we just throw an exception
    EXCEPTION("Current timestep not found in result file");
  }else{
    step = stepIter->first;
  }
  return step;
}

template<class DATA_TYPE>
void CoefFunctionGridBase<DATA_TYPE>::ReadXmlNode(PtrParamNode configNode){

  std::string solString;
  configNode->Get("input")->GetValue("inputId", inputId_);
  configNode->Get("input")->GetValue("gridId", gridId_);
  configNode->Get("input")->GetValue("quantity", solString);
  solType_ = SolutionTypeEnum.Parse(solString );

  PtrParamNode regionListNode = configNode->Get("regionList", ParamNode::PASS );

  ParamNodeList regionNodes;

  if( regionListNode)
    regionNodes = regionListNode->GetList("region");

  UInt numRegions = regionNodes.GetSize();
  // iterate over all regions
  for(UInt i = 0; i < numRegions; ++i ) {
    // get data from node
    std::string region;
    regionNodes[i]->GetValue( "srcRegionName", region );
    srcRegions_.insert(region);
  }

  UInt aSeqStep = domain->GetDriver()->GetActSequenceStep();
  //obtain availResults
  StdVector<shared_ptr<ResultInfo> > results;
  domain->GetResultHandler()->GetResultTypes(inputId_,aSeqStep,results,false);
  //now we search for the appropriate result
  for(UInt i = 0;i<results.GetSize();i++){
	if( results[i]->resultType == solType_ ) {
      resultInfo_ = results[i];
	  break;
	}
  }
  if(!resultInfo_.get())
	  EXCEPTION("Could not find result " << solString);

  ResultHandler* resHandler = domain->GetResultHandler();
  resHandler->GetStepValues(inputId_,aSeqStep,resultInfo_,stepValueMap_,false);
  feFunct_ = resHandler->GetFeFunction<DATA_TYPE>(inputId_,aSeqStep,stepValueMap_.begin()->first,solType_,srcRegions_);
}


// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionGridBase<Double>;
  template class CoefFunctionGridBase<Complex>;
#endif

}
