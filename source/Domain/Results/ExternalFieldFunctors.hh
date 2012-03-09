// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ExternalFieldFunctors.hh
 *       \brief    ENTER BRIEF DESCRIPTION HERE
 *
 *       \date     Feb 11, 2012
 *       \author   ahueppe
 */ 
//================================================================================================

#ifndef EXTERNALFIELDFUNCTORS_HH_
#define EXTERNALFIELDFUNCTORS_HH_

#include "ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionGrid.hh"
#include <boost/tr1/type_traits.hpp>

/*! \class ExternalFieldFunctor
 *    \brief Functor for handling user given fields
 *    \date     Feb 11, 2012
 *    \author   ahueppe
 *
 *    This class can handle multiple coefficient functions and can be passed to integrators
 *    to evaluate the field at a given point. It still needs an FeFunction representing the
 *    given Field. This can be seen as an interface because this way an FeFunction can be
 *    defined by a mixture of given coefficient functions disregarding if they are given by
 *    an analytic expression or defined on an external grid.
 *    It provides also an EvalResult function to compute the field values for postprocessing
 */

template<class B_OP, class DATA_TYPE>
class ExternalFieldFunctor : public FieldInterpolFunctor<B_OP,DATA_TYPE>{
public:
  ExternalFieldFunctor(shared_ptr<BaseFeFunction> feFct,
                        shared_ptr<ResultInfo> inf)
      : FieldInterpolFunctor<B_OP,DATA_TYPE>( feFct, inf){

  }

  virtual ~ExternalFieldFunctor(){

  }

  //! Evaluate field for complete result
  virtual void EvalResult( shared_ptr<BaseResult> res ){
    Grid* curGrid = res->GetEntityList()->GetGrid();

    // call to appropriate coefFunction
    ResultInfo& resInfo = *(res->GetResultInfo() );
    UInt numDofs = resInfo.dofNames.GetSize();

    shared_ptr<EntityList> list = res->GetEntityList();
    Vector<DATA_TYPE> & actSol = dynamic_cast<Result<DATA_TYPE>&>(*res).GetVector();
    actSol.Resize( list->GetSize() * numDofs );

    EntityIterator it = list->GetIterator();
    actSol.Init();

    StdVector<Integer> eqnNums;
    UInt pos = 0;
    std::string regName = res->GetEntityList()->GetName();
    RegionIdType reId = curGrid->GetRegion().Parse(regName);
    if(regionCoefFunctions_.find(reId) == regionCoefFunctions_.end())
      return;

    for ( it.Begin(); !it.IsEnd(); it++ ) {
      Vector<Double> globCoord;
      Vector<Double> locCoord;
      Vector<DATA_TYPE> entResult;
      if(it.GetType()!= EntityList::NODE_LIST){
    	  EXCEPTION("Interpoation for extract result not implemented for the Element case");
      }

      UInt curNodeNum = it.GetNode();
      StdVector<UInt> nodes(1);
      StdVector<RegionIdType> regs(1);
      nodes[0] = curNodeNum;
      regs[0] = reId;

      curGrid->GetNodeCoordinate(globCoord,curNodeNum,true);
      //Obtain intersecting element
      StdVector<Elem*> myElems;
      Elem* myElem = NULL;
      curGrid->GetElemsNextToNodes(myElems,nodes,regs);
      //now we find the first element which belongs to the current region
      for(UInt i=0;i<myElems.GetSize();i++){
        if(myElems[i]->regionId == reId){
          myElem = myElems[i];
          break;
        }
      }

      if(!myElem){
        for(UInt iDim = 0; iDim < numDofs; iDim++ ) {
          actSol[pos++] = 0.0;
        }
        continue;
      }

      shared_ptr<ElemShapeMap> esm = curGrid->GetElemShapeMap( myElem, true );
      esm->Global2Local(locCoord,globCoord);
      LocPoint myLp = locCoord;
      LocPointMapped lpm;
      lpm.Set(myLp,esm);
      //try to dermine the correct coefficient function
      regionCoefFunctions_[reId]->GetVector(entResult,lpm);

      for(UInt iDim = 0; iDim < entResult.GetSize(); iDim++ ) {
        actSol[pos++] = entResult[iDim];
      }
    }
  }

  virtual void GetVector(Vector<DATA_TYPE>& vec, const LocPointMapped& lpm ){
    RegionIdType curRegion = lpm.ptEl->regionId;
    if(regionFlags_.find(curRegion) == regionFlags_.end()){
      //this region seems to be processed already we simply return
      EXCEPTION("Region " << curRegion << " is not contained in functor");
      return;
    }
    regionCoefFunctions_[curRegion]->GetVector(vec,lpm);
  }

  //! Add src regions to a tagret region by passing a parameternode
  virtual void AddRegion(RegionIdType region, PtrParamNode config){

    if(regionFlags_.find(region) != regionFlags_.end()){
      //this region seems to be processed already we simply return
      WARN("Region " << region << " already contained in functor");
      return;
    }

    bool complex = std::tr1::is_same<DATA_TYPE,Complex>::value;

    this->feFct_->GetPDE()->ReadUserFieldValues(this->ptGrid_->GetRegion().ToString(region),
                                                config, this->resultInfo_->dofNames,
                                                this->resultInfo_->entryType,complex,
                                                regionCoefFunctions_[region]);

    // create new entity list
    shared_ptr<ElemList> actSDList( new ElemList(this->ptGrid_ ) );
    actSDList->SetRegion( region );
    regionCoefFunctions_[region]->AddEntities(actSDList);

    if(config->Has("grid")){
      regionFlags_[region] = GRID_BASED;
    }else if(config->Has("values")){
      regionFlags_[region] = ANALYTICAL;
    }else if(config->Has("comp")){
      regionFlags_[region] = ANALYTICAL;
    }else{
      EXCEPTION("Faild to process xml information");
    }

  }

protected:


  typedef enum{
    GRID_BASED,ANALYTICAL
  }SourceFieldType;

  //! Map storing the targetRegion id with the flag if its Analytical or GridBased
  std::map<RegionIdType, SourceFieldType> regionFlags_;

  //! Map storing coefFunction of the analytical regions
  std::map<RegionIdType,shared_ptr<CoefFunction> > regionCoefFunctions_;

};


#endif /* EXTERNALFIELDFUNCTORS_HH_ */
