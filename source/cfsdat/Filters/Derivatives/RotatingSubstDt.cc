// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     RotatingSubstDt.cc
 *       \brief    <Description>
 *
 *       \date     Nov 1, 2015
 *       \author   ahueppe, sschoder
 */
//================================================================================================

#include "RotatingSubstDt.hh"
#include "Utils/AnalyticalFields.hh"
#include "Utils/EqnNumberingSimple.hh"
#include <vector>
#include <algorithm>

namespace CFSDat{

ResultIdList RotatingSubstDt::SetUpstreamResults(){
  ResultIdList generated;
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  for(;aIt!=filterResIds.End();++aIt){
    std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;

    //first a single gradient input
    std::string gradInResult = params_->Get("presGrad")->Get("resultName")->As<std::string>();
    gradId_ = resultManager_->AddResult(gradInResult,this->filterTag_);

    //second a single meanflow input
    hasMeanFlow_ = false;
    if(params_->Has("meanFlow")){
    	hasMeanFlow_ = true;
        std::string meanFlowInResult = params_->Get("meanFlow")->Get("resultName")->As<std::string>();
        meanFlowId_ = resultManager_->AddResult(meanFlowInResult,this->filterTag_);

        //copy the timeline from input result
        resultManager_->SetTimeLine(meanFlowId_,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
        generated.Push_back(meanFlowId_);
    }


    //now a time cached time result
    std::string timeInResult = params_->Get("pressure")->Get("resultName")->As<std::string>();
    //hardcode a fith order stencil
    //formular 2(p(1)-p(-1)) + p(2) - p(-2) / 8dt
    CF::StdVector<Integer> timeLine(4);
    timeLine[0] = -2;
    timeLine[1] = -1;
    timeLine[2] = 1;
    timeLine[3] = 2;
    timeId_ = resultManager_->AddResult(timeInResult,this->filterTag_,timeLine);

    //copy the timeline from input result
    resultManager_->SetTimeLine(timeId_,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
    resultManager_->SetTimeLine(gradId_,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));

    generated.Push_back(timeId_);
    generated.Push_back(gradId_);

  }
  return generated;
}


void RotatingSubstDt::AdaptFilterResults(){

  for(UInt aRes = 0; aRes < filterResIds.GetSize(); aRes++){
    ResultManager::ConstInfoPtr aInfo = resultManager_->GetExtInfo(filterResIds[aRes]);
    std::string aFiltResName = aInfo->resultName;

    ResultManager::ConstInfoPtr gradInfo = resultManager_->GetExtInfo(gradId_);
    ResultManager::ConstInfoPtr timeInfo = resultManager_->GetExtInfo(timeId_);

    //got the upstream result validated?
    if(!gradInfo->isValid){
      EXCEPTION("Problem in filter pipeline detected. Pressure gradient input result \"" <<  gradInfo->resultName << "\" could not be provided.")
    }
    if(!timeInfo->isValid){
      EXCEPTION("Problem in filter pipeline detected. Time derivative input result \"" <<  timeInfo->resultName << "\" could not be provided.")
    }

    CF::StdVector<Double>& curTime = *timeInfo->timeLine.get();
    if(curTime.GetSize() == 0){
      EXCEPTION("No upstream filter provided time information for the result "  << aFiltResName);
    }


    //TODO: Known bug mesh result is not consitently defined
    //check if everything is a mesh result
   // if(!gradInfo->isMeshResult || !timeInfo->isMeshResult){
   //   EXCEPTION("Rotating time derivative filter need input results defined on mesh data!")
   // }

    //check if both inputs are defined on Nodes or ELements
    if(gradInfo->definedOn != timeInfo->definedOn){
      EXCEPTION("Rotating time derivative filter needs both inputs on nodes or both inputs on elements.")
    }

    if(hasMeanFlow_){
    	ResultManager::ConstInfoPtr meanFlowInfo = resultManager_->GetExtInfo(meanFlowId_);
        if(!meanFlowInfo->isValid){
          EXCEPTION("Problem in filter pipeline detected. Mean flow input result \"" <<  gradInfo->resultName << "\" could not be provided.")
        }
        if(gradInfo->definedOn != timeInfo->definedOn ||
        		gradInfo->definedOn != meanFlowInfo->definedOn ||
				meanFlowInfo->definedOn != timeInfo->definedOn ){
          EXCEPTION("Rotating time derivative filter either needs all inputs on nodes or on elements.")
        }
        if(gradInfo->dofNames.GetSize() != meanFlowInfo->dofNames.GetSize()){
        	EXCEPTION("Input of the pressure gradient (Dimension: "<<gradInfo->dofNames.GetSize()<<
        			" and the mean flow velocity (Dimension: "<< meanFlowInfo->dofNames.GetSize() <<
					" have a different dimension.")
        }
    }


    //determine the timestep by subtracting first and second entry in timeline
    //test if we have enough steps
    if(curTime.GetSize()<2)
      dt_ = 1;
    else
      dt_ = std::abs( curTime[1] - curTime[0]);

    //we can almost copy everything from time input (also scalar, etc.)
    resultManager_->CopyResultData(timeId_,filterResIds[aRes]);
    resultManager_->SetValid(filterResIds[aRes]);
  }
}

void RotatingSubstDt::FinishInit(){
  std::cout << "\t ---> Rotating Time Derivative Filter prepares geometric information" << std::endl;
  //the only thing left to do is the identification
  //of rotating entity numbers if requested
  if(params_->Has("rotatingDomain")){
    this->rpm_ = params_->Get("rotatingDomain")->Get("rpm")->As<Double>();
    if(params_->Get("rotatingDomain")->Has("cylinder")){
      if(resultManager_->GetExtInfo(gradId_)->ptGrid->GetDim() == 2){
        ExtractCylinderVelocities<2>(params_->Get("rotatingDomain")->Get("cylinder"));
      }else{
        ExtractCylinderVelocities<3>(params_->Get("rotatingDomain")->Get("cylinder"));
      }
    }
    std::cout << "\t\tFound " << rotEnts_.GetSize() << " cells/nodes in rotating region" <<   std::endl;
  }
  gradDim_ = resultManager_->GetExtInfo(gradId_)->dofNames.GetSize();
}


bool RotatingSubstDt::Run(){
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();
  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1)
      continue;
    //we now set the time for the associated results
    Double aTF = resultManager_->GetStepValue(*aIter);
    std::string resName = resultManager_->GetExtInfo(*aIter)->resultName;

    resultManager_->SetTimeValue(gradId_,aTF);
    resultManager_->SetTimeValue(timeId_,aTF);
    resultManager_->ActivateResult(gradId_);
    resultManager_->ActivateResult(timeId_);
    if(hasMeanFlow_){
    	resultManager_->SetTimeValue(meanFlowId_,aTF);
    	resultManager_->ActivateResult(meanFlowId_);
    }
    resultManager_->DeactivateResult(*aIter);
  }
  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }
  aIter = activeResults.begin();
  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1)
      continue;
    std::string resName = resultManager_->GetExtInfo(*aIter)->resultName;

    //obtain again active data

    CF::StdVector<UInt> eqnNums; //they should be equal...
    Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(*aIter,eqnNums);

    returnVec.Init();

    // computation of the actual derivative 5 point stencil
    Vector<Double>& rM2 = resultManager_->GetResultVector<Double>(timeId_,eqnNums,-2);
    Vector<Double>& rM1 = resultManager_->GetResultVector<Double>(timeId_,eqnNums,-1);
    Vector<Double>& rP1 = resultManager_->GetResultVector<Double>(timeId_,eqnNums,1);
    Vector<Double>& rP2 = resultManager_->GetResultVector<Double>(timeId_,eqnNums,2);

    UInt last = (eqnNums.GetSize() == 0)? returnVec.GetSize() : eqnNums.GetSize();

    // Smoothed noise robust derivatives
    // http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
    // scheme error of O(h^3)
    for(UInt i=0;i<last;++i){
      returnVec[i] = ((2.0*(rP1[i]-rM1[i])+rP2[i]-rM2[i])/(8.0*dt_));
    }

    Vector<Double>& gradient = resultManager_->GetResultVector<Double>(gradId_,eqnNums);
    //now we add the substantial part
    UInt gIdx = 0;

    Vector<Double> meanFlow;
    if(hasMeanFlow_){
      meanFlow.Resize(gradient.GetSize());
      meanFlow = resultManager_->GetResultVector<Double>(meanFlowId_,eqnNums);
    }



    EqnMapSimple& mapping = *resultManager_->GetResultAdapter(*aIter)->mapping.get();
    StdVector<UInt> aEqn(1);
    for(UInt aEnt = 0; aEnt <rotEnts_.GetSize();aEnt++){

      mapping.GetEquation(aEqn,rotEnts_[aEnt],resultManager_->GetExtInfo(*aIter)->definedOn);

      //this may be a little hackisch but knowing how eqnationMapSimple works
      //its a little faster
      //returnVec[aEqn[0]] = rotField_[aEnt][1];
      for(UInt d =0;d<gradDim_;++d){
        gIdx = aEqn[0]*gradDim_+d;
        returnVec[aEqn[0]] += rotField_[aEnt][d]*gradient[gIdx];
        if(hasMeanFlow_){
        	returnVec[aEqn[0]] += meanFlow[gIdx]*gradient[gIdx];
        }
      }
    }
    //now we loop over the entity equations
    resultManager_->ActivateResult(*aIter);
  }
  return true;
}


template<unsigned int D>
void RotatingSubstDt::ExtractCylinderVelocities(CF::PtrParamNode cylNode){
  Double r;

  CF::Vector<Double> p1(D);
  CF::Vector<Double> p2(D);
  p1[0]   = cylNode->Get("p1")->Get("x")->As<Double>();
  p1[1]   = cylNode->Get("p1")->Get("y")->As<Double>();
  p1[2]   = cylNode->Get("p1")->Get("z")->As<Double>();
  p2[0]   = cylNode->Get("p2")->Get("x")->As<Double>();
  p2[1]   = cylNode->Get("p2")->Get("y")->As<Double>();
  p2[2]   = cylNode->Get("p2")->Get("z")->As<Double>();
  r       = cylNode->Get("radius")->Get("value")->As<Double>();

  ResultManager::ConstInfoPtr gradInfo = resultManager_->GetExtInfo(gradId_);
  ResultManager::ConstResPtr gradRes = resultManager_->GetResultAdapter(gradId_);
  EqnMapSimple& mapping = *resultManager_->GetResultAdapter(timeId_,1)->mapping.get();
  StdVector<UInt>& eNums = *gradInfo->entityNumbers.get();
  bool hasExtraction = eNums.GetSize() != 0;

  Vector<Double> entCoord(D);
  Vector<Double> velocity(D);
  StdVector<UInt> entEqn(D);
  //! create counterrotating field
  CylinderVortex<D> vortexComp(p1,p2,r,rpm_*-1.0);

  Grid* pGrid = gradInfo->ptGrid;

  //obtain vector of Centroids or Nodes
  if(gradInfo->definedOn == ExtendedResultInfo::NODE){
    if(hasExtraction){
      for(UInt aNum = 0; aNum < eNums.GetSize(); ++aNum){
        pGrid->GetNodeCoordinate(entCoord,eNums[aNum],true);
        vortexComp.ComputeVortexVelocity(velocity,entCoord);
        if(velocity.NormL2() > 0){
          mapping.GetEquation(entEqn,eNums[aNum],CF::ResultInfo::NODE);
          rotEnts_.Push_back(eNums[aNum]);
          rotField_.Push_back(velocity);
        }
      }
    }else{
      std::set<std::string>::iterator regIter = gradInfo->regNames->begin();
      std::set<std::string>::iterator endIter = gradInfo->regNames->end();
      StdVector<UInt> regNodes;
      for(; regIter != endIter; ++regIter){
        RegionIdType rId = gradInfo->ptGrid->GetRegion().Parse(*regIter);
        gradInfo->ptGrid->GetNodesByRegion(regNodes,rId);
        for(UInt aNum = 0;aNum<regNodes.GetSize();++aNum){
          pGrid->GetNodeCoordinate(entCoord,regNodes[aNum]+1,true);
          vortexComp.ComputeVortexVelocity(velocity,entCoord);
          if(velocity.NormL2() > 0){
            mapping.GetEquation(entEqn,regNodes[aNum],CF::ResultInfo::NODE);
            rotEnts_.Push_back(regNodes[aNum]);
            rotField_.Push_back(velocity);
          }
        }
      }
    }
  }else if(gradInfo->definedOn == ExtendedResultInfo::ELEMENT){
    if(hasExtraction){
      for(UInt aNum = 0; aNum < eNums.GetSize(); ++aNum){
        pGrid->GetElemCentroid(entCoord,eNums[aNum],true);
        vortexComp.ComputeVortexVelocity(velocity,entCoord);
        if(velocity.NormL2() > 0){
          mapping.GetEquation(entEqn,eNums[aNum],CF::ResultInfo::ELEMENT);
          rotEnts_.Push_back(eNums[aNum]);
          rotField_.Push_back(velocity);
        }
      }
    }else{
      std::set<std::string>::iterator regIter = gradInfo->regNames->begin();
      std::set<std::string>::iterator endIter = gradInfo->regNames->end();
      StdVector<Elem*> regElems;

      for(; regIter != endIter; ++regIter){
        std::cout << "\t\tComputing rotational velocities for region " << *regIter << std::endl;
        RegionIdType rId = pGrid->GetRegion().Parse(*regIter);
        regElems.Resize(pGrid->GetNumElems(rId));
        pGrid->GetElems(regElems,rId);
        for(UInt aNum = 0;aNum<regElems.GetSize();++aNum){
          pGrid->GetElemCentroid(entCoord,regElems[aNum]->elemNum,true);
          vortexComp.ComputeVortexVelocity(velocity,entCoord);
          if(velocity.NormL2() > 0){
            mapping.GetEquation(entEqn,regElems[aNum]->elemNum,CF::ResultInfo::ELEMENT);
            rotEnts_.Push_back(regElems[aNum]->elemNum);
            rotField_.Push_back(velocity);
          }
        }
      }
    }
  }else{
    EXCEPTION("Rotating time derivative filter can only handle mesh results on nodes or elements.")
  }
}
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template void RotatingSubstDt::ExtractCylinderVelocities<2>(CF::PtrParamNode);
  template void RotatingSubstDt::ExtractCylinderVelocities<3>(CF::PtrParamNode);
#endif

}
