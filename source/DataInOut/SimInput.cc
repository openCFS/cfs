// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;



#include "FeBasis/FeSpace.hh"

#include "SimInput.hh"
#include "Domain/Results/BaseResults.hh"

namespace CoupledField 
{
  
  SimInput::GetResultArguments::GetResultArguments(UInt sequenceStep, SolutionType solType, std::string regionName):
    sequenceStep_(sequenceStep),
    solType_(solType),
    regionName_(regionName) {
  }
  
   bool SimInput::GetResultArguments::operator<(const GetResultArguments &r) const {
    return (sequenceStep_ < r.sequenceStep_) || (solType_ < r.solType_) || (regionName_ < r.regionName_);
  }

  shared_ptr<BaseResult> 
  SimInput::GetResult( UInt sequenceStep,
                       UInt stepValue,
                       SolutionType solType,
                       const std::string& regionName )
  {
    GetResultArguments getResultArguments(sequenceStep, solType, regionName);
    
    if (rawBaseResults_.find(getResultArguments) == rawBaseResults_.end()) {
      // get all defined result types
      StdVector<shared_ptr<ResultInfo> > infos;
      bool isHistory = false;
      
      GetResultTypes( sequenceStep, infos, isHistory );
  
      // find correct one; if multiple are present -> Exception
      bool found = false;
      shared_ptr<ResultInfo> actInfo;
      for( UInt i = 0; i < infos.GetSize(); i++ ) {
        if( infos[i]->resultType == solType ) {
          // check, if result was already found
          if(found) { 
            EXCEPTION( "A result of type '" << SolutionTypeEnum.ToString(solType) << "' was already "
                        << "found in sequence Step " << sequenceStep );
          }
          actInfo = infos[i];
          found = true;
        }
      }
      
      // check if any result at all was found
      if( !actInfo ) {
      EXCEPTION( "Result was not found in file '"
                  << fileName_ << "'" );
      }
  
      // get all regions for given resulinfo object
      StdVector<shared_ptr<EntityList> > entList;
      GetResultEntities( sequenceStep, actInfo, entList, false );
  
      // find correct one; if none is found -> Exception
      shared_ptr<EntityList> actList;
      for( UInt i = 0; i < entList.GetSize(); i++ ) {
        if( entList[i]->GetName() == regionName )
          actList = entList[i];
      }
      
      // check if any region at all was  found
      if( !actList) {
        EXCEPTION( "No entitylist found for result '"
                  << solType << "' on region '" << regionName << "'" ); 
      }
        
  
      // determine analysistype of current multi sequence step
      std::map<UInt, BasePDE::AnalysisType> analysis;
      std::map<UInt, UInt> numSteps;
      GetNumMultiSequenceSteps( analysis, numSteps, false );
      
      // create new raw result object and put it to the rawBaseResults_ map
      shared_ptr<BaseResult> rawResult;
      if( analysis[sequenceStep] != BasePDE::HARMONIC &&
          analysis[sequenceStep] != BasePDE::EIGENFREQUENCY ) {
        rawResult = shared_ptr<BaseResult>(new Result<Double>() );
      } else {
        rawResult = shared_ptr<BaseResult>(new Result<Complex>() );
      }
      rawResult->SetResultInfo( actInfo );
      rawResult->SetEntityList( actList );
      rawBaseResults_[getResultArguments] = rawResult;
    }
    
    // clone result object from rawBaseResults_ map, fill it and return it
    shared_ptr<BaseResult> result = rawBaseResults_[getResultArguments]->Clone();
    GetResult( sequenceStep, stepValue, result);
    return result;
  } 
  
  template<typename TYPE>
  shared_ptr<FeFunction<TYPE> >
  SimInput::GetFeFunction( UInt sequenceStep,
                           UInt stepValue,
                           SolutionType solType,
                           std::set<std::string> & regionNames )
  {
    // get grid
    shared_ptr<BaseResult> Bres =   GetResult( sequenceStep,stepValue, solType, *regionNames.begin() );
    Grid * ptGrid = dynamic_cast<Result<TYPE>&>(*Bres).GetEntityList()->GetGrid();
    
    shared_ptr<ResultInfo> actInfo = dynamic_cast<Result<TYPE>&>(*Bres).GetResultInfo();

    // Create FeSpace (Only Nodal Right now)
    //For higher order results and such we need some additional functionality in
    // the fe space here....
    shared_ptr<FeSpace> nodalSpace = 
        FeSpace::CreateInstance(myParam_,myInfo_->Get("ResultHandler"),FeSpace::H1, ptGrid);
    nodalSpace->SetDefaultRegionApproximation();

    //create FeFunction
    shared_ptr<FeFunction<TYPE> > myFunc(new FeFunction<TYPE>(NULL));
    
    myFunc->SetGrid(ptGrid);
    // Take care of cyclic dependency between FeFunction and FeSpace
    myFunc->SetFeSpace(nodalSpace);
    nodalSpace->AddFeFunction(myFunc);
    myFunc->SetResultInfo(actInfo);

    actInfo->SetFeFunction(myFunc);
    nodalSpace->AddFeFunction(myFunc);
    // iterate over all regionNames
    StdVector<shared_ptr<BaseResult > > results;
    results.Resize( regionNames.size() );
    std::set<std::string>::iterator regIter = regionNames.begin();
    UInt pos = 0;
    for( ; regIter != regionNames.end(); ++regIter) {

      // obtain result object
      results[pos++] = GetResult(sequenceStep, stepValue, solType, *regIter );

      // pass it ot eqnMap
      shared_ptr<EntityList> entList = ptGrid->GetEntityList(EntityList::ELEM_LIST,*regIter);

      myFunc->AddEntityList( entList );
    }
    
    // finalize feSpace
    myFunc->SetFctId(PSEUDO_FCT_ID);

    nodalSpace->Finalize();
    myFunc->Finalize();

    //if(progOpts->DoListMapping())
    //  eqnMap->ToInfo(info->Get(ParamNode::HEADER)->Get("mappings", ParamNode::APPEND));

    //now we map the result vector to the coefficient vector
    
    // iterate over all regions
    //Double max = 0.0;
    regIter = regionNames.begin();
    for( UInt i = 0; regIter != regionNames.end(); ++i,++regIter) {

      // get result and entitylist
      shared_ptr<EntityList> regionList = results[i]->GetEntityList();

      // get related nodelist
      shared_ptr<EntityList> nodeList
      = ptGrid->GetEntityList(EntityList::NODE_LIST,
                              regionList->GetName());

      EntityIterator it= nodeList->GetIterator();
      Vector<TYPE> & resVec = dynamic_cast<Result<TYPE>&>(*results[i]).GetVector();
      //get grip of singlevector
      SingleVector* coefVec = myFunc->GetSingleVector();

      UInt pos = 0;
      StdVector<Integer> eqns;
      for( it.Begin(); !it.IsEnd(); it++ ) {

        // fetch equations
        nodalSpace->GetEqns( eqns, it );

        // iterate over all equations
        for( UInt iEqn = 0; iEqn < eqns.GetSize(); iEqn++ ) {
          if(eqns[iEqn] > 0)
            coefVec->SetEntry(eqns[iEqn]-1,resVec[pos++]);
        } 
      }
    }
    return myFunc;
  } 

// Explicit template instantiation
  template shared_ptr<FeFunction<Double> >
  SimInput::GetFeFunction( UInt sequenceStep,
                           UInt stepValue,
                           SolutionType solType,
                           std::set<std::string> & regionNames );

  template shared_ptr<FeFunction<Complex> >
  SimInput::GetFeFunction( UInt sequenceStep,
                           UInt stepValue,
                           SolutionType solType,
                           std::set<std::string> & regionNames );

}
