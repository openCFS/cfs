// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_hdf5.hh>

#include <boost/filesystem/fstream.hpp>

#include "DataInOut/SimOutput.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "Domain/CoordinateSystems/DefaultCoordSystem.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeSpace.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/CurlOperator.hh"


#ifdef USE_HDF5
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#endif

#include "ParamsInit.hh"
#include "HelperFuncs.hh"
#include "WVT.hh"

using namespace CoupledField;

namespace CFSTool {

  // Definition of WVT input file type mappings
  static EnumTuple inputFileTypeTuples[] = {
    EnumTuple(WVT::PRIMARY_MODE,   "primary mode"), 
    EnumTuple(WVT::SECONDARY_MODE, "secondary mode"),
    EnumTuple(WVT::MEAN_FLOW,      "mean flow")
  };
  Enum<WVT::InputFileType> WVT::inputFileType =                       \
       Enum<WVT::InputFileType>("WVT input file types",
           sizeof(inputFileTypeTuples) / sizeof(EnumTuple),
           inputFileTypeTuples);


  WVT::WVT( const PtrParamNode& param,
            const PtrParamNode& info ) :
    param_(param),
    info_(info)
  {
    // get filenames from parameter
    priModeFile_ = param_->Get("file1")->As<std::string>();
    secModeFile_ = param_->Get("file2")->As<std::string>();
    meanFlowFile_ = param_->Get("file3")->As<std::string>();
    outFile_ = param_->Get("file4")->As<std::string>();

    writeOutputFile_ = outFile_ != "";

    bool printGridOnly = false;

    // Linear correction  factor for scaling  the actual mean  flow integrated
    // over the fluid domain to the desired mean flow.
    Double correct = 0.0;
       

    inputs_[PRIMARY_MODE] = GetInputReader( priModeFile_, param_, info_ );
    inputs_[SECONDARY_MODE] = GetInputReader( secModeFile_, param_, info_ );
    inputs_[MEAN_FLOW] = GetInputReader( meanFlowFile_, param_, info_ );

    shared_ptr<MathParser> mp(new MathParser());       

    // check capabilities of input class
    InputsType::iterator iit, iend;
    iit = inputs_.begin();
    iend = inputs_.end();    
    for( ; iit != iend; iit++) 
    {
      std::string readerCapsResult;
      if(!CheckReaderCapabilities(inputFileType.ToString(iit->first),
                                  iit->second,
                                  readerCapsResult))
      {
        EXCEPTION(readerCapsResult);
      }
    }

       // read in mesh of input1
       inputs_[PRIMARY_MODE]->InitModule();
       UInt dim = inputs_[PRIMARY_MODE]->GetDim();
       Grid * ptGrid1 = new GridCFS(dim, param_, info_);
       inputs_[PRIMARY_MODE]->ReadMesh(ptGrid1);
       ptGrid1->FinishInit();

       // std::cout << "##############" << xmlFile << std::endl
       //    << "##############" << matFile << std::endl;
       PtrParamNode pNode = GetParamNodeFromHDF5(inputs_[PRIMARY_MODE], "ParameterFile");

       ParamNodeList list;
       
       PtrParamNode multiSeqStep1 = pNode->GetByVal("sequenceStep", "index", 1);

       // list = multiSeqStep1->GetChildren();
       // for(UInt i=0; i<list.GetSize(); i++) {
       //  std::cout << "############## " << list[i]->GetName() << std::endl;
       // }

       // PtrParamNode multiSeqStep2 = pNode->GetByVal("sequenceStep", "index", 2);

       bool dirCoupled = false;
       
       bool ms1HasMech = multiSeqStep1->Get("pdeList")->Has("mechanic");
       bool ms1HasFlow = multiSeqStep1->Get("pdeList")->Has("fluidMech");

       PtrParamNode mechNode, fluidNode;

       StdVector<std::string> cplSurfList;
       
       if( ms1HasMech && ms1HasFlow && multiSeqStep1->Has("couplingList"))
       {
         // For directly coupled WVT, mechanic and fluid analyses need to be
         // in the first MS step. There needs to be also a coupling section.
         // Normally you do not need to apply WVT to directly coupled system
         // since the phase difference between the sensors can be readily
         // obtained by comparing the phases of the mech. displacements.
         // But we want to cross-check our directly and iteratively coupled
         // results.

         PtrParamNode cplNode = multiSeqStep1->Get("couplingList");
         
         if( cplNode->Get("direct")->Has("fluidMechDirect") ) 
         {  
           dirCoupled = true;
           std::cout << "We have a <couplingList> tag." << std::endl;

           mechNode = multiSeqStep1->Get("pdeList")->Get("mechanic");
           fluidNode = multiSeqStep1->Get("pdeList")->Get("fluidMech");         

           list = cplNode->Get("direct")->Get("fluidMechDirect")->Get("surfRegionList")->GetChildren();
           for(UInt i=0; i<list.GetSize(); i++) {
             std::cout << "############## " << list[i]->GetName() << ": " << list[i]->Get("name")->As<std::string>() << std::endl;
             cplSurfList.Push_back(list[i]->Get("name")->As<std::string>());
           }
         }
         else
         {
           EXCEPTION("There exist mechanic and fluidMech PDEs in the " <<
                     "first multi-sequence step,\nbut no coupling node." <<
                     "This case is not yet handled by WVT.");
         }
       }
       else
       {
         // For iteratively coupled WVT, we want to have the mechanic analyis
         // in the first MS step and the fluid analysis in the second MS step 
         // in the future. Since this is not implemented at the moment, we
         // require fluid in the first MS and get the mechanic velocities at
         // sensor positions in lateral mode from somewhere else.

         if(ms1HasMech) 
         {
           EXCEPTION("There exists a mechanic PDE in the first " << 
                     "multi-sequencestep.\n" <<
                     "This case is not yet handled by WVT.");
         }

         if(ms1HasFlow) 
         {
           fluidNode = multiSeqStep1->Get("pdeList")->Get("fluidMech");

           list = fluidNode->Get("bcsAndLoads")->GetList("velocity");
           for(UInt i=0; i<list.GetSize(); i++) {
             std::cout << "############## " << list[i]->GetName() << ": " << list[i]->Get("name")->As<std::string>() << std::endl;
             cplSurfList.Push_back(list[i]->Get("name")->As<std::string>());
           }

         }
       }
       
       std::cout << dirCoupled << std::endl;
       

       // list = fluidNode->GetChildren();
       // for(UInt i=0; i<list.GetSize(); i++) {
       //   std::cout << "############## " << list[i]->GetName() << std::endl;
       // }

       std::map<std::string, std::string> regionMatMap;
       std::map<std::string, std::string>::const_iterator regMatIt, regMatEnd;

       list = fluidNode->Get("regionList")->GetChildren();
       ParamNodeList domainRegionList;
       domainRegionList = pNode->Get("domain")->Get("regionList")->GetChildren();

       for(UInt i=0; i<list.GetSize(); i++) {
         std::string regionName = list[i]->Get("name")->As<std::string>();
         
         std::cout << "############## " << list[i]->GetName() << std::endl;
         std::cout << "-------------> " << regionName << std::endl;
         for(UInt j=0; j<domainRegionList.GetSize(); j++) {
           std::string domRegionName = domainRegionList[j]->Get("name")->As<std::string>();

           if(regionName == domRegionName) 
           {
             std::string domMatName = domainRegionList[j]->Get("material")->As<std::string>();
             // std::cout << "=============> " << domMatName << std::endl;
             regionMatMap[regionName] = domMatName;
           }
         }
       }

       pNode = GetParamNodeFromHDF5(inputs_[PRIMARY_MODE], "MaterialFile");

       std::map<std::string, Double> regionDensityMap;
       std::map<std::string, Double> regionViscosityMap;

       regMatIt = regionMatMap.begin();
       regMatEnd = regionMatMap.end();
       for( ; regMatIt != regMatEnd; regMatIt++ ){
         std::string matName = regMatIt->second;
         std::string regionName = regMatIt->first;
         PtrParamNode mNode = pNode->GetByVal("material",
                                              "name",
                                              matName);

         PtrParamNode flowNode = mNode->Get("flow");

         Double density = flowNode->Get("density")->As<Double>();
         Double viscosity = 0.0;
         if(flowNode->Has("dynamicViscosity")) 
         {
           viscosity = flowNode->Get("dynamicViscosity")->As<Double>();
         }
         else
         {
           viscosity = flowNode->Get("kinematicViscosity")->As<Double>();
           viscosity *= density;
         }
         
         std::cout << "density " << density << std::endl;         
         std::cout << "viscosity " << viscosity << std::endl;         
         regionDensityMap[regionName] = density;
         regionViscosityMap[regionName] = viscosity;
       }

       // read in mesh of input2
       inputs_[SECONDARY_MODE]->InitModule();
       Grid * ptGrid2 = new GridCFS(dim, param_, info_);
       inputs_[SECONDARY_MODE]->ReadMesh(ptGrid2);
       ptGrid2->FinishInit();

       // read in mesh of input3
       inputs_[MEAN_FLOW]->InitModule();
       Grid * ptGrid3 = new GridCFS(dim, param_, info_);
       inputs_[MEAN_FLOW]->ReadMesh(ptGrid3);
       ptGrid3->FinishInit();

       // obtain output writer
       shared_ptr<SimOutput> output;
       if( writeOutputFile_ ) {
         output = GetOutputWriter( outFile_, param_, info_ );
         output->Init( ptGrid1, printGridOnly);
       }

       // obtain number of Sequence Steps and get analysis types
       std::map<UInt, BasePDE::AnalysisType> types;
       std::map<UInt, UInt> numSteps;
       inputs_[PRIMARY_MODE]->GetNumMultiSequenceSteps( types, numSteps, false );

       std::cout << "\nFound " << types.size() << " sequence step(s) in '" << inputs_[PRIMARY_MODE]->GetFileName() << "'\n";
       std::map<UInt, BasePDE::AnalysisType> types2;
       std::map<UInt, UInt> numSteps2;
       inputs_[SECONDARY_MODE]->GetNumMultiSequenceSteps( types2, numSteps2, false );
       std::cout << "\nFound " << types2.size() << " sequence step(s) in '" << inputs_[SECONDARY_MODE]->GetFileName() << "'\n";
       
       if(types.size() != types2.size()){
         std::cout << "'" << inputs_[PRIMARY_MODE]->GetFileName() << "' and '" << inputs_[SECONDARY_MODE]->GetFileName()
            << "' have different number of sequence steps!\n";
         exit(EXIT_FAILURE);
       }

       if(types.size() != 1) {
         std::cout << "WVT input files should have only one sequence step!\n";
         exit(EXIT_FAILURE);
       }

       std::cout << "types.size() " << types.size() << " types2.size() " << types2.size() << "\n";
       std::cout << "types[1] " << types[1] << " types2[1] " << types2[1] << "\n";
       
       if(  types[1] != BasePDE::HARMONIC ||
           types2[1] != BasePDE::HARMONIC ) {
         std::cout << "WVT is only available for harmonic analyses at the moment!\n";
         exit(EXIT_FAILURE);
       }

       UInt actMsStep = 1;
       std::cout << " Computing WVT for sequence step " << actMsStep << std::endl
                 << "-------------------------------------\n\n";
       
       // get resulttypes
       StdVector<shared_ptr<ResultInfo> > infos, infos2, infos3, infos_mean_flow;
       inputs_[PRIMARY_MODE]->GetResultTypes( actMsStep, infos, false );
       inputs_[SECONDARY_MODE]->GetResultTypes( actMsStep, infos2, false );
       inputs_[MEAN_FLOW]->GetResultTypes( actMsStep, infos3, false );
       
       StdVector<shared_ptr<BaseResult> > inResults1, inResults2, inResults_mean_flow,
         outResults, outResults2, outResults3, outResults4, outResults5;
       // stepnumbers, for which at least one result is defined
       std::map<UInt, Double> stepVals, stepVals2, stepVals_mean_flow;
       // contains the stepnumbers/-values in which the particular result is
       // defined in
       std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
       std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps2;
       std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps_mean_flow;
       
       std::map<std::string, bool> surfMap;

       // First determine number of results in current multi sequence step in mean flow file
       // since it can be different than in the lateral and coriolis mode files.
       
       for( UInt iRes = 0; iRes < infos3.GetSize(); iRes++) {
         shared_ptr<ResultInfo> actRes3 = infos3[iRes];
         
         if(actRes3->resultType != FLUIDMECH_VELOCITY) 
         {
           continue;
         }
         
         // get stepvalues of mean flow file
         inputs_[MEAN_FLOW]->GetStepValues( actMsStep, actRes3,
                                resultSteps_mean_flow[actRes3], false);
         stepVals_mean_flow.insert( resultSteps_mean_flow[actRes3].begin(),
                                    resultSteps_mean_flow[actRes3].end() );
         infos_mean_flow.Push_back(actRes3);
       }
       
       
       if( infos.GetSize() > 0 ){
         std::cout << "Computing WVT from the following results:\n";
       }
       
       // iterate over all result types of input1
       for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {
         
         if(infos[iRes]->resultType != FLUIDMECH_VELOCITY)
           continue;
         
         std::cout << "\t" << infos[iRes]->resultName << "\n";
         
         // get stepvalues of reference file
         shared_ptr<ResultInfo> actRes = infos[iRes];
         shared_ptr<ResultInfo> actRes2 = infos2[iRes];
         inputs_[PRIMARY_MODE]->GetStepValues( actMsStep, actRes,
                                resultSteps[actRes], false);
         stepVals.insert( resultSteps[actRes].begin(),
                          resultSteps[actRes].end() );
         
         // get stepvalues of second file
         inputs_[SECONDARY_MODE]->GetStepValues( actMsStep, actRes2,
                                resultSteps2[actRes2], false);
         stepVals2.insert( resultSteps2[actRes2].begin(),
                           resultSteps2[actRes2].end() );
         
         // Loop over all step values in both sets and compare them. Thus we can see
         // differences e.g. in eigenfrequency analysis
         std::map<UInt, Double>::const_iterator svIt1, svIt2;
         svIt1 = stepVals.begin();
         svIt2 = stepVals2.begin();
         for( ; svIt1 != stepVals.end(); ++svIt1, ++svIt2 ) {
           if( svIt1->first != svIt2->first ) {
             EXCEPTION( "Encountered different result steps for result " << 
                        infos[iRes]->resultName );
           } else {
             
             Double val1 = svIt1->second;
             Double val2 = svIt2->second;
             Double relDiff = std::abs(std::abs(val1-val2))/std::abs(val1);
             if ( relDiff > 1e-4 ) {
               EXCEPTION("Time / Frequency values of step " << svIt1->first << " differ by " 
                         << relDiff*100.0 << " %:\n"
                         << "\treference value: " << val1 <<" Hz\n\tcompared value:  " << val2 
                         << " Hz" << std::endl );
             }
           }
         }


         // iterate over all regions
         StdVector< shared_ptr<EntityList> > regions, resRegions;
         inputs_[PRIMARY_MODE]->GetResultEntities( actMsStep, infos[iRes],
                                        resRegions, false );

         std::map<std::string, Double>::iterator it, end;
         it = regionDensityMap.begin();
         end = regionDensityMap.end();
         for( ; it != end; it++ ) {
           std::cout << "region: " << it->first << " density: " << it->second << std::endl;
           std::string regionName = it->first;

           for( UInt iRegion = 0, nReg = resRegions.GetSize(); iRegion < nReg; iRegion++ ) {
             std::string resRegionName = resRegions[iRegion]->GetName();

             if(regionName == resRegionName) 
             {
               regions.Push_back(resRegions[iRegion]);
               surfMap[regionName] = false;
             }
           }           
         }

         for( UInt iRegion = 0, nReg = cplSurfList.GetSize(); iRegion < nReg; iRegion++ ) {
           std::string cplRegionName = cplSurfList[iRegion];
           
           Enum<RegionIdType> regionEnum = ptGrid1->GetRegion();
           RegionIdType regionId = regionEnum.Parse(cplRegionName);

           ElemList* cplElemEntities = new ElemList(ptGrid1);
           cplElemEntities->SetRegion(regionId);
           shared_ptr<EntityList> cplElems(cplElemEntities);             

           regions.Push_back(cplElems);
           surfMap[cplRegionName] = true;
         }           

         for( UInt iRegion = 0, nReg = regions.GetSize(); iRegion < nReg; iRegion++ ) {
           // generate new result object and add it to output writer
           shared_ptr<BaseResult > inResult1, inResult2, inResult3, outResult,
             outResult2, outResult3, outResult4, outResult5;
           
           std::string regionName = regions[iRegion]->GetName();

           //           if(infos[iRes]->resultType != FLUIDMECH_VELOCITY) 
           //           {
           //             continue;
           //           }

           //           if(regionDensityMap.find( regionName ) == regionDensityMap.end() &&
           //              std::find(cplSurfList.Begin(), cplSurfList.End(), regionName) == cplSurfList.End() )
           //           {
           //             continue;
           //           }

           std::cout << "regions[iRegion]->GetName() " << regionName << std::endl;
           //           EXCEPTION("HALT");
           
           inResult1 = shared_ptr<BaseResult>( new Result<Complex>() );
           inResult2 = shared_ptr<BaseResult>( new Result<Complex>() );
           inResult3 = shared_ptr<BaseResult>( new Result<Double>() );
           outResult = shared_ptr<BaseResult>( new Result<Complex>() );
           outResult2 = shared_ptr<BaseResult>( new Result<Complex>() );
           outResult3 = shared_ptr<BaseResult>( new Result<Complex>() );
           outResult4 = shared_ptr<BaseResult>( new Result<Complex>() );
           outResult5 = shared_ptr<BaseResult>( new Result<Complex>() );
           
           inResult1->SetEntityList( regions[iRegion] );
           inResult2->SetEntityList( regions[iRegion] );
           inResult3->SetEntityList( regions[iRegion] );
           
           // Since derivatives need to be taken the output entities need to be elements
           Enum<RegionIdType> regionEnum = ptGrid1->GetRegion();
           RegionIdType regionId = regionEnum.Parse(regions[iRegion]->GetName());
           
           ElemList* outElemEntities = new ElemList(ptGrid1);
           outElemEntities->SetRegion(regionId);
           shared_ptr<EntityList> outElems(outElemEntities);             
           outResult->SetEntityList( outElems );
           outResult2->SetEntityList( outElems );
           outResult3->SetEntityList( outElems );
           outResult4->SetEntityList( outElems );
           outResult5->SetEntityList( outElems );
           
           inResult1->SetResultInfo( infos[iRes] );
           inResult2->SetResultInfo( infos[iRes] );
           inResult3->SetResultInfo( infos_mean_flow[0] );
           
           shared_ptr<ResultInfo> outInfo(new ResultInfo());
           outInfo->resultType = FLUIDMECH_WEIGHT_VECTOR;
           outInfo->resultName = SolutionTypeEnum.ToString(outInfo->resultType);
           outInfo->dofNames = infos[iRes]->dofNames;
           outInfo->unit = MapSolTypeToUnit(outInfo->resultType);
           outInfo->entryType = ResultInfo::VECTOR;
           outInfo->definedOn = ResultInfo::ELEMENT;
           
           outResult->SetResultInfo( outInfo );

           shared_ptr<ResultInfo> outInfo2(new ResultInfo());
           outInfo2->resultType = MEAN_FLUIDMECH_VELOCITY;
           outInfo2->resultName = SolutionTypeEnum.ToString(outInfo2->resultType);
           outInfo2->dofNames = infos[iRes]->dofNames;
           outInfo2->unit = MapSolTypeToUnit(outInfo2->resultType);
           outInfo2->entryType = ResultInfo::VECTOR;
           outInfo2->definedOn = ResultInfo::ELEMENT;

           outResult2->SetResultInfo( outInfo2 );

           shared_ptr<ResultInfo> outInfo3(new ResultInfo());
           outInfo3->resultType = FLUIDMECH_WEIGHT_DENSITY;
           outInfo3->resultName = SolutionTypeEnum.ToString(outInfo3->resultType);
           outInfo3->dofNames = "";
           outInfo3->unit = MapSolTypeToUnit(outInfo3->resultType);
           outInfo3->entryType = ResultInfo::SCALAR;
           outInfo3->definedOn = ResultInfo::ELEMENT;

           outResult3->SetResultInfo( outInfo3 );

           shared_ptr<ResultInfo> outInfo4(new ResultInfo());
           outInfo4->resultType = FLUIDMECH_WEIGHT_VECTOR_PHI;
           outInfo4->resultName = SolutionTypeEnum.ToString(outInfo4->resultType);
           outInfo4->dofNames = infos[iRes]->dofNames;
           outInfo4->unit = MapSolTypeToUnit(outInfo4->resultType);
           outInfo4->entryType = ResultInfo::VECTOR;
           outInfo4->definedOn = ResultInfo::ELEMENT;
           
           outResult4->SetResultInfo( outInfo4 );

           shared_ptr<ResultInfo> outInfo5(new ResultInfo());
           outInfo5->resultType = FLUIDMECH_WEIGHT_DENSITY_PHI;
           outInfo5->resultName = SolutionTypeEnum.ToString(outInfo5->resultType);
           outInfo5->dofNames = "";
           outInfo5->unit = MapSolTypeToUnit(outInfo5->resultType);
           outInfo5->entryType = ResultInfo::SCALAR;
           outInfo5->definedOn = ResultInfo::ELEMENT;

           outResult5->SetResultInfo( outInfo5 );
           
           inResults1.Push_back( inResult1 );
           inResults2.Push_back( inResult2 );
           inResults_mean_flow.Push_back( inResult3 );
           outResults.Push_back( outResult );
           outResults2.Push_back( outResult2 );
           outResults3.Push_back( outResult3 );
           outResults4.Push_back( outResult4 );
           outResults5.Push_back( outResult5 );
           if( output) {
             // Hardcoded: set output format to AMPL_PHASE
             //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
             outResult->GetResultInfo()->complexFormat = REAL_IMAG;
             
             // CAUTION: begin, inc, end are hardcoded and noch checked for each result
             output->RegisterResult( outResult, 1, 1,
                                     resultSteps[actRes].size(),
                                     false );

             // Hardcoded: set output format to AMPL_PHASE
             //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
             outResult2->GetResultInfo()->complexFormat = REAL_IMAG;
             
             // CAUTION: begin, inc, end are hardcoded and noch checked for each result
             output->RegisterResult( outResult2, 1, 1,
                                     resultSteps[actRes].size(),
                                     false );

             // Hardcoded: set output format to AMPL_PHASE
             //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
             outResult3->GetResultInfo()->complexFormat = REAL_IMAG;
             
             // CAUTION: begin, inc, end are hardcoded and noch checked for each result
             output->RegisterResult( outResult3, 1, 1,
                                     resultSteps[actRes].size(),
                                     false );

             // Hardcoded: set output format to AMPL_PHASE
             //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
             outResult4->GetResultInfo()->complexFormat = REAL_IMAG;
             
             // CAUTION: begin, inc, end are hardcoded and noch checked for each result
             output->RegisterResult( outResult4, 1, 1,
                                     resultSteps[actRes].size(),
                                     false );
             // Hardcoded: set output format to AMPL_PHASE
             //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
             outResult5->GetResultInfo()->complexFormat = REAL_IMAG;
             
             // CAUTION: begin, inc, end are hardcoded and noch checked for each result
             output->RegisterResult( outResult5, 1, 1,
                                     resultSteps[actRes].size(),
                                     false );
           }
         }
       }
       
       // notify writer
       if( output) {
         output->BeginMultiSequenceStep( actMsStep, types[actMsStep],
                                         numSteps[actMsStep] );
       }
       
       // iterate over all frequency steps
       for( UInt iStep = 0; iStep < numSteps[actMsStep]; iStep++ ) {
         
         // check, if current step contains any results
         if( stepVals.find( iStep+1 ) == stepVals.end() )
           continue;
         UInt actStepNum = iStep+1;
         Double actStepVal = stepVals[iStep+1];
         Double freq = actStepVal;
         
         if( output) {
           output->BeginStep( actStepNum, actStepVal );
         }
         std::cout << "\n\t=============================================\n";
         std::cout << "\t  Treating step " << actStepNum << ", " << actStepVal
                   << " Hz\n";
         std::cout << "\t=============================================\n";
         
         // iterate over all results
         for( UInt iRes = 0; iRes < inResults1.GetSize(); iRes++) {
           
           std::string regionName = inResults1[iRes]->GetEntityList()->GetName();
           bool isSurf = surfMap[regionName];

           std::cout << "\n\t-- Computing " << (isSurf ? "surface" : "volume") << " weight vector for " <<
             inResults1[iRes]->GetResultInfo()->resultName << " on " 
                     << regionName << " --\n";
           
           // check if current result is defined within this step
           if( resultSteps[inResults1[iRes]->GetResultInfo()].find(actStepNum)
               == resultSteps[inResults1[iRes]->GetResultInfo()].end() ) {
             continue;
           }
           
           // obtain both result objects for current step
           inputs_[PRIMARY_MODE]->GetResult( actMsStep, actStepNum, inResults1[iRes], false );
           inputs_[SECONDARY_MODE]->GetResult( actMsStep, actStepNum, inResults2[iRes], false );
           inputs_[MEAN_FLOW]->GetResult( actMsStep, actStepNum, inResults_mean_flow[iRes], false );
           
           
           std::set<std::string> volRegionNames;
           std::set<RegionIdType> volRegions;
           // volRegionNames.insert(inResults1[iRes]->GetEntityList()->GetName());
           std::map<std::string, Double>::iterator it, end;
           it = regionDensityMap.begin();
           end = regionDensityMap.end();
           for( ; it != end; it++ ) {
             const std::string& regionName = it->first;
             std::cout << "region: " << it->first << " density: " << it->second << std::endl;
             volRegionNames.insert( regionName );

             for( UInt rd=0, rnd=ptGrid1->regionData.GetSize(); rd < rnd; rd++ ) {
               if(ptGrid1->regionData[rd].name == regionName) 
               {
                 volRegions.insert(ptGrid1->regionData[rd].id);
               }
             }             
           }

           // get number of dofs of result
           //             UInt numDofs = inResults1[iRes]->GetResultInfo()->dofNames.GetSize();
           
           if( types[actMsStep] == BasePDE::HARMONIC ) { 
          #if 0
             Vector<Complex> & inVec1 =
               dynamic_cast<Result<Complex>& >(*inResults1[iRes]).GetVector();
             Vector<Complex> & inVec2 =
               dynamic_cast<Result<Complex>& >(*inResults2[iRes]).GetVector();
          #endif
             
             Vector<Complex> & outVec =
               dynamic_cast<Result<Complex>& >(*outResults[iRes]).GetVector();
             Vector<Complex> & outVec2 =
               dynamic_cast<Result<Complex>& >(*outResults2[iRes]).GetVector();
             Vector<Complex> & outVec3 =
               dynamic_cast<Result<Complex>& >(*outResults3[iRes]).GetVector();
             Vector<Complex> & outVec4 =
               dynamic_cast<Result<Complex>& >(*outResults4[iRes]).GetVector();
             Vector<Complex> & outVec5 =
               dynamic_cast<Result<Complex>& >(*outResults5[iRes]).GetVector();
             UInt numElems = outResults[iRes]->GetEntityList()->GetSize();
             UInt numDofs = outResults[iRes]->GetResultInfo()->dofNames.GetSize();
             UInt dim = ptGrid1->GetDim();

             outVec.Resize( numElems * numDofs );
             outVec2.Resize( numElems * numDofs );
             outVec3.Resize( numElems );
             outVec4.Resize( numElems * numDofs );
             outVec5.Resize( numElems );
             
             shared_ptr<BaseFeFunction> u1FeFunc = inputs_[PRIMARY_MODE]->GetFeFunction<Complex>( actMsStep,
                                                                             actStepNum,
                                                                             inResults1[iRes]->GetResultInfo()->resultType,
                                                                             volRegionNames);
             shared_ptr<BaseFeFunction> u2FeFunc = inputs_[SECONDARY_MODE]->GetFeFunction<Complex>( actMsStep,
                                                                             actStepNum,
                                                                             inResults2[iRes]->GetResultInfo()->resultType,
                                                                             volRegionNames);

             shared_ptr<BaseFeFunction> VFeFunc = inputs_[MEAN_FLOW]->GetFeFunction<Double>( actMsStep,
                                                                              actStepNum,
                                                                              inResults_mean_flow[iRes]->GetResultInfo()->resultType,
                                                                              volRegionNames );


             shared_ptr<BaseBOperator> gradOp;
             shared_ptr<BaseBOperator> curlOp;
             switch(dim) {
             case 2:
               gradOp.reset( new GradientOperator<FeH1,2> );
               curlOp.reset( new CurlOperator<FeH1,2,Double> );
               break;
             case 3:
               gradOp.reset( new GradientOperator<FeH1,3> );
               curlOp.reset( new CurlOperator<FeH1,3,Double> );
               break;               
             } 


             EntityIterator eIt = outResults[iRes]->GetEntityList()->GetIterator();
             
             Complex u_p_prime = Complex(0.0, 0.0);
             Double deltaPhi = 0.0;
             Double vol = 0.0;
             Double durchfluss = 0.0;
             Double u_p_real = 1.0;
             Double u_p_imag = 0;
             Vector<Complex> W_sum(3);
             W_sum.Init();
             Complex u_p = Complex(1.0, 0);
             
             PtrParamNode wvtNode = param->Get("wvt");
             if(wvtNode->Has("u_p")) {
               PtrParamNode upNode = wvtNode->Get("u_p");
               
               if(upNode->Has("real")) {
                 u_p_real = upNode->Get("real")->As<Double>();
               }
               if(upNode->Has("imag")) {
                 u_p_imag = upNode->Get("imag")->As<Double>();
               }
               u_p = Complex(u_p_real, u_p_imag);
               u_p *= Complex(0,1);
               u_p *= 2*PI*freq;
             }
             if(wvtNode->Has("v_p")) {
               PtrParamNode vpNode = wvtNode->Get("v_p");
               
               if(vpNode->Has("real")) {
                 u_p_real = vpNode->Get("real")->As<Double>();
               }
               if(vpNode->Has("imag")) {
                 u_p_imag = vpNode->Get("imag")->As<Double>();
               }
             }

             StdVector<std::string> realVal(3);
             realVal[0] = "0";
             realVal[1] = "0";
             realVal[2] = "1";
             Double meanVel = 1;
             if(wvtNode->Has("V")) {
               PtrParamNode vNode = wvtNode->Get("V");
               
               if(vNode->Has("meanVel")) {
                 meanVel = vNode->Get("meanVel")->As<Double>();
               }

               if(vNode->Has("profile")) {
                 PtrParamNode profileNode = vNode->Get("profile");
                 
                 // realVal[2] = "1-(sqrt(x^2+y^2)/13.15e-3)^10";
                 realVal[2] = profileNode->As<std::string>();
               }
             }

             PtrCoefFct coefV = CoefFunction::Generate(mp.get(), Global::REAL, realVal);
             shared_ptr<DefaultCoordSystem> coordSys( new DefaultCoordSystem(ptGrid1) );
             coefV->SetCoordinateSystem(coordSys.get());
             

             FeFunction<Double>* feFctPtr = new FeFunction<Double>(mp.get());
             shared_ptr<BaseFeFunction> feFct(feFctPtr);
             feFct->SetFeSpace( VFeFunc->GetFeSpace() );
             feFct->AddEntityList( inResults_mean_flow[0]->GetEntityList() );
             feFct->SetGrid(ptGrid3);
             feFct->SetResultInfo(inResults_mean_flow[iRes]->GetResultInfo());
             feFct->AddExternalDataSource( coefV,
                                           inResults_mean_flow[0]->GetEntityList() ); // <--- Hier muss unbedingt die zugehörige Volumenregion zur Oberflächenregion hinein!
             feFct->SetFctId(PSEUDO_FCT_ID);
             feFct->Finalize();
             feFct->ApplyExternalData();

             //CoupledField::PtrCoefFct&, boost::shared_ptr<CoupledField::EntityList>
             //CoupledField::PtrCoefFct, boost::shared_ptr<CoupledField::EntityList>&
           

             Double rho_0 = regionDensityMap[inResults1[iRes]->GetEntityList()->GetName()];


             StdVector< Vector<Double> > globMidPointCoords;

             for(UInt elIdx=0 ; !eIt.IsEnd(); eIt++, elIdx++ ) 
             {
               const Elem* el1 = eIt.GetElem();
               
               // Obtain FE element from feSpace and integration scheme
               LocPoint lp = Elem::shapes[el1->type].midPointCoord;
               shared_ptr<ElemShapeMap> esm1 = ptGrid1->GetElemShapeMap( el1, true );
        
               Vector<Double> p(dim);
               esm1->Local2Global(p, lp);

               globMidPointCoords.Push_back(p);

               //               if(el1->elemNum == 2771) 
               //               {
               //                 std::cout << "elIdx: " << elIdx << std::endl;
               //               }
             }

             StdVector< LocPoint > localCoords;
             StdVector< const Elem* > elems;
             std::set<RegionIdType> srcRegions;
             ptGrid2->GetElemsAtGlobalCoords(globMidPointCoords, localCoords, elems, srcRegions);

             eIt = outResults[iRes]->GetEntityList()->GetIterator();
             for(UInt elIdx=0 ; !eIt.IsEnd(); eIt++, elIdx++ ) 
             {
               const Elem* el1 = eIt.GetElem();
               const Elem* el2 = ptGrid2->GetElem(el1->elemNum);
               const Elem* el3 = ptGrid3->GetElem(el1->elemNum);
               
               // Obtain FE element from feSpace and integration scheme
               IntegOrder order;
               IntScheme::IntegMethod method;
               BaseFE* ptFe1 = u1FeFunc->GetFeSpace()->GetFe( eIt, method, order );
               BaseFE* ptFe2 = u2FeFunc->GetFeSpace()->GetFe( el2->elemNum );
               BaseFE* ptFe3 = NULL;

               // Get integration points
               shared_ptr<IntScheme> intScheme = u1FeFunc->GetFeSpace()->GetIntScheme();
               
               StdVector<LocPoint> intPoints;
               StdVector<Double> weights;
               intScheme->GetIntPoints( Elem::GetShapeType(el1->type), method, order,
                                        intPoints, weights );


               // std::cout << "Integration order " << order.GetIsoOrder() << " num. int points: " << intPoints.GetSize() << " weight: " << weights[0] << std::endl;

               for( UInt i=0; i < numDofs; i++ )
               {
                 // Init weight vector
                 outVec[elIdx*numDofs+i] = 0.0;
                 
                 // Init mean velocity
                 outVec2[elIdx*numDofs+i] = 0.0;
                 
                 // Init weight vector W_Phi as given in Hemp1994 (19).
                 outVec4[elIdx*numDofs+i] = 0.0;
               }
               
               // Init weight vector density result
               outVec3[elIdx] = 0.0;

               // Init weight  vector density  result  as inner  product
               // between V and W_Phi as given in Hemp1994 (18).
               outVec5[elIdx] = 0.0;
               
               Double elementVolume = 0.0;


               shared_ptr<LocPointMapped> lpm1, lpm2, lpm3;
               LocPointMapped lpm1Surf, lpm2Surf, lpm3Surf;
               // BaseFE* ptFe1Surf = NULL;
               // BaseFE* ptFe2Surf = NULL;
               // BaseFE* ptFe3Surf = NULL;
               // const Elem* el1Surf = NULL;
               // const Elem* el2Surf = NULL;
               // const Elem* el3Surf = NULL;
               shared_ptr<ElemShapeMap> esm1 = ptGrid1->GetElemShapeMap( el1, true );
               shared_ptr<ElemShapeMap> esm2 = ptGrid2->GetElemShapeMap( el2, true );
               shared_ptr<ElemShapeMap> esm3 = ptGrid3->GetElemShapeMap( el3, true );
               
               if( isSurf ) 
               {
                 // el1Surf = el1;
                 // ptFe1Surf = ptFe1;
                 
                 // el2Surf = el2;
                 // ptFe2Surf = ptFe2;

                 // el3Surf = el3;
                 // ptFe3Surf = ptFe3;
               }
               
               for(UInt ip=0, nip=intPoints.GetSize(); ip < nip; ip++) 
               {

                 // BaseFE* ptFe3 = VFeFunc->GetFeSpace()->GetFe( el3->elemNum );
                 //LocPoint lp = Elem::shapes[el1->type].midPointCoord;
                 LocPoint lp = intPoints[ip];

                 if( isSurf ) 
                 {
                   lpm1Surf.Set( lp, esm1, volRegions, 0.0 );
                   lpm1 = lpm1Surf.lpmVol;
                   lpm2Surf.Set( lp, esm2, volRegions, 0.0 );
                   lpm2 = lpm2Surf.lpmVol;
                   lpm3Surf.Set( lp, esm3, volRegions, 0.0 );
                   lpm3 = lpm3Surf.lpmVol;

                   el1 = lpm1->ptEl;
                   ptFe1 = u1FeFunc->GetFeSpace()->GetFe( lpm1->ptEl->elemNum );
                   
                   el2 = lpm2->ptEl;
                   ptFe2 = u2FeFunc->GetFeSpace()->GetFe( lpm2->ptEl->elemNum );
                   
                   el3 = lpm3->ptEl;
                   ptFe3 = feFct->GetFeSpace()->GetFe( lpm3->ptEl->elemNum );
                 } else 
                 {
                   lpm1.reset(new LocPointMapped());
                   lpm1->Set( lp, esm1, 0.0 );
                   lpm2.reset(new LocPointMapped());
                   lpm2->Set( lp, esm2, 0.0 );
                   lpm3.reset(new LocPointMapped());
                   lpm3->Set( lp, esm3, 0.0 );

                   ptFe3 = feFct->GetFeSpace()->GetFe( lpm3->ptEl->elemNum );
                 }
               
               

                 // const Elem* el2_new = elems[elIdx];
                 //               std::cout << "el2: " << el2->elemNum <<  " el2_new: " << el2_new << " " << std::endl;
                 Vector<Double> p(dim);
                 esm1->Local2Global(p, lp);

                 // Matrix<Double> bMat1, bMat2, bMat3;
                 Vector<Complex> eSol1, eSol2;
                 Vector<Double> eSol3;
                 StdVector< Vector<Complex> > eSol1Comp(numDofs), eSol2Comp(numDofs);
                 StdVector< Vector<Double> > eSol3Comp(numDofs);
                 Vector<Complex> eSol2_x, eSol2_y, esol2_z;
                 Vector<Complex> u1, u2;
                 Vector<Double> V;
                 Vector<Complex> W;
                 Vector<Complex> X;

                 // gradOp->CalcOpMat( bMat1, *lpm1, ptFe1 );
                 // gradOp->CalcOpMat( bMat2, *lpm2, ptFe2 );
                 //curlOp->CalcOpMat( bMat3, *lpm1, ptFe1 );
                 // gradOp->CalcOpMat( bMat3, *lpm1, ptFe1 );
                 
                 u1FeFunc->GetVector( u1, *lpm1 );
                 u2FeFunc->GetVector( u2, *lpm2 );
                 // Get V from HDF5 file
                 // VFeFunc->GetVector( V, *lpm3 );
                 // Get analytical V
                 // coefV->GetVector(V, *lpm1);
                 feFct->GetVector(V, *lpm3);

#if 0
                 Double radius = sqrt(p[0]*p[0] + p[1]*p[1]);
                 Double R = 0.01315;
                 V[0] = 0;
                 V[1] = 0;
                 V[2] = 1-(radius*radius/(R*R));
#endif

                 dynamic_cast<FeFunction<Complex>*>(u1FeFunc.get())->GetElemSolution( eSol1, el1 );
                 dynamic_cast<FeFunction<Complex>*>(u2FeFunc.get())->GetElemSolution( eSol2, el2 );
                 dynamic_cast<FeFunction<Double>*>(feFct.get())->GetElemSolution( eSol3, el3 );
                 
                 StdVector< Vector<Complex> > u1Derivs(numDofs), u2Derivs(numDofs);
                 Vector<Double> curlV;
                 StdVector< Vector<Complex> > u2Tau(numDofs);
                 Double viscosity = regionViscosityMap[ptGrid1->regionData[el1->regionId].name];
                 
                 for(UInt dof=0; dof < numDofs; dof++) 
                 {
                   UInt eN=eSol1.GetSize()/numDofs;
                   
                   eSol1Comp[dof].Resize(eN);
                   eSol2Comp[dof].Resize(eN);
                   eSol3Comp[dof].Resize(eN);
                   
                   for(UInt eI=0; eI < eN; eI++) 
                   {
                     eSol1Comp[dof][eI] = eSol1[eI*numDofs+dof];
                     eSol2Comp[dof][eI] = eSol2[eI*numDofs+dof];
                     eSol3Comp[dof][eI] = eSol3[eI*numDofs+dof];
                   }
                   
                   gradOp->ApplyOp( u1Derivs[dof], *lpm1, ptFe1, eSol1Comp[dof] );
                   gradOp->ApplyOp( u2Derivs[dof], *lpm2, ptFe2, eSol2Comp[dof] );
                   //                   gradOp->ApplyOp( curlV, *lpm3, ptFe3, eSol3Comp[dof] );
                 }
                 
                 curlOp->ApplyOp( curlV, *lpm3, ptFe3, eSol3 );

                 //                std::cout << "elemNum " << el1->elemNum << " numDofs " << numDofs << " dim " << dim << std::endl;
                 Vector<Complex> Xl(numDofs), Xr(numDofs);
                 Vector<Double> n(numDofs);

                 if(isSurf) 
                 {
                   n = -lpm3Surf.normal;
                 } else 
                 {
                   n.Init();
                 }                 

                 Xl.Resize(numDofs);
                 Xr.Resize(numDofs);
                 X.Resize(numDofs);
                 Xl.Init();
                 Xr.Init();
                 X.Init();

#if 0
                 Vector<Double> t(numDofs);
                 t[0] = n[1];
                 t[1] = -n[0];
                 t[2] = 0;
                 t *= 1/t.NormL2();
                 
                 curlV = 2*radius/(R*R) * t;
#endif

                 for(UInt dof=0; dof < numDofs; dof++) 
                 {
                   UInt eN=u2Derivs[dof].GetSize();
                   u2Tau[dof].Resize(eN);

                   for(UInt eI=0; eI < eN; eI++) 
                   {
                     if(dof == eI) 
                     {
                       u2Tau[dof][eI] = 0.0;
                     }
                     else
                     {
                       u2Tau[dof][eI] = u2Derivs[dof][eI]+u2Derivs[eI][dof];
                     }
                     
                     if(isSurf) 
                     {
                       Xl[dof] += u2Tau[dof][eI] * n[eI];
                     }                     
                   }
                 }
                 
                 if(isSurf) 
                 {
                   Complex dp;
                   dp = u1[0] * n[0] + u1[1] * n[1] + u1[2] * n[2];
                   Xr[0] = n[0] * dp;
                   Xr[1] = n[1] * dp;
                   Xr[2] = n[2] * dp;
                   //                 CrossProd(Xl, Xr, X);
                   
                   Complex complexFreq(0, 2*PI*freq);
                   
                   X[0] =  Xl[1] * Xr[2] - Xl[2] * Xr[1];
                   X[1] = -Xl[0] * Xr[2] + Xl[2] * Xr[0];
                   X[2] =  Xl[0] * Xr[1] - Xl[1] * Xr[0];

                   X *= viscosity / complexFreq;
                 }                 

                 // std::cout << "elemNum " << el3->elemNum << " curl(V) " << curlV[2] << " V " << V << std::endl;
                 // std::cout << "eSol3 " << eSol3 << V << std::endl;

                 W.Resize(u1.GetSize());
                 W.Init();
                 
                 for( UInt i=0; i < numDofs; i++ )
                 {
                   for( UInt j=0; j < dim; j++ )
                   {
                  #if 0
                     W[i] += -rho_0 * (u2[j] * u1Derivs[i][j] - u1[j] * u2Derivs[j][i]);
                  #endif               
                     
                     if(isSurf) 
                     {
                       // W[i] += lpm3Surf.normal[i];
                       // W[i] += curlV[2][i];
                       W[i] = V[i];
                     }
                     else 
                     {
                       // Compute components of weight vector as given in Hemp1994 (15).
                       W[i] += -rho_0 * (u2[j] * u1Derivs[j][i] - u1[j] * u2Derivs[i][j]);
                     }
                     
                   }
                 }
                 
                 Double volume = 0.0;
                 
                 if(isSurf) 
                 {
                   volume = lpm1Surf.jacDet * weights[ip];
                 }
                 else
                 {
                   volume = (*lpm1).jacDet * weights[ip];
                 }
                 
                 elementVolume += volume;

                 if (isSurf) 
                 {
                   durchfluss = 1.0;
                   
                   for( UInt i=0; i < numDofs; i++ )
                   {
                     // Prepare weight vector
                     outVec[elIdx*numDofs+i] += X[i] * volume;
                     
                     // Prepare mean velocity
                     outVec2[elIdx*numDofs+i] += V[i] * volume;
                     
                     // Prepare weight vector density result
                     outVec3[elIdx] += X[i]*curlV[i] * volume;
                     
                     
                     // Prepare weight vector W_Phi as given in Hemp1994 (19).
                     outVec4[elIdx*numDofs+i] = curlV[i] * volume;
                     
                     // Prepare  weight  vector density  result  as inner  product
                     // between V and W_Phi as given in Hemp1994 (18).
                     //                     outVec5[elIdx] += W[i]/u_p * V[i];

                     outVec5[elIdx] += X[i]/u_p * curlV[i] * volume;
                   }

                 }
                 else 
                 {         
                   for( UInt i=0; i < numDofs; i++ )
                   {
                     // Prepare weight vector
                     outVec[elIdx*numDofs+i] += W[i] * volume;
                     
                     W_sum[i] += W[i] * volume;
                     
                     // Prepare mean velocity
                     //                   V[i] *= meanVel;                   
                     
                     durchfluss += V[i] * volume;
                     
                     outVec2[elIdx*numDofs+i] += V[i] * volume;
                     
                     // Prepare weight vector density result
                     outVec3[elIdx] += W[i]*V[i]*volume;
                     
                     
                     // Prepare weight vector W_Phi as given in Hemp1994 (19).
                     outVec4[elIdx*numDofs+i] += W[i]/u_p * volume;
                     
                     // Prepare  weight  vector density  result  as inner  product
                     // between V and W_Phi as given in Hemp1994 (18).
                     outVec5[elIdx] += W[i]/u_p * V[i] * volume;
                   }
                 }
                 
                 // u_p_prime (0.000359091,-2.8713e-05)  
                 // Compute u_p_prime as given in Hemp1994 (13).
                 u_p_prime += outVec3[elIdx];

                 // Compute deltaPhi as given in Hemp1994 (18) and (19).
                 deltaPhi += outVec5[elIdx].imag();

               } // loop over integration points

               for( UInt i=0; i < numDofs; i++ )
               {
                 outVec[elIdx*numDofs+i] /= elementVolume;
                 outVec2[elIdx*numDofs+i] /= elementVolume;
                 outVec4[elIdx*numDofs+i] /= elementVolume;
               }               
               outVec3[elIdx] /= elementVolume;
               outVec5[elIdx] /= elementVolume;

               vol += elementVolume;
             }
             W_sum /= vol;

             if(!isSurf) 
             {
               correct = (meanVel*vol/durchfluss);
             }
             
             std::cout << "\nProfilkorrekturfaktor " << correct << std::endl;
             std::cout << "Volumen " << vol << std::endl;
             outVec2 *= correct;
             
             deltaPhi *= correct;
             
             boost::filesystem::ofstream csv;

             if(isSurf)
             {
               csv.open("output.csv", std::ios_base::binary | std::ios_base::app);
             }
             else
             {
               csv.open("output.csv", std::ios_base::binary);
               csv << "freq,u_p_prime_real,u_p_prime_imag,u_p_prime_ampl,u_p_prime_phase,deltaPhi[rad],deltaPhi[deg],volume,real(W_x),imag(W_x),real(W_y),imag(W_y),real(W_z),imag(W_z)" << std::endl;
             }
             
             std::cout << "\nu_p_prime " << u_p_prime << std::endl;
             std::cout << "\ndeltaPhi [rad]: " << deltaPhi << std::endl;
             std::cout << "\ndeltaPhi [deg]: " << deltaPhi/3.14159*180 << std::endl;

             csv << freq << ","
                 << u_p_prime.real() << ","
                 << u_p_prime.imag() << ","
                 << abs(u_p_prime) << ","
                 << arg(u_p_prime) << ","
                 << deltaPhi << ","
                 << (deltaPhi/PI*180) << ","
                 << vol << ","
                 << W_sum[0].real() << ","
                 << W_sum[0].imag() << ","
                 << W_sum[1].real() << ","
                 << W_sum[1].imag() << ","
                 << W_sum[2].real() << ","
                 << W_sum[2].imag() << std::endl;

             csv.close();         
           } // switch: Analysis type
           // add result to output file
           if (output ) {
             output->AddResult( outResults[iRes] );
             output->AddResult( outResults2[iRes] );
             output->AddResult( outResults3[iRes] );
             output->AddResult( outResults4[iRes] );
             output->AddResult( outResults5[iRes] );
           }
           
         } // loop over results
         
         if( output )
           output->FinishStep();
       }
       if( output )
         output->FinishMultiSequenceStep();
       
       if( output ) {

         std::stringstream sstr;
         param->ToXML(sstr, 2);
         
         WriteStringToHDF5UserData(output, "cfstoolConfig", sstr.str());

         output->Finalize();
         std::cout << "\nOutput successfully written to " << outFile_ << std::endl;
       }
       delete ptGrid1;
       delete ptGrid2;
       
  } //WVT

}

