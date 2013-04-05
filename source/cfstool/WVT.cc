// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_hdf5.hh>

#include "DataInOut/SimOutput.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeSpace.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Forms/Operators/GradientOperator.hh"

#ifdef USE_HDF5
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#endif

#include "ParamsInit.hh"
#include "HelperFuncs.hh"
#include "WVT.hh"

using namespace CoupledField;

namespace CFSTool {

  WVT::WVT( const std::string& lateral_mode_file,
            const std::string& coriolis_mode_file,
            const std::string& mean_flow_file,
            const std::string& outFile,
            const PtrParamNode& param,
            const PtrParamNode& info) :
    param_(param),
    info_(info),
    writeOutputFile_(outFile != ""),
    inputs_(3)
  {

       bool printGridOnly = false;

       inputs_[0] = GetInputReader( lateral_mode_file, param_, info_ );
       inputs_[1] = GetInputReader( coriolis_mode_file, param_, info_ );
       inputs_[2] = GetInputReader( mean_flow_file, param_, info_ );

       // check capabilities of input class
       StdVector<bool> readerCaps(3);
       StdVector<std::string> readerDescriptions(3);
       readerDescriptions[0] = "lateral mode";
       readerDescriptions[1] = "coriolis mode";
       readerDescriptions[2] = "mean flow";
       StdVector<std::string> readerCapsResults(3);
       for(UInt i=0; i<inputs_.GetSize(); i++) 
       {
         readerCaps[i] = CheckReaderCapabilities(readerDescriptions[i],
                                                 inputs_[i],
                                                 readerCapsResults[i]);
       }
       
       if( !readerCaps[0] || !readerCaps[1] || !readerCaps[2] ) {
         std::cerr << "Some input files are only capable of handling meshes, not results!\n\n";

         for(UInt i=0; i<inputs_.GetSize(); i++) 
         {
           std::cerr << readerCapsResults[i] << std::endl;
         }

         exit(EXIT_FAILURE);
       }

       // read in mesh of input1
       inputs_[0]->InitModule();
       UInt dim = inputs_[0]->GetDim();
       Grid * ptGrid1 = new GridCFS(dim, param_, info_);
       inputs_[0]->ReadMesh(ptGrid1);
       ptGrid1->FinishInit();

       // std::cout << "##############" << xmlFile << std::endl
       //    << "##############" << matFile << std::endl;
       PtrParamNode pNode = GetParamNodeFromHDF5(inputs_[0], "ParameterFile");

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
       ParamNodeList domainRegionList = pNode->Get("domain")->Get("regionList")->GetChildren();

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

       pNode = GetParamNodeFromHDF5(inputs_[0], "MaterialFile");

       std::map<std::string, Double> regionDensityMap;

       regMatIt = regionMatMap.begin();
       regMatEnd = regionMatMap.end();
       for( ; regMatIt != regMatEnd; regMatIt++ ){
         std::string matName = regMatIt->second;
         PtrParamNode mNode = pNode->GetByVal("material",
                                              "name",
                                              matName);

         Double density = mNode->Get("flow")->Get("density")->As<Double>();
         std::cout << "density " << density << std::endl;         
         regionDensityMap[matName] = density;
       }

       // read in mesh of input2
       inputs_[1]->InitModule();
       Grid * ptGrid2 = new GridCFS(dim, param_, info_);
       inputs_[1]->ReadMesh(ptGrid2);
       ptGrid2->FinishInit();

       // read in mesh of input3
       inputs_[2]->InitModule();
       Grid * ptGrid3 = new GridCFS(dim, param_, info_);
       inputs_[2]->ReadMesh(ptGrid3);
       ptGrid3->FinishInit();

       // obtain output writer
       shared_ptr<SimOutput> output;
       if( outFile != "" ) {
         output = GetOutputWriter( outFile, param_, info_ );
         output->Init( ptGrid1, printGridOnly);
       }

       // obtain number of Sequence Steps and get analysis types
       std::map<UInt, BasePDE::AnalysisType> types;
       std::map<UInt, UInt> numSteps;
       inputs_[0]->GetNumMultiSequenceSteps( types, numSteps, false );

       std::cout << "\nFound " << types.size() << " sequence step(s) in '" << inputs_[0]->GetFileName() << "'\n";
       std::map<UInt, BasePDE::AnalysisType> types2;
       std::map<UInt, UInt> numSteps2;
       inputs_[1]->GetNumMultiSequenceSteps( types2, numSteps2, false );
       std::cout << "\nFound " << types2.size() << " sequence step(s) in '" << inputs_[1]->GetFileName() << "'\n";
       
       if(types.size() != types2.size()){
         std::cout << "'" << inputs_[0]->GetFileName() << "' and '" << inputs_[1]->GetFileName()
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
       inputs_[0]->GetResultTypes( actMsStep, infos, false );
       inputs_[1]->GetResultTypes( actMsStep, infos2, false );
       inputs_[2]->GetResultTypes( actMsStep, infos3, false );
       
       StdVector<shared_ptr<BaseResult> > inResults1, inResults2, inResults_mean_flow, outResults, outResults2, outResults3;
       // stepnumbers, for which at least one result is defined
       std::map<UInt, Double> stepVals, stepVals2, stepVals_mean_flow;
       // contains the stepnumbers/-values in which the particular result is
       // defined in
       std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
       std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps2;
       std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps_mean_flow;
       
       // First determine number of results in current multi sequence step in mean flow file
       // since it can be different than in the lateral and coriolis mode files.
       
       for( UInt iRes = 0; iRes < infos3.GetSize(); iRes++) {
         shared_ptr<ResultInfo> actRes3 = infos3[iRes];
         
         if(actRes3->resultType != FLUIDMECH_VELOCITY) 
         {
           continue;
         }
         
         // get stepvalues of mean flow file
         inputs_[2]->GetStepValues( actMsStep, actRes3,
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
         
         if(infos[iRes]->resultName != "fluidMechVelocity")
           continue;
         
         std::cout << "\t" << infos[iRes]->resultName << "\n";
         
         // get stepvalues of reference file
         shared_ptr<ResultInfo> actRes = infos[iRes];
         shared_ptr<ResultInfo> actRes2 = infos2[iRes];
         inputs_[0]->GetStepValues( actMsStep, actRes,
                                resultSteps[actRes], false);
         stepVals.insert( resultSteps[actRes].begin(),
                          resultSteps[actRes].end() );
         
         // get stepvalues of second file
         inputs_[1]->GetStepValues( actMsStep, actRes2,
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
         StdVector<shared_ptr<EntityList> > regions;
         inputs_[0]->GetResultEntities( actMsStep, infos[iRes],
                                    regions, false );
         for( UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++ ) {
           // generate new result object and add it to output writer
           shared_ptr<BaseResult > inResult1, inResult2, inResult3, outResult, outResult2, outResult3;
           
           if(infos[iRes]->resultType != FLUIDMECH_VELOCITY) 
           {
             continue;
           }
           
           
           inResult1 = shared_ptr<BaseResult>( new Result<Complex>() );
           inResult2 = shared_ptr<BaseResult>( new Result<Complex>() );
           inResult3 = shared_ptr<BaseResult>( new Result<Double>() );
           outResult = shared_ptr<BaseResult>( new Result<Complex>() );
           outResult2 = shared_ptr<BaseResult>( new Result<Complex>() );
           outResult3 = shared_ptr<BaseResult>( new Result<Complex>() );
           
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
           
           inResults1.Push_back( inResult1 );
           inResults2.Push_back( inResult2 );
           inResults_mean_flow.Push_back( inResult3 );
           outResults.Push_back( outResult );
           outResults2.Push_back( outResult2 );
           outResults3.Push_back( outResult3 );
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
         
         if( output) {
           output->BeginStep( actStepNum, actStepVal );
         }
         std::cout << "\n\t=============================================\n";
         std::cout << "\t  Treating step " << actStepNum << ", " << actStepVal
                   << " Hz\n";
         std::cout << "\t=============================================\n";
         
         // iterate over all results
         for( UInt iRes = 0; iRes < inResults1.GetSize(); iRes++) {
           
           if(inResults1[iRes]->GetResultInfo()->resultName != "fluidMechVelocity")
             continue;
           
           // check if current result is defined within this step
           if( resultSteps[inResults1[iRes]->GetResultInfo()].find(actStepNum)
               == resultSteps[inResults1[iRes]->GetResultInfo()].end() ) {
             continue;
           }
           
           // obtain both result objects for current step
           inputs_[0]->GetResult( actMsStep, actStepNum, inResults1[iRes], false );
           inputs_[1]->GetResult( actMsStep, actStepNum, inResults2[iRes], false );
           inputs_[2]->GetResult( actMsStep, actStepNum, inResults_mean_flow[iRes], false );
           
           std::cout << "\n\t-- Comparing result " <<
             inResults1[iRes]->GetResultInfo()->resultName << " on " 
                     << inResults1[iRes]->GetEntityList()->GetName() << " --\n";
           
           std::set<std::string> regionNames;
           regionNames.insert(inResults1[iRes]->GetEntityList()->GetName());
           
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
             UInt numElems = outResults[iRes]->GetEntityList()->GetSize();
             UInt numDofs = outResults[iRes]->GetResultInfo()->dofNames.GetSize();
             UInt dim = ptGrid1->GetDim();

             outVec.Resize( numElems * numDofs );
             outVec2.Resize( numElems * numDofs );
             outVec3.Resize( numElems );
             
             shared_ptr<BaseFeFunction> u1FeFunc = inputs_[0]->GetFeFunction<Complex>( actMsStep,
                                                                             actStepNum,
                                                                             inResults1[iRes]->GetResultInfo()->resultType,
                                                                             regionNames);
             shared_ptr<BaseFeFunction> u2FeFunc = inputs_[1]->GetFeFunction<Complex>( actMsStep,
                                                                             actStepNum,
                                                                             inResults2[iRes]->GetResultInfo()->resultType,
                                                                             regionNames);

             shared_ptr<BaseFeFunction> VFeFunc = inputs_[2]->GetFeFunction<Double>( actMsStep,
                                                                              actStepNum,
                                                                              inResults_mean_flow[iRes]->GetResultInfo()->resultType,
                                                                              regionNames );


             GradientOperator<FeH1,3> gradOp;
             
             EntityIterator eIt = outResults[iRes]->GetEntityList()->GetIterator();
             
             Complex u_p_prime = Complex(0.0, 0.0);
             Double rho_0 = regionDensityMap[inResults1[iRes]->GetEntityList()->GetName()];


             StdVector< Vector<Double> > globMidPointCoords;

             for(UInt elIdx=0 ; !eIt.IsEnd(); eIt++, elIdx++ ) 
             {
               const Elem* el1 = eIt.GetElem();
               
               // Obtain FE element from feSpace and integration scheme
               LocPoint lp = Elem::shapes[el1->type].midPointCoord;
               shared_ptr<ElemShapeMap> esm1 = ptGrid1->GetElemShapeMap( el1, true );
        
               Vector<Double> p(3);
               esm1->Local2Global(p, lp);

               globMidPointCoords.Push_back(p);

               if(el1->elemNum == 2771) 
               {
                 std::cout << "elIdx: " << elIdx << std::endl;
               }
               
             }

             StdVector< LocPoint > localCoords;
             StdVector< const Elem* > elems;
             std::set<RegionIdType> srcRegions;
             ptGrid2->GetElemsAtGlobalCoords(globMidPointCoords, localCoords, elems, srcRegions);

             eIt = outResults[iRes]->GetEntityList()->GetIterator();
             for(UInt elIdx=0 ; !eIt.IsEnd(); eIt++, elIdx++ ) 
             {
               if(inResults1[iRes]->GetEntityList()->GetName() == "pipe_inner") continue;
               
               const Elem* el1 = eIt.GetElem();
               const Elem* el2 = ptGrid2->GetElem(el1->elemNum);
               const Elem* el3 = ptGrid3->GetElem(el1->elemNum);
               
               // Obtain FE element from feSpace and integration scheme
               IntegOrder order;
               IntScheme::IntegMethod method;
               BaseFE* ptFe1 = u1FeFunc->GetFeSpace()->GetFe( eIt, method, order );
               BaseFE* ptFe2 = u2FeFunc->GetFeSpace()->GetFe( el2->elemNum );
               // BaseFE* ptFe3 = VFeFunc->GetFeSpace()->GetFe( el3->elemNum );
               LocPoint lp = Elem::shapes[el1->type].midPointCoord;
               LocPointMapped lpm1, lpm2, lpm3;
               shared_ptr<ElemShapeMap> esm1 = ptGrid1->GetElemShapeMap( el1, true );
               shared_ptr<ElemShapeMap> esm2 = ptGrid2->GetElemShapeMap( el2, true );
               shared_ptr<ElemShapeMap> esm3 = ptGrid3->GetElemShapeMap( el3, true );
               lpm1.Set( lp, esm1, 0.0 );
               lpm2.Set( lp, esm2, 0.0 );
               lpm3.Set( lp, esm3, 0.0 );
        
               const Elem* el2_new = elems[elIdx];
               std::cout << "el2: " << el2->elemNum <<  " el2_new: " << el2_new << " " << std::endl;
               

               Matrix<Double> bMat1, bMat2;
               Vector<Complex> eSol1, eSol2;
               StdVector< Vector<Complex> > eSol1Comp(numDofs), eSol2Comp(numDofs);
               Vector<Complex> eSol2_x, eSol2_y, esol2_z;
               Vector<Complex> u1, u2;
               Vector<Double> V;
               Vector<Complex> W;
               gradOp.CalcOpMat( bMat1, lpm1, ptFe1 );
               gradOp.CalcOpMat( bMat2, lpm2, ptFe2 );
               
               
               u1FeFunc->GetVector( u1, lpm1 );
               u2FeFunc->GetVector( u2, lpm2 );
               VFeFunc->GetVector( V, lpm3 );

               dynamic_cast<FeFunction<Complex>*>(u1FeFunc.get())->GetElemSolution( eSol1, el1 );
               dynamic_cast<FeFunction<Complex>*>(u2FeFunc.get())->GetElemSolution( eSol2, el2 );

               StdVector< Vector<Complex> > u1Derivs(numDofs), u2Derivs(numDofs);

               for(UInt dof=0; dof < numDofs; dof++) 
               {
                 UInt eN=eSol1.GetSize()/numDofs;
                 
                 eSol1Comp[dof].Resize(eN);
                 eSol2Comp[dof].Resize(eN);

                 for(UInt eI=0; eI < eN; eI++) 
                 {
                   eSol1Comp[dof][eI] = eSol1[eI*numDofs+dof];
                   eSol2Comp[dof][eI] = eSol2[eI*numDofs+dof];
                 }

                 gradOp.ApplyOp( u1Derivs[dof], lpm1, ptFe1, eSol1Comp[dof] );
                 gradOp.ApplyOp( u2Derivs[dof], lpm2, ptFe2, eSol2Comp[dof] );
               }

#if 0
               Matrix<Complex> u1Derivs(3,3);
               Matrix<Complex> u2Derivs(3,3);
               u1Derivs.Init();
               u2Derivs.Init();

               UInt numFncs = ptFe1->GetNumFncs();
               
               for( UInt dof=0; dof < numDofs; dof++ ) 
               {
                 for( UInt f=0; f < numFncs; f++ ) 
                 {
                   u1Derivs[dof][0] += eSol1[f*numDofs+0]*bMat1[dof][f*numDofs+0];
                   u1Derivs[dof][1] += eSol1[f*numDofs+1]*bMat1[dof][f*numDofs+1];
                   u1Derivs[dof][2] += eSol1[f*numDofs+2]*bMat1[dof][f*numDofs+2];
                   
                   u2Derivs[dof][0] += eSol2[f*numDofs+0]*bMat2[dof][f*numDofs+0];
                   u2Derivs[dof][1] += eSol2[f*numDofs+1]*bMat2[dof][f*numDofs+1];
                   u2Derivs[dof][2] += eSol2[f*numDofs+2]*bMat2[dof][f*numDofs+2];
                 }                   
               }
#endif               
               
               std::cout << "elemNum " << el1->elemNum << std::endl;
               W.Resize(u1.GetSize());
               W.Init();
               
               for( UInt i=0; i < numDofs; i++ )
               {
                 for( UInt j=0; j < dim; j++ )
                 {
#if 0
                   W[i] += -rho_0 * (u2[j] * u1Derivs[i][j] - u1[j] * u2Derivs[j][i]);
#endif               

                   W[i] += -rho_0 * (u2[j] * u1Derivs[j][i] - u1[j] * u2Derivs[i][j]);
                 }
               }
               
               outVec[elIdx*numDofs+0] = W[0];
               outVec[elIdx*numDofs+1] = W[1];
               outVec[elIdx*numDofs+2] = W[2];

#if 0               
               outVec[elIdx*numDofs+0] = u1Derivs[2][0];
               outVec[elIdx*numDofs+1] = u1Derivs[2][1];
               outVec[elIdx*numDofs+2] = u1Derivs[2][2];
#endif

               outVec2[elIdx*numDofs+0] = V[0];
               outVec2[elIdx*numDofs+1] = V[1];
               outVec2[elIdx*numDofs+2] = V[2];

               outVec3[elIdx] = W[0]*V[0] + W[1]*V[1] + W[2]*V[2];

               // Get integration points
               shared_ptr<IntScheme> intScheme = u1FeFunc->GetFeSpace()->GetIntScheme();
               
               StdVector<LocPoint> intPoints;
               StdVector<Double> weights;
               intScheme->GetIntPoints( Elem::GetShapeType(el1->type), method, 0,
                                         intPoints, weights );

               u_p_prime += outVec3[elIdx] * lpm1.jacDet * weights[0];
               
               
             }

             std::cout << "\nu_p_prime " << u_p_prime << std::endl;
         
           } // switch: Analysitype
           // add result to output file
           if (output )
             output->AddResult( outResults[iRes] );
             output->AddResult( outResults2[iRes] );
             output->AddResult( outResults3[iRes] );
           
         } // loop over results
         
         if( output )
           output->FinishStep();
       }
       if( output )
         output->FinishMultiSequenceStep();
       
       if( output ) {
         output->Finalize();
         std::cout << "\nOutput successfully written to " << outFile << std::endl;
       }
       delete ptGrid1;
       delete ptGrid2;
       
  } //WVT

}

