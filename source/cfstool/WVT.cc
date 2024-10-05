// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Domain/CoefFunction/CoefFunctionScatteredData.hh"
#include "Utils/mathParser/mathParser.hh"
#include "DataInOut/SimOutput.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "Domain/CoordinateSystems/DefaultCoordSystem.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeSpace.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/CurlOperator.hh"

#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"

#include "ParamsInit.hh"
#include "HelperFuncs.hh"
#include "MatFile.hh"
#include "WVT.hh"

using namespace CoupledField;

namespace CFSTool {

  // Definition of WVT input and output file type mappings
  static EnumTuple fileTypeTuples[] = {
    EnumTuple(WVT::PRIMARY_MODE,   "primary mode"), 
    EnumTuple(WVT::SECONDARY_MODE, "secondary mode"),
    EnumTuple(WVT::MEAN_FLOW,      "mean flow"),
    EnumTuple(WVT::EVAL_GRID,      "evaluation grid"),
    EnumTuple(WVT::OUTPUT,         "output")
  };
  Enum<WVT::FileType> WVT::fileType =                       \
       Enum<WVT::FileType>("WVT file types",
           sizeof(fileTypeTuples) / sizeof(EnumTuple),
           fileTypeTuples);

  // Definition of WVT mean flow data type mappings
  static EnumTuple meanFlowDataTypeTuples[] = {
    EnumTuple(WVT::MF_GRID_DATA,   "mean flow data on grid points"), 
    EnumTuple(WVT::MF_SCATTERED_DATA, "mean flow data on point cloud"), 
    EnumTuple(WVT::MF_ANALYTIC_EXP, "mean flow data given as analytic expression with respect to grid")
  };
  Enum<WVT::MeanFlowDataType> WVT::meanFlowDataType =                       \
       Enum<WVT::MeanFlowDataType>("WVT mean flow data types",
           sizeof(meanFlowDataTypeTuples) / sizeof(EnumTuple),
           meanFlowDataTypeTuples);

  WVT::WVT( const PtrParamNode& param,
            const PtrParamNode& info ) :
    param_(param),
    info_(info),
    dim_(0),
    numDofs_(0),
    writeOutputFile_ (false),
    sensorNodeName_("s_1"),
    integOrder_(1),
    dirCoupled_(false),
    u_p_(1.0, 0.0),
    nodalResults_(true),
    meanFlowType_(MF_GRID_DATA),
    writeIntPointsFile_(false),
    intPointsFileName_("intpoints.vtk"),
    intPointsFilePos_(0),
    numIntPoints_(0)
  {
  }


  void WVT::Initialize() 
  {
    PtrParamNode wvtNode = param_->Get("wvt");
    PtrParamNode inputNode = param_->Get("fileFormats")->Get("input");
    PtrParamNode outputNode = param_->Get("fileFormats")->Get("output");

    std::string evalGridId = "evalGrid";
    std::string priModeType = "grid";
    std::string priModeId = "u1";
    std::string secModeType = "grid";
    std::string secModeId = "u2";
    std::string meanFlowType = "grid";
    std::string meanFlowId = "V";

    priModeType = wvtNode->Get("primaryMode")->Get("type")->As<std::string>();
    if(priModeType == "grid") 
    {
      priModeId = wvtNode->Get("primaryMode")->Get("grid")
                  ->Get("id")->As<std::string>();

      // TODO: Hier auch andere reader als HDF5 zulassen
      priModeFile_ = inputNode->GetByVal("hdf5", "id", priModeId)
                     ->Get("fileName")->As<std::string>();
      std::cout << "Primary mode file: " << priModeFile_ << std::endl;

      // Generate input (HDF5) readers for primary and secondary modes.
      inputs_[PRIMARY_MODE] = GetInputReader( priModeFile_, param_, info_ );
    }
    else
    {
      EXCEPTION("Primary mode type '" << priModeType
                << "' cannot be handled yet!");      
    }

    secModeType = wvtNode->Get("secondaryMode")->Get("type")->As<std::string>();
    if(secModeType == "grid") 
    {
      secModeId = wvtNode->Get("secondaryMode")->Get("grid")
                  ->Get("id")->As<std::string>();

      // TODO: Hier auch andere reader als HDF5 zulassen
      secModeFile_ = inputNode->GetByVal("hdf5", "id", secModeId)
                     ->Get("fileName")->As<std::string>();
      std::cout << "Secondary mode file: " << secModeFile_ << std::endl;

      inputs_[SECONDARY_MODE] = GetInputReader( secModeFile_, param_, info_ );
    }
    else
    {
      EXCEPTION("Secondary mode type '" << secModeType
                << "' cannot be handled yet!");      
    }

    meanFlowType = wvtNode->Get("meanFlow")->Get("type")->As<std::string>();
    if(meanFlowType == "grid") 
    {
      meanFlowId = wvtNode->Get("meanFlow")->Get("grid")
                  ->Get("id")->As<std::string>();

      // TODO: Hier auch andere reader als HDF5 zulassen
      meanFlowFile_ = inputNode->GetByVal("hdf5", "id", meanFlowId)
                     ->Get("fileName")->As<std::string>();
      std::cout << "Mean flow file: " << meanFlowFile_ << std::endl;
      meanFlowType_ = MF_GRID_DATA;
    }
    else if(meanFlowType == "scatteredData") 
    {
      meanFlowFile_ = "";
      meanFlowType_ = MF_SCATTERED_DATA;
    }
    else if(meanFlowType == "analytical") 
    {
      meanFlowFile_ = "";
      meanFlowType_ = MF_ANALYTIC_EXP;
    }
    else
    {
      EXCEPTION("Mean flow type '" << meanFlowType
                << "' cannot be handled yet!");      
    }

    if(wvtNode->Has("evaluationGrid"))
    {
      // TODO: Hier noch die Regionen angeben, auf denen ausgwertet werden soll.
      evalGridId = wvtNode->Get("evaluationGrid")->Get("id")->As<std::string>();
    }
    else 
    {
      evalGridId = priModeId;
    }
    evalGridFile_ = inputNode->GetByVal("hdf5", "id", evalGridId)
                    ->Get("fileName")->As<std::string>();
    std::cout << "Evaluation grid file: " << evalGridFile_ << std::endl;    


    if(wvtNode->Has("output"))
    {
      PtrParamNode wvtOutputNode = wvtNode->Get("output");
      
      if(wvtOutputNode->Has("gridResults"))
      {
        std::string id;
        id = wvtOutputNode->Get("gridResults")->Get("id")->As<std::string>();
        outFile_ = outputNode->GetByVal("hdf5", "id", id)
                   ->Get("fileName")->As<std::string>();
        writeOutputFile_ = true;        

        std::cout << "Output file for grid results: " << outFile_ << std::endl;
      }

      // TODO: Hier noch abchecken, ob .csv oder .mat geschrieben werden soll
    }
    
    if(wvtNode->Has("integOrder")) {
      PtrParamNode integOrderNode = wvtNode->Get("integOrder");
      
      integOrder_ = integOrderNode->As<UInt>();      
    }

    if(wvtNode->Has("intPointsFile")) {
      PtrParamNode intPointsFileNode = wvtNode->Get("intPointsFile");
      
      std::string fn  = intPointsFileNode->As<std::string>();
      if(fn != "") 
      {
        intPointsFileName_ = fn;
      }
      writeIntPointsFile_ = true;
    }

    if(wvtNode->Has("sensorNodeName")) {
      PtrParamNode sensorNodeNameNode = wvtNode->Get("sensorNodeName");
      
      sensorNodeName_  = sensorNodeNameNode->As<std::string>();
    }

    wvtNode->GetValue("nodalResults", nodalResults_, ParamNode::PASS );

  }

  WVT::~WVT()
  {
    CloseIntPointsFile();
  }
  
  void WVT::PostProcess() 
  {
    PtrParamNode wvtNode = param_->Get("wvt");

    Initialize();    

    OpenIntPointsFile();

    std::cout << "Integration order: " << integOrder_ << std::endl;
    std::cout << "Sensor node name (cf. P in Fig. 3 [Hemp1]): "
              << sensorNodeName_ << std::endl;
    
    // Linear correction  factor for scaling  the actual mean  flow integrated
    // over the fluid domain to the desired mean flow. Valid only for straight pipe!
    Double correct = 0.0;


    mp_.reset(new MathParser());

    // check capabilities of input class
    IOType::iterator iit, iend;
    iit = inputs_.begin();
    iend = inputs_.end();    
    for( ; iit != iend; iit++) 
    {
      std::string readerCapsResult;
      if(!CheckReaderCapabilities(fileType.ToString(iit->first),
                                  iit->second,
                                  readerCapsResult))
      {
        EXCEPTION(readerCapsResult);
      }
    }

    // read in mesh of input1
    inputs_[PRIMARY_MODE]->InitModule();
    dim_ = inputs_[PRIMARY_MODE]->GetDim();
    Grid * ptGrid1 = new GridCFS(dim_, param_, info_);
    inputs_[PRIMARY_MODE]->ReadMesh(ptGrid1);
    ptGrid1->FinishInit();

    // Init dofNames array.
    switch(dim_) {
    case 2:
      dofNames_ = "x", "y";
      break;
    case 3:
      dofNames_ = "x", "y", "z";
      break;               
    } 

    numDofs_ = dofNames_.GetSize();
    
    // Read the XML parameter file from the HDF5 result file of the primary/drive mode.
    PtrParamNode pNode = GetParamNodeFromHDF5(inputs_[PRIMARY_MODE],
                                              "ParameterFile");
    
    ParamNodeList list;
    
    PtrParamNode multiSeqStep1;
    list = pNode->GetList("sequenceStep");
    if(list.GetSize() > 1) 
    {
      multiSeqStep1 = pNode->GetByVal("sequenceStep", "index", 1);
    }
    else
    {
      multiSeqStep1 = pNode->Get("sequenceStep");
    }

    bool ms1HasMech = multiSeqStep1->Get("pdeList")->Has("mechanic");
    bool ms1HasFlow = multiSeqStep1->Get("pdeList")->Has("fluidMech");
    bool ms1HasAcou = multiSeqStep1->Get("pdeList")->Has("acoustic");
    bool ms1HasCplLst = multiSeqStep1->Has("couplingList");
    
    PtrParamNode mechNode, fluidNode;
    
    StdVector<std::string> cplSurfList;
    
    if( ms1HasCplLst )
    {
      PtrParamNode cplNode = multiSeqStep1->Get("couplingList");
      std::cout << "We have a <couplingList> tag." << std::endl;

      if( ms1HasMech ) 
      {
        if( ms1HasFlow)
        {
          // For directly coupled WVT, mechanic and fluid analyses need to be
          // in the first MS step. There needs to be also a coupling section.
          // Normally you do not need to apply WVT to directly coupled system
          // since the phase difference between the sensors can be readily
          // obtained by comparing the phases of the mech. displacements.
          // But we want to cross-check our directly and iteratively coupled
          // results.
          
          if( cplNode->Get("direct")->Has("fluidMechDirect") ) 
          {  
            dirCoupled_ = true;
            
            mechNode = multiSeqStep1->Get("pdeList")->Get("mechanic");
            fluidNode = multiSeqStep1->Get("pdeList")->Get("fluidMech");         

            std::cout << "Coupling type: fluidMechDirect" << std::endl;
            
            list = cplNode->Get("direct")->Get("fluidMechDirect")
                   ->Get("surfRegionList")->GetChildren();
            for(UInt i=0; i<list.GetSize(); i++) {
              std::cout << "Coupling " << list[i]->GetName() << ": "
                        << list[i]->Get("name")->As<std::string>() << std::endl;
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

        if( ms1HasAcou )
        {
          if( cplNode->Get("direct")->Has("acouMechDirect") ) 
          {  
            dirCoupled_ = true;
            
            mechNode = multiSeqStep1->Get("pdeList")->Get("mechanic");
            fluidNode = multiSeqStep1->Get("pdeList")->Get("acoustic");
            
            std::cout << "Coupling type: acouMechDirect" << std::endl;

            // Since the acoustic particle velocity is the negative gradient
            // of the acoustic potential, we require acouPotential formulation.
            if( !fluidNode->Has("formulation") ) 
            {
              EXCEPTION("The acoustic PDE does not have a formulation " <<
                        "attribute!\nWVT requires acouPotential " <<
                        "formulation at the moment.");
            }

            // Now check if the acouPotential result is available in the HDF5.
            try 
            {
              fluidNode->Get("storeResults")->GetByVal("nodeResult",
                                                       "type",
                                                       "acouPotential");
            } catch( Exception& ex ) 
            {
              RETHROW_EXCEPTION(ex, "No nodal result 'acouPotential' " <<
                                "available in HDF5 file!");
            }            
            
            // Obtain list of coupling surfaces.
            list = cplNode->Get("direct")->Get("acouMechDirect")
                   ->Get("surfRegionList")->GetChildren();
            for(UInt i=0; i<list.GetSize(); i++) {
              std::cout << "Coupling " << list[i]->GetName() << ": "
                        << list[i]->Get("name")->As<std::string>() << std::endl;
              cplSurfList.Push_back(list[i]->Get("name")->As<std::string>());
            }
          }
          else
          {
            EXCEPTION("There exist mechanic and acoustic PDEs in the " <<
                      "first multi-sequence step,\nbut no coupling node." <<
                      "This case is not yet handled by WVT.");
          }
        }

        // Now read mech displacement at sensor position 1 of lateral mode.
        SimInputHDF5* hdf5Reader = dynamic_cast<SimInputHDF5*>(
          inputs_[PRIMARY_MODE].get()
          );
        StdVector<Complex> result;
        hdf5Reader->GetNamedNodeResult(sensorNodeName_, "mechDisplacement", result);
        
        u_p_ = result[0];
      }
    }
    else
    {
      // For iteratively coupled WVT, we want to have the mechanic analyis
      // in the first MS step and the fluid analysis in the second MS step 
      // in the future. Since this is not implemented at the moment, we
      // require fluid in the first MS and get the mechanic velocities at
      // sensor positions in lateral mode from somewhere else (namely the
      // wvtArgs.json file)
      
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
          std::cout << "Volume " << list[i]->GetName() << ": "
                    << list[i]->Get("name")->As<std::string>() << std::endl;
          cplSurfList.Push_back(list[i]->Get("name")->As<std::string>());
        }
        
      }
    }
    
    std::cout << "Direct coupled: " << dirCoupled_ << std::endl;
    
    // Obtain volume flow regions along with the according material names.
    std::map<std::string, std::string> regionMatMap;
    std::map<std::string, std::string>::const_iterator regMatIt, regMatEnd;
    
    list = fluidNode->Get("regionList")->GetChildren();
    ParamNodeList domainRegionList;
    domainRegionList = pNode->Get("domain")->Get("regionList")->GetChildren();
    
    for(UInt i=0; i<list.GetSize(); i++) {
      std::string regionName = list[i]->Get("name")->As<std::string>();
      
      std::cout << "Volume " << list[i]->GetName()
                << ": " << regionName << std::endl;
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
    
    // Read the XML material file from the HDF5 result file of the primary/drive mode.
    pNode = GetParamNodeFromHDF5(inputs_[PRIMARY_MODE], "MaterialFile");

    // Obtain densities and viscosities for volume flow regions.
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
      
      std::cout << "Density for region " << regionName << ": "
                << density << std::endl;
      std::cout << "Dynamic Viscosity for region " << regionName << ": "
                << viscosity << std::endl;
      regionDensityMap[regionName] = density;
      regionViscosityMap[regionName] = viscosity;
    }
    
    // read in mesh of input2
    inputs_[SECONDARY_MODE]->InitModule();
    Grid * ptGrid2 = new GridCFS(dim_, param_, info_);
    inputs_[SECONDARY_MODE]->ReadMesh(ptGrid2);
    ptGrid2->FinishInit();

    // Only read third grid if we really need it.
    Grid * ptGrid3 = NULL;
    // Generate input readers for mean flow data
    switch(meanFlowType_)
    {
    case MF_GRID_DATA:
      ptGrid3 = new GridCFS(dim_, param_, info_);

      inputs_[MEAN_FLOW] = GetInputReader( meanFlowFile_, param_, info_ );
      inputs_[MEAN_FLOW]->InitModule();
      // read in mesh of input3
      inputs_[MEAN_FLOW]->ReadMesh(ptGrid3);
      ptGrid3->FinishInit();
      break;
      
    case MF_ANALYTIC_EXP:
      ptGrid3 = new GridCFS(dim_, param_, info_);

      if(evalGridFile_ == priModeFile_) 
      {
        inputs_[MEAN_FLOW] = inputs_[PRIMARY_MODE];
      }
      else
      {
        inputs_[MEAN_FLOW] = GetInputReader( evalGridFile_, param_, info_ );
        inputs_[MEAN_FLOW]->InitModule();
      }
      
      // read in mesh of input3
      inputs_[MEAN_FLOW]->ReadMesh(ptGrid3);
      ptGrid3->FinishInit();

      std::cout << "Reading " << meanFlowDataType.ToString(meanFlowType_)
                << " from '" << meanFlowFile_ << "'"
                << "..." << std::endl;
      break;
      
    case MF_SCATTERED_DATA:
      PtrParamNode scatteredDataNode = wvtNode->Get("meanFlow")->Get("scatteredData");

      meanFlowCoefScattered_.reset(new CoefFunctionScatteredData<Double, 3>(
                                     scatteredDataNode)
        );
      break;
    }

    // obtain output writer
    shared_ptr<SimOutput> output;
    if( writeOutputFile_ ) {
      output = GetOutputWriter( outFile_, param_, info_ );
      output->Init(ptGrid1, false);
    }
    
    // obtain number of Sequence Steps and get analysis types
    std::map<UInt, BasePDE::AnalysisType> types;
    std::map<UInt, UInt> numSteps;
    inputs_[PRIMARY_MODE]->GetNumMultiSequenceSteps( types, numSteps, false );
    
    std::cout << "\nFound " << types.size() << " sequence step(s) in '"
              << inputs_[PRIMARY_MODE]->GetFileName() << "'\n";
    std::map<UInt, BasePDE::AnalysisType> types2;
    std::map<UInt, UInt> numSteps2;
    inputs_[SECONDARY_MODE]->GetNumMultiSequenceSteps( types2, numSteps2, false );
    std::cout << "\nFound " << types2.size() << " sequence step(s) in '"
              << inputs_[SECONDARY_MODE]->GetFileName() << "'\n";
    
    if(types.size() != types2.size()){
      std::cout << "'" << inputs_[PRIMARY_MODE]->GetFileName() << "' and '"
                << inputs_[SECONDARY_MODE]->GetFileName()
                << "' have different number of sequence steps!\n";
      exit(EXIT_FAILURE);
    }
    
    // Check for single sequence step.
    if(types.size() != 1) {
      std::cout << "WVT input files should have only one sequence step!\n";
      exit(EXIT_FAILURE);
    }
    
    // Check if primary and secondary mode input files are harmonic analyses.
    if( types[1] != BasePDE::HARMONIC ||
        types2[1] != BasePDE::HARMONIC ) {
      std::cout << "WVT is only available for harmonic analyses!\n";
      exit(EXIT_FAILURE);
    }
    
    UInt actMsStep = 1;
    std::cout << " Computing WVT for sequence step " << actMsStep
              << std::endl
              << "-------------------------------------\n\n";
    
    // get resulttypes
    StdVector<shared_ptr<ResultInfo> > infos, infos2, infos3, infos_mean_flow;
    inputs_[PRIMARY_MODE]->GetResultTypes( actMsStep, infos, false );
    inputs_[SECONDARY_MODE]->GetResultTypes( actMsStep, infos2, false );
    if(meanFlowType_ != MF_SCATTERED_DATA) 
    {
      inputs_[MEAN_FLOW]->GetResultTypes( actMsStep, infos3, false );
    }
    
    StdVector<shared_ptr<BaseResult> > inResults1, inResults2, inResults_mean_flow,
      outResults, outResults2, outResults3, outResults4, outResults5,
      outResults6, outResults7, outResults8;

    // stepnumbers, for which at least one result is defined
    std::map<UInt, Double> stepVals, stepVals2, stepVals_mean_flow;
    // contains the stepnumbers/-values in which the particular result is
    // defined in
    std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
    std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps2;
    std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps_mean_flow;
    
    std::map<std::string, bool> surfMap;
    
    if(meanFlowType_ == MF_GRID_DATA) 
    {
      // First determine number of results in current multi sequence step in mean flow file
      // since it can be different than in the lateral and coriolis mode files.    
      for( UInt iRes = 0; iRes < infos3.GetSize(); iRes++) {
        shared_ptr<ResultInfo> actRes3 = infos3[iRes];
        
        if(actRes3->resultType != FLUIDMECH_VELOCITY)
        {
          if(actRes3->resultType != MEAN_FLUIDMECH_VELOCITY) 
          {
          continue;
          }
        }
        
        // get stepvalues of mean flow file
        inputs_[MEAN_FLOW]->GetStepValues( actMsStep, actRes3,
                                           resultSteps_mean_flow[actRes3], false);
        stepVals_mean_flow.insert( resultSteps_mean_flow[actRes3].begin(),
                                   resultSteps_mean_flow[actRes3].end() );
        infos_mean_flow.Push_back(actRes3);
      }
    }
    
    if( infos.GetSize() > 0 ){
      std::cout << "Computing WVT from the following results:\n";
    }
    
    // iterate over all result types of input1
    for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {
      
      if(ms1HasCplLst && ms1HasAcou)
      {
        if(infos[iRes]->resultType != ACOU_POTENTIAL)
          continue;
      }
      else 
      {
        if(infos[iRes]->resultType != FLUIDMECH_VELOCITY)
          continue;
      }
      
      
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
      
      
      // iterate over all regions of primary and secondary modes
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
      
      
      StdVector< shared_ptr<EntityList> > meanFlowRegions;
      if(meanFlowType_ != MF_SCATTERED_DATA) 
      {      
        // iterate over all regions of mean flow grid
        inputs_[MEAN_FLOW]->GetResultEntities( actMsStep, infos[iRes],
                                               resRegions, false );
      
        for( UInt iRegion = 0, nReg = resRegions.GetSize(); iRegion < nReg; iRegion++ ) {
          std::string resRegionName = resRegions[iRegion]->GetName();
          
          //           if(resRegionName == "fluid__________________________________________________________________________") 
          //           {
          meanFlowRegions.Push_back(resRegions[iRegion]);
          //           }
        }           
      }
      
      
      
      for( UInt iRegion = 0, nReg = regions.GetSize(); iRegion < nReg; iRegion++ ) {
        // generate new result object and add it to output writer
        shared_ptr<BaseResult >
          inResult1, inResult2, inResult3, outResult,
          outResult2, outResult3, outResult4, outResult5,
          outResult6, outResult7, outResult8;
        
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
        
        std::cout << "Defining results for region " << regionName << "..."
                  << std::endl;
        //           EXCEPTION("HALT");
        
        inResult1 = shared_ptr<BaseResult>( new Result<Complex>() );
        inResult2 = shared_ptr<BaseResult>( new Result<Complex>() );
        outResult = shared_ptr<BaseResult>( new Result<Complex>() );
        outResult2 = shared_ptr<BaseResult>( new Result<Complex>() );
        outResult3 = shared_ptr<BaseResult>( new Result<Complex>() );
        outResult4 = shared_ptr<BaseResult>( new Result<Complex>() );
        outResult5 = shared_ptr<BaseResult>( new Result<Complex>() );
        outResult6 = shared_ptr<BaseResult>( new Result<Complex>() );
        outResult7 = shared_ptr<BaseResult>( new Result<Complex>() );
        outResult8 = shared_ptr<BaseResult>( new Result<Complex>() );
        
        inResult1->SetEntityList( regions[iRegion] );
        inResult2->SetEntityList( regions[iRegion] );

        // Since derivatives need to be taken the output entities need to be elements
        Enum<RegionIdType> regionEnum = ptGrid1->GetRegion();
        RegionIdType regionId = regionEnum.Parse(regions[iRegion]->GetName());
        ResultInfo::EntityUnknownType definedOn;
        definedOn = nodalResults_ ? ResultInfo::NODE : ResultInfo::ELEMENT;

        ElemList* outElemEntities = new ElemList(ptGrid1);
        outElemEntities->SetRegion(regionId);
        shared_ptr<EntityList> outElems(outElemEntities);             
        outResult->SetEntityList( outElems );
        outResult2->SetEntityList( outElems );
        outResult3->SetEntityList( outElems );
        outResult4->SetEntityList( outElems );
        outResult5->SetEntityList( outElems );
        outResult6->SetEntityList( outElems );
        outResult7->SetEntityList( outElems );
        outResult8->SetEntityList( outElems );
        
        inResult1->SetResultInfo( infos[iRes] );
        inResult2->SetResultInfo( infos[iRes] );
        
        shared_ptr<ResultInfo> outInfo(new ResultInfo());
        outInfo->resultType = FLUIDMECH_WVT;
        outInfo->resultName = SolutionTypeEnum.ToString(outInfo->resultType);
        outInfo->dofNames = dofNames_;
        outInfo->unit = MapSolTypeToUnit(outInfo->resultType);
        outInfo->entryType = ResultInfo::VECTOR;
        outInfo->definedOn = definedOn;
        
        outResult->SetResultInfo( outInfo );
        
        shared_ptr<ResultInfo> outInfo2(new ResultInfo());
        outInfo2->resultType = MEAN_FLUIDMECH_VELOCITY;
        outInfo2->resultName = SolutionTypeEnum.ToString(outInfo2->resultType);
        outInfo2->dofNames = dofNames_;
        outInfo2->unit = MapSolTypeToUnit(outInfo2->resultType);
        outInfo2->entryType = ResultInfo::VECTOR;
        outInfo2->definedOn = definedOn;
        
        outResult2->SetResultInfo( outInfo2 );
        
        shared_ptr<ResultInfo> outInfo3(new ResultInfo());
        outInfo3->resultType = FLUIDMECH_WVT_DENSITY;
        outInfo3->resultName = SolutionTypeEnum.ToString(outInfo3->resultType);
        outInfo3->dofNames = "";
        outInfo3->unit = MapSolTypeToUnit(outInfo3->resultType);
        outInfo3->entryType = ResultInfo::SCALAR;
        outInfo3->definedOn = definedOn;
        
        outResult3->SetResultInfo( outInfo3 );
        
        shared_ptr<ResultInfo> outInfo4(new ResultInfo());
        outInfo4->resultType = FLUIDMECH_WVT_PHI;
        outInfo4->resultName = SolutionTypeEnum.ToString(outInfo4->resultType);
        outInfo4->dofNames = dofNames_;
        outInfo4->unit = MapSolTypeToUnit(outInfo4->resultType);
        outInfo4->entryType = ResultInfo::VECTOR;
        outInfo4->definedOn = definedOn;
        
        outResult4->SetResultInfo( outInfo4 );
        
        shared_ptr<ResultInfo> outInfo5(new ResultInfo());
        outInfo5->resultType = FLUIDMECH_WVT_DENSITY_PHI;
        outInfo5->resultName = SolutionTypeEnum.ToString(outInfo5->resultType);
        outInfo5->dofNames = "";
        outInfo5->unit = MapSolTypeToUnit(outInfo5->resultType);
        outInfo5->entryType = ResultInfo::SCALAR;
        outInfo5->definedOn = definedOn;
        
        outResult5->SetResultInfo( outInfo5 );

        shared_ptr<ResultInfo> outInfo6(new ResultInfo());
        outInfo6->resultType = FLUIDMECH_WVT_U1;
        outInfo6->resultName = SolutionTypeEnum.ToString(outInfo6->resultType);
        outInfo6->dofNames = dofNames_;
        outInfo6->unit = MapSolTypeToUnit(outInfo6->resultType);
        outInfo6->entryType = ResultInfo::VECTOR;
        outInfo6->definedOn = definedOn;
        
        outResult6->SetResultInfo( outInfo6 );

        shared_ptr<ResultInfo> outInfo7(new ResultInfo());
        outInfo7->resultType = FLUIDMECH_WVT_U2;
        outInfo7->resultName = SolutionTypeEnum.ToString(outInfo7->resultType);
        outInfo7->dofNames = dofNames_;
        outInfo7->unit = MapSolTypeToUnit(outInfo7->resultType);
        outInfo7->entryType = ResultInfo::VECTOR;
        outInfo7->definedOn = definedOn;
        
        outResult7->SetResultInfo( outInfo7 );

        shared_ptr<ResultInfo> outInfo8(new ResultInfo());
        outInfo8->resultType = FLUIDMECH_WVT_F;
        outInfo8->resultName = SolutionTypeEnum.ToString(outInfo8->resultType);
        outInfo8->dofNames = dofNames_;
        outInfo8->unit = MapSolTypeToUnit(outInfo8->resultType);
        outInfo8->entryType = ResultInfo::VECTOR;
        outInfo8->definedOn = definedOn;
        
        outResult8->SetResultInfo( outInfo8 );
        
        inResults1.Push_back( inResult1 );
        inResults2.Push_back( inResult2 );
        outResults.Push_back( outResult );
        outResults2.Push_back( outResult2 );
        outResults3.Push_back( outResult3 );
        outResults4.Push_back( outResult4 );
        outResults5.Push_back( outResult5 );
        outResults6.Push_back( outResult6 );
        outResults7.Push_back( outResult7 );
        outResults8.Push_back( outResult8 );

        if(meanFlowType_ == MF_GRID_DATA) 
        {
          inResult3 = shared_ptr<BaseResult>( new Result<Double>() );          
          inResult3->SetEntityList( meanFlowRegions[iRegion] );

          inResult3->SetResultInfo( infos_mean_flow[0] );

          inResults_mean_flow.Push_back( inResult3 );
        }


        if( output) {
          outResult->GetResultInfo()->complexFormat = REAL_IMAG;
          // CAUTION: begin, inc, end are hardcoded and noch checked for each result
          output->RegisterResult( outResult, 1, 1,
                                  resultSteps[actRes].size(),
                                  false );
          
          outResult2->GetResultInfo()->complexFormat = REAL_IMAG;
          output->RegisterResult( outResult2, 1, 1,
                                  resultSteps[actRes].size(),
                                  false );
          
          outResult3->GetResultInfo()->complexFormat = REAL_IMAG;
          output->RegisterResult( outResult3, 1, 1,
                                  resultSteps[actRes].size(),
                                  false );
          
          outResult4->GetResultInfo()->complexFormat = REAL_IMAG;
          output->RegisterResult( outResult4, 1, 1,
                                  resultSteps[actRes].size(),
                                  false );

          outResult5->GetResultInfo()->complexFormat = REAL_IMAG;
          output->RegisterResult( outResult5, 1, 1,
                                  resultSteps[actRes].size(),
                                  false );

          outResult6->GetResultInfo()->complexFormat = REAL_IMAG;
          output->RegisterResult( outResult6, 1, 1,
                                  resultSteps[actRes].size(),
                                  false );

          outResult7->GetResultInfo()->complexFormat = REAL_IMAG;
          output->RegisterResult( outResult7, 1, 1,
                                  resultSteps[actRes].size(),
                                  false );

          outResult8->GetResultInfo()->complexFormat = REAL_IMAG;
          output->RegisterResult( outResult8, 1, 1,
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
        
        // We are only interested in the volume WVT from Hemp1 at the moment.
        if(isSurf)
          continue;
        
        
        std::cout << "\n\t-- Computing " << (isSurf ? "surface" : "volume")
                  << " weight vector for " 
                  << inResults1[iRes]->GetResultInfo()->resultName << " on "
                  << regionName << " --\n";
        
        // check if current result is defined within this step
        if( resultSteps[inResults1[iRes]->GetResultInfo()].find(actStepNum)
            == resultSteps[inResults1[iRes]->GetResultInfo()].end() ) {
          continue;
        }
        
        // obtain both result objects for current step
        std::set<std::string> volRegionNames;
        std::set<std::string> meanFlowVolRegionNames;
        inputs_[PRIMARY_MODE]->GetResult( actMsStep, actStepNum,
                                          inResults1[iRes], false );
        inputs_[SECONDARY_MODE]->GetResult( actMsStep, actStepNum,
                                            inResults2[iRes], false );
        if(meanFlowType_ == MF_GRID_DATA) 
        {
          inputs_[MEAN_FLOW]->GetResult( actMsStep, actStepNum,
                                         inResults_mean_flow[iRes],
                                         false );
          meanFlowVolRegionNames.insert( regionName );
        }
        
        std::set<RegionIdType> volRegions;
        // volRegionNames.insert(inResults1[iRes]->GetEntityList()->GetName());
        std::map<std::string, Double>::iterator it, end;
        it = regionDensityMap.begin();
        end = regionDensityMap.end();
        for( ; it != end; it++ ) {
          const std::string& regionName = it->first;
          std::cout << "Density of region " << it->first << ": "
                    << it->second << std::endl;
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
          Vector<Complex> & outVec6 =
            dynamic_cast<Result<Complex>& >(*outResults6[iRes]).GetVector();
          Vector<Complex> & outVec7 =
            dynamic_cast<Result<Complex>& >(*outResults7[iRes]).GetVector();
          Vector<Complex> & outVec8 =
            dynamic_cast<Result<Complex>& >(*outResults8[iRes]).GetVector();
          
          std::string entityListName = outResults[iRes]->GetEntityList()->GetName();
          RegionIdType regionId = ptGrid1->GetRegion().Parse(entityListName);

          shared_ptr<EntityList> entityList;
          NodeList* nl = new NodeList(ptGrid1);
          nl->SetNodesOfRegion(regionId);
          entityList.reset(nl);  
          const StdVector<UInt> & nodes = nl->GetNodes();

          UInt numEntities = nodalResults_ ? nodes.GetSize() : 
                             ptGrid1->GetNumElems(regionId);

          outVec.Resize( numEntities * numDofs_ );  outVec.Init();
          outVec2.Resize( numEntities * numDofs_ ); outVec2.Init();
          outVec3.Resize( numEntities );            outVec3.Init();
          outVec4.Resize( numEntities * numDofs_ ); outVec4.Init();
          outVec5.Resize( numEntities );            outVec5.Init();
          outVec6.Resize( numEntities * numDofs_ ); outVec6.Init();
          outVec7.Resize( numEntities * numDofs_ ); outVec7.Init();
          outVec8.Resize( numEntities * numDofs_ ); outVec8.Init();

          Vector<Double> nodalWeights;
          if(nodalResults_) 
          {
            nodalWeights.Resize(numEntities); nodalWeights.Init();
          }

          shared_ptr<BaseFeFunction> u1FeFunc =
            inputs_[PRIMARY_MODE]->GetFeFunction<Complex>( actMsStep,
                                                           actStepNum,
                                                           inResults1[iRes]->GetResultInfo()->resultType,
                                                           volRegionNames);
          shared_ptr<BaseFeFunction> u2FeFunc = 
            inputs_[SECONDARY_MODE]->GetFeFunction<Complex>( actMsStep,
                                                             actStepNum,
                                                             inResults2[iRes]->GetResultInfo()->resultType,
                                                             volRegionNames);
          
          
          shared_ptr<BaseFeFunction> VFeFunc;
          if(meanFlowType_ == MF_GRID_DATA) 
          {
            VFeFunc = 
            inputs_[MEAN_FLOW]->GetFeFunction<Double>( actMsStep,
                                                       actStepNum,
                                                       inResults_mean_flow[iRes]->GetResultInfo()->resultType,
                                                       meanFlowVolRegionNames );
          }
          
          
          shared_ptr<BaseBOperator> gradOp;
          shared_ptr<BaseBOperator> curlOp;
          switch(dim_) {
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
          
          if(dirCoupled_) 
          {
            u_p_ *= Complex(0,1);
            u_p_ *= 2*M_PI*freq;
          }
          else
          {
            if(wvtNode->Has("u_p")) {
              PtrParamNode upNode = wvtNode->Get("u_p");
              
              if(upNode->Has("real")) {
                u_p_real = upNode->Get("real")->As<Double>();
              }
              if(upNode->Has("imag")) {
                u_p_imag = upNode->Get("imag")->As<Double>();
              }
              u_p_ = Complex(u_p_real, u_p_imag);
              u_p_ *= Complex(0,1);
              u_p_ *= 2*M_PI*freq;
            }
            if(wvtNode->Has("v_p")) {
              PtrParamNode vpNode = wvtNode->Get("v_p");
              
              if(vpNode->Has("real")) {
                u_p_real = vpNode->Get("real")->As<Double>();
              }
              if(vpNode->Has("imag")) {
                u_p_imag = vpNode->Get("imag")->As<Double>();
              }
              u_p_ = Complex(u_p_real, u_p_imag);
            }
          }
          
          StdVector<std::string> realVal(3);
          realVal[0] = "0";
          realVal[1] = "0";
          realVal[2] = "1";
          StdVector<std::string> imagVal(3);
          realVal[0] = "0";
          realVal[1] = "0";
          realVal[2] = "0";
          Double meanVel = 1;
          // FeFunction pointer for analytical V-field.
          FeFunction<Double>* feFctVPtr = NULL;
          PtrCoefFct coefFctV;
          shared_ptr<DefaultCoordSystem> coordSys( new DefaultCoordSystem(ptGrid1) );
          shared_ptr<BaseFeFunction> meanFlowFeFct;

          if(meanFlowType_ == MF_ANALYTIC_EXP) 
          {          
            PtrParamNode meanFlowNode = wvtNode->Get("meanFlow");
            if(meanFlowNode->Has("analytical")) {
              PtrParamNode vNode = meanFlowNode->Get("analytical");
              
              if(vNode->Has("meanVel")) {
                meanVel = vNode->Get("meanVel")->As<Double>();
              }
              
              if(vNode->Has("profile")) {
                PtrParamNode profileNode = vNode->Get("profile");
                
                // realVal[2] = "1-(sqrt(x^2+y^2)/13.15e-3)^10";
                realVal[2] = profileNode->As<std::string>();
                
                coefFctV = CoefFunction::Generate(mp_.get(), Global::REAL, realVal);
                coefFctV->SetCoordinateSystem(coordSys.get());
                
                feFctVPtr = new FeFunction<Double>(mp_.get());
                
              }
            }
          }
          
          if(meanFlowType_ == MF_GRID_DATA) 
          {          
            // If we have an analytical mean flow field, we use it. Otherwise we use
            // the one read from the mean flow grid file.
            if(!feFctVPtr) 
            {
              meanFlowFeFct = VFeFunc;
            }
          }
          
          shared_ptr<BaseFeFunction> priFeFct = u1FeFunc;
          shared_ptr<BaseFeFunction> secFeFct = u2FeFunc;
          
          if(ms1HasCplLst && ms1HasAcou) 
          {
            AcouPot2AcouVelFeFct(u1FeFunc, gradOp, priFeFct);
            AcouPot2AcouVelFeFct(u2FeFunc, gradOp, secFeFct);
          }          

          if(wvtNode->Has("u1")) {
            PtrParamNode u1Node = wvtNode->Get("u1");
            
            realVal[0] = realVal[1] = realVal[2] = "0.0";
            imagVal[0] = imagVal[1] = imagVal[2] = "0.0";
            if(u1Node->Has("x")) {
              PtrParamNode u1xNode = u1Node->Get("x");
              
              imagVal[0] = u1xNode->As<std::string>();
            }
            
            PtrCoefFct coefu1 = CoefFunction::Generate(mp_.get(), Global::COMPLEX, realVal, imagVal);
            shared_ptr<DefaultCoordSystem> coordSys( new DefaultCoordSystem(ptGrid1) );
            coefu1->SetCoordinateSystem(coordSys.get());
               
            
            FeFunction<Complex>* feFctPtr = new FeFunction<Complex>(mp_.get());
            shared_ptr<BaseFeFunction> feFct(feFctPtr);
            feFct->SetFeSpace( u1FeFunc->GetFeSpace() );
            feFct->AddEntityList( inResults1[iRes]->GetEntityList() );
            feFct->SetGrid(ptGrid1);
            feFct->SetResultInfo(inResults1[iRes]->GetResultInfo());
            feFct->AddExternalDataSource( coefu1,
                                          inResults1[iRes]->GetEntityList() ); // <--- Hier muss unbedingt die zugehoerige Volumenregion zur Oberflaechenregion hinein!
            feFct->SetFctId(PSEUDO_FCT_ID);
            feFct->Finalize();
            feFct->ApplyExternalData();

            priFeFct = feFct;
          }
          
          if(wvtNode->Has("u2")) {
            PtrParamNode u2Node = wvtNode->Get("u2");
            
            realVal[0] = realVal[1] = realVal[2] = "0.0";
            imagVal[0] = imagVal[1] = imagVal[2] = "0.0";
            if(u2Node->Has("x")) {
              PtrParamNode u2xNode = u2Node->Get("x");
              
              imagVal[0] = u2xNode->As<std::string>();
            }
            
            PtrCoefFct coefu2 = CoefFunction::Generate(mp_.get(), Global::COMPLEX, realVal, imagVal);
            shared_ptr<DefaultCoordSystem> coordSys( new DefaultCoordSystem(ptGrid1) );
            coefu2->SetCoordinateSystem(coordSys.get());
            
            
            FeFunction<Complex>* feFctPtr = new FeFunction<Complex>(mp_.get());
            shared_ptr<BaseFeFunction> feFct(feFctPtr);
            feFct->SetFeSpace( u2FeFunc->GetFeSpace() );
            feFct->AddEntityList( inResults1[iRes]->GetEntityList() );
            feFct->SetGrid(ptGrid1);
            feFct->SetResultInfo(inResults1[iRes]->GetResultInfo());
            feFct->AddExternalDataSource( coefu2,
                                          inResults1[iRes]->GetEntityList() ); // <--- Hier muss unbedingt die zugehoerige Volumenregion zur Oberflaechenregion hinein!
            feFct->SetFctId(PSEUDO_FCT_ID);
            feFct->Finalize();
            feFct->ApplyExternalData();
            secFeFct = feFct;
          }
          
          
          
          Double rho_0 = regionDensityMap[inResults1[iRes]->GetEntityList()->GetName()];
          
          // fluid________________________ region
          
          eIt = outResults[iRes]->GetEntityList()->GetIterator();

          for(UInt elIdx=0 ; !eIt.IsEnd(); eIt++, elIdx++ ) 
          {
            const Elem* el1 = eIt.GetElem();
            const Elem* el2 = ptGrid2->GetElem(el1->elemNum);
            const Elem* el3 = NULL;
            if(ptGrid3) { 
              el3 = ptGrid3->GetElem(el1->elemNum);
            }
            
            
            // Obtain FE element from feSpace and integration scheme
            IntegOrder order;
            IntScheme::IntegMethod method;
            BaseFE* ptFe1 = priFeFct->GetFeSpace()->GetFe( eIt, method, order );
            BaseFE* ptFe2 = secFeFct->GetFeSpace()->GetFe( el2->elemNum );
            BaseFE* ptFe3 = NULL;
            
            // Get integration points
            shared_ptr<IntScheme> intScheme = priFeFct->GetFeSpace()->GetIntScheme();
            
            StdVector<LocPoint> intPoints;
            StdVector<Double> weights;
            order.SetIsoOrder(integOrder_);
            
            intScheme->GetIntPoints( Elem::GetShapeType(el1->type), method, order,
                                     intPoints, weights );
            
            
            //               std::cout << "Integration order " << order.GetIsoOrder() << " num. int points: " << intPoints.GetSize() << " weight: " << weights[0] << std::endl;
            
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
            shared_ptr<ElemShapeMap> esm3;
            if(ptGrid3) 
            {
              esm3 = ptGrid3->GetElemShapeMap( el3, true );
            }
            
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
              Vector<Double> p(dim_);
              LocPoint lp3;
              
              if( isSurf ) 
              {
                lpm1Surf.SetWithSurface( lp, esm1, volRegions, 0.0 );
                lpm1 = lpm1Surf.lpmVol;
                lpm2Surf.SetWithSurface( lp, esm2, volRegions, 0.0 );
                lpm2 = lpm2Surf.lpmVol;
                
                el1 = lpm1->ptEl;
                ptFe1 = priFeFct->GetFeSpace()->GetFe( lpm1->ptEl->elemNum );
                
                el2 = lpm2->ptEl;
                ptFe2 = secFeFct->GetFeSpace()->GetFe( lpm2->ptEl->elemNum );

                if(ptGrid3) 
                {
                  lpm3Surf.SetWithSurface( lp, esm3, volRegions, 0.0 );
                  lpm3 = lpm3Surf.lpmVol;
                  el3 = lpm3->ptEl;
                  ptFe3 = meanFlowFeFct->GetFeSpace()->GetFe( lpm3->ptEl->elemNum );
                }                
              } else 
              {
                lpm1.reset(new LocPointMapped());
                lpm1->Set( lp, esm1, 0.0 );
                lpm2.reset(new LocPointMapped());
                lpm2->Set( lp, esm2, 0.0 );

                if(ptGrid3) 
                {
                  lpm3.reset(new LocPointMapped());
                  /*
                    esm1->Local2Global(p, lp);
                    
                    Vector<Double> locPoint3;
                    try 
                    {
                    esm3->Global2Local(locPoint3, p);
                    } catch (Exception& ex) 
                    {
                    locPoint3 = lp.coord;
                    locPoint3.Init();
                    }                   
                    
                    LocPoint lp3_2(locPoint3);
                    lp3 = lp3_2;
                  */
                  
                  LocPoint lp3 = Elem::GetShape( Elem::GetShapeType( el3->type) ).midPointCoord;
                  
                  lpm3->Set( lp3, esm3, 0.0 );
                  
                  if(meanFlowType_ == MF_GRID_DATA) 
                  {          
                    ptFe3 = meanFlowFeFct->GetFeSpace()->GetFe( lpm3->ptEl->elemNum );                  
                  }
                  
                }
              }
              
              
              if(writeIntPointsFile_) 
              {
                esm1->Local2Global(p, lp);
                intPointsFile_ << p[0] << " " << p[1] << " " << p[2] << std::endl;
                numIntPoints_++;
              }
              
              // Matrix<Double> bMat1, bMat2, bMat3;
              Vector<Complex> eSol1, eSol2;
              Vector<Double> eSol3;
              StdVector< Vector<Complex> > eSol1Comp(numDofs_), eSol2Comp(numDofs_);
              StdVector< Vector<Double> > eSol3Comp(numDofs_);
              Vector<Complex> eSol2_x, eSol2_y, esol2_z;
              Vector<Complex> u1, u2;
              Vector<Double> V;
              Vector<Complex> W;
              Vector<Complex> f(3);
              Vector<Complex> X;
              
              // gradOp->CalcOpMat( bMat1, *lpm1, ptFe1 );
              // gradOp->CalcOpMat( bMat2, *lpm2, ptFe2 );
              //curlOp->CalcOpMat( bMat3, *lpm1, ptFe1 );
              // gradOp->CalcOpMat( bMat3, *lpm1, ptFe1 );
              
              // Get primary and secondary perturbed velocity fields.
              priFeFct->GetVector( u1, *lpm1 );
              secFeFct->GetVector( u2, *lpm2 );

              // Get mean velocity V either from file or from analytical expression.
              switch(meanFlowType_)
              {
              case MF_GRID_DATA:
                meanFlowFeFct->GetVector(V, *lpm3);
                break;
                
              case MF_ANALYTIC_EXP:
                coefFctV->GetVector(V, *lpm3);
                break;
                
              case MF_SCATTERED_DATA:
                meanFlowCoefScattered_->GetVector(V, *lpm1);
                break;
              }
                 
              dynamic_cast<FeFunction<Complex>*>(priFeFct.get())->GetElemSolution( eSol1, el1 );
              dynamic_cast<FeFunction<Complex>*>(secFeFct.get())->GetElemSolution( eSol2, el2 );

              if(ptGrid3) 
              {
                  if(meanFlowType_ == MF_GRID_DATA) 
                  {          
                    dynamic_cast<FeFunction<Double>*>(meanFlowFeFct.get())->GetElemSolution( eSol3, el3 );
                  }
              }
              
              StdVector< Vector<Complex> > u1Derivs(numDofs_), u2Derivs(numDofs_);
              StdVector< Vector<Double> > VDerivs(numDofs_);
              Vector<Double> curlV;
              StdVector< Vector<Complex> > u2Tau(numDofs_);
              Double viscosity = regionViscosityMap[ptGrid1->regionData[el1->regionId].name];
                 
              for(UInt dof=0; dof < numDofs_; dof++) 
              {
                UInt eN=eSol1.GetSize()/numDofs_;
                   
                eSol1Comp[dof].Resize(eN);
                eSol2Comp[dof].Resize(eN);
                if(ptGrid3) 
                {
                  eSol3Comp[dof].Resize(eN);
                }
                
                for(UInt eI=0; eI < eN; eI++) 
                {
                  eSol1Comp[dof][eI] = eSol1[eI*numDofs_+dof];
                  eSol2Comp[dof][eI] = eSol2[eI*numDofs_+dof];
                  if(ptGrid3) 
                  {
                    if(meanFlowType_ == MF_GRID_DATA) 
                    {          
                      eSol3Comp[dof][eI] = eSol3[eI*numDofs_+dof];
                    }
                  }
                }
                   
                gradOp->ApplyOp( u1Derivs[dof], *lpm1, ptFe1, eSol1Comp[dof] );
                gradOp->ApplyOp( u2Derivs[dof], *lpm2, ptFe2, eSol2Comp[dof] );
                // gradOp->ApplyOp( VDerivs[dof], *lpm3, ptFe3, eSol3Comp[dof] );
                VDerivs[dof].Resize(3);
                VDerivs[dof].Init();
              }
                 
              if(ptGrid3) 
              {
                if(meanFlowType_ == MF_GRID_DATA) 
                {
                  curlOp->ApplyOp( curlV, *lpm3, ptFe3, eSol3 );
                }
              }
              
              //                std::cout << "elemNum " << el1->elemNum << " numDofs_ " << numDofs_ << " dim " << dim << std::endl;
              Vector<Complex> Xl(numDofs_), Xr(numDofs_);
              Vector<Double> n(numDofs_);

              if(isSurf) 
              {
                n = -lpm3Surf.normal;
              } else 
              {
                n.Init();
              }                 

              Xl.Resize(numDofs_);
              Xr.Resize(numDofs_);
              X.Resize(numDofs_);
              Xl.Init();
              Xr.Init();
              X.Init();

           #if 0
              Vector<Double> t(numDofs_);
              t[0] = n[1];
              t[1] = -n[0];
              t[2] = 0;
              t *= 1/t.NormL2();
                 
              curlV = 2*radius/(R*R) * t;
           #endif

              for(UInt dof=0; dof < numDofs_; dof++) 
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
                   
                Complex complexFreq(0, 2*M_PI*freq);
                   
                X[0] =  Xl[1] * Xr[2] - Xl[2] * Xr[1];
                X[1] = -Xl[0] * Xr[2] + Xl[2] * Xr[0];
                X[2] =  Xl[0] * Xr[1] - Xl[1] * Xr[0];

                X *= viscosity / complexFreq;
              }                 

              // std::cout << "elemNum " << el3->elemNum << " curl(V) " << curlV[2] << " V " << V << std::endl;
              // std::cout << "eSol3 " << eSol3 << V << std::endl;

              W.Resize(u1.GetSize());
              W.Init();
                 
              for( UInt i=0; i < numDofs_; i++ )
              {
                for( UInt j=0; j < dim_; j++ )
                {
                     
                  if(isSurf) 
                  {
                    // W[i] += lpm3Surf.normal[i];
                    // W[i] += curlV[2][i];
                    W[i] = V[i];
                  }
                  else 
                  {
                    // Compute components of weight vector as given in Hemp1994 (15).
                    // Just take 2*C1 also with WVT!
                    W[i] += -rho_0 * 2 * (u2[j] * u1Derivs[j][i]);
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
                   
                for( UInt i=0; i < numDofs_; i++ )
                {
                  if(!nodalResults_) 
                  {
                    // Prepare weight vector
                    outVec[elIdx*numDofs_+i] += X[i] * volume;
                    
                    // Prepare mean velocity
                    outVec2[elIdx*numDofs_+i] += Complex(V[i], 0) * volume;
                    
                    // Prepare weight vector density result
                    outVec3[elIdx] += X[i]*curlV[i] * volume;
                    
                    
                    // Prepare weight vector W_Phi as given in Hemp1994 (19).
                    outVec4[elIdx*numDofs_+i] = curlV[i] * volume;
                    
                    // Prepare  weight  vector density  result  as inner  product
                    // between V and W_Phi as given in Hemp1994 (18).
                    //                     outVec5[elIdx] += W[i]/u_p * V[i];
                    
                    outVec5[elIdx] += X[i]/u_p_ * curlV[i] * volume;
                  }
                }

              }
              else 
              {
                  Vector<Complex> out(numDofs_);  out.Init();
                  Vector<Complex> out2(numDofs_); out2.Init();
                  Complex out3 = Complex(0,0);
                  Vector<Complex> out4(numDofs_); out4.Init();
                  Complex out5 = Complex(0,0);
                  Vector<Complex> out6(numDofs_); out6.Init();
                  Vector<Complex> out7(numDofs_); out7.Init();
                  Vector<Complex> out8(numDofs_); out8.Init();

                  for( UInt i=0; i < numDofs_; i++ )
                  {
                    // Prepare weight vector
                    out[i] += W[i] * volume;

                    durchfluss += V[i] * volume;
                      
                    // Prepare mean velocity output vector
                    out2[i] += Complex(V[i], 0) * volume;
                      
                    // Prepare weight vector density result
                    Complex wvtDensity = W[i]*V[i]*volume;
                    out3 += wvtDensity;
                      
                    // Compute u_p_prime as given in Hemp1994 (13).
                    u_p_prime += wvtDensity;
                      
                    f[i] = - rho_0 * (
                      u1Derivs[i][0] * V[0] + u1Derivs[i][1] * V[1] + u1Derivs[i][2] * V[2] // + // C1
                      // + VDerivs[i][0] * u1[0] + VDerivs[i][1] * u1[1] + VDerivs[i][2] * u1[2] // C2
                      );
                      
                    // Prepare weight vector W_Phi as given in Hemp1994 (19).
                    out4[i] += W[i]/u_p_ * volume;
                      
                    // Prepare  weight  vector density  result  as inner  product
                    // between V and W_Phi as given in Hemp1994 (18).
                    Complex wvtDensityPhi(0,0);
                      
                    wvtDensityPhi = W[i]/u_p_ * V[i] * volume;
                      
                    // Compute "weight vector density" result from force and u2
                    // Just take 2*C1 also with WVT!
                    // wvtDensityPhi = 2*(f[i] * u2[i])/u_p_ * volume;
                      
                    out5 += wvtDensityPhi;
                      
                    // Prepare u1 velocity output vector
                    out6[i] += u1[i] * volume;
                      
                    // Prepare u2 velocity output vector
                    out7[i] += u2[i] * volume;
                      
                    // Prepare force output vector
                    out8[i] += f[i] * volume;
                      
                    // Compute deltaPhi as given in Hemp1994 (18) and (19).
                    deltaPhi += wvtDensityPhi.imag();
                  }

                if(!nodalResults_) 
                {
                  // Assemble values to element vectors                  
                  for( UInt i=0; i < numDofs_; i++ )
                  {
                    outVec[elIdx*numDofs_+i] += out[i];
                    outVec2[elIdx*numDofs_+i] += out2[i];
                    outVec3[elIdx] += out3;
                    outVec4[elIdx*numDofs_+i] += out4[i];
                    outVec5[elIdx] += out5;
                    outVec6[elIdx*numDofs_+i] += out6[i];
                    outVec7[elIdx*numDofs_+i] += out7[i];
                    outVec8[elIdx*numDofs_+i] += out8[i];
                  }
                }
                else 
                {
                  // Assemble values to nodal vectors
                  for(UInt nIdx=0, numNodes=el1->connect.GetSize(); nIdx < numNodes; nIdx++)
                  {
                    StdVector<UInt>::const_iterator nodeIt = std::find(nodes.Begin(), nodes.End(), el1->connect[nIdx]);
                    UInt nodeIndex = std::distance(nodes.Begin(), nodeIt);

                    nodalWeights[nodeIndex] += volume;

                    for( UInt i=0; i < numDofs_; i++ )
                    {
                      outVec[nodeIndex*numDofs_+i] += out[i];
                      outVec2[nodeIndex*numDofs_+i] += out2[i];
                      outVec3[nodeIndex] += out3;
                      outVec4[nodeIndex*numDofs_+i] += out4[i];
                      outVec5[nodeIndex] += out5;
                      outVec6[nodeIndex*numDofs_+i] += out6[i];
                      outVec7[nodeIndex*numDofs_+i] += out7[i];
                      outVec8[nodeIndex*numDofs_+i] += out8[i];
                    }
                  }
                }
              }
            } // loop over integration points

            // Normalize element results by element volume
            if(!nodalResults_) 
            {
              for( UInt i=0; i < numDofs_; i++ )
              {
                outVec[elIdx*numDofs_+i] /= elementVolume;
                outVec2[elIdx*numDofs_+i] /= elementVolume;
                outVec4[elIdx*numDofs_+i] /= elementVolume;
                outVec6[elIdx*numDofs_+i] /= elementVolume;
                outVec7[elIdx*numDofs_+i] /= elementVolume;
                outVec8[elIdx*numDofs_+i] /= elementVolume;
              }               
              outVec3[elIdx] /= elementVolume;
              outVec5[elIdx] /= elementVolume;
            }
            
            vol += elementVolume;
          }

          if(!isSurf) 
          {
            correct = (meanVel*vol/durchfluss);
          }
             
          std::cout << std::endl 
                    << "Volume of fluid domain (vol):               "
                    << vol
                    << " m^3" << std::endl << std::endl;
          std::cout << "The following values are only valid for a "
                    << "straight pipe:" << std::endl;
          std::cout << "  Flow rate (R = \\int_{\\Omega} V d\\Omega):  "
                    << (durchfluss) << " m/s m^3" << std::endl;
          std::cout << "  Mean Velocity (V_ = R/vol):               "
                    << (durchfluss/vol)  << " m/s" << std::endl;
          std::cout << "  Profile correction factor (meanVel/V_):   "
                    << correct << std::endl;
          //             outVec2 *= correct;
             
          // deltaPhi *= correct;

          // Normalize node results by nodal weights
          if(nodalResults_) 
          {
            for(UInt nodeIndex=0, n=nodalWeights.GetSize(); nodeIndex < n; nodeIndex++)
            {
              Double w = nodalWeights[nodeIndex];
              
              for( UInt i=0; i < numDofs_; i++ )
              {
                outVec[nodeIndex*numDofs_+i] /= w;
                outVec2[nodeIndex*numDofs_+i] /= w;
                outVec4[nodeIndex*numDofs_+i] /= w;
                outVec6[nodeIndex*numDofs_+i] /= w;
                outVec7[nodeIndex*numDofs_+i] /= w;
                outVec8[nodeIndex*numDofs_+i] /= w;
              }               
              outVec3[nodeIndex] /= w;
              outVec5[nodeIndex] /= w;
            }
          }
             
          WriteResultsToCSV(freq, u_p_prime, deltaPhi, 0,
                            vol, meanVel, correct);
          WriteResultsToMatlab(freq, u_p_prime, deltaPhi, 0,
                               vol, meanVel, correct);

        } // switch: Analysis type


        // add result to output file
        if (output ) {
          output->AddResult( outResults[iRes] );
          output->AddResult( outResults2[iRes] );
          output->AddResult( outResults3[iRes] );
          output->AddResult( outResults4[iRes] );
          output->AddResult( outResults5[iRes] );
          output->AddResult( outResults6[iRes] );
          output->AddResult( outResults7[iRes] );
          output->AddResult( outResults8[iRes] );
        }
           
      } // loop over results
         
      if( output )
        output->FinishStep();
    }
    if( output )
      output->FinishMultiSequenceStep();
       
    if( output ) {

      std::stringstream sstr;
      param_->ToXML(sstr, 2);
         
      WriteStringToHDF5UserData(output, "cfstoolConfig", sstr.str());

      output->Finalize();
      std::cout << "\nOutput successfully written to " << outFile_ << std::endl;
    }
    delete ptGrid1;
    delete ptGrid2;
  } //WVT

  void WVT::OpenIntPointsFile() 
  {
    if(writeIntPointsFile_) {
      std::cout << "Opening integration point output file '" << intPointsFileName_
                << "'.." << std::endl;

      intPointsFile_.open(intPointsFileName_, std::ofstream::out | std::ofstream::trunc);
      
      intPointsFile_ << "# vtk DataFile Version 2.0" << std::endl
                     << "Integration points from openCFS for STAR-CCM+ Arbitrary Probe" << std::endl
                     << "ASCII" << std::endl
                     << "DATASET POLYDATA" << std::endl;
      intPointsFilePos_ = intPointsFile_.tellp();
      numIntPoints_ = 0;
      intPointsFile_ << "POINTS"
                     << "                                     "
                     << "                                     " << std::endl;
      intPointsFile_.setf( std::ios::fixed, std::ios::floatfield );
      intPointsFile_.precision(12);
    }
  }
  
  void WVT::CloseIntPointsFile()
  {
    if(writeIntPointsFile_) {
      std::cout << "Number of integration points written: " << numIntPoints_ 
                << std::endl;
      std::cout << "Closing integration point file '" << intPointsFileName_
                << "'.." << std::endl;

      std::stringstream sstr;
      sstr << "POINTS " << numIntPoints_ << " float";
      intPointsFile_.seekp(intPointsFilePos_);
      intPointsFile_.write(sstr.str().c_str(), sstr.str().length());
      intPointsFile_.close();
    }
  }

  void WVT::AcouPot2AcouVelFeFct(const shared_ptr<BaseFeFunction>& acouPotFeFct,
                                 const shared_ptr<BaseBOperator>& gradOp,
                                 shared_ptr<BaseFeFunction>& acouVelFeFct) 
  {
    // Create new CoefFunctionBOp for acoustic particle velocity
    PtrCoefFct partVelCoefFunc;

    shared_ptr<ResultInfo> partVelInfoElem(new ResultInfo());
    partVelInfoElem->resultType = ACOU_VELOCITY;
    partVelInfoElem->resultName = SolutionTypeEnum.ToString(partVelInfoElem->resultType);
    partVelInfoElem->dofNames = dofNames_;
    partVelInfoElem->unit = MapSolTypeToUnit(partVelInfoElem->resultType);
    partVelInfoElem->entryType = ResultInfo::VECTOR;
    partVelInfoElem->definedOn = ResultInfo::ELEMENT;

    CoefFunctionBOp<Complex>* pVCF = NULL;
    pVCF = new CoefFunctionBOp<Complex>(acouPotFeFct, partVelInfoElem, -1);
            
    const std::set<RegionIdType> regions = acouPotFeFct->GetRegions();
    std::set<RegionIdType>::const_iterator regIt, regEnd;
    regIt = regions.begin();
    regEnd = regions.end();
    
    for( ; regIt != regEnd; regIt++) 
    {
      RegionIdType id = *regIt;
      pVCF->AddBOperator(gradOp.get(), id);
    }
            
    partVelCoefFunc.reset(pVCF);

    // Create new FeFunction for acoustic particle velocity
    shared_ptr<ResultInfo> partVelInfoNode(new ResultInfo());
    partVelInfoNode->resultType = ACOU_VELOCITY;
    partVelInfoNode->resultName = SolutionTypeEnum.ToString(partVelInfoNode->resultType);
    partVelInfoNode->dofNames = dofNames_;
    partVelInfoNode->unit = MapSolTypeToUnit(partVelInfoNode->resultType);
    partVelInfoNode->entryType = ResultInfo::VECTOR;
    partVelInfoNode->definedOn = ResultInfo::NODE;

    // Create a dummy ParamNode for the FeSpace
    PtrParamNode aNode = PtrParamNode(new ParamNode());
    aNode->SetName("feSpaceNode");

    // Create a dummy info node for the FeSpace
    std::stringstream infoNodeStr;
    infoNodeStr << "{\"infoNode\":{" << std::endl
                << "\"h1Nodal\":{\"regionList\":{\"default\":{\"order\":\"\"}}}"
                << "}}" << std::endl;
    PtrParamNode infoNode = ParamNodeFromPropertyTree( infoNodeStr.str(),
                                                       "infoNode" );

    shared_ptr<FeSpace> nodalSpace = 
      FeSpace::CreateInstance(aNode, infoNode,
                              FeSpace::H1, acouPotFeFct->GetGrid());
    nodalSpace->SetDefaultRegionApproximation();

    acouVelFeFct.reset(new FeFunction<Complex>(mp_.get()));
    acouVelFeFct->SetGrid(acouPotFeFct->GetGrid());
    // Take care of cyclic dependency between FeFunction and FeSpaces
    acouVelFeFct->SetFeSpace(nodalSpace);
    nodalSpace->AddFeFunction(acouVelFeFct);
    acouVelFeFct->SetResultInfo(partVelInfoNode);
    partVelInfoNode->SetFeFunction(acouVelFeFct);
    nodalSpace->AddFeFunction(acouVelFeFct);

    StdVector< shared_ptr<EntityList> > eList = acouPotFeFct->GetEntityList();            
    for( UInt i=0,n=eList.GetSize() ; i<n; i++) 
    {
      acouVelFeFct->AddEntityList( eList[i] );
    }

    acouVelFeFct->AddExternalDataSource( partVelCoefFunc,
                                  eList );
    acouVelFeFct->SetFctId(PSEUDO_FCT_ID);
    nodalSpace->Finalize();
    acouVelFeFct->Finalize();
    acouVelFeFct->ApplyExternalData();
  }

  void WVT::WriteResultsToCSV(Double freq,
                              Complex u_p_prime, 
                              Double deltaPhiVol,
                              Double deltaPhiSurf,
                              Double vol,
                              Double meanVel,
                              Double meanVelCorrectionFactor)
  {
    boost::filesystem::ofstream csv;
    
    csv.open("output.csv", std::ios_base::binary);
    csv << "freq,u_p_prime_real,u_p_prime_imag,u_p_prime_ampl,"
        << "u_p_prime_phase,deltaPhi[rad],deltaPhi[deg],fluid_volume,mean_velocity,mean_velocity_correction_factor"
        << std::endl;

    std::cout << "\n(u')_P (cf. [Hemp1] (13)):                  "
              << u_p_prime << std::endl;
    std::cout << "\nDeltaPhi (cf. [Hemp1] (18)):                "
              << deltaPhiVol/M_PI*180 << " deg" << std::endl;
    
    csv << freq << ","
        << u_p_prime.real() << ","
        << u_p_prime.imag() << ","
        << abs(u_p_prime) << ","
        << arg(u_p_prime) << ","
        << deltaPhiVol << ","
        << (deltaPhiVol/M_PI*180) << ","
        << vol << ","
        << meanVel << ","
        << meanVelCorrectionFactor
        << std::endl;
    csv.close();         
  }

  void WVT::WriteResultsToMatlab(Double freq,
                                 Complex u_p_prime, 
                                 Double deltaPhiVol,
                                 Double deltaPhiSurf,
                                 Double vol,
                                 Double meanVel,
                                 Double meanVelCorrectionFactor)
  {
    // Create new MatFile object.
    CoupledField::MatFile mf("wvtout.mat");

    Double val = 0.0;

    mf.WriteVector("freq", 1, &freq);
    val = u_p_prime.real();
    mf.WriteVector("u_p_prime_real", 1, &val);
    val = u_p_prime.imag();
    mf.WriteVector("u_p_prime_imag", 1, &val);
    val = abs(u_p_prime);
    mf.WriteVector("u_p_prime_ampl", 1, &val);
    val = arg(u_p_prime);
    mf.WriteVector("u_p_prime_phase", 1, &val);
    mf.WriteVector("deltaPhi_rad", 1, &deltaPhiVol);
    val = (deltaPhiVol/M_PI*180);
    mf.WriteVector("deltaPhi_deg", 1, &val);
    mf.WriteVector("fluid_volume", 1, &vol);
    mf.WriteVector("mean_velocity", 1, &meanVel);
    mf.WriteVector("mean_velocity_correction_factor", 1, &meanVelCorrectionFactor);
    
    // Close .mat file and automatically prepend a header.
    mf.Close();
  }


}

