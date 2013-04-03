// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// Include headers which define what types
// of in/output files CFS++ should support
#include <boost/lexical_cast.hpp>

#include <def_cfs_stats.hh>
#include <def_use_mesh.hh>
#include <def_use_gidpost.hh>
#include <def_use_hdf5.hh>
#include <def_use_gmsh.hh>
#include <def_use_gmv.hh>
#include <def_use_unv.hh>
#include <def_use_ansysrst.hh>
#include <def_use_comsol.hh>

#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "ParamsInit.hh"
#include "General/Environment.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/SimInput.hh"
#include "DataInOut/SimOutput.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "Forms/Operators/GradientOperator.hh"

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/SimInputMESH.hh"
#endif

#ifdef USE_GMSH
#include "DataInOut/SimInOut/gmsh/SimInputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/SimOutputGmsh.hh"
#endif

#ifdef USE_GMV
#include "DataInOut/SimInOut/gmv/SimInputGMV.hh"
#include "DataInOut/SimInOut/gmv/SimOutGMV.hh"
#endif

#ifdef USE_HDF5
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"

#include "DataInOut/SimInOut/xdmf/SimOutputXDMF.hh"
#endif

#include "DataInOut/SimInOut/RefElems/SimInputRefElems.hh"

#ifdef USE_GIDPOST
#include "DataInOut/SimInOut/GiD/SimOutGiD.hh"
#endif

#ifdef USE_UNV
#include "DataInOut/SimInOut/Unverg/SimInputUnv.hh"
#include "DataInOut/SimInOut/Unverg/SimOutputUnv.hh"
#endif

#ifdef USE_ANSYSRST
#include "DataInOut/SimInOut/AnsysRST/SimOutputRST.hh"
#endif

#ifdef USE_COMSOL
#include "DataInOut/SimInOut/COMSOL/SimInputMPHTXT.hh"
#endif

#include "DataInOut/SimInOut/TextOutput/TextSimOutput.hh"


using namespace CoupledField;

//! Global parameter node and info node instance
PtrParamNode param, info;

namespace CFSTool {

  


  void setFreeCoord(std::string coordSysId="default",
      std::string node_name="averageDomain");

  inline void setParamNode(PtrParamNode paramNode, std::string name, std::string value,
      ParamNodeList* children = NULL);

  // Here we define our hard-coded node and element limits for testsuite
  // tests. These limits come into effect, if the reference file of a
  // Diff() operation has the extension '.h5ref' or if the --checkLimits
  // command line option is set to true.
  const UInt SOFT_NODE_LIMIT = 800;
  const UInt HARD_NODE_LIMIT = 1000;
  const UInt SOFT_ELEM_LIMIT = 1000;
  const UInt HARD_ELEM_LIMIT = 6000;
  
  // Here we print out the hard-coded node and element limit in response
  // to the --printLimits command line option.
  void PrintLimits(const std::string& type) 
  {
    bool typeCMake = (type == "cmake");

    std::cout << std::endl;
    
    if(typeCMake) { std::cout << "SET("; }
    std::cout << "SOFT_NODE_LIMIT " << SOFT_NODE_LIMIT;
    if(typeCMake) { std::cout << ")"; }
    std::cout << std::endl;

    if(typeCMake) { std::cout << "SET("; }
    std::cout << "HARD_NODE_LIMIT " << HARD_NODE_LIMIT;
    if(typeCMake) { std::cout << ")"; }
    std::cout << std::endl;

    if(typeCMake) { std::cout << "SET("; }
    std::cout << "SOFT_ELEM_LIMIT " << SOFT_ELEM_LIMIT;
    if(typeCMake) { std::cout << ")"; }
    std::cout << std::endl;

    if(typeCMake) { std::cout << "SET("; }
    std::cout << "HARD_ELEM_LIMIT " << HARD_ELEM_LIMIT;
    if(typeCMake) { std::cout << ")"; }
    std::cout << std::endl;

    std::cout << std::endl;
  }

  shared_ptr<SimInput> GetInputReader(const std::string& fileName)
  {
    // determine suffix of fileName
    shared_ptr<SimInput> reader;

    if (fileName.find( ".refelem") == std::string::npos &&
        !fs::exists( fileName ))
      EXCEPTION( "\nFile '" << fileName << "' does not exist!");

    if( fileName.find( ".mesh") != std::string::npos ) {
#ifdef USE_MESH
      reader = shared_ptr<SimInput>(new SimInputMESH( fileName, PtrParamNode(), info ) );
#else
      EXCEPTION( "No support for MESH input file format." );
#endif
    } else if( fileName.find( ".h5") != std::string::npos ) {
#ifdef USE_HDF5
      reader = shared_ptr<SimInput>(new SimInputHDF5(fileName, param, info));
#else
      EXCEPTION( "No support for HDF5 input file format." );
#endif
    } else if( fileName.find( ".refelem") != std::string::npos ) {
      reader = shared_ptr<SimInput>(new SimInputRefElems(fileName, param, info));
    } else if( fileName.find( ".msh") != std::string::npos ) {
#ifdef USE_GMSH
      PtrParamNode gmshNode(new ParamNode (ParamNode::EX, ParamNode::ELEMENT));
      reader = shared_ptr<SimInput>(new SimInputGmsh(fileName, gmshNode, info) );
#else  
      EXCEPTION( "No support for Gmsh input file format." );
#endif
    } else if( fileName.find( ".mphtxt") != std::string::npos ) {
#ifdef USE_COMSOL
      PtrParamNode comsolNode(new ParamNode (ParamNode::EX, ParamNode::ELEMENT));
      reader = shared_ptr<SimInput>(new SimInputMPHTXT(fileName, comsolNode, info) );
#else  
      EXCEPTION( "No support for Comsol .mphtxt input file format." );
#endif
    } else if( fileName.find( ".gmv") != std::string::npos ) {
#ifdef USE_GMV
      reader = shared_ptr<SimInput>(new SimInputGMV(fileName, param, info));
#else
      EXCEPTION( "No support for GMV input file format." );
#endif
    } else if( fileName.find( ".unv") != std::string::npos ||
        fileName.find( ".unverg") != std::string::npos ||
        fileName.find( ".unvref") != std::string::npos ) {
#ifdef USE_UNV
      reader = shared_ptr<SimInput>(new SimInputUnv( fileName, param, info ));
#else
      EXCEPTION( "No support for UNV input file format." );
#endif
    } else {
      EXCEPTION( "Found not suitable reader for file '" << fileName
          << "'" );
    }

    return reader;
  }

  shared_ptr<SimOutput> GetOutputWriter( const std::string& fileName ) {
    // determine suffix for fileName
    shared_ptr<SimOutput> writer;
    std::string baseName;

    if( fileName.find( ".post") != std::string::npos ) {
#ifdef USE_GIDPOST
      baseName = std::string(fileName, 0, fileName.find(".post"));
      PtrParamNode gidNode( new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode binary (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      binary->SetName( "binaryFormat");
      if( fileName.find( ".bin") != std::string::npos ) {
        binary->SetValue( "yes");
      } else {
        binary->SetValue( "false");
      }
      gidNode->AddChildNode( binary);
      writer = shared_ptr<SimOutput>( new SimOutputGiD( baseName, gidNode, info ) );
#else
      EXCEPTION( "No support for GiD output file format." );
#endif
    } else if( fileName.find( ".gmv") != std::string::npos ) {
#ifdef USE_GMV
      baseName = std::string(fileName, 0, fileName.find(".gmv"));
      PtrParamNode gmvNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode binary (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));

      binary->SetName("binaryFormat");
      binary->SetValue( "yes" );
      PtrParamNode fixedGrid (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      fixedGrid->SetName("fixedGrid");
      fixedGrid->SetValue( "yes" );
      gmvNode->AddChildNode( binary);
      gmvNode->AddChildNode( fixedGrid );
      writer = shared_ptr<SimOutput>( new SimOutputGMV( baseName, gmvNode, info ) );
#else
      EXCEPTION( "No support for GMV output file format." );
#endif
    } else if( fileName.find( ".msh") != std::string::npos ) {
#ifdef USE_GMSH
      baseName = std::string(fileName, 0, fileName.find(".msh"));
      PtrParamNode gmshNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode binary (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      binary->SetName("binaryFormat");
      binary->SetValue( "yes" );
      PtrParamNode bigEndian (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      bigEndian->SetName("endianness");
      bigEndian->SetValue( "big" );
      gmshNode->AddChildNode(binary);
      gmshNode->AddChildNode(bigEndian);
      writer = shared_ptr<SimOutput>( new SimOutputGmsh( baseName, gmshNode, info ) );
#else 
      EXCEPTION( "No support for GMsh output file format." );
#endif
    } else if(fileName.find( ".h5") != std::string::npos) {
#ifdef USE_HDF5
      baseName = std::string(fileName, 0, fileName.find(".h5"));
      PtrParamNode h5Node (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      eFiles->SetName("externalFiles");
      eFiles->SetValue( "false" );
      h5Node->AddChildNode(eFiles);
      writer =  shared_ptr<SimOutput>( new SimOutputHDF5( baseName, h5Node, info ) );
#else
      EXCEPTION( "No support for HDF5 output file format." );
#endif
    } else if(fileName.find( ".xmf") != std::string::npos) {
#ifdef USE_HDF5
      baseName = std::string(fileName, 0, fileName.find(".xmf"));
      PtrParamNode h5Node (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      eFiles->SetName("externalFiles");
      eFiles->SetValue( "false" );
      h5Node->AddChildNode(eFiles);
      writer =  shared_ptr<SimOutput>( new SimOutputXDMF( baseName, h5Node, info ) );
#else
      EXCEPTION( "No support for HDF5 output file format. Cannot write XDMF files." );
#endif
    } else if(fileName.find( ".rst") != std::string::npos) {
#ifdef USE_ANSYSRST
      baseName = std::string(fileName, 0, fileName.find(".rst"));
      PtrParamNode rstNode (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      writer =  shared_ptr<SimOutput>( new SimOutputRST( baseName, rstNode ) );
#else
      EXCEPTION( "No support for ANSYS .rst output file format." );
#endif
    } else if(fileName.find( ".unv") != std::string::npos) {
#ifdef USE_UNV
      baseName = std::string(fileName, 0, fileName.find(".unv"));
      PtrParamNode unvNode (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      if(fileName.find( ".unverg") != std::string::npos) 
      {
        PtrParamNode flavor (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
        flavor->SetName("flavor");
        flavor->SetValue( "CAPA" );
        unvNode->AddChildNode(flavor);
      }
      writer =  shared_ptr<SimOutput>( new SimOutputUnv( baseName, unvNode, info ) );
#else
      EXCEPTION( "No support for IDEAS universal output file format." );
#endif
    } else {
      EXCEPTION( "Output format not supported!" );
    }
    // create output writer
    return writer;
  }

  void Convert( const std::string& inFile, const std::string& outFile ) {

    // obtain input reader for inFile
    shared_ptr<SimInput> input = GetInputReader( inFile );

    // read in mesh
    input->InitModule();
    UInt dim = input->GetDim();
    Grid * ptGrid = new GridCFS(dim, param, info);
    input->ReadMesh(ptGrid);
    ptGrid->FinishInit();



    // obtain output writer
    shared_ptr<SimOutput> output = GetOutputWriter( outFile );

    // obtain number of multiSequenceSteps and get analysis types
    std::map<UInt, BasePDE::AnalysisType> types;
    std::map<UInt, UInt> numSteps;
    input->GetNumMultiSequenceSteps( types, numSteps );
    std::cout << "\nFound " << types.size() << " sequence step(s)\n";

    // check if the input reader has results
    bool printGridOnly = false;
    if( types.size() == 0) {
      printGridOnly = true;
      std::cerr << "Printing only grid\n";
    }

     output->Init( ptGrid, printGridOnly);

     // only iterate over results, if not only the mesh is converted
     if( !printGridOnly ) {


       // iterate over all multiSequenceSteps
       std::map<UInt,UInt>::iterator it;
       for( it = numSteps.begin(); it != numSteps.end(); it++ ) {

         UInt actMsStep = it->first;
         std::cout << "\n----------------------------\n"
                   << " Converting sequence step " << actMsStep << std::endl
                   << "----------------------------\n\n";
         // get resulttypes
         StdVector<shared_ptr<ResultInfo> > infos;
         input->GetResultTypes( actMsStep, infos );

         StdVector<shared_ptr<BaseResult> > results;
         std::map<UInt, Double> stepVals;
         std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;

         if( infos.GetSize() > 0 ){
             std::cout << "Converting the following results:\n";
         }
         // iterate over all result types
         for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {

           shared_ptr<ResultInfo> actRes = infos[iRes];
           std::cout << "\t" << actRes->resultName << std::endl;

           // get stepvalues
           input->GetStepValues( actMsStep, actRes, resultSteps[actRes] );
           stepVals.insert( resultSteps[actRes].begin(),
                            resultSteps[actRes].end() );

           // iterate over all regions
           StdVector<shared_ptr<EntityList> > resEntities;
           input->GetResultEntities( actMsStep, actRes, resEntities );
           for( UInt iEntity = 0; iEntity < resEntities.GetSize(); iEntity++ ) {
             // generate new result object and add it to output writer
             shared_ptr<BaseResult > result;
             if( types[actMsStep] == BasePDE::STATIC || 
                 types[actMsStep] == BasePDE::TRANSIENT ) {
               result = shared_ptr<BaseResult>( new Result<Double>() );
             } else {
               result = shared_ptr<BaseResult>( new Result<Complex>() );
             }
             result->SetEntityList( resEntities[iEntity] );
             result->SetResultInfo( actRes );
             results.Push_back( result );
             // Note: as the real values of saveBegin, saveInc and saveEnd are almost
             // nevert queried within an output format. we simply set saveBegin = 1,
             // saveInc = 1, saveEnd = number of result steps.
             output->RegisterResult( result, 1, 1, resultSteps[actRes].size(), false );
           }
         }

         // notify writer
         output->BeginMultiSequenceStep( actMsStep, types[actMsStep], numSteps[actMsStep] );

         // iterate over all stepvalues of this multisequence step
         for( UInt iStep = 0; iStep < numSteps[actMsStep]; iStep++ ) {

           // check, if current step contains any results
           if( stepVals.find( iStep+1) == stepVals.end() )
             continue;

           UInt actStepNum = iStep+1;
           Double actStepVal = stepVals[iStep+1];

           output->BeginStep( actStepNum, actStepVal );

           // iterate over all results
           for( UInt iRes = 0; iRes < results.GetSize(); iRes++) {

             // check if current result is defined within this step
             if( resultSteps[results[iRes]->GetResultInfo()].find(actStepNum)
                 == resultSteps[results[iRes]->GetResultInfo()].end() ) {
               continue;
             }

             try {
               input->GetResult( actMsStep, actStepNum, results[iRes] );
               output->AddResult( results[iRes] );
             } catch (Exception& ex ) {
             std::cerr <<  "\nResult '" << results[iRes]->GetResultInfo()->resultName
                       << "' in MsStep" << actMsStep << ", step " << actStepNum
                       << " could not be converted:\n\n";
             std::cerr << ex.what() << std::endl;
             }
           }
           output->FinishStep();
         }
         output->FinishMultiSequenceStep();

       } // loop over multisequence steps
     } // printGridOnly
     output->Finalize();
     delete ptGrid;
     std::cout << "\nOutput successfully written to " << outFile << std::endl;
  } //Convert()

  double RadPhase( const Complex& c ) {
    return std::atan2(c.imag() ,c.real() );
  }

  bool CheckReaderCapabilities(const std::string& readerDescription,
                               shared_ptr<SimInput>& reader,
                               std::string& result) 
  {
    std::stringstream sstr;
    bool ret = false;
    std::string fileName = reader->GetFileName();
    
    ret = std::find( reader->GetCapabilities().begin(),
                      reader->GetCapabilities().end(),
                      SimInput::MESH_RESULTS )
      != reader->GetCapabilities().end();

    sstr << "Reader for " << readerDescription << " with input file '"
         << fileName << "' " <<
      (ret ? "can provide results." :  "CANNOT provide results!");

    result = sstr.str();

    return ret;
  }

  Double Diff( const std::string& inFile_ref,
               const std::string& inFile_fut,
               const std::string& outFile,
               bool normedtomax,
               bool isHistory,
               std::string& maxDiffResultName) {

    // obtain input reader for inFiles
    StdVector< shared_ptr<SimInput> > inputs(2);
    inputs[0] = GetInputReader( inFile_ref );
    inputs[1] = GetInputReader( inFile_fut );

    // check capabilities of input class
    StdVector<bool> readerCaps(2);
    StdVector<std::string> readerDescriptions(2);
    readerDescriptions[0] = "reference results";
    readerDescriptions[1] = "file under test";
    StdVector<std::string> readerCapsResults(2);
    for(UInt i=0; i<inputs.GetSize(); i++) 
    {
      readerCaps[i] = CheckReaderCapabilities(readerDescriptions[i],
                                              inputs[i],
                                              readerCapsResults[i]);
    }

    if( !readerCaps[0] || !readerCaps[1] ) {
      std::cerr << "Some input files are only capable of handling meshes, not results!\n\n";

      for(UInt i=0; i<inputs.GetSize(); i++) 
      {
        std::cerr << readerCapsResults[i] << std::endl;
      }

      exit(EXIT_FAILURE);
    }

    // obtain input reader for inFiles
    shared_ptr<SimInput>& input_ref = inputs[0];
    shared_ptr<SimInput>& input_fut = inputs[1];

    // read in mesh of reference file
    input_ref->InitModule();
    UInt dim = input_ref->GetDim();
    Grid * ptGrid_ref = new GridCFS(dim, param, info);
    input_ref->ReadMesh(ptGrid_ref);
    ptGrid_ref->FinishInit();

    // Obtain settings for node and element limit check
    bool checkLimits = false;
    if (param->Has("checkLimits"))
    {
      checkLimits = param->Get("checkLimits")->As<bool>();         
    }
    else
    {
      if(inFile_ref.find( ".h5ref") != std::string::npos) 
      {
        checkLimits = true;
      }         
    }

    // Do the actual limits check
    if(checkLimits) 
    {
      UInt numNodes = ptGrid_ref->GetNumNodes();
      UInt numElems = ptGrid_ref->GetNumElems();
      std::string msg =
          "\nIn the interests of fast test execution and the ability to postprocess\n" \
          "using certain postprocessors with limitations (e.g. GiD without a license),\n" \
          "please keep the node and element count below the limits.\n";

      if(numNodes > SOFT_NODE_LIMIT) 
      {
        WARN("Number of nodes " << numNodes << 
             " exceeds soft limit " << SOFT_NODE_LIMIT << ".");
      }

      if(numElems > SOFT_ELEM_LIMIT) 
      {
        WARN("Number of elements " << numElems << 
             " exceeds soft limit " << SOFT_ELEM_LIMIT << ".");
      }

      if(numNodes > HARD_NODE_LIMIT) 
      {
        delete ptGrid_ref;
        EXCEPTION("Number of nodes " << numNodes << 
                  " exceeds hard limit " << HARD_NODE_LIMIT << "." <<
                  msg);
      }

      if(numElems > HARD_ELEM_LIMIT) 
      {
        delete ptGrid_ref;
        EXCEPTION("Number of elements " << numElems << 
                  " exceeds hard limit " << HARD_ELEM_LIMIT << "." << 
                  msg);
      }

    }


    // read in mesh of file under test
    input_fut->InitModule();
    Grid * ptGrid_fut = new GridCFS(dim, param, info);
    input_fut->ReadMesh(ptGrid_fut);
    ptGrid_fut->FinishInit();

    // obtain output writer
    bool printGridOnly = false;
    shared_ptr<SimOutput> output;
    if( outFile != "" ) {
      output = GetOutputWriter( outFile );
      output->Init( ptGrid_fut, printGridOnly);
    }

    // obtain number of Sequence Steps and get analysis types
    std::map<UInt, BasePDE::AnalysisType> types;
    std::map<UInt, UInt> numSteps;
    input_fut->GetNumMultiSequenceSteps( types, numSteps, isHistory );

    std::cout << "\nFound " << types.size() << " sequence step(s) in '" << inFile_fut << "'\n";
    std::map<UInt, BasePDE::AnalysisType> types_ref;
    std::map<UInt, UInt> numSteps_ref;
    input_ref->GetNumMultiSequenceSteps( types_ref, numSteps_ref, isHistory );
    std::cout << "\nFound " << types_ref.size() << " sequence step(s) in '" << inFile_ref << "'\n";

    if(types.size() != types_ref.size()){
      std::cout << "Reference file '" << inFile_ref << "' and file under test '" << inFile_fut
          << "' have different number of sequence steps!\n";
      exit(EXIT_FAILURE);
    }


    // iterate over all Sequence Steps
    Double maxDiff = 0.0;
    std::map<UInt,UInt>::iterator it;
    for( it = numSteps_ref.begin(); it != numSteps_ref.end(); it++ ) {

      UInt actMsStep = it->first;
      std::cout << " Diffing sequence step " << actMsStep << std::endl
          << "-------------------------\n\n";

      // get resulttypes
      StdVector<shared_ptr<ResultInfo> > infos, infos_ref;
      input_ref->GetResultTypes( actMsStep, infos_ref, isHistory );
      input_fut->GetResultTypes( actMsStep, infos, isHistory );

      StdVector<shared_ptr<BaseResult> > inResults_fut, inResults_ref, outResults;
      // stepnumbers, for which at least one result is defined
      std::map<UInt, Double> stepVals, stepVals_ref;
      // contains the stepnumbers/-values in which the particular result is
      // defined in
      std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
      std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps_ref;

      if( infos.GetSize() > 0 ){
        std::cout << "Performing diff on the following results:\n";
      }
      // iterate over all result types of input_ref
      for( UInt iRes = 0; iRes < infos_ref.GetSize(); iRes++) {

        std::cout << "\t" << infos_ref[iRes]->resultName << "\n";

        // get stepvalues of reference file
        shared_ptr<ResultInfo> actRes_ref = infos_ref[iRes];
        input_ref->GetStepValues( actMsStep, actRes_ref,
                                  resultSteps_ref[actRes_ref], isHistory);
        stepVals_ref.insert( resultSteps_ref[actRes_ref].begin(),
                             resultSteps_ref[actRes_ref].end() );

        // get stepvalues of file under test
        shared_ptr<ResultInfo> actRes = infos[iRes];
        input_fut->GetStepValues( actMsStep, actRes,
                                  resultSteps[actRes], isHistory);
        stepVals.insert( resultSteps[actRes].begin(),
                         resultSteps[actRes].end() );

        // Loop over all step values in both sets and compare them. Thus we can see
        // differences e.g. in eigenfrequency analysis
        std::map<UInt, Double>::const_iterator svIt_fut, svIt_ref;
        svIt_fut = stepVals.begin();
        svIt_ref = stepVals_ref.begin();
        for( ; svIt_ref != stepVals_ref.end(); ++svIt_ref, ++svIt_fut ) {
          if( svIt_fut->first != svIt_ref->first ) {
            EXCEPTION( "Encountered different result steps for result " << 
                       infos[iRes]->resultName );
          } else {

            Double val_fut = svIt_fut->second;
            Double val_ref = svIt_ref->second;
            Double relDiff = std::abs(std::abs(val_fut-val_ref))/std::abs(val_fut);
            if ( relDiff > 1e-4 ) {
              EXCEPTION("Time / Frequency values of step " << svIt_ref->first << " differ by " 
                        << relDiff*100.0 << " %:\n"
                        << "\treference value: " << val_ref <<" s / Hz\n\tcompared value:  " << val_fut
                        << " s / Hz" << std::endl );
            }
          }
        }


        // iterate over all regions
        StdVector<shared_ptr<EntityList> > regions;
        input_ref->GetResultEntities( actMsStep, infos_ref[iRes],
                                      regions, isHistory );
        for( UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++ ) {
          // generate new result object and add it to output writer
          shared_ptr<BaseResult > inResult_fut, inResult_ref, outResult;
          if( types[actMsStep] == BasePDE::TRANSIENT ||
              types[actMsStep] == BasePDE::STATIC ) {
            inResult_fut = shared_ptr<BaseResult>( new Result<Double>() );
            inResult_ref = shared_ptr<BaseResult>( new Result<Double>() );
            outResult = shared_ptr<BaseResult>( new Result<Double>() );
          } else {
            inResult_fut = shared_ptr<BaseResult>( new Result<Complex>() );
            inResult_ref = shared_ptr<BaseResult>( new Result<Complex>() );
            outResult = shared_ptr<BaseResult>( new Result<Complex>() );
          }
          inResult_fut->SetEntityList( regions[iRegion] );
          inResult_ref->SetEntityList( regions[iRegion] );
          outResult->SetEntityList( regions[iRegion] );

          inResult_fut->SetResultInfo( infos[iRes] );

          inResult_ref->SetResultInfo( infos_ref[iRes] );
          outResult->SetResultInfo( infos[iRes] );

          inResults_fut.Push_back( inResult_fut );
          inResults_ref.Push_back( inResult_ref );
          outResults.Push_back( outResult );
          if( output) {
            // Hardcoded: set output format to AMPL_PHASE
            //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
            outResult->GetResultInfo()->complexFormat = REAL_IMAG;

            // CAUTION: begin, inc, end are hardcoded and noch checked for each result
            output->RegisterResult( outResult, 1, 1,
                                    resultSteps[actRes].size(),
                                    isHistory );
          }
        }
      }

      Vector<Double> maxResVec_ref;
      maxResVec_ref.Resize( inResults_ref.GetSize() );

      // For transient simulation find maximum amplitude over all timesteps
      if( types[actMsStep] == BasePDE::STATIC ||
          types[actMsStep] == BasePDE::TRANSIENT ) {

        // iterate over all results
        for( UInt iRes = 0; iRes < inResults_ref.GetSize(); iRes++) {

          if(numSteps[actMsStep] != numSteps_ref[actMsStep]){
            std::cout << "Reference file '" << inFile_ref << "' has " << numSteps_ref[actMsStep] << " and '" << inFile_fut
                << "' has " << numSteps[actMsStep] << " time steps!\n";
            exit(EXIT_FAILURE);
          }

          maxResVec_ref[iRes] = 0.0;
          // iterate over all time steps
          for( UInt iStep = 0; iStep < numSteps_ref[actMsStep]; iStep++ ) {

            UInt actStepNum = iStep+1;
            // check if current result is defined within this step
            if( resultSteps[inResults_ref[iRes]->GetResultInfo()].find(actStepNum)
                == resultSteps[inResults_ref[iRes]->GetResultInfo()].end() ) {
              continue;
            }

            input_ref->GetResult( actMsStep, actStepNum, inResults_ref[iRes], isHistory );
            Vector<Double> & inVec_ref =
                dynamic_cast<Result<Double>& >(*inResults_ref[iRes]).GetVector();

            for( UInt i = 0; i<inVec_ref.GetSize(); i++ ) {
              if( std::abs(inVec_ref[i]) > maxResVec_ref[iRes] )
                maxResVec_ref[iRes] = std::abs(inVec_ref[i]);
            }
          }
          std::cout << "For result '" << inResults_ref[iRes]->GetResultInfo()->resultName
              << "' maximum amplitude is: " << maxResVec_ref[iRes] << "\n";
        }
      }


      // notify writer
      if( output) {
        output->BeginMultiSequenceStep( actMsStep, types_ref[actMsStep],
                                        numSteps_ref[actMsStep] );
      }

      // iterate over all time/frequency steps
      for( UInt iStep = 0; iStep < numSteps_ref[actMsStep]; iStep++ ) {

        // check, if current step contains any results
        if( stepVals.find( iStep+1) == stepVals.end() )
          continue;
        UInt actStepNum = iStep+1;
        Double actStepVal = stepVals[iStep+1];

        if( output) {
          output->BeginStep( actStepNum, actStepVal );
        }
        std::cout << "\n\t=============================================\n";
        std::cout << "\t  Treating step " << actStepNum << ", " << actStepVal
            << "s / Hz\n";
        std::cout << "\t=============================================\n";

        // iterate over all results
        for( UInt iRes = 0; iRes < inResults_fut.GetSize(); iRes++) {
          // check if current result is defined within this step
          if( resultSteps[inResults_fut[iRes]->GetResultInfo()].find(actStepNum)
              == resultSteps[inResults_fut[iRes]->GetResultInfo()].end() ) {
            continue;
          }

          // obtain both result objects for current step
          input_fut->GetResult( actMsStep, actStepNum, inResults_fut[iRes], isHistory );
          input_ref->GetResult( actMsStep, actStepNum, inResults_ref[iRes], isHistory );

          std::cout << "\n\t-- Comparing result " <<
              inResults_fut[iRes]->GetResultInfo()->resultName << " on " 
              << inResults_fut[iRes]->GetEntityList()->GetName() << " --\n";

          // get number of dofs of result
          UInt numDofs = inResults_fut[iRes]->GetResultInfo()->dofNames.GetSize();

          // cast result objects, get vector and calculate difference vector
          if( types[actMsStep] == BasePDE::STATIC || 
              types[actMsStep] == BasePDE::TRANSIENT ) {
            Vector<Double> & inVec_fut =
                dynamic_cast<Result<Double>& >(*inResults_fut[iRes]).GetVector();
            Vector<Double> & inVec_ref =
                dynamic_cast<Result<Double>& >(*inResults_ref[iRes]).GetVector();
            Vector<Double> & outVec =
                dynamic_cast<Result<Double>& >(*outResults[iRes]).GetVector();
            outVec.Resize( inVec_fut.GetSize() );

            // find maximum amplitude of inResult_ref
            for( UInt i = 0; i<inVec_ref.GetSize(); i++ ) {
              if( std::abs(inVec_ref[i]) > maxResVec_ref[iRes])
                maxResVec_ref[iRes] = std::abs(inVec_ref[i]);
            }

            // calculate difference entrywise
            outVec = inVec_fut - inVec_ref;
            if (normedtomax == true) {
              outVec /= maxResVec_ref[iRes];
            }

            // find maximum entry in difference vector
            for( UInt i = 0; i < outVec.GetSize(); i++ ) {
              if( std::abs(outVec[i]) > maxDiff) {
                maxDiff = std::abs(outVec[i]) ;
                maxDiffResultName = inResults_fut[iRes]->GetResultInfo()->resultName;
              }
            }

          } else if( types[actMsStep] == BasePDE::EIGENFREQUENCY ) {
            // in the eigenfrequency case we only compare the absolute value,
            // as the sign is not defined uniquely
            Vector<Complex> & inVec_fut =
                dynamic_cast<Result<Complex>& >(*inResults_fut[iRes]).GetVector();
            Vector<Complex> & inVec_ref =
                dynamic_cast<Result<Complex>& >(*inResults_ref[iRes]).GetVector();
            Vector<Complex> & outVec =
                dynamic_cast<Result<Complex>& >(*outResults[iRes]).GetVector();
            outVec.Resize( inVec_fut.GetSize() );

            // find maximum amplitude of inResult_ref in every frequency step
            Double maxRes_ref = 0.0;
            for( UInt i = 0; i<inVec_ref.GetSize(); i++ ) {
              if( std::abs(inVec_ref[i]) > maxRes_ref)
                maxRes_ref = std::abs(inVec_ref[i]);
            }

            Double aDiff, pDiff;
            // iterate over all dofs
            for (UInt dof = 0; dof<numDofs ; dof++) {
              // iterate over number of entities
              for( UInt i = 0, nEntities=UInt(inVec_ref.GetSize()/numDofs); i<nEntities; i++ ) {

                // index to access entity 'i' of dof 'dof'
                UInt actIndex = i * numDofs + dof;

                // amplitude difference
                if (normedtomax == true)
                  aDiff = std::abs(( std::abs(inVec_fut[actIndex]) - std::abs(inVec_ref[actIndex]) )/maxRes_ref);
                else 
                  aDiff = std::abs(std::abs(inVec_fut[actIndex]) - std::abs(inVec_ref[actIndex]));


                // phase difference in multiples of pi (just used for
                // grid diffs
                pDiff = RadPhase(inVec_fut[actIndex]) - RadPhase(inVec_ref[actIndex]);

                // correct 2*pi-offset if phase angles have different signs
                if ( (std::abs(pDiff)>PI) && (pDiff<0) )
                  pDiff+= 2*PI;
                if ( (std::abs(pDiff)>PI) && (pDiff>0) )
                  pDiff-= 2*PI;

                // Dirty hack! Write differences in real_imag format.
                outVec[actIndex] = Complex( aDiff, pDiff*180/PI );

                // return maximum of amplitude difference
                if  (aDiff > maxDiff)
                  maxDiff = aDiff;
              }


            } // loop dof
            std::cout << "\n\tMaximum overall rel. difference = " << maxDiff << "\n\n";

          }  else if( types[actMsStep] == BasePDE::HARMONIC ) { 
            Vector<Complex> & inVec_fut =
                dynamic_cast<Result<Complex>& >(*inResults_fut[iRes]).GetVector();
            Vector<Complex> & inVec_ref =
                dynamic_cast<Result<Complex>& >(*inResults_ref[iRes]).GetVector();
            Vector<Complex> & outVec =
                dynamic_cast<Result<Complex>& >(*outResults[iRes]).GetVector();
            outVec.Resize( inVec_fut.GetSize() );

            // find maximum amplitude of inResult_ref in every frequency step
            Double maxRes_ref = 0.0;
            for( UInt i = 0; i<inVec_ref.GetSize(); i++ ) {
              if( std::abs(inVec_ref[i]) > maxRes_ref)
                maxRes_ref = std::abs(inVec_ref[i]);
            }

            Double aDiff, pDiff, aMax=0.0, aMin=0.0, pMax=0.0, pMin=0.0;
            Double rDiff, iDiff, rMax=0.0, iMax=0.0;

            // iterate over all dofs
            for (UInt dof = 0; dof<numDofs ; dof++) {
              // iterate over number of entities
              for( UInt i = 0, nEntities=UInt(inVec_ref.GetSize()/numDofs); i<nEntities; i++ ) {

                // index to access entity 'i' of dof 'dof'
                UInt actIndex = i * numDofs + dof;

                // amplitude difference
                if (normedtomax == true)
                  aDiff = ( std::abs(inVec_fut[actIndex]) - std::abs(inVec_ref[actIndex]) )/maxRes_ref;
                else
                  aDiff = std::abs(inVec_fut[actIndex]) - std::abs(inVec_ref[actIndex]);

                // phase difference in multiples of pi
                pDiff = RadPhase(inVec_fut[actIndex]) - RadPhase(inVec_ref[actIndex]);

                // correct 2*pi-offset if phase angles have different signs
                if ( (std::abs(pDiff)>PI) && (pDiff<0) )
                  pDiff+= 2*PI;
                if ( (std::abs(pDiff)>PI) && (pDiff>0) )
                  pDiff-= 2*PI;

                // Dirty hack! Write differences in real_imag format.
                outVec[actIndex] = Complex( aDiff, pDiff*180/PI );

                // maximum and minimum values
                if( pDiff > pMax )
                  pMax = pDiff;
                if( pDiff < pMin )
                  pMin = pDiff;
                if( aDiff > aMax )
                  aMax = aDiff;
                if( aDiff < aMin )
                  aMin = aDiff;

                // maximum difference in real and imaginary part
                rDiff = std::abs( inVec_fut[actIndex].real() - inVec_ref[actIndex].real() )/maxRes_ref;
                iDiff = std::abs( inVec_fut[actIndex].imag() - inVec_ref[actIndex].imag() )/maxRes_ref;
                if ( rDiff > rMax )
                  rMax = rDiff;
                if ( iDiff > iMax)
                  iMax = iDiff;
              }
              if( normedtomax == true)
                std::cout << "\n\tMaximum rel. + amplitude difference:  " << aMax*100 << " %\n"
                << "\tMaximum rel. - amplitude difference: " << aMin*100 << " %\n";
              else
                std::cout << "\n\tMaximum + amplitude difference:  " << aMax <<  "\n"
                << "\tMaximum - amplitude difference: " << aMin <<  "\n";

              std::cout << "\tMaximum + phase difference:      " << pMax*180/PI <<  " deg\n"
                  << "\tMaximum - phase difference:     " << pMin*180/PI <<  " deg\n";

              // return maxDiff for differences in real and imaginary part
              if ( (rMax > iMax) && (rMax > maxDiff) )
                maxDiff = rMax;
              else if ( (iMax > rMax) && (iMax > maxDiff) )
                maxDiff = iMax;
            } // loop dof
            std::cout << "\n\tMaximum overall rel. difference = " << maxDiff << "\n\n";

          } // switch: Analysitype
          // add result to output file
          if (output )
            output->AddResult( outResults[iRes] );
        } // loop over results
        if( output )
          output->FinishStep();
      }
      if( output )
        output->FinishMultiSequenceStep();

    }
    if( output ) {
      output->Finalize();
      std::cout << "\nOutput successfully written to " << outFile << std::endl;
    }
    delete ptGrid_fut;
    delete ptGrid_ref;

    return maxDiff;

  } //Diff

  void WVT( const std::string& lateral_mode_file,
            const std::string& coriolis_mode_file,
            const std::string& mean_flow_file,
            const std::string& outFile ) {

       bool printGridOnly = false;

       // obtain input reader for inFiles
       StdVector< shared_ptr<SimInput> > inputs(3);
       inputs[0] = GetInputReader( lateral_mode_file );
       inputs[1] = GetInputReader( coriolis_mode_file );
       inputs[2] = GetInputReader( mean_flow_file );

       // check capabilities of input class
       StdVector<bool> readerCaps(3);
       StdVector<std::string> readerDescriptions(3);
       readerDescriptions[0] = "lateral mode";
       readerDescriptions[1] = "coriolis mode";
       readerDescriptions[2] = "mean flow";
       StdVector<std::string> readerCapsResults(3);
       for(UInt i=0; i<inputs.GetSize(); i++) 
       {
         readerCaps[i] = CheckReaderCapabilities(readerDescriptions[i],
                                                 inputs[i],
                                                 readerCapsResults[i]);
       }
       
       if( !readerCaps[0] || !readerCaps[1] || !readerCaps[2] ) {
         std::cerr << "Some input files are only capable of handling meshes, not results!\n\n";

         for(UInt i=0; i<inputs.GetSize(); i++) 
         {
           std::cerr << readerCapsResults[i] << std::endl;
         }

         exit(EXIT_FAILURE);
       }

       // read in mesh of input1
       inputs[0]->InitModule();
       UInt dim = inputs[0]->GetDim();
       Grid * ptGrid1 = new GridCFS(dim, param, info);
       inputs[0]->ReadMesh(ptGrid1);
       ptGrid1->FinishInit();

       // read XML and material file from lateral mode file
       SimInputHDF5* hdf5Reader = dynamic_cast<SimInputHDF5*>(inputs[0].get());
       std::string xmlFile, matFile;
       hdf5Reader->ReadStringFromUserData("ParameterFile", xmlFile);
       hdf5Reader->ReadStringFromUserData("MaterialFile", matFile);

       // std::cout << "##############" << xmlFile << std::endl
       //    << "##############" << matFile << std::endl;
       PtrParamNode pNode;

       {
         CoupledField::Xerces xerces;
         xerces.SetString(xmlFile);
         pNode = xerces.CreateParamNodeInstance();
       }       

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

       {
         CoupledField::Xerces matFileParser;
         matFileParser.SetString(matFile);
         pNode = matFileParser.CreateParamNodeInstance();
       }
       
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
       inputs[1]->InitModule();
       Grid * ptGrid2 = new GridCFS(dim, param, info);
       inputs[1]->ReadMesh(ptGrid2);
       ptGrid2->FinishInit();

       // read in mesh of input3
       inputs[2]->InitModule();
       Grid * ptGrid3 = new GridCFS(dim, param, info);
       inputs[2]->ReadMesh(ptGrid3);
       ptGrid3->FinishInit();

       // obtain output writer
       shared_ptr<SimOutput> output;
       if( outFile != "" ) {
         output = GetOutputWriter( outFile );
         output->Init( ptGrid1, printGridOnly);
       }

       // obtain number of Sequence Steps and get analysis types
       std::map<UInt, BasePDE::AnalysisType> types;
       std::map<UInt, UInt> numSteps;
       inputs[0]->GetNumMultiSequenceSteps( types, numSteps, false );

       std::cout << "\nFound " << types.size() << " sequence step(s) in '" << lateral_mode_file << "'\n";
       std::map<UInt, BasePDE::AnalysisType> types2;
       std::map<UInt, UInt> numSteps2;
       inputs[1]->GetNumMultiSequenceSteps( types2, numSteps2, false );
       std::cout << "\nFound " << types2.size() << " sequence step(s) in '" << coriolis_mode_file << "'\n";
       
       if(types.size() != types2.size()){
         std::cout << "'" << lateral_mode_file << "' and '" << coriolis_mode_file
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
       inputs[0]->GetResultTypes( actMsStep, infos, false );
       inputs[1]->GetResultTypes( actMsStep, infos2, false );
       inputs[2]->GetResultTypes( actMsStep, infos3, false );
       
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
         inputs[2]->GetStepValues( actMsStep, actRes3,
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
         inputs[0]->GetStepValues( actMsStep, actRes,
                                resultSteps[actRes], false);
         stepVals.insert( resultSteps[actRes].begin(),
                          resultSteps[actRes].end() );
         
         // get stepvalues of second file
         inputs[1]->GetStepValues( actMsStep, actRes2,
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
         inputs[0]->GetResultEntities( actMsStep, infos[iRes],
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
           inputs[0]->GetResult( actMsStep, actStepNum, inResults1[iRes], false );
           inputs[1]->GetResult( actMsStep, actStepNum, inResults2[iRes], false );
           inputs[2]->GetResult( actMsStep, actStepNum, inResults_mean_flow[iRes], false );
           
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
             
             shared_ptr<BaseFeFunction> u1FeFunc = inputs[0]->GetFeFunction<Complex>( actMsStep,
                                                                             actStepNum,
                                                                             inResults1[iRes]->GetResultInfo()->resultType,
                                                                             regionNames);
             shared_ptr<BaseFeFunction> u2FeFunc = inputs[1]->GetFeFunction<Complex>( actMsStep,
                                                                             actStepNum,
                                                                             inResults2[iRes]->GetResultInfo()->resultType,
                                                                             regionNames);

             shared_ptr<BaseFeFunction> VFeFunc = inputs[2]->GetFeFunction<Double>( actMsStep,
                                                                              actStepNum,
                                                                              inResults_mean_flow[iRes]->GetResultInfo()->resultType,
                                                                              regionNames );


             GradientOperator<FeH1,3> gradOp;
             
             EntityIterator eIt = outResults[iRes]->GetEntityList()->GetIterator();
             
             Complex u_p_prime = Complex(0.0, 0.0);
             Double rho_0 = regionDensityMap[inResults1[iRes]->GetEntityList()->GetName()];

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
  
  /**
   * calcAverage calclulates the average of a physical field at each timestep. It
   * stores the results in outFile
   * @param inFile The file which carries the data (e.g. acoustic pressure , mechanical
   * displacement)
   * @param meshFile The mesh file to associate each value with a point in space
   * Import not all points should be included into the averagin.
   * @return void
   */

  void calcAverage( const std::string& inFile,
              const std::string& outFile)
  {
    std::string isFreecoord = param->Get("freeCoord")->As<std::string>();
    if (isFreecoord != "")
    {
      std::cerr << "selection of region not implemented!\n";
      exit(EXIT_FAILURE);
      // TODO: this is the call then for setting everything up
      //setFreeCoord();
    }

    // obtain input reader for inFiles
    shared_ptr<SimInput> input = GetInputReader(inFile);

    // check capabilities of input class
    bool printGridOnly = false;
    if (std::find( input->GetCapabilities().begin(),
          input->GetCapabilities().end(),
          SimInput::MESH_RESULTS )
        == input->GetCapabilities().end())
    {
      std::cerr << "input file is only capable of handling mesh, not results!\n";
      exit(EXIT_FAILURE);
    }

    // read in mesh of input
    input->InitModule();
    UInt dim = input->GetDim();
    Grid * ptGrid = new GridCFS(dim, param, info);
    input->ReadMesh(ptGrid);
    ptGrid->FinishInit();

    // obtain output writer
    shared_ptr<SimOutput> output;
    output = GetOutputWriter( outFile );
    output->Init( ptGrid, printGridOnly);

    // obtain number of Sequence Steps and get analysis types
    std::map<UInt, BasePDE::AnalysisType> types;
    std::map<UInt, UInt> numSteps;
    input->GetNumMultiSequenceSteps( types, numSteps, false );

    std::cout << "\nFound " << types.size() << " sequence step(s)\n";

    // iterate over all Sequence Steps
    std::map<UInt,UInt>::iterator it;
    for (it = numSteps.begin(); it != numSteps.end(); it++)
    {
      UInt actMsStep = it->first;
      std::cout << " averaging step " << actMsStep << std::endl
        << "-------------------------\n\n";

      // get resulttypes
      StdVector<shared_ptr<ResultInfo> > infos;
      input->GetResultTypes( actMsStep, infos, false );

      StdVector<shared_ptr<BaseResult> > inResults, outResults;
      // stepnumbers, for which at least one result is defined
      std::map<UInt, Double> stepVals;
      // contains the stepnumbers/-values in which the particular result is
      // defined in
      std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;

      if (infos.GetSize() > 0)
      {
        std::cout << "Performing average on the following results:\n";
      }
      // iterate over all result types of input
      for (UInt iRes = 0; iRes < infos.GetSize(); iRes++)
      {
        std::cout << "\t" << infos[iRes]->resultName << "\n\n";

        // get stepvalues
        shared_ptr<ResultInfo> actRes = infos[iRes];
        input->GetStepValues( actMsStep, actRes,
            resultSteps[actRes], false);
        stepVals.insert( resultSteps[actRes].begin(),
            resultSteps[actRes].end() );

        // iterate over all regions
        StdVector<shared_ptr<EntityList> > regions;
        input->GetResultEntities( actMsStep, infos[iRes],
            regions, false );
        for (UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++)
        {
          // generate new result object and add it to output writer
          shared_ptr<BaseResult > inResult, outResult;
          if (types[actMsStep] != BasePDE::HARMONIC)
          {
            inResult  = shared_ptr<BaseResult>( new Result<Double>() );
            outResult = shared_ptr<BaseResult>( new Result<Double>() );
          } else {
            std::cerr << "Averaging over harmonic results does not make sense\n";
            exit(EXIT_FAILURE);
          }
          inResult->SetEntityList( regions[iRegion] );
          outResult->SetEntityList( regions[iRegion] );

          inResult->SetResultInfo( infos[iRes] );
          outResult->SetResultInfo( infos[iRes] );

          inResults.Push_back( inResult );
          outResults.Push_back( outResult );
          if (output)
          {
            // Hardcoded: set output format to AMPL_PHASE
            //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
            outResult->GetResultInfo()->complexFormat = REAL_IMAG;

            // CAUTION: begin, inc, end are hardcoded and noch checked for each result
            output->RegisterResult( outResult, 1, 1,
                                    resultSteps[actRes].size(),
                                    false );
          }
        }
      }

      // notify writer
      if (output)
      {
        output->BeginMultiSequenceStep( actMsStep, types[actMsStep],
            numSteps[actMsStep] );
      }

      // iterate over all time/frequency steps
      for (UInt iStep = 1; iStep <= numSteps[actMsStep]; ++iStep)
      {
        // check, if current step contains any results
        if (stepVals.find(iStep) == stepVals.end())
        {
          continue;
        }
        Double actStepVal = stepVals[iStep];

        if (output)
        {
          output->BeginStep(iStep, actStepVal);
        }
        // iterate over all results
        for (UInt iRes = 0; iRes < inResults.GetSize(); iRes++)
        {
          // check if current result is defined within this step
          if (resultSteps[inResults[iRes]->GetResultInfo()].find(iStep)
              == resultSteps[inResults[iRes]->GetResultInfo()].end())
          {
            continue;
          }

          // obtain result objects for current step
          input->GetResult(actMsStep, iStep, inResults[iRes], false);

          // get number of dofs of result
#if 0 // TODO check if needed
          UInt numDofs = inResults[iRes]->GetResultInfo()->dofNames.GetSize();
#endif

          // cast result objects, get vector and calculate difference vector
          if (types[actMsStep] != BasePDE::HARMONIC)
          {
            Double meanVal = 0.0;
            Vector<Double> & inVec =
              dynamic_cast<Result<Double>& >(*inResults[iRes]).GetVector();
            Vector<Double> & outVec =
              dynamic_cast<Result<Double>& >(*outResults[iRes]).GetVector();
            UInt inVec_size = inVec.GetSize();

            // sum up and divide by number of entries <- averaging
            for (UInt i = 0; i < inVec_size; ++i)
            {
              meanVal += inVec[i];
            }
            meanVal /= (Double)inVec_size;
            outVec.Resize(inVec_size);
            for (UInt i = 0; i < inVec_size; ++i)
            {
              outVec[i] = meanVal;
            }
          } else {
            std::cerr << "Averaging over harmonic results does not make sense\n";
            exit(EXIT_FAILURE);
          }
          // add result to output file
          if (output)
          {
              output->AddResult(outResults[iRes]);
          }
        }
        if (output)
        {
            output->FinishStep();
        }
      }
      if (output)
      {
          output->FinishMultiSequenceStep();
      }

    }
    if (output)
    {
      output->Finalize();
      std::cout << "\nOutput successfully written to " << outFile << std::endl;
    }
    delete ptGrid;
  } //calcAverage

  /** Initialize static Enums.
   * todo: do better once - Fabian */
  void InitEnums()
  {
    SetEnvironmentEnums();
    BasePDE::SetEnums();
    EntityList::SetEnums();
  }

  void setFreeCoord(std::string coordSysId,
      std::string node_name)
  {

    std::vector<std::string> freeCoord;
    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");

    // Initialize vector with output fields
    Tok tokenizer(param->Get("freeCoord")->As<std::string>(), sep);
    std::copy(tokenizer.begin(), tokenizer.end(),
              std::back_inserter(freeCoord));
    if (freeCoord.size() != 4)
    {
      EXCEPTION("Not enought arguments " << freeCoord.size() \
          << ". Need 4 arguments (comp, start, stop, inc')");
    }

    PtrParamNode parent;

    ParamNodeList childVec;
    childVec.Push_back(PtrParamNode(new ParamNode())); // comp
    childVec.Push_back(PtrParamNode(new ParamNode())); // start
    childVec.Push_back(PtrParamNode(new ParamNode())); // stop
    childVec.Push_back(PtrParamNode(new ParamNode())); // inc

    // create a tree (Paramnode):
    // create leafes (comp, start, stop, inc)
    setParamNode(childVec[0], "comp",  freeCoord[0]);
    setParamNode(childVec[1], "start", freeCoord[1]);
    setParamNode(childVec[2], "stop",  freeCoord[2]);
    setParamNode(childVec[3], "inc",   freeCoord[3]);
    

    // create node (freeCoord) with comp,start,stop,inc as children
    parent = PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
    setParamNode(parent, "freeCoord", "", &childVec);
    
    // create node (list) with freeCoord and coordSysId as children
    childVec.Clear();
    childVec.Push_back(PtrParamNode( new ParamNode())); // coordSysId
    childVec.Push_back(parent);
    setParamNode(childVec[0], "coordSysId", coordSysId);
    childVec.Push_back(PtrParamNode( new ParamNode())); // gridId
    childVec.Push_back(parent);
    setParamNode(childVec[1], "gridId", "default");
    parent = PtrParamNode( new ParamNode()); //list
    setParamNode(parent, "list", "", &childVec);

    // create leaf (name)
    childVec.Clear();
    childVec.Push_back(parent);
    parent = PtrParamNode( new ParamNode()); //nodes name=""
    setParamNode(parent, "name", node_name);
    childVec.Push_back(parent);


    // create node (nodes) with name and list as childer
    parent = PtrParamNode( new ParamNode()); //nodes
    setParamNode(parent, "nodes", "", &childVec);

    // create node (nodeList) with nodes as child
    childVec.Clear();
    childVec.Push_back(parent);
    parent = PtrParamNode( new ParamNode()); //nodeList
    setParamNode(parent, "nodeList", "", &childVec);

    // create node (domain) with nodeList as child
    childVec.Clear();
    childVec.Push_back(parent); //domain
    //parent = new ParamNode(false); //nodeList
    //setParamNode(param, "domain", "", &childVec);
    //childVec.Clear();

    // put the childs into param
    param->Get("domain",ParamNode::INSERT )->GetChildren().Push_back(parent);
  }

  inline void setParamNode(PtrParamNode paramNode, std::string name, std::string value,
      ParamNodeList* children  )
  {
    paramNode->SetName(name);
    if (name != "")
    {
      paramNode->SetValue(value);
    }
    if (children != NULL)
    {
      ParamNodeList &childTmp = paramNode->GetChildren();
      childTmp = *children;
    }
  }



} //Namespace CFSTool


int main(int argc, char** argv)
{

  // todo: do better once! - Fabian
  CFSTool::InitEnums(); 
  ElemShape::Initialize();

  domain = NULL;
  
  std::string maxDiffResultName; // that's chaos man!
  std::string infoFileName;
  try
  {
    param = PtrParamNode(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));

    ParamsInit(argc, argv, param);

    // Switch this flag tc true for debugging
    if (param->Get("forceSegFault")->As<bool>())
    {
      Exception::segfault_ = true;
    } else {
      Exception::segfault_ = false;
    }

    // Print out hard-coded node and element limits and exit.
    if (param->Has("printLimits"))
    {
      CFSTool::PrintLimits(param->Get("printLimits")->As<std::string>());
      return 0;
    }

    std::cout << std::endl
              << "============================================================"
              << "===========" << std::endl;
    std::cout << " CFSTOOL - File Conversions/Comparisons for CFS++" << std::endl << std::endl
              << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'"
              << " (rev " << CFS_WC_REVISION << ")" << std::endl
              << " compiled " << __DATE__
              << " as " << CMAKE_BUILD_TYPE << std::endl;
    std::cout << "============================================================"
              << "==========="
              << std::endl << std::endl;

    std::string param_mode = param->Get("mode")->As<std::string>();

    // get filenames from parameter
    std::string file1 = param->Get("file1")->As<std::string>();
    std::string file2 = param->Get("file2")->As<std::string>();
    std::string file3 = param->Get("file3")->As<std::string>();
    std::string file4 = param->Get("file4")->As<std::string>();

    // Initialize vector with output fields
    UInt num_files = 0;
    if (file1 != "")
    {
      ++num_files;
      if (file2 != "")
      {
        ++num_files;
        if (file3 != "")
        {
          ++num_files;
          if (file4 != "")
          {
            ++num_files;
          }
        }
      }
    }

    if (file1 != "")
    {
      // This is necessary to run, but I do not know what it is for
      info = PtrParamNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT));
      infoFileName = param_mode + "_" + file1 + ".info.xml";
      info->SetName("cfsInfo"); 
    }

    if (param_mode == "calcAverage")
    {
      if (file3 != "")
      {
        EXCEPTION( "Two many arguments, please only provide two files. (in- and output file)" );
      }
      file3 = file2;
      if (num_files != 2)
      {
        EXCEPTION( "Please provide a reference file and output File" );
      }
      CFSTool::calcAverage(file1, file3);
    } else if (param_mode == "convert") {
      if (file3 != "")
      {
        EXCEPTION( "Two many arguments, please only provide two files. (in- and output file)" );
      }
      file3 = file2;
      if (num_files != 2)
      {
        EXCEPTION( "Please provide 'input_file' and 'output_file'" );
      }
      CFSTool::Convert( file1, file3 );
    } else if (param_mode == "scalardiff") {
      Double tolerance = param->Get("eps")->As<Double>();
      if (num_files != 2)
      {
        EXCEPTION( "Please provide 'reference_file' and 'file_under_test'" );
      }
      Double maxDiffMesh = 0.0, maxDiffHist = 0.0;
      std::cout << "Checking for mesh results:\n"
        << "==========================\n";
      maxDiffMesh = CFSTool::Diff( file1, file2, "", \
                                  true, false, maxDiffResultName);
      std::cout << "Checking for history results:\n"
        << "=============================\n";
      maxDiffHist = CFSTool::Diff( file1, file2, "", \
                                  true, true, maxDiffResultName );
      Double maxDiff = std::max( maxDiffMesh, maxDiffHist );
      if( maxDiff > tolerance ) {
        std::cout << "'" << file1 << "' and '" << file2
          << "' have maximum difference " << maxDiff
          << " at '" << maxDiffResultName << "'\n";
        exit(EXIT_FAILURE);
      } else {
        std::cout << "  No differences larger than tolerance found.\n";
        exit(EXIT_SUCCESS);
      }
    } else if (param_mode == "meshdiff") {
      if (num_files != 3)
      {
        EXCEPTION( "Please provide 'reference_file', 'file_under_test' and 'out_file'" );
      }
      CFSTool::Diff( file1, file2, file3, \
                     false, false, maxDiffResultName);
    } else if (param_mode == "meshdiffnormed") {
      if (num_files != 3)
      {
        EXCEPTION( "Please provide 'reference_file', 'file_under_test' and 'out_file'" );
      }
      CFSTool::Diff( file1, file2, file3, \
                     true, false, maxDiffResultName);
    } else if (param_mode == "wvt") {
      if (num_files != 4)
      {
        EXCEPTION( "Please provide 'lateral_mode_file', 'coriolis_mode_file',\n 'mean_velocity_file' and 'out_file'" );
      }

      CFSTool::WVT( file1, file2, file3, file4 );

    } else {
      EXCEPTION( "No such mode: " << param_mode <<". See help for available modes." );
      return EXIT_FAILURE;
    }
  } catch(std::exception& ex) {
    std::cerr << "The following error occured during program execution:\n\n" << ex.what();

    if (info != NULL)
    {
      info->Get(ParamNode::PN_ERROR)->SetValue(ex.what());
      info->ToFile(infoFileName);
    }

    return -1;
  }
  info->ToFile(infoFileName);
  return 0;
}
