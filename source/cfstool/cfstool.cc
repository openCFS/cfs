// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <complex>
#include <exception>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/path.hpp"
// Include headers which define what types
// of in/output files CFS++ should support
#include "boost/lexical_cast.hpp"
#include "boost/tokenizer.hpp"
#include "def_cfs_stats.hh"
#include "def_use_ansysrst.hh"
#include "def_use_gidpost.hh"
#include "def_use_gmsh.hh"
#include "def_use_gmv.hh"
#include "def_use_hdf5.hh"
#include "def_use_mesh.hh"
#include "def_use_unv.hh"
#include "def_use_fftw.hh"

namespace fs = boost::filesystem;

#include "DataInOut/simInput.hh"
#include "DataInOut/simOutput.hh"
#include "Domain/GridCFS/grid_cfs.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/vector.hh"
#include "PDE/basePDE.hh"
#include "ParamsInit.hh"
#include "Utils/StdVector.hh"
#include "Utils/result.hh"

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/simInputMESH.hh"
#endif

#ifdef USE_GMSH
#include "DataInOut/SimInOut/gmsh/simInputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/simOutputGmsh.hh"
#endif

#ifdef USE_GMV_INPUT
#include "DataInOut/SimInOut/gmv/simInputGMV.hh"
#endif

#ifdef USE_GMV_OUTPUT
#include "DataInOut/SimInOut/gmv/simOutGMV.hh"
#endif

#ifdef USE_HDF5
#include "DataInOut/SimInOut/hdf5/simInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/simOutputHDF5.hh"
#include "DataInOut/SimInOut/xdmf/simOutputXDMF.hh"
#endif

#ifdef USE_GIDPOST
#include "DataInOut/SimInOut/GiD/simOutGiD.hh"
#endif

#ifdef USE_UNV
#include "DataInOut/SimInOut/Unverg/SimInputUnv.hh"
#include "DataInOut/SimInOut/Unverg/SimOutputUnv.hh"
#endif

#ifdef USE_ANSYSRST
#include "DataInOut/SimInOut/AnsysRST/simOutputRST.hh"
#endif


#ifdef CFSTOOL_FFTW
#include <fftw3.h>
#include "fft.hh"

//include boost stuff for file system access for chunced fft
#include <unistd.h>
#include <fstream>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/progress.hpp>
#endif

#include "DataInOut/ParamHandling/ParamNode.hh"

using namespace CoupledField;

namespace CFSTool {

#ifdef CFSTOOL_FFTW
  void prepareChunkResults(std::vector< std::vector< std::vector<Complex> > >& fftVals,
                              UInt size, UInt nChunk,StdVector<std::string> fNames,
                              UInt numRes);

  void updateChunkResult(std::vector< std::vector< std::vector<Complex> > >& fftVals,
                            UInt freqStep,const StdVector<std::string> & fNames,
                            UInt chunkSize,UInt numRes);

  void updateChunkResultone(std::vector< std::vector< std::vector<Double> > >& fftVals,
                            UInt freqStep,const StdVector<std::string> & fNames,
                            UInt chunkSize,UInt numRes);
#endif

  void setFreeCoord(std::string coordSysId="default",
      std::string node_name="averageDomain");
  inline void setParamNode(PtrParamNode paramNode, std::string name, std::string value,
      ParamNodeList* children = NULL);

  shared_ptr<SimInput> GetInputReader(const std::string& fileName)
  {
    // determine suffix of fileName
    shared_ptr<SimInput> reader;

    if (!fs::exists( fileName ))
      EXCEPTION( "\nFile '" << fileName << "' does not exist!");

    if( fileName.find( ".mesh") != std::string::npos ) {
#ifdef USE_MESH
      reader = shared_ptr<SimInput>(new SimInputMESH( fileName, PtrParamNode() ) );
#else
      EXCEPTION( "No support for MESH input file format." );
#endif
    } else if( fileName.find( ".h5") != std::string::npos ) {
#ifdef USE_HDF5
      reader = shared_ptr<SimInput>(new SimInputHDF5(fileName, param));
#else
      EXCEPTION( "No support for HDF5 input file format." );
#endif
    } else if( fileName.find( ".msh") != std::string::npos ) {
#ifdef USE_GMSH
      PtrParamNode gmshNode(new ParamNode (ParamNode::EX, ParamNode::ELEMENT));
      reader = shared_ptr<SimInput>(new SimInputGmsh(fileName, gmshNode) );
#else  
      EXCEPTION( "No support for Gmsh input file format." );
#endif
    } else if( fileName.find( ".gmv") != std::string::npos ) {
#ifdef USE_GMV_INPUT
      reader = shared_ptr<SimInput>(new SimInputGMV(fileName, param));
#else
      EXCEPTION( "No support for GMV input file format." );
#endif
    } else if( fileName.find( ".unv") != std::string::npos ||
        fileName.find( ".unverg") != std::string::npos ||
        fileName.find( ".unvref") != std::string::npos ) {
#ifdef USE_UNV
      reader = shared_ptr<SimInput>(new SimInputUnv( fileName, param ));
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
      writer = shared_ptr<SimOutput>( new SimOutputGiD( baseName, gidNode ) );
#else
      EXCEPTION( "No support for GiD output file format." );
#endif
    } else if( fileName.find( ".gmv") != std::string::npos ) {
#ifdef USE_GMV_OUTPUT
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
      writer = shared_ptr<SimOutput>( new SimOutputGMV( baseName, gmvNode ) );
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
      writer = shared_ptr<SimOutput>( new SimOutputGmsh( baseName, gmshNode ) );
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
      writer =  shared_ptr<SimOutput>( new SimOutputHDF5( baseName, h5Node ) );
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
      writer =  shared_ptr<SimOutput>( new SimOutputXDMF( baseName, h5Node ) );
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
      writer =  shared_ptr<SimOutput>( new SimOutputUnv( baseName, unvNode ) );
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
    Grid * ptGrid = new GridCFS(dim);
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
             if( types[actMsStep] != BasePDE::HARMONIC ) {
               result = shared_ptr<BaseResult>( new Result<Double>() );
             } else {
               result = shared_ptr<BaseResult>( new Result<Complex>() );
             }
             result->SetEntityList( resEntities[iEntity] );
             result->SetResultInfo( actRes );
             results.Push_back( result );
             // Note: as the real values of saveBegin, saveInc and saveEnd are almost
             // never queried within an output format. we simply set saveBegin = 1,
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

  Double Diff( const std::string& inFile1,
               const std::string& inFile2,
               const std::string& outFile,
               bool normedtomax,
               bool isHistory,
               std::string& maxDiffResultName) {

       // obtain input reader for inFiles
       shared_ptr<SimInput> input1 = GetInputReader( inFile1 );
       shared_ptr<SimInput> input2 = GetInputReader( inFile2 );

       // check capabilities of input class
       bool printGridOnly = false;
       if( std::find( input1->GetCapabilities().begin(),
                      input1->GetCapabilities().end(),
                      SimInput::MESH_RESULTS )
           == input1->GetCapabilities().end() ) {
         std::cerr << "input files are only capable of handling meshes, not results!\n";
         exit(EXIT_FAILURE);
       }

       // read in mesh of input1
       input1->InitModule();
       UInt dim = input1->GetDim();
       Grid * ptGrid1 = new GridCFS(dim);
       input1->ReadMesh(ptGrid1);
       ptGrid1->FinishInit();

       // read in mesh of input2
       input2->InitModule();
       Grid * ptGrid2 = new GridCFS(dim);
       input2->ReadMesh(ptGrid2);
       ptGrid2->FinishInit();

       // obtain output writer
       shared_ptr<SimOutput> output;
       if( outFile != "" ) {
         output = GetOutputWriter( outFile );
         output->Init( ptGrid1, printGridOnly);
       }

       // obtain number of Sequence Steps and get analysis types
       std::map<UInt, BasePDE::AnalysisType> types;
       std::map<UInt, UInt> numSteps;
       input1->GetNumMultiSequenceSteps( types, numSteps, isHistory );

       std::cout << "\nFound " << types.size() << " sequence step(s) in '" << inFile1 << "'\n";
       std::map<UInt, BasePDE::AnalysisType> types2;
       std::map<UInt, UInt> numSteps2;
       input2->GetNumMultiSequenceSteps( types2, numSteps2, isHistory );
       std::cout << "\nFound " << types2.size() << " sequence step(s) in '" << inFile2 << "'\n";
       
       if(types.size() != types2.size()){
         std::cout << "'" << inFile1 << "' and '" << inFile2
            << "' have different number of sequence steps!\n";
         exit(EXIT_FAILURE);
       }

       // iterate over all Sequence Steps
       Double maxDiff = 0.0;
       std::map<UInt,UInt>::iterator it;
       for( it = numSteps.begin(); it != numSteps.end(); it++ ) {

         UInt actMsStep = it->first;
         std::cout << " Diffing sequence step " << actMsStep << std::endl
                   << "-------------------------\n\n";

         // get resulttypes
         StdVector<shared_ptr<ResultInfo> > infos;
         input1->GetResultTypes( actMsStep, infos, isHistory );

         StdVector<shared_ptr<BaseResult> > inResults1, inResults2, outResults;
         // stepnumbers, for which at least one result is defined
         std::map<UInt, Double> stepVals;
         // contains the stepnumbers/-values in which the particular result is
         // defined in
         std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;

         if( infos.GetSize() > 0 ){
           std::cout << "Performing diff on the following results:\n";
         }
         // iterate over all result types of input1
         for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {

           std::cout << "\t" << infos[iRes]->resultName << "\n";

           // get stepvalues
           shared_ptr<ResultInfo> actRes = infos[iRes];
           input1->GetStepValues( actMsStep, actRes,
                                  resultSteps[actRes], isHistory);
           stepVals.insert( resultSteps[actRes].begin(),
                            resultSteps[actRes].end() );

           // iterate over all regions
           StdVector<shared_ptr<EntityList> > regions;
           input1->GetResultEntities( actMsStep, infos[iRes],
                                      regions, isHistory );
           for( UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++ ) {
             // generate new result object and add it to output writer
             shared_ptr<BaseResult > inResult1, inResult2, outResult;
             if( types[actMsStep] != BasePDE::HARMONIC ) {
               inResult1 = shared_ptr<BaseResult>( new Result<Double>() );
               inResult2 = shared_ptr<BaseResult>( new Result<Double>() );
               outResult = shared_ptr<BaseResult>( new Result<Double>() );
             } else {
               inResult1 = shared_ptr<BaseResult>( new Result<Complex>() );
               inResult2 = shared_ptr<BaseResult>( new Result<Complex>() );
               outResult = shared_ptr<BaseResult>( new Result<Complex>() );
             }
             inResult1->SetEntityList( regions[iRegion] );
             inResult2->SetEntityList( regions[iRegion] );
             outResult->SetEntityList( regions[iRegion] );

             inResult1->SetResultInfo( infos[iRes] );

             inResult2->SetResultInfo( infos[iRes] );
             outResult->SetResultInfo( infos[iRes] );

             inResults1.Push_back( inResult1 );
             inResults2.Push_back( inResult2 );
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

         Vector<Double> maxResVec2;
         maxResVec2.Resize( inResults2.GetSize() );

         // For transient simulation find maximum amplitude over all timesteps
         if( types[actMsStep] != BasePDE::HARMONIC ) {

           // iterate over all results
           for( UInt iRes = 0; iRes < inResults2.GetSize(); iRes++) {
             
             if(numSteps[actMsStep] != numSteps2[actMsStep]){
               std::cout << "'" << inFile1 << "' has " << numSteps[actMsStep] << " and '" << inFile2
                  << "' has " << numSteps2[actMsStep] << " time steps!\n";
               exit(EXIT_FAILURE);
             }

             maxResVec2[iRes] = 0.0;
             // iterate over all time steps
             for( UInt iStep = 0; iStep < numSteps[actMsStep]; iStep++ ) {

               UInt actStepNum = iStep+1;
               // check if current result is defined within this step
               if( resultSteps[inResults2[iRes]->GetResultInfo()].find(actStepNum)
                   == resultSteps[inResults2[iRes]->GetResultInfo()].end() ) {
                 continue;
               }

               input2->GetResult( actMsStep, actStepNum, inResults2[iRes], isHistory );
               Vector<Double> & inVec2 =
                 dynamic_cast<Result<Double>& >(*inResults2[iRes]).GetVector();

               for( UInt i = 0; i<inVec2.GetSize(); i++ ) {
                 if( std::abs(inVec2[i]) > maxResVec2[iRes] )
                     maxResVec2[iRes] = std::abs(inVec2[i]);
               }
             }
             std::cout << "For result '" << inResults2[iRes]->GetResultInfo()->resultName
                       << "' maximum amplitude is: " << maxResVec2[iRes] << "\n";
           }
         }


         // notify writer
         if( output) {
           output->BeginMultiSequenceStep( actMsStep, types[actMsStep],
                                           numSteps[actMsStep] );
         }

         // iterate over all time/frequency steps
         for( UInt iStep = 0; iStep < numSteps[actMsStep]; iStep++ ) {

           // check, if current step contains any results
           if( stepVals.find( iStep+1) == stepVals.end() )
             continue;
           UInt actStepNum = iStep+1;
           Double actStepVal = stepVals[iStep+1];

           if( output) {
             output->BeginStep( actStepNum, actStepVal );
           }

           // iterate over all results
           for( UInt iRes = 0; iRes < inResults1.GetSize(); iRes++) {
             // check if current result is defined within this step
             if( resultSteps[inResults1[iRes]->GetResultInfo()].find(actStepNum)
                 == resultSteps[inResults1[iRes]->GetResultInfo()].end() ) {
               continue;
             }

             // obtain both result objects for current step
             input1->GetResult( actMsStep, actStepNum, inResults1[iRes], isHistory );
             input2->GetResult( actMsStep, actStepNum, inResults2[iRes], isHistory );

             // get number of dofs of result
             UInt numDofs = inResults1[iRes]->GetResultInfo()->dofNames.GetSize();

             // cast result objects, get vector and calculate difference vector
             if( types[actMsStep] != BasePDE::HARMONIC ) {
               Vector<Double> & inVec1 =
                 dynamic_cast<Result<Double>& >(*inResults1[iRes]).GetVector();
               Vector<Double> & inVec2 =
                  dynamic_cast<Result<Double>& >(*inResults2[iRes]).GetVector();
               Vector<Double> & outVec =
                 dynamic_cast<Result<Double>& >(*outResults[iRes]).GetVector();
               outVec.Resize( inVec1.GetSize() );

               // find maximum amplitude of inResult2
               for( UInt i = 0; i<inVec2.GetSize(); i++ ) {
                 if( std::abs(inVec2[i]) > maxResVec2[iRes])
                   maxResVec2[iRes] = std::abs(inVec2[i]);
               }

               // calculate difference entrywise
               outVec = inVec1 - inVec2;
               if (normedtomax == true) {
                 outVec /= maxResVec2[iRes];
               }

               // find maximum entry in difference vector
               for( UInt i = 0; i < outVec.GetSize(); i++ ) {
                 if( std::abs(outVec[i]) > maxDiff) {
                   maxDiff = std::abs(outVec[i]) ;
                   maxDiffResultName = inResults1[iRes]->GetResultInfo()->resultName;
                 }
               }

             } else {
               Vector<Complex> & inVec1 =
                 dynamic_cast<Result<Complex>& >(*inResults1[iRes]).GetVector();
               Vector<Complex> & inVec2 =
                 dynamic_cast<Result<Complex>& >(*inResults2[iRes]).GetVector();
               Vector<Complex> & outVec =
                 dynamic_cast<Result<Complex>& >(*outResults[iRes]).GetVector();
               outVec.Resize( inVec1.GetSize() );

               // find maximum amplitude of inResult2 in every frequency step
               Double maxRes2 = 0.0;
               for( UInt i = 0; i<inVec2.GetSize(); i++ ) {
                 if( std::abs(inVec2[i]) > maxRes2)
                   maxRes2 = std::abs(inVec2[i]);
               }

               Double aDiff, pDiff, aMax=0.0, aMin=0.0, pMax=0.0, pMin=0.0;
               Double rDiff, iDiff, rMax=0.0, iMax=0.0;

               // iterate over all dofs
               for (UInt dof = 0; dof<numDofs ; dof++) {
                 // iterate over number of entities
                 for( UInt i = 0; i<UInt(inVec2.GetSize()/numDofs); i++ ) {

                   // index to access entity 'i' of dof 'dof'
                   UInt actIndex = i * numDofs + dof;

                   // amplitude difference
                   if (normedtomax == true)
                     aDiff = ( std::abs(inVec1[actIndex]) - std::abs(inVec2[actIndex]) )/maxRes2;
                   else
                     aDiff = std::abs(inVec1[actIndex]) - std::abs(inVec2[actIndex]);

                   // phase difference in multiples of pi
                   pDiff = RadPhase(inVec1[actIndex]) - RadPhase(inVec2[actIndex]);

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
                   rDiff = std::abs( inVec1[actIndex].real() - inVec2[actIndex].real() )/maxRes2;
                   iDiff = std::abs( inVec1[actIndex].imag() - inVec2[actIndex].imag() )/maxRes2;
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
               }
               std::cout << "\n\tMaximum overall rel. difference = " << maxDiff << "\n\n";
             }

             // add result to output file
             if (output )
               output->AddResult( outResults[iRes] );
           }
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
       delete ptGrid1;
       delete ptGrid2;

       return maxDiff;

  } //Diff

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
    Grid * ptGrid = new GridCFS(dim);
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


  //==========================================================================================================
  // FFT SECTION
  //==========================================================================================================
#ifdef CFSTOOL_FFTW
  
  void updateChunkResult(std::vector< std::vector< std::vector<Complex> > >& fftVals,
                              UInt freqStep,const StdVector<std::string> & fNames,
                              UInt chunkSize,UInt numRes){
      std::cout << "Reading data from chunk files...";
      std::cout.flush();
      std::vector< std::vector<Complex> > inData;
      fftVals.clear();
      fftVals.resize(numRes);
      //now lets read in each chunk file, for each result
      StdVector<std::string>::const_iterator fIter = fNames.Begin();
      for(;fIter!=fNames.End();fIter++){
        std::fstream curFile((*fIter).c_str(),std::ios_base::in|std::ios_base::binary);
        boost::archive::binary_iarchive ia(curFile);
        ia >> inData;
        //analyse which result/chunk we are dealing with
        boost::char_separator<char> sep("_");
        boost::tokenizer< boost::char_separator<char> > tokenizer(*fIter, sep);
        std::vector<std::string> fTok;
        std::copy(tokenizer.begin(), tokenizer.end(),
                        std::back_inserter(fTok));
        UInt resNum = boost::lexical_cast<UInt>(fTok[1]);
        //UInt cNum = boost::lexical_cast<UInt>(*(tokenizer+2));
        for(UInt i=0;i<inData.size();i++){
          fftVals[resNum].push_back(std::vector<Complex>(chunkSize));
          //alternatively perform just a copy command...
          for(UInt j = 0;j<chunkSize && j+freqStep < inData[i].size();j++){
            fftVals[resNum][fftVals[resNum].size()-1][j] = inData[i][freqStep+j];
          }
        }
      }
      //open the chunkfile and just assign it to fftVals
      std::cout << "done" << std::endl;
    }

  void updateChunkResultone(std::vector< std::vector< std::vector<Double> > >& fftVals,
                              UInt freqStep,const StdVector<std::string> & fNames,
                              UInt chunkSize,UInt numRes){
      std::vector< std::vector<Double> > inData;
      fftVals.clear();
      fftVals.resize(numRes);
      //now lets read in each chunk file, for each result
      StdVector<std::string>::const_iterator fIter = fNames.Begin();
      for(;fIter!=fNames.End();fIter++){
        std::fstream curFile((*fIter).c_str(),std::ios_base::in|std::ios_base::binary);
        boost::archive::binary_iarchive ia(curFile);
        ia >> inData;
        //analyse which result/chunk we are dealing with
        boost::char_separator<char> sep("_");
        boost::tokenizer< boost::char_separator<char> > tokenizer(*fIter, sep);
        std::vector<std::string> fTok;
        std::copy(tokenizer.begin(), tokenizer.end(),
                        std::back_inserter(fTok));
        UInt resNum = boost::lexical_cast<UInt>(fTok[1]);
        //UInt cNum = boost::lexical_cast<UInt>(*(tokenizer+2));
        for(UInt i=0;i<inData.size();i++){
          fftVals[resNum].push_back(std::vector<Double>(chunkSize));
          //alternatively perform just a copy command...
          for(UInt j = 0;j<chunkSize && j+freqStep < inData[i].size();j++){
            fftVals[resNum][fftVals[resNum].size()-1][j] = inData[i][freqStep+j];
          }
        }
      }
      //open the chunkfile and just assign it to fftVals
    }

  /**
   * calcFFT calclulates the frequency transformation of a given input file
   * @param inFile The file which carries the data (e.g. acoustic pressure , mechanical
   * displacement)
   * @param meshFile The mesh file to associate each value with a point in space
   * Import not all points should be included into the averagin.
   * @param regs various programm options
   * @param iFFT performs an inversrse fourier transform at the end to obtain transient data again
   * @param window performs an additional windowing of time data before computing the fft
   * @return void
   */
  void PerformFFT(const std::string& inFile,
                    const std::string& outFile,
                    std::vector<std::string> regs,
                    bool iFFT,bool window){


    Double maxMem = 0;
    UInt numNodes = 0;

    std::vector<std::string> chunkFileNames;

    PtrParamNode maxMemNode = param->Get("maxMemory",ParamNode::PASS);
    if(maxMemNode){
      try{
        maxMem = boost::lexical_cast<Double>(maxMemNode->As<std::string>());
      } catch(bad_lexical_cast &) {
      EXCEPTION("Error converting maxMem parameter");
      }
    }else{
      std::cerr << "No maxMemory parameter specified. Trying to guess from system and use 20% of it..." << std::endl;
      long pages = sysconf(_SC_PHYS_PAGES);
      long page_size = sysconf(_SC_PAGE_SIZE);
      long mem = pages * page_size;
      //compute in MB
      maxMem = mem / (1024*1024) * 0.2;
      std::cerr << "\t Maximum available memory in MB for FFT: " << maxMem << std::endl;
    }

    //================================================================================
    // INITIALIZATION
    //================================================================================
    shared_ptr<SimInput> input = GetInputReader( inFile );
    input->InitModule();
    UInt dim = input->GetDim();
    Grid * ptGrid = new GridCFS(dim);
    input->ReadMesh(ptGrid);

    ptGrid->FinishInit();

    // obtain output writer
    if( outFile == "" ) {
      Exception("Outputfile is needed for FFT");
    }
    shared_ptr<SimOutput> output;
    output = GetOutputWriter( outFile );
    output->Init( ptGrid, false);


    //=================================================================================
    //GATHER INFORMATION ABOUT THE INPUT RESULT FILE
    //=================================================================================

    std::map<UInt, BasePDE::AnalysisType> types;
    std::map<UInt, UInt> numSteps;
    std::map<UInt, UInt>::iterator it;

    input->GetNumMultiSequenceSteps( types, numSteps, false );

    std::cout << "\nFound " << types.size() << " sequence step(s) in '" << inFile << "'\n";

    for( it = numSteps.begin(); it != numSteps.end(); it++ ) {
      StdVector< std::pair<UInt,UInt> > nodeChunks;
      UInt freqChunkSize;
      UInt numFreqChunks = 1;
      UInt timeChunkSize;
      UInt numTimeChunks = 1;

      StdVector<shared_ptr<BaseResult> > inResults, outResults;
      std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
      StdVector<shared_ptr<ResultInfo> > infos;
      UInt actMsStep = it->first;
     // std::string param_analysisList = param->Get("analysisList")->As<std::string>();

      //if(param_analysisList == "TRANSIENT"){
      //std::cout<< BasePDE::HARMONIC<<endl;
      /*if( types[actMsStep] == BasePDE::HARMONIC) {
        std::cerr << "FFT does not make sense on anything but transient data.... going to continue" << std::endl;
        continue;
        cout<< "what the hell is this"<<endl;
      }*/
      //}
      std::cout << " Performing FFT of Sequence step " << actMsStep << std::endl
                << "-------------------------\n\n";

      // get resulttypes
      input->GetResultTypes( actMsStep, infos, false );

      if( infos.GetSize() > 0 ){
        std::cout << "Performing FFT on the following results:\n";
      }else{
        std::cout << "Found no results for Sequence Step " << actMsStep << std::endl;
        continue;
      }

      // iterate over all result types of input
      for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {
        StdVector<shared_ptr<EntityList> > regions;
        shared_ptr<ResultInfo> actRes = infos[iRes];

        input->GetResultEntities( actMsStep, infos[iRes],regions, false );
        input->GetStepValues( actMsStep, actRes,resultSteps[actRes], false);
        std::cout << "\t" << infos[iRes]->resultName << "\n";

        //iterate over each region the result is definied on
        for( UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++ ) {
          std::string regName = regions[iRegion]->GetName();
          //chekc if we need this region
          if(std::find(regs.begin(),regs.end(), regName) == regs.end() &&
             std::find(regs.begin(),regs.end(), "allRegions") == regs.end()){
            continue;
          }
          UInt regNodes = 0;
          if(regions[iRegion]->GetType() == EntityList::NODE_LIST){
            regNodes = regions[iRegion]->GetSize();
            numNodes += regNodes;
          }else{
            std::cerr << "Got entity list which is not a node list. Cannot cope with that" << std::endl;
            continue;
          }
          //=====================================
          //Determine Nodal Chunks Information
          //=====================================
          //calc memory requirements in MB
          Double memRequired = regNodes * sizeof(Double) * resultSteps[actRes].size()/1048576.0;
//maxMem=10;
          //now determine the number of chunks
          if(memRequired > maxMem){
            UInt numChunks =  std::ceil(memRequired/maxMem);
            UInt chunksize = std::ceil(regNodes/numChunks);
            nodeChunks.Push_back(std::pair<UInt,UInt>(numChunks,chunksize));
            //std::cout << numChunks << endl;
          }else{
            UInt numChunks =  1;
            UInt chunksize = regNodes;
            nodeChunks.Push_back(std::pair<UInt,UInt>(numChunks,chunksize));
          }

          std::cout << "\t\t Adding region " << regName  << " with " << regNodes << " nodes -> "\
              "process in " << nodeChunks[iRes].first << " chunks" <<std::endl;


          if(iFFT==true){

          shared_ptr<BaseResult > inResult = shared_ptr<BaseResult>( new Result<Double>() );
          shared_ptr<BaseResult > outResult = shared_ptr<BaseResult>( new Result<Double>() );

          //define input regions and results
          inResult->SetResultInfo( actRes );
          inResult->SetEntityList( regions[iRegion] );
          inResults.Push_back( inResult );

          //define the output
          outResult->SetEntityList( regions[iRegion] );
          outResult->SetResultInfo( infos[iRes] );
          outResult->GetResultInfo()->complexFormat = REAL_IMAG;
          outResults.Push_back( outResult );

          }else{

            shared_ptr<BaseResult > inResult = shared_ptr<BaseResult>( new Result<Double>() );
            shared_ptr<BaseResult > outResult = shared_ptr<BaseResult>( new Result<Complex>() );

              //define input regions and results
              inResult->SetResultInfo( actRes );
              inResult->SetEntityList( regions[iRegion] );
              inResults.Push_back( inResult );

              //define the output
              outResult->SetEntityList( regions[iRegion] );
              outResult->SetResultInfo( infos[iRes] );
              outResult->GetResultInfo()->complexFormat = REAL_IMAG;
              outResults.Push_back( outResult );

          }

        }
      }

      StdVector<std::string> fileNames;

      //==============================================================================
      //READ IN TIMEVALUES AND PERFORM FFT (CHUNKED)
      //==============================================================================
      UInt numFrequencies = 0;
      UInt numTimes = 0;

      for(UInt aRes = 0; aRes < inResults.GetSize(); aRes++){
       // std::cout << "----> Processing results...." << std::endl;
      //std::cout<<inResults.GetSize()<<endl;
        UInt actNumChunks = nodeChunks[aRes].first;
        //std::cout<<actNumChunks<<endl;
        UInt actChunkSize = nodeChunks[aRes].second;
        //std::cout<< actChunkSize<<endl;
        UInt numRegNodes = inResults[aRes]->GetEntityList()->GetSize();
        UInt remainder = (numRegNodes - (actNumChunks-1)*actChunkSize);
        UInt offset = 0;
        //READ IN TIMESTEPS
        UInt numSteps = resultSteps[inResults[aRes]->GetResultInfo()].size();

        //=======================================================================
        // Guess or just read in cut on/off freqeuncies and determine timestep
        //=======================================================================
        Double ts1 = resultSteps[inResults[aRes]->GetResultInfo()][1];
        Double ts2 = resultSteps[inResults[aRes]->GetResultInfo()][2];
        Double tStep =abs(ts2-ts1);
        //total runtime
        Double fs = 1/tStep;
        Double df = fs/numSteps;
        std::cout << std::endl <<  "Time freqeuncy information for Result " << inResults[aRes]->GetResultInfo()->resultName << std::endl;
        std::cout << "\t Found time step size to be :" << tStep << std::endl;
        std::cout << "\t Frequency resolution :" << df << std::endl;

        std::vector<std::string> fRange;
        // Initialize vector with output fields
        Double fMin=0;
        Double fMax=0;

        typedef boost::tokenizer< boost::char_separator<char> > Tok;
        boost::char_separator<char> sep(";| ");

        //lets read in frequency range the user specified
        PtrParamNode freqRangeNode = param->Get("frequencyRange",ParamNode::PASS);
        if(freqRangeNode){
          std::string freq = freqRangeNode->As<std::string>();
          Tok tokenizer(freq, sep);
          std::copy(tokenizer.begin(), tokenizer.end(),
              std::back_inserter(fRange));
          if(fRange.size() != 2){
              EXCEPTION("frequency range parameter invalid "  \
                  << "specify 'min max'");
          }

          try {
            fMin = boost::lexical_cast<Double>(fRange[0]);
            fMax = boost::lexical_cast<Double>(fRange[1]);
          } catch(bad_lexical_cast &) {
             EXCEPTION("Error converting freqeuncy range to double parameter");
          }

        }else{
          std::cerr << "\n\tNo frequency range specified. Trying to guess from time data..." << std::endl;
          //compute ampling freqeuncy and divide by two
          Double fs = 1/tStep;
          fMax = fs/2;
        }
        fMin = 0;
        fMax = 10000;

        std::cout << "\t Cut on freqeuncy= " << fMin << " Hz" <<  std::endl;
        std::cout << "\t Cut off freqeuncy = " << fMax << " Hz" << std::endl << std::endl;



        for(UInt aChunk = 0; aChunk < actNumChunks;aChunk++){
          std::cout << "Processing '" << inResults[aRes]->GetResultInfo()->resultName \
                    << "' on region '" << inResults[aRes]->GetEntityList()->GetName() << "' chunk " << aChunk+1 << " :" << std::endl;
          std::cout.flush();


          UInt curNumNodes = (aChunk==actNumChunks-1)? remainder : actChunkSize;
          Matrix<Double> tVals(numSteps,curNumNodes );
          boost::progress_display show_progress(numSteps);
          for( UInt iStep = 0; iStep < numSteps; iStep++ ) {
            UInt actStepNum = iStep+1;
            input->GetResult( actMsStep, actStepNum, inResults[aRes], false );
            Vector<Double> & inVec = dynamic_cast<Result<Double>& >(*inResults[aRes]).GetVector();
            for(UInt i=0; i< curNumNodes;i++){
              tVals[iStep][i]=inVec[offset+i];
            }
            ++show_progress;
          }

          offset +=curNumNodes;
          //FFT::WindowType WINDOW;



          FFT::WindowType windows;
          PtrParamNode windowNode = param->Get("windowsList",ParamNode::PASS);
          std::string param_windowsList;
          if(windowNode){
            param_windowsList = windowNode->As<std::string>();
          }else{
            std::cout << "no Window specified... using RECTANGULAR." << std::endl;
            param_windowsList = "RECTANGULAR";
          }
          //PERFORM FFT
          FFT ff;
          std::vector< std::vector<Complex> > fVals;
          Matrix<Double> wtVals;
          std::vector< std::vector<Complex> > wfVals;
         // Vector<Complex>   winfftVals;
          StdVector<Double> freqs;
          std::vector< std::vector<Double> > ifft;
          StdVector<Double> times;
          //std::vector< std::vector<Complex> > wtempVals;


          ff.FT(tVals,fs,fMin,fMax,fVals,freqs);

          //ff.FT(tVals,ts1,fVals,freqs);

       /*   for(UInt i=0; i<curNumNodes; i++){
            for( UInt j = 0; j < 1000; j++ ) {
               if(freqs[j]>=fMin && freqs[j]<fMax){
            cout<<freqs[j]<<"\n"<<fVals[i][j]<<endl;


               }
            }

          }*/

          if(param_windowsList == "HANNING"){
            ff.WIN(windows = FFT::HANNING,tVals,wtVals);
          }else if(param_windowsList == "RECTANGULAR"){
            ff.WIN(windows = FFT::RECTANGULAR,tVals,wtVals);
          }else if(param_windowsList == "BLACKMAN"){
            ff.WIN(windows = FFT::BLACKMAN,tVals,wtVals);
          }else if(param_windowsList == "BLACKMANHARRIS"){
            ff.WIN(windows = FFT::BLACKMANHARRIS,tVals,wtVals);
          }else if(param_windowsList == "NUTTALL"){
            ff.WIN(windows = FFT::BLACKMANHARRIS,tVals,wtVals);
          }else if(param_windowsList == "BLACKMANNUTTALL"){
            ff.WIN(windows = FFT::BLACKMANHARRIS,tVals,wtVals);
          }


/*
          if(param_windowsList == "HANNING"){
          ff.WIN(windows = FFT::HANNING,tVals,winfftVals);
          }else if(param_windowsList == "RECTANGULAR"){
           ff.WIN(windows = FFT::RECTANGULAR,tVals,winfftVals);
           }else if(param_windowsList == "BLACKMAN"){
           ff.WIN(windows = FFT::BLACKMAN,tVals,winfftVals);
           }else if(param_windowsList == "HAMMING"){
           ff.WIN(windows = FFT::HAMMING,tVals,winfftVals);
           }*/


           //ff.WINFFT(tVals,wtVals,fMin,fMax,wfVals,freqs);


           //ff.WINMUL(tVals,fVals,wfVals,freqs,fMin,fMax,ts1,wtempVals,winfftVals);

          //ff.IFT(tVals,wfVals,ifft,tStep,freqs,times);

          //ff.IFFT(tVals,wfVals,ifft,tStep,freqs,times);

          //WRITE FREQUENCY SCALE FOR FIRST RESULT ONLY
          //IF WE HAVE DIFFERENT TIMESTEPVALUES IN THE RESULTS WE THROW AN ERROR



          if(iFFT==true){

            std::fstream scaleFile("scale.ifft", std::ios_base::out | std::ios_base::trunc);
            scaleFile << times.GetSize() << std::endl;
              scaleFile << times << std::endl;
              scaleFile.close();
            numTimes = times.GetSize();

          }else{
          if(aRes == 0){
            std::fstream scaleFile("scaleone.fft", std::ios_base::out | std::ios_base::trunc);
            scaleFile << freqs.GetSize() << std::endl;
            scaleFile << freqs << std::endl;
            scaleFile.close();
            numFrequencies = freqs.GetSize();
          }else{
            if(numFrequencies != freqs.GetSize()){
              EXCEPTION("The FFTs from different results have different frequency scalings. This cannot be handled right now");
            }
          }
          }
          //WRITE THE CHUNK TO A FILE
          std::stringstream fName;
          if(iFFT==true){
          fName << "values_" << aRes << "_" << aChunk << ".ifft";

          }
          else{
            fName << "values_" << aRes << "_" << aChunk << ".fft";
          }
          std::fstream chunkFile(fName.str().c_str(),std::ios_base::out|std::ios_base::binary|std::ios_base::trunc);
          if ( chunkFile.good() ) {
            try{
              if(iFFT == true){
              boost::archive::binary_oarchive oa(chunkFile);
              oa << ((const std::vector< std::vector<Double> > &) ifft);
              std::cout<< "ifft is successful" << std::endl;
              }

              if(window == true){
              boost::archive::binary_oarchive oa(chunkFile);
              oa << ((const std::vector< std::vector<Complex> > &) wfVals);
              std::cout<< "window is successful"<<std::endl;

              }
              else{
               std::cout << "Buffering result...";
               std::cout.flush();
                    boost::archive::binary_oarchive oa(chunkFile);
                    oa << ((const std::vector< std::vector<Complex> > &) fVals);
                    std::cout << "done."<<std::endl;
                    //cout<<"no window"<<endl;
              }
            }catch(std::exception &ex){
              EXCEPTION("The following problem occurred while trying to "\
                        << "write conservative interpolation weights to '"\
                        << fName.str() << "': " << ex.what());
            }
            chunkFile.close();
            chunkFileNames.push_back(fName.str());
          }else{
          WARN("An error occured while writing the chunk file ");
          }
          std::cout.flush();
        }


       if(iFFT == true){
        output->RegisterResult( outResults[aRes], 1, 1,
                numTimes,
                                false );
        }
        else{
            output->RegisterResult( outResults[aRes], 1, 1,
                                    numFrequencies,
                                    false );
        }

      }

      //std::cout<<"done"<< std::endl;

      if(iFFT==true){
      //==============================================================================================
      //WRITE OUT IFFT CHUNKED ASSUPTION
      //==============================================================================================
      //create datastructure for IFFT values
        std::vector< std::vector< std::vector<Double> > > curIFFTVals;

      //Read in Frequency Scale
      std::fstream sFile("scale.ifft",std::ios_base::in);
      std::string curline = "";
      getline(sFile,curline);
      UInt numTsteps = boost::lexical_cast<UInt>(curline);

      //Determine memory requirement for one IFFT step in MB
      Double memPerTime = sizeof(Double) * numNodes / 1048576.0;
      if(memPerTime>maxMem){
        EXCEPTION("Not even a single frequency step fits into memory. Incease the maximum Memory to at least"\
                 << memPerTime << "MBytes" );
      }
      if(memPerTime*numTsteps > maxMem){
        numTimeChunks = std::ceil(memPerTime*(Double)numTsteps/maxMem);
        timeChunkSize = std::ceil((Double)numTsteps/(Double)numTimeChunks);
      }else{
        numTimeChunks = 1;
        timeChunkSize = numTsteps;
      }

      //prepareChunkResults(curFFTVals,freqChunkSize,numFreqChunks,outResults.GetSize());

      output->BeginMultiSequenceStep( actMsStep, BasePDE::TRANSIENT,
                                      numSteps[actMsStep] );

      UInt fCounter = 0;
      updateChunkResultone(curIFFTVals,0,chunkFileNames, timeChunkSize,outResults.GetSize());

      for(UInt i = 0 ; i<numTsteps ; i++){
        if(fCounter >= curIFFTVals[0][0].size()){
          updateChunkResultone(curIFFTVals,i,chunkFileNames, timeChunkSize,outResults.GetSize());
          fCounter = 0;
        }
        Double curStep = 0;
        getline(sFile,curline);
        boost::erase_all(curline, " ");
        curStep = boost::lexical_cast<Double>(curline);
        //std::cout << "writing frequency step " << i << "->" << curStep << "Hz" << std::endl;
        //std::cout.flush();

        if(curStep>0){
        output->BeginStep( i, curStep );
        std::cout << "writing time step " << i << "->" << curStep << "Sec" << std::endl;
        std::cout.flush();
        for(UInt aRes = 0; aRes < outResults.GetSize(); aRes++){

          Vector<Double> & outVec =  dynamic_cast<Result<Double>& >(*outResults[aRes]).GetVector();

          UInt actNumNodes = outResults[aRes]->GetEntityList()->GetSize();
          outVec.Resize(actNumNodes);
          for(UInt actN = 0;actN< actNumNodes ;actN++){
              outVec[actN] = curIFFTVals[aRes][actN][fCounter];
          }
          output->AddResult( outResults[aRes] );
        }
        fCounter++;
        output->FinishStep();
      }
       }
      } else{
      //==============================================================================================
      //WRITE OUT FFT CHUNKED. ASSUPTION: at least one frequency fits completely into memory
      //==============================================================================================
      //create datastructure for FFT values
      std::vector< std::vector< std::vector<Complex> > > curFFTVals;

      //Read in Frequency Scale
      std::fstream sFile("scaleone.fft",std::ios_base::in);
      std::string curline = "";
      getline(sFile,curline);
      UInt numFsteps = boost::lexical_cast<UInt>(curline);



      //Determine memory requirement for one FFT step in MB
      Double memPerFreq = sizeof(Complex) * numNodes / 1048576.0;
      if(memPerFreq>maxMem){
        EXCEPTION("Not even a single frequency step fits into memory. Incease the maximum Memory to at least"\
                 << memPerFreq << "MBytes" );
      }
      if(memPerFreq*numFsteps > maxMem){
        numFreqChunks = std::ceil(memPerFreq*(Double)numFsteps/maxMem);
        freqChunkSize = std::ceil((Double)numFsteps/(Double)numFreqChunks);
      }else{
        numFreqChunks = 1;
        freqChunkSize = numFsteps;
      }


      //prepareChunkResults(curFFTVals,freqChunkSize,numFreqChunks,outResults.GetSize());

      output->BeginMultiSequenceStep( actMsStep, BasePDE::HARMONIC,
                                      numSteps[actMsStep] );

      UInt fCounter = 0;
      updateChunkResult(curFFTVals,0,chunkFileNames,freqChunkSize,outResults.GetSize());
      for(UInt i = 0 ; i<numFsteps; i++){
        //cout<<curFFTVals[0][0].size()<<endl;
        //cout<<fCounter<<endl;
        if(fCounter >= curFFTVals[0][0].size()){
          updateChunkResult(curFFTVals,i,chunkFileNames,freqChunkSize,outResults.GetSize());
          fCounter = 0;
        }
        Double curStep = 0;
        getline(sFile,curline);
        boost::erase_all(curline, " ");
        curStep = boost::lexical_cast<Double>(curline);


        output->BeginStep( i, curStep );
        std::cout << "writing frequency step " << i << "->" << curStep << "Hz" << std::endl;
        std::cout.flush();

      //

        for(UInt aRes = 0; aRes < outResults.GetSize(); aRes++){

          Vector<Complex> & outVec =  dynamic_cast<Result<Complex>& >(*outResults[aRes]).GetVector();

          UInt actNumNodes = outResults[aRes]->GetEntityList()->GetSize();
          outVec.Resize(actNumNodes);
          //cout<<"actno.nodes"<<actNumNodes<<endl;
          for(UInt actN = 0;actN< actNumNodes ;actN++){
              outVec[actN] = curFFTVals[aRes][actN][fCounter];
              //cout<< outVec[actN]<<endl;

           }
          output->AddResult( outResults[aRes] );
        }
        fCounter++;
        output->FinishStep();

       }
      }
       output->FinishMultiSequenceStep();
     }
     output->Finalize();
     std::cout << "\nOutput successfully written to " << outFile << std::endl;
     //================================================================================
     // DELETE ALL TEMPORARY FILES
     //================================================================================
     if(iFFT==true){
     chunkFileNames.push_back("scale.ifft");
     }
     else{
       chunkFileNames.push_back("scaleone.fft");
     }
     for(UInt i=0;i<chunkFileNames.size();i++){
       if(fs::exists(chunkFileNames[i])){
         fs::remove(chunkFileNames[i]);
       }
     }
   }

#else //ifdef CFSTOOL_FFTW
  void PerformFFT(const std::string& inFile,
                    const std::string& outFile,
                    std::vector<std::string> regs,
                    bool iFFT,bool window){
    EXCEPTION("This executable does not support FFT. "
        << "If you need this feature, recompile with CFSTOOL_FFTW=ON.");
  }
#endif
  //===================================================================================
  // END OF FFT SECTION
  //===================================================================================



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
    // domain->nodeList->nodes->(name, list (coordSysId, freeCoord->(comp, start, stop, inc)))
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

  std::cout << std::endl
            << "============================================================"
            << "===========" << std::endl;
  std::cout << " CFSTOOL - File Conversions/Comparisons for CFS++" << std::endl << std::endl
            << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'"
            << " (rev " << CFS_SUBVERSION_REV << ")" << std::endl
            << " compiled " << __DATE__
            << " as " << CMAKE_BUILD_TYPE << std::endl;
  std::cout << "============================================================"
            << "==========="
            << std::endl << std::endl;

  std::string maxDiffResultName; // that's chaos man!
  std::string infoFileName;
  try
  {
    param = PtrParamNode(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));

    ParamsInit(argc, argv);

    // Switch this flag tc true for debugging
    if (param->Get("forceSegFault")->As<bool>())
    {
      Exception::segfault_ = true;
    } else {
      Exception::segfault_ = false;
    }

      Exception::segfault_ = true;

    std::string param_mode = param->Get("mode")->As<std::string>();

    // get filenames from parameter
    std::string inputFile = param->Get("inputFile")->As<std::string>();
    std::string compareFile = param->Get("compareFile")->As<std::string>();
    std::string outputFile = param->Get("outputFile")->As<std::string>();

    std::cout << "inputFile " << inputFile << std::endl;
    std::cout << "compareFile " << compareFile << std::endl;
    std::cout << "outputFile " << outputFile << std::endl;

    // Initialize vector with output fields
    UInt num_files = 0;
    if (inputFile != "")
    {
      ++num_files;
      if (compareFile != "")
      {
        ++num_files;
        if (outputFile != "")
        {
          ++num_files;
        }
      }
    }
    if (inputFile != "")
    {
      // This is necessary to run, but I do not know what it is for
      info = PtrParamNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT));
      infoFileName = param_mode + "_" + inputFile + ".info.xml";
      info->SetName("cfsInfo"); 
    }

    if (param_mode == "calcAverage")
    {
      if (outputFile != "")
      {
        EXCEPTION( "Two many arguments, please only provide two files. (in- and output file)" );
      }
      outputFile = compareFile;
      if (num_files != 2)
      {
        EXCEPTION( "Please provide a reference file and output File" );
      }
      CFSTool::calcAverage(inputFile, outputFile);
    } else if (param_mode == "fft" || param_mode == "ifft" || param_mode == "window") {
      if( outputFile == ""){
        std::cerr << "no output-file specified for FFT mode. Trying to guess from compare file..." << std::endl;
        outputFile = compareFile;
        if( outputFile == ""){
          EXCEPTION( "FAILED: Please provide <inFile> AND <outFile>" );
        }else{
          std::cout << "Output file is: " << outputFile << std::endl;
        }
      }



      //std::vector<std::string> fRange;
      std::vector<std::string> regList;
      typedef boost::tokenizer< boost::char_separator<char> > Tok;
      boost::char_separator<char> sep(";| ");
      // Initialize vector with output fields

   /* Double fMin;
      Double fMax;
      std::string freq = param->Get("frequencyRange")->As<std::string>();
      Tok tokenizer(freq, sep);
      std::copy(tokenizer.begin(), tokenizer.end(),
                std::back_inserter(fRange));
      if(fRange.size() != 2){
        EXCEPTION("frequency range parameter invalid "  \
            << "specify 'min max'");
      }


      try {
         fMin = boost::lexical_cast<Double>(fRange[0]);
         fMax = boost::lexical_cast<Double>(fRange[1]);
      } catch(bad_lexical_cast &) {
        EXCEPTION("Error converting freqeuncy range to double parameter");
      }*/

      PtrParamNode regNode = param->Get("regionList",ParamNode::PASS);
      std::string regions;
      if(regNode){
        regions = regNode->As<std::string>();
        if(regions == ""){
          regions = "allRegions";
        }
      }else{
        regions = "allRegions";
      }

      Tok tokenizer1(regions, sep);
      std::copy(tokenizer1.begin(), tokenizer1.end(),
                std::back_inserter(regList));

      if(param_mode == "fft"){
        CFSTool::PerformFFT( inputFile, outputFile, \
                            regList, false,false);
      }else if(param_mode == "ifft") {
        CFSTool::PerformFFT( inputFile, outputFile, \
                             regList,true,false);
       }
      else if(param_mode == "window") {
              CFSTool::PerformFFT( inputFile, outputFile, \
                                   regList, false,true);
      }
    } else if (param_mode == "convert") {
      if (outputFile != "")
      {
        EXCEPTION( "Two many arguments, please only provide two files. (in- and output file)" );
      }
      outputFile = compareFile;
      if (num_files != 2)
      {
        EXCEPTION( "Please provide <infFile> and <outFile>" );
      }
      CFSTool::Convert( inputFile, outputFile );
    } else if (param_mode == "scalardiff") {
      Double tolerance = param->Get("eps")->As<Double>();
      if (num_files != 2)
      {
        EXCEPTION( "Please provide <inFile1> and <inFile2>" );
      }
      Double maxDiffMesh = 0.0, maxDiffHist = 0.0;
      std::cout << "Checking for mesh results:\n"
        << "==========================";
      maxDiffMesh = CFSTool::Diff( inputFile, compareFile, "", \
                                  true, false, maxDiffResultName);
      std::cout << "Checking for history results:\n"
        << "=============================";
      maxDiffHist = CFSTool::Diff( inputFile, compareFile, "", \
                                  true, true, maxDiffResultName );
      Double maxDiff = std::max( maxDiffMesh, maxDiffHist );
      if( maxDiff > tolerance ) {
        std::cout << "'" << inputFile << "' and '" << compareFile
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
        EXCEPTION( "Please provide <inFile1>, <inFile2> and <outFile>" );
      }
      // Double maxDiff = 0.0; // TODO: Unused variable maxDiff
      // maxDiff = CFSTool::Diff( inputFile, compareFile, outputFile, 

      CFSTool::Diff( inputFile, compareFile, outputFile, \
                     false, false, maxDiffResultName);
    } else if (param_mode == "meshdiffnormed") {
      if (num_files != 3)
      {
        EXCEPTION( "Please provide <inFile1>, <inFile2> and <outFile>" );
      }
      // Double maxDiff = 0.0; // TODO: Unused variable maxDiff
      // maxDiff = CFSTool::Diff( inputFile, compareFile, outputFile, 

      CFSTool::Diff( inputFile, compareFile, outputFile, 
                     true, false, maxDiffResultName);
    } else {
      EXCEPTION( "No such mode: " << param_mode <<". See help for available modes" );
      return EXIT_FAILURE;
    }
  } catch(std::exception& ex) {
    std::cerr << "The following error occured:\n" << ex.what();
    if (info != NULL)
    {
      info->Get(ParamNode::ERROR)->SetValue(ex.what());
      info->ToFile(infoFileName);
    }
    std::cerr << "The following error occured during program execution:\n\n" << ex.what();
    return -1;
  }
  info->ToFile(infoFileName);
  return 0;
}
