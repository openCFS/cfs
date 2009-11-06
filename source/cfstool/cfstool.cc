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
#include <def_use_gmv.hh>
#include <def_use_unv.hh>
#include <def_use_ansysrst.hh>

#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include "ParamsInit.hh"
#include "General/environment.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/simInput.hh"
#include "DataInOut/simOutput.hh"
#include "Domain/GridCFS/grid_cfs.hh"

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/simInputMESH.hh"
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
#endif

#ifdef USE_GIDPOST
#include "DataInOut/SimInOut/GiD/simOutGiD.hh"
#endif

#ifdef USE_UNV
#include "DataInOut/SimInOut/Unverg/simInputUnv.hh"
#include "DataInOut/SimInOut/Unverg/simOutUnv.hh"
#endif

#ifdef USE_ANSYSRST
#include "DataInOut/SimInOut/AnsysRST/simOutputRST.hh"
#endif

#include "DataInOut/SimInOut/TextOutput/textSimOutput.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"

using namespace CoupledField;

namespace CFSTool {

  void setFreeCoord(std::string coordSysId="default",
      std::string node_name="averageDomain");
  inline void setParamNode(ParamNode *paramNode, std::string name, std::string value,
      StdVector<ParamNode*>* children=NULL);

  shared_ptr<SimInput> GetInputReader(const std::string& fileName)
  {
    // determine suffix of fileName
    shared_ptr<SimInput> reader;

    if (!fs::exists( fileName ))
      EXCEPTION( "\nFile '" << fileName << "' does not exist!");

    if( fileName.find( ".mesh") != std::string::npos ) {
#ifdef USE_MESH
      reader = shared_ptr<SimInput>(new SimInputMESH( fileName, NULL ) );
#else
      EXCEPTION( "No support for MESH input file format." );
#endif
    } else if( fileName.find( ".h5") != std::string::npos ) {
#ifdef USE_HDF5
      reader = shared_ptr<SimInput>(new SimInputHDF5(fileName, param));
#else
      EXCEPTION( "No support for HDF5 input file format." );
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
      ParamNode * gidNode = new ParamNode(false);
      ParamNode * binary = new ParamNode(false);
      binary->SetName( "binaryFormat");
      if( fileName.find( ".bin") != std::string::npos ) {
        binary->SetValue( "yes");
      } else {
        binary->SetValue( "false");
      }
      gidNode->GetChildren().Push_back(binary);
      writer = shared_ptr<SimOutput>( new SimOutputGiD( baseName, gidNode ) );
#else
      EXCEPTION( "No support for GiD output file format." );
#endif
    } else if( fileName.find( ".gmv") != std::string::npos ) {
#ifdef USE_GMV_OUTPUT
      baseName = std::string(fileName, 0, fileName.find(".gmv"));
      ParamNode * gmvNode = new ParamNode(false);
      ParamNode * binary = new ParamNode(false);
      binary->SetName("binaryFormat");
      binary->SetValue( "yes" );
      ParamNode * fixedGrid = new ParamNode(false);
      fixedGrid->SetName("fixedGrid");
      fixedGrid->SetValue( "yes" );
      gmvNode->GetChildren().Push_back(binary);
      gmvNode->GetChildren().Push_back(fixedGrid);
      writer = shared_ptr<SimOutput>( new SimOutputGMV( baseName, gmvNode ) );
#else
      EXCEPTION( "No support for GMV output file format." );
#endif
    } else if(fileName.find( ".h5") != std::string::npos) {
#ifdef USE_HDF5
      baseName = std::string(fileName, 0, fileName.find(".h5"));
      ParamNode * h5Node = new ParamNode(false);
      ParamNode * eFiles = new ParamNode(false);
      eFiles->SetName("externalFiles");
      eFiles->SetValue( "false" );
      h5Node->GetChildren().Push_back(eFiles);
      writer =  shared_ptr<SimOutput>( new SimOutputHDF5( baseName, h5Node ) );
#else
      EXCEPTION( "No support for HDF5 output file format." );
#endif
    } else if(fileName.find( ".rst") != std::string::npos) {
#ifdef USE_ANSYSRST
      baseName = std::string(fileName, 0, fileName.find(".rst"));
      ParamNode * rstNode = new ParamNode(false);
      writer =  shared_ptr<SimOutput>( new SimOutputRST( baseName, rstNode ) );
#else
      EXCEPTION( "No support for ANSYS .rst output file format." );
#endif
    } else if(fileName.find( ".unv") != std::string::npos) {
#ifdef USE_UNV
      baseName = std::string(fileName, 0, fileName.find(".unv"));
      ParamNode * unvNode = new ParamNode(false);
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

       std::cout << "\nFound " << types.size() << " sequence step(s)\n";

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
    std::string isFreecoord = param->Get("freeCoord")->AsString();
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
    Tok tokenizer(param->Get("freeCoord")->AsString(), sep);
    std::copy(tokenizer.begin(), tokenizer.end(),
              std::back_inserter(freeCoord));
    if (freeCoord.size() != 4)
    {
      EXCEPTION("Not enought arguments " << freeCoord.size() \
          << ". Need 4 arguments (comp, start, stop, inc')");
    }

    ParamNode* parent;

    StdVector<ParamNode*> childVec;
    childVec.Push_back(new ParamNode(false));  //comp
    childVec.Push_back(new ParamNode(false));  //start
    childVec.Push_back(new ParamNode(false));  //stop
    childVec.Push_back(new ParamNode(false));  //inc

    // create a tree (Paramnode):
    // domain->nodeList->nodes->(name, list (coordSysId, freeCoord->(comp, start, stop, inc)))
    // create leafes (comp, start, stop, inc)
    setParamNode(childVec[0], "comp",  freeCoord[0]);
    setParamNode(childVec[1], "start", freeCoord[1]);
    setParamNode(childVec[2], "stop",  freeCoord[2]);
    setParamNode(childVec[3], "inc",   freeCoord[3]);

    // create node (freeCoord) with comp,start,stop,inc as children
    parent = new ParamNode(false); //freeCoord
    setParamNode(parent, "freeCoord", "", &childVec);

    // create node (list) with freeCoord and coordSysId as children
    childVec.Clear();
    childVec.Push_back(new ParamNode(false)); // coordSysId
    childVec.Push_back(parent);
    setParamNode(childVec[0], "coordSysId", coordSysId);
    childVec.Push_back(new ParamNode(false)); // gridId
    childVec.Push_back(parent);
    setParamNode(childVec[1], "gridId", "default");
    parent = new ParamNode(false); //list
    setParamNode(parent, "list", "", &childVec);

    // create leaf (name)
    childVec.Clear();
    childVec.Push_back(parent);
    parent = new ParamNode(false); //nodes name=""
    setParamNode(parent, "name", node_name);
    childVec.Push_back(parent);

    // create node (nodes) with name and list as childer
    parent = new ParamNode(false); //nodes
    setParamNode(parent, "nodes", "", &childVec);

    // create node (nodeList) with nodes as child
    childVec.Clear();
    childVec.Push_back(parent);
    parent = new ParamNode(false); //nodeList
    setParamNode(parent, "nodeList", "", &childVec);

    // create node (domain) with nodeList as child
    childVec.Clear();
    childVec.Push_back(parent); //domain
    //parent = new ParamNode(false); //nodeList
    //setParamNode(param, "domain", "", &childVec);
    //childVec.Clear();

    // put the childs into param
    param->Get("domain")->GetChildren().Push_back(parent);
  }

  inline void setParamNode(ParamNode *paramNode, std::string name, std::string value,
      StdVector<ParamNode*>* children)
  {
    paramNode->SetName(name);
    if (name != "")
    {
      paramNode->SetValue(value);
    }
    if (children != NULL)
    {
      StdVector<ParamNode*> &childTmp = paramNode->GetChildren();
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

  try
  {
    param = new ParamNode(false);
    ParamsInit(argc, argv);

    // Switch this flag tc true for debugging
    if (param->Get("forceSegFault")->AsBool())
    {
      Exception::segfault_ = true;
    } else {
      Exception::segfault_ = false;
    }

    std::string param_mode = param->Get("mode")->AsString();

    // get filenames from parameter
    std::vector<std::string> param_files;
    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");

    // Initialize vector with output fields
    Tok tokenizer(param->Get("files")->AsString(), sep);
    std::copy(tokenizer.begin(), tokenizer.end(),
              std::back_inserter(param_files));
    UInt num_files = param_files.size();
    if (param_files.size() >= 1)
    {
      // This is necessary to run, but I do not know what it is for
      info = new InfoNode(param_mode + "_" + param_files[0] + ".info.xml", "<?xml version=\"1.0\"?>"); 
      info->SetName("cfsInfo"); 
    }

    if (param_mode == "calcAverage")
    {
      if (num_files != 2)
      {
        EXCEPTION( "Please provide <infFile> and <outFile>" );
      }
      CFSTool::calcAverage(param_files[0], param_files[1]);
    } else if (param_mode == "convert") {
      if (num_files != 2)
      {
        EXCEPTION( "Please provide <infFile> and <outFile>" );
      }
      CFSTool::Convert( param_files[0], param_files[1] );
    } else if (param_mode == "scalardiff") {
      Double tolerance = param->Get("eps")->AsDouble();
      if (num_files != 2)
      {
        EXCEPTION( "Please provide <infFile1> and <inFile2>" );
      }
      Double maxDiffMesh = 0.0, maxDiffHist = 0.0;
      std::cout << "Checking for mesh results:\n"
        << "==========================";
      maxDiffMesh = CFSTool::Diff( param_files[0], param_files[1], "", \
                                  true, false, maxDiffResultName);
      std::cout << "Checking for history results:\n"
        << "=============================";
      maxDiffHist = CFSTool::Diff( param_files[0], param_files[1], "", \
                                  true, true, maxDiffResultName );
      Double maxDiff = std::max( maxDiffMesh, maxDiffHist );
      if( maxDiff > tolerance ) {
        std::cout << "'" << param_files[0] << "' and '" << param_files[1]
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
        EXCEPTION( "Please provide <infFile1>, <inFile2> and <outFile>" );
      }
      Double maxDiff = 0.0;
      maxDiff = CFSTool::Diff( param_files[0], param_files[1], param_files[2], \
                                false, false, maxDiffResultName);
    } else if (param_mode == "meshdiffnormed") {
      if (num_files != 3)
      {
        EXCEPTION( "Please provide <infFile1>, <inFile2> and <outFile>" );
      }
      Double maxDiff = 0.0;
      maxDiff = CFSTool::Diff( param_files[0], param_files[1], param_files[2], \
                                true, false, maxDiffResultName);
    } else {
      EXCEPTION( "No such mode: " << param_mode <<". See help for available modes" );
      return EXIT_FAILURE;
    }
  } catch(std::exception& ex) {
    std::cerr << "The following error occured:\n" << ex.what();
    if (info != NULL)
    {
      info->Get(InfoNode::ERROR)->SetValue(ex.what());
      info->ToFile();
    }
    std::cerr << "The following error occured during program execution:\n\n" << ex.what();
    return -1;
  }
  info->ToFile();
  delete info;

  return 0;
}
