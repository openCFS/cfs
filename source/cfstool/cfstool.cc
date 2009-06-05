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
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

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
  
 
  
  shared_ptr<SimInput> GetInputReader( const std::string& fileName ) {
    
    // determine suffix of fileName
    shared_ptr<SimInput> reader;

    if ( !fs::exists( fileName ) )
      EXCEPTION( "\nFile '" << fileName << "' does not exist!");

    if( fileName.find( ".mesh") != std::string::npos ) {
#ifdef USE_MESH
      reader = shared_ptr<SimInput>(new SimInputMESH( fileName, NULL ) );
#else
      EXCEPTION( "No support for MESH input file format." );
#endif
    } else if( fileName.find( ".h5") != std::string::npos ) {
#ifdef USE_HDF5
      ParamNode * hdf5Node = new ParamNode(false);
      reader = shared_ptr<SimInput>(new SimInputHDF5(fileName, hdf5Node) );
#else  
      EXCEPTION( "No support for HDF5 input file format." );
#endif
    } else if( fileName.find( ".gmv") != std::string::npos ) {
#ifdef USE_GMV_INPUT
      ParamNode * gmvNode = new ParamNode(false);
      reader = shared_ptr<SimInput>(new SimInputGMV(fileName, gmvNode) );
#else  
      EXCEPTION( "No support for GMV input file format." );
#endif
    } else if( fileName.find( ".unv") != std::string::npos ||
        fileName.find( ".unverg") != std::string::npos ||
        fileName.find( ".unvref") != std::string::npos ) {
#ifdef USE_UNV
      ParamNode * unvNode = new ParamNode(false);
      reader = shared_ptr<SimInput>(new SimInputUnv( fileName, unvNode ) );
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

         // notify writer
         output->BeginMultiSequenceStep( actMsStep, types[actMsStep], numSteps[actMsStep] );
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
               bool isHistory ) {
       
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

           std::cout << "\t" << infos[iRes]->resultName << "\n\n"; 

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
             if ( output )
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
  
  /** Initialize static Enums. 
   * todo: do better once - Fabian */
  void InitEnums()
  {
    SetEnvironmentEnums();
    BasePDE::SetEnums();
  }
  
  
} //Namespace CFSTool


int main(int argc, char** argv) {
  
  // Switch this flag to true for debugging
#ifndef NDEBUG
  Exception::segfault_ = true;
#else
  Exception::segfault_ = false;
#endif

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
 
  try {
    if( argc < 2) {
      std::cout << "CFS TOOL 1.0 \n\n"
                << "Usage:\tcfstool modus [args]\n\n"
                << "The following modi are vailable:\n\n"
                << "\tconvert <infile> <outfile>\n"
                << "\tscalardiff <reference_file> <compare_file> <tolerance>\n"
                << "\t\tIs max(in_1 - in_2) / max(in_2) < tolerance?\n"
                << "\tmeshdiff <reference_file> <compare_file> <outfile>\n"
                << "\t\ti.e. out = in_1 -in_2\n"
                << "\tmeshdiffnormed <reference_file> <compare_file> <outfile>\n"
                << "\t\ti.e. out = (in_1 - in_2) / max(in_2)\n"
                << "\tversion"
                << "\n\n"
                << "List of supported formats:\n"
                << "\tinput: .mesh .h5 .unv[erg|ref] .gmv\n"
                << "\toutput: .h5, .post.res, .post.bin, .gmv .rst\n"
                << "\n\n"
                << "Please note that the input/output format is chosen "
                << "depending on the file suffix\n\n";
      return EXIT_FAILURE;
    }
    // check command line arguments
    // 4 modi:
    // 1) cfstool convert <infile> <outfile>
    // 2) cfstool scalardiff <infile1> <infile2>
    // 3) cfstool meshdiff <infile1> <infile2> <outfile>
    // 4) cfstool meshdiffnormed <infile1> <infile2> <outfile>
    std::string modus = argv[1];
    
    // we have to instantiate the global InfoNode object as some classes 
    // assumes it for its logging. It might also be helpfull in providing more details 
    // do this only if we have a change of valid command line 
    if(argc >= 4) { 
      info = new InfoNode(modus + "_" + argv[2] + ".info.xml", "<?xml version=\"1.0\"?>"); 
      info->SetName("cfsInfo"); 
    } else { 
      info = NULL; 
    } 
    
    
    if( modus == "convert" ) {
      if( argc != 4 ) {
        EXCEPTION( "Please provide <infFile> and <outFile>" );
      }
      CFSTool::Convert( argv[2], argv[3] );

    } else if( modus == "scalardiff" ) {
      if( argc != 5 ) {
        EXCEPTION( "Please provide <infFile1>, <inFile2> and <tolerance>" );
      }
      Double tolerance = 0.0;
      try {
        tolerance = boost::lexical_cast<Double>(argv[4]);
      } catch (std::exception& ex ) {
        EXCEPTION( "Could not convert '" << argv[4] << "' to double value");
      }
      Double maxDiffMesh = 0.0, maxDiffHist = 0.0;
      std::cout << "Checking for mesh results:\n"
                << "==========================";
      maxDiffMesh = CFSTool::Diff( argv[2], argv[3], "", true, false );
      std::cout << "Checking for history results:\n"
                << "=============================";
      maxDiffHist = CFSTool::Diff( argv[2], argv[3], "", true, true );
      Double maxDiff = std::max( maxDiffMesh, maxDiffHist );
      if( maxDiff > tolerance ) {
        std::cout << "  Files '" << argv[2] << "' and '" << argv[3]
                  << "' differ with maximum difference of "
                  << maxDiff << "\n";
        exit(EXIT_FAILURE);
      } else {
        std::cout << "  No differences larger than tolerance found.\n";
        exit(EXIT_SUCCESS);
      }

    } else if( modus == "meshdiff" ) {
      if( argc != 5 ) {
        EXCEPTION( "Please provide <infFile1>, <inFile2> and <outFile>" );
      }
      Double maxDiff = 0.0;
      maxDiff = CFSTool::Diff( argv[2], argv[3], argv[4], false, false );

    } else if( modus == "meshdiffnormed" ) {
      if( argc != 5 ) {
        EXCEPTION( "Please provide <infFile1>, <inFile2> and <outFile>" );
      }
      Double maxDiff = 0.0;
      maxDiff = CFSTool::Diff( argv[2], argv[3], argv[4], true, false );
    } else if( modus == "version" ) {
      ProgramOptions::GetVersionString(std::cout, false);
      return EXIT_SUCCESS;
    } else {
      std::cerr << "mode '" << modus << "' not known\n";
      return EXIT_FAILURE;
    }


  }  
  catch(std::exception& ex) 
  { 
    std::cerr << "The following error occured:\n" << ex.what(); 
    if(info != NULL)  { 
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
