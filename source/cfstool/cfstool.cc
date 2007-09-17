// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// Include headers which define what types
// of in/output files CFS++ should support
#include <boost/lexical_cast.hpp>

#include <def_use_mesh.hh>
#include <def_use_gidpost.hh>
#include <def_use_hdf5.hh>
#include <def_use_gmv.hh>
#include <def_use_unv.hh>
#include <def_use_ansysrst.hh>

#include <iostream>
#include "General/environment.hh"
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

using namespace CoupledField;
namespace CFSTool {
  
 
  
  shared_ptr<SimInput> GetInputReader( const std::string& fileName ) {
    
    // determine suffix of fileName
    shared_ptr<SimInput> reader;

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
//    } else if( fileName.find( ".unv") != std::string::npos ||
//        fileName.find( ".unverg") != std::string::npos ||
//        fileName.find( ".unvref") != std::string::npos ) {
//#ifdef USE_UNV
//      reader = shared_ptr<SimInput>(new SimInputUnv( fileName, NULL ) );
//#else
//      EXCEPTION( "No support for UNV input file format." );
//#endif  
    } else {
      EXCEPTION( "Found not suitalbe reader for file '" << fileName
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
      writer = 
        shared_ptr<SimOutput>( new SimOutputGiD( baseName, gidNode ) );
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
      writer = 
              shared_ptr<SimOutput>( new SimOutputGMV( baseName, gmvNode ) );
#else 
      EXCEPTION( "No support for GMV output file format." );
#endif
    }
    else {
      EXCEPTION( "No support for GiD format" );
    }
    // create output writer
    return writer;
  }

  void Convert( const std::string& inFile, const std::string& outFile ) {
     
     // obtain input reader for inFile
     shared_ptr<SimInput> input = GetInputReader( inFile );
     
     // check capabilities of input class
     bool printGridOnly = false;
     if( std::find( input->GetCapabilities().begin(),
                    input->GetCapabilities().end(),
                    SimInput::MESH_RESULTS ) 
         == input->GetCapabilities().end() ) {
       printGridOnly = true;
     }
                     
     // read in mesh
     input->InitModule();
     UInt dim = input->GetDim();
     Grid * ptGrid = new GridCFS(dim);
     input->ReadMesh(ptGrid);
     ptGrid->FinishInit();
     
     // obtain output writer
     shared_ptr<SimOutput> output = GetOutputWriter( outFile );
     output->Init( ptGrid, printGridOnly);  
     
     // only iterate over results, if not only the mesh is converted
     if( !printGridOnly ) {

       // obtain number of multiSequenceSteps and get analysis types
       StdVector<AnalysisType> types;
       input->GetNumMultiSequenceSteps( types );
       std::cout << "Found " << types.GetSize() << " multisequence steps\n\n";

       // iterate over all multiSequenceSteps
       UInt stepOffset = 0;
       for( UInt iMsStep = 0; iMsStep < types.GetSize(); iMsStep++ ) {
         
         std::cout << "Converting multi sequence step " << iMsStep+1 << std::endl;
         // get resulttypes
         StdVector<shared_ptr<ResultInfo> > infos;
         input->GetResultTypes( iMsStep+1, infos );

         // get stepvalues
         StdVector<Double> stepVals;
         input->GetStepValues( iMsStep+1, stepVals );
         
         // notify writer
         output->BeginMultiSequenceStep( iMsStep+1, types[iMsStep], stepVals.GetSize() );

         StdVector<shared_ptr<BaseResult> > results;
         
         if( infos.GetSize() > 0 ){
             std::cout << "Found the following results:\n";
         }
         // iterate over all result types
         for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {

           std::string resultName;
           Enum2String( infos[iRes]->resultType, resultName);
           std::cout << "\t" << resultName << std::endl;
           // iterate over all regions
           StdVector<shared_ptr<EntityList> > regions;
           input->GetResultEntities( iMsStep+1, infos[iRes], regions );
           for( UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++ ) {
             // generate new result object and add it to output writer
             shared_ptr<BaseResult > result;
             if( types[iMsStep] != HARMONIC ) {
               result = shared_ptr<BaseResult>( new Result<Double>() );
             } else {
               result = shared_ptr<BaseResult>( new Result<Complex>() );
             }
             result->SetEntityList( regions[iRegion] );
             result->SetResultInfo( infos[iRes] );
             results.Push_back( result );
             // CAUTION: begin, inc, end are hardcoded and noch checked for each result!
             output->RegisterResult( result, 1, 1, stepVals.GetSize() ); 
           }
           
         }
         // iterate over all stepvalues
         for( UInt iStep = 0; iStep < stepVals.GetSize(); iStep++ ) {
           Double actStepVal = stepVals[iStep];

           output->BeginStep( iStep+1+stepOffset, actStepVal );
           // iterate over all results
           for( UInt iRes = 0; iRes < results.GetSize(); iRes++) {
             input->GetResult( iMsStep+1, iStep+1+stepOffset, results[iRes] );
             output->AddResult( results[iRes] );
           }
           output->FinishStep();

         }
         output->FinishMultiSequenceStep();
         
         // add number of steps to stepoffset for subsequente 
         // multisequence steps
         stepOffset += stepVals.GetSize();
       }
     } // printGridOnly
     output->Finalize();
     delete ptGrid;
     std::cout << "\nOutput successfully written to " << outFile << std::endl;
  } //Convert()
  
  
  Double Diff( const std::string& inFile1, 
               const std::string& inFile2,
               const std::string& outFile ) {
       
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
       
       // obtain number of multiSequenceSteps and get analysis types
       StdVector<AnalysisType> types;
       input1->GetNumMultiSequenceSteps( types );
       if( output )
         std::cout << "Found " << types.GetSize() << " multisequence steps\n\n";

       // iterate over all multiSequenceSteps
       UInt stepOffset = 0;
       Double maxDiff = 0.0;
       for( UInt iMsStep = 0; iMsStep < types.GetSize(); iMsStep++ ) {

         if( output )
           std::cout << "Performing diff on multi sequence step " << iMsStep+1 << std::endl;
         
         // get resulttypes
         StdVector<shared_ptr<ResultInfo> > infos;
         input1->GetResultTypes( iMsStep+1, infos );

         // get stepvalues
         StdVector<Double> stepVals;
         input1->GetStepValues( iMsStep+1, stepVals );

         // notify writer
         if( output) {
           output->BeginMultiSequenceStep( iMsStep+1, types[iMsStep], stepVals.GetSize() );
         }

         StdVector<shared_ptr<BaseResult> > inResults1, inResults2, outResults;

         // iterate over all result types of input1
         for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {

           std::string resultName;
           Enum2String( infos[iRes]->resultType, resultName);
           if( output )
             std::cout << "\t" << resultName << std::endl;
           
           // iterate over all regions
           StdVector<shared_ptr<EntityList> > regions;
           input1->GetResultEntities( iMsStep+1, infos[iRes], regions );
           for( UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++ ) {
             // generate new result object and add it to output writer
             shared_ptr<BaseResult > inResult1, inResult2, outResult;
             if( types[iMsStep] != HARMONIC ) {
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
               outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
               
               // CAUTION: begin, inc, end are hardcoded and noch checked for each result
               output->RegisterResult( outResult, 1, 1, stepVals.GetSize() );
             }
           }
         }

         // iterate over all stepvalues
         for( UInt iStep = 0; iStep < stepVals.GetSize(); iStep++ ) {
           Double actStepVal = stepVals[iStep];
           
           if( output) {
             output->BeginStep( iStep+1+stepOffset, actStepVal );
           }
           
           // iterate over all results
           for( UInt iRes = 0; iRes < inResults1.GetSize(); iRes++) {
             // obtain both result objects for current step
             input1->GetResult( iMsStep+1, iStep+1+stepOffset, inResults1[iRes] );
             input2->GetResult( iMsStep+1, iStep+1+stepOffset, inResults2[iRes] );
             if( types[iMsStep] != HARMONIC ) {
               Vector<Double> & inVec1 = dynamic_cast<Result<Double>& >(*inResults1[iRes]).GetVector();
               Vector<Double> & inVec2 = dynamic_cast<Result<Double>& >(*inResults2[iRes]).GetVector();
               Vector<Double> & outVec = dynamic_cast<Result<Double>& >(*outResults[iRes]).GetVector();
               
               outVec.Resize( inVec1.GetSize() );
               Double max = 0.0;
               for( UInt i = 0; i<inVec1.GetSize(); i++ ) {
                 if( std::abs(inVec1[i]) > max) 
                   max = std::abs(inVec1[i]);
               }
               outVec = inVec1 - inVec2;
               outVec /= max;
               for( UInt i = 0; i < outVec.GetSize(); i++ ) {
                 if( std::abs(outVec[i]) > maxDiff) {
                   maxDiff = std::abs(outVec[i]) ;
                 }
               }
             
             } else {
               Vector<Complex> & inVec1 = dynamic_cast<Result<Complex>& >(*inResults1[iRes]).GetVector();
               Vector<Complex> & inVec2 = dynamic_cast<Result<Complex>& >(*inResults2[iRes]).GetVector();
               Vector<Complex> & outVec = dynamic_cast<Result<Complex>& >(*outResults[iRes]).GetVector();

               outVec.Resize( inVec1.GetSize() );
               Double max = 0.0;
               for( UInt i = 0; i<inVec1.GetSize(); i++ ) {
                 if( std::abs(inVec1[i]) > max) 
                   max = std::abs(inVec1[i]);
               }
               outVec = inVec1 - inVec2;
               outVec /= max;
               
               for( UInt i = 0; i < outVec.GetSize(); i++ ) {
                 if( std::abs(outVec[i]) > maxDiff) {
                   maxDiff = std::abs(outVec[i]) ;
                 }
               }
               
             }
             
             // cast result objects, get vector and calculate difference vector
             // iterate over all entities and dofs
             // for vector results: 
             // 1) calculate vector norm of reference result
             // 2) 
             if (output) {
               output->AddResult( outResults[iRes] ); 
             }
           }
           if( output )
             output->FinishStep();

         }
         if( output )
           output->FinishMultiSequenceStep();

         // add number of steps to stepoffset for subsequente 
         // multisequence steps
         stepOffset += stepVals.GetSize();
       }
       if( output ) {
         output->Finalize();
         std::cout << "\nOutput successfully written to " << outFile << std::endl;
       }
       delete ptGrid1;
       delete ptGrid2;
       
       return maxDiff;
       
  } //MeshDiff
}


int main(int argc, char** argv) {
  
  // 
  Exception::segfault_ = false;
  
  try {
    if( argc < 2) {
      std::cout << "CFS TOOL 1.0 \n\n"
                << "Usage:\tcfstool modus [args]\n\n"
                << "The following modi are vailable:\n\n"
                << "\tconvert <infile> <outfile>\n"
                << "\tdiff <infile1> <infile2> <tolerance>\n"
                << "\tmeshdiff <infile1> <infile2> <outfile>\n"
                << "\n\n"
                << "List of supported formats:\n"
                << "\tinput: .mesh .h5\n"
                << "\toutput: .h5, .post.res, .post.bin, .gmv\n"
                << "\n\n"
                << "Please not that the input / output format is chosen "
                << "depending on the file suffix\n\n";
      return EXIT_FAILURE;
    }
    // check command line arguments
    // 3 modi:
    // 1) cfstool convert <infile> <outfile>
    // 2) cfstool meshdiff <infile1> <infile2> <outfile>
    // 3) cfstool scalardiff <infile1> <infile2>
    std::string modus = argv[1];

    if( modus == "convert" ) {
      if( argc != 4 ) {
        EXCEPTION( "Please provide <infFile> and <outFile>" );
      }
      CFSTool::Convert( argv[2], argv[3] );

    } else if( modus == "scalardiff" ) {
      if( argc != 5 ) {
        EXCEPTION( "Please provide <infFile1>, <inFile2> and <tolerance>" );
      }
      Double maxDiff = 0.0;
      Double tolerance = 0.0;
      try {
        tolerance = boost::lexical_cast<Double>(argv[4]);
      } catch (std::exception& ex ) {
        EXCEPTION( "Could not convert '" << argv[4] << "' to double value");
      }
      maxDiff = CFSTool::Diff( argv[2], argv[3], "" );
      if( maxDiff > tolerance ) {
        std::cerr << "Files '" << argv[2] << "' and '"
        << argv[3] << "' differ with maximum difference of "
        << maxDiff << "\n";
        exit(EXIT_FAILURE);
      } else {
        exit(EXIT_SUCCESS);
      }

    } else if( modus == "meshdiff" ) {
      if( argc != 5 ) {
        EXCEPTION( "Please provide <infFile1>, <inFile2> and <outFile>" );
      }
      Double maxDiff = 0.0;
      maxDiff = CFSTool::Diff( argv[2], argv[3], argv[4] );

    } else {
      std::cerr << "modus '" << modus << "' not known\n";
      return EXIT_FAILURE;
    }


  }  catch(std::exception& ex) {
    std::cerr << "The following error occured during program execution:\n\n" 
    << ex.what ();
  }

  return 0;
}
