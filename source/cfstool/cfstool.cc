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
    } else {
      EXCEPTION( "Output format not supported!" );
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
       std::cout << "\nFound " << types.GetSize() << " sequence steps:\n";

       // iterate over all multiSequenceSteps
       UInt stepOffset = 0;
       for( UInt iMsStep = 0; iMsStep < types.GetSize(); iMsStep++ ) {
         
         std::cout << "\tConverting step no. " << iMsStep+1 << std::endl;
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
             std::cout << "\tfor results ";
         }
         // iterate over all result types
         for( UInt iRes = 0; iRes < infos.GetSize(); iRes++) {

           std::cout << "\t" << infos[iRes]->resultName << std::endl;

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
             try {
             input->GetResult( iMsStep+1, iStep+1+stepOffset, results[iRes] );
             output->AddResult( results[iRes] );
             } catch (Exception& ex ) {
               // do nothing
               // Note: This is currently the "hard-coded" way of handling results,
               // which have a saveBegin/End/Increment which differs from the 
               // default values
             }
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
 
  double RadPhase( const Complex& c ) {
    return std::atan2(c.imag() ,c.real() ); 
  }
  
  Double Diff( const std::string& inFile1, 
               const std::string& inFile2,
               const std::string& outFile,
               const bool normedtomax) {
       
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
       StdVector<AnalysisType> types;
       input1->GetNumMultiSequenceSteps( types );
       if( output )
         std::cout << "\nFound " << types.GetSize() << " sequence steps:\n";

       // iterate over all Sequence Steps
       UInt stepOffset = 0;
       Double maxDiff = 0.0;
       for( UInt iMsStep = 0; iMsStep < types.GetSize(); iMsStep++ ) {

         if( output )
           std::cout << "\tPerforming diff on step no. " << iMsStep+1 << std::endl;
         
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

           if( output )
             std::cout << "\tfor result '" << infos[iRes]->resultName << "'" << std::endl;

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
               //outResult->GetResultInfo()->complexFormat = AMPLITUDE_PHASE;
               outResult->GetResultInfo()->complexFormat = REAL_IMAG;
               
               // CAUTION: begin, inc, end are hardcoded and noch checked for each result
               output->RegisterResult( outResult, 1, 1, stepVals.GetSize() );
             }
           }
         }

         // iterate over all time/frequency steps
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
             
             // cast result objects, get vector and calculate difference vector
             if( types[iMsStep] != HARMONIC ) {
               Vector<Double> & inVec1 = dynamic_cast<Result<Double>& >(*inResults1[iRes]).GetVector();
               Vector<Double> & inVec2 = dynamic_cast<Result<Double>& >(*inResults2[iRes]).GetVector();
               Vector<Double> & outVec = dynamic_cast<Result<Double>& >(*outResults[iRes]).GetVector();
               outVec.Resize( inVec1.GetSize() );

               // normalize to maximum value of inResult2
               Double maxRes2 = 0.0;
               for( UInt i = 0; i<inVec2.GetSize(); i++ ) {
                 if( std::abs(inVec2[i]) > maxRes2) 
                   maxRes2 = std::abs(inVec2[i]);
               }
               outVec = inVec1 - inVec2;
               if (normedtomax == true) {
                 outVec /= maxRes2;
               }
               
               // find maximum entry in difference vector
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

               // find maximum amplitude of inResult2
               Double maxRes2 = 0.0;
               for( UInt i = 0; i<inVec2.GetSize(); i++ ) {
                 if( std::abs(inVec2[i]) > maxRes2) 
                   maxRes2 = std::abs(inVec2[i]);
               }

               Double aDiff, pDiff, aMax=0.0, aMin=0.0, pMax=0.0, pMin=0.0;
               for( UInt i = 0; i<inVec2.GetSize(); i++ ) {
                 aDiff = std::abs(inVec1[i]) - std::abs(inVec2[i]);
                 // phase difference in multiples of pi
                 pDiff = RadPhase(inVec1[i]) - RadPhase(inVec2[i]);
                 
                 // correct 2*pi-offset if phase angles have different signs
                 if ( (std::abs(pDiff)>PI) && (pDiff<0) )
                   pDiff+= 2*PI;
                 if ( (std::abs(pDiff)>PI) && (pDiff>0) )
                   pDiff-= 2*PI;

                 
                 if ( normedtomax == true ) {
                   aDiff /= maxRes2; // norm by maximum amplitude of inResult2
                   pDiff /= 2*PI; // norm relative phase difference by 2 pi
                   
                   if( std::abs(aDiff) > std::abs(pDiff) ) {
                     if( std::abs(aDiff) > maxDiff)
                       maxDiff = std::abs(aDiff);
                   } else {
                     if( std::abs(pDiff) > maxDiff)
                       maxDiff = std::abs(pDiff);
                   }
                   
                 } else {
                   
                   // if no relative differences are calculated return maximum amplitude difference
                   if( std::abs(aDiff) > maxDiff)
                     maxDiff = std::abs(aDiff);
                 }
                 
                 // maximum values in positive and negative direction                                  
                 if( pDiff > pMax )
                   pMax = pDiff;
                 if( pDiff < pMin )
                   pMin = pDiff;
                 if( aDiff > aMax )
                   aMax = aDiff;
                 if( aDiff < aMin )
                   aMin = aDiff;
                 
                 // ------------------------------------------------------------------------------------
                 // sign of amplitude difference gets lost, otherwise
                 //  negative sign of aDiff would mean phase shift by +/- pi
                 //outVec[i] = Complex( std::abs(aDiff)*std::cos(pDiff), 
                 //                     std::abs(aDiff)*std::sin(pDiff) );
                 
                 // Dirty hack! Write differences in real_imag format.
                 outVec[i] = Complex( aDiff, pDiff*180/PI );
                 // ------------------------------------------------------------------------------------
               }
               
               if( normedtomax == true)
                 std::cerr << "\n\tMaximum rel. + amplitude difference:  " << aMax*100 << " %\n"
                           << "\tMaximum rel. - amplitude difference: " << aMin*100 << " %\n"
                           << "\tMaximum rel. + phase difference:      " << pMax*100 << " %\n"
                           << "\tMaximum rel. - phase difference:     " << pMin*100 << " %\n";
               else
                 std::cerr << "\n\tMaximum + amplitude difference:  " << aMax <<  "\n"
                           << "\tMaximum - amplitude difference: " << aMin <<  "\n"
                           << "\tMaximum + phase difference:      " << pMax*180/PI <<  "°\n"
                           << "\tMaximum - phase difference:     " << pMin*180/PI <<  "°\n";
               
               
               std::cerr << "\n\tmaxDiff = " << maxDiff << std::endl;
             }
             
             
             if ( output ) {
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
       
  } //Diff
  
} //Namespace CFSTool


int main(int argc, char** argv) {
  
  // 
  Exception::segfault_ = true;
  
  try {
    if( argc < 2) {
      std::cout << "CFS TOOL 1.0 \n\n"
                << "Usage:\tcfstool modus [args]\n\n"
                << "The following modi are vailable:\n\n"
                << "\tconvert <infile> <outfile>\n"
                << "\tscalardiff <infile1> <infile2> <tolerance>\n"
                << "\tmeshdiff <infile1> <infile2> <outfile>\n"
                << "\tmeshdiffnormed <infile1> <infile2> <outfile>\n"
                << "\n\n"
                << "List of supported formats:\n"
                << "\tinput: .mesh .h5\n"
                << "\toutput: .h5, .post.res, .post.bin, .gmv\n"
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
      Double maxDiff = 0.0;
      maxDiff = CFSTool::Diff( argv[2], argv[3], "", true );
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
      maxDiff = CFSTool::Diff( argv[2], argv[3], argv[4], false );

    } else if( modus == "meshdiffnormed" ) {
      if( argc != 5 ) {
        EXCEPTION( "Please provide <infFile1>, <inFile2> and <outFile>" );
      }
      Double maxDiff = 0.0;
      maxDiff = CFSTool::Diff( argv[2], argv[3], argv[4], true );

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
