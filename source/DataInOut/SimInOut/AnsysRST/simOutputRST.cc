// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include <boost/algorithm/string/trim.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Domain/elem.hh"

#include "simOutputRST.hh"

#include "Utils/DynamicLibrary.hh"
#include "Utils/DynamicObject.hh"
#include "AnsysBinlibIface.hh"

#include <def_use_ansysrst.hh>
#include <def_cfs_stats.hh>

namespace CoupledField {
  
  SimOutputRST::SimOutputRST( const std::string& fileName,
                              ParamNode * outputNode )
    : SimOutput( fileName, outputNode ) {

    // Initialize variables
    formatName_ = "rst";
    fileName_ = fileName;
    dirName_ = "results_" + formatName_;
    outputNode->Get("directory", dirName_, false );
    outputNode_ = outputNode;

    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

  }


  SimOutputRST::~SimOutputRST() {
    // Get rid of the binlibIface_ object
    binlibIface_->deleteSelf();
    binlibIface_ = NULL;
    
    // Close the dynamic library
    delete dynLibrary_;
  }

  void SimOutputRST::Init( Grid* ptGrid, bool printGridOnly ) {
    std::stringstream sstr;

    std::map<std::string, std::string> revisionDLLMap;
    std::map<std::string, std::string> archMachIdMap;

    try 
    {
      sysPathSep_ = fs::path("/").native_directory_string();
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    revisionDLLMap["10.0"] = "libCFSAnsysBinlibIfaceRev100.so";
    revisionDLLMap["11.0"] = "libCFSAnsysBinlibIfaceRev110.so";
    
    archMachIdMap["I386"] = "linia32";
    archMachIdMap["X86_64"] = "linop64";

    sstr.clear(); sstr.str("");
    sstr << outputNode_->Get("ansysRevision")->AsString();
    sstr >> ansysBinlibRev_;
    sstr.clear(); sstr.str("");
    sstr.precision(1);
    sstr << std::fixed << ansysBinlibRev_;
    std::string requestedRev;
    sstr >> requestedRev;
    sstr.clear(); sstr.str("");

    if(revisionDLLMap.find(requestedRev) == revisionDLLMap.end()) 
    {
      std::map<std::string, std::string>::const_iterator it, end;
      
      sstr << "Supported revisions include: ";
      it = revisionDLLMap.begin();
      end = revisionDLLMap.end();
      for( ; it != end; it++ ) 
      {
        sstr << std::fixed << it->first << " ";
      }
      sstr << "\n";
           
      EXCEPTION("Writing .rst files for ANSYS revision " << requestedRev 
                << " is not supported.\n"
                << sstr.str());
    }
    
    sstr.clear(); sstr.str("");
    sstr << CFS_ANSYS_LD_PATH << sysPathSep_ << revisionDLLMap[requestedRev];
    std::string machineId = archMachIdMap[CFS_ARCH];

    // Load the dynamic library with the high-pass filter class
    dynLibrary_ = 
      DynamicLoader::loadObjectFile(sstr.str().c_str(), RTLD_NOW);
    if(dynLibrary_ == NULL) {
      EXCEPTION("Couldn't load the dynamic library '" << sstr.str() 
                << "'.\n"
                << "Please make sure you have the necessary paths:\n"
                << "'" << CFS_ANSYS_LD_PATH << "'\n and"
                << "'ansys/vXXX/ansys/syslib/" << machineId << "' and\n"
                << "'ansys/vXXX/ansys/customize/misc/" << machineId << "'\n"
                << "in your LD_LIBRARY_PATH variable!");
    }
    // Load a HighPassBinlibIface_ object
    binlibIface_ = 
      dynamic_cast<AnsysBinlibIface*>(dynLibrary_->newObject("AnsysBinlibIfaceGeneral", 0, NULL));
    if(binlibIface_ == NULL) {
      std::cerr << "Couldn't create binlibIface_ object." << std::endl;
    }

    msStep_ = 0;
    
    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }
    
    // Variables for construction of file name
    std::string filename;
    std::stringstream fnstream;
    
    // Generate basename for output file
    // If just the grid is written the file will end with
    // "_grid.rst". If we are in multistep N the file will end
    // with "_msN.rst".
    fnstream << dirName_ << sysPathSep_ << fileName_;
    fnstream << "_grid";
    fnstream << ".rst";
    
    filename = fnstream.str();

    RemoveFile(filename, "Error while removing previous .rst file." );

    GetAnsysRuntimeInfos();

    if(ansysBinlibRev_ != compiledAnsysRev_)
    {
      EXCEPTION("You have requested to write .rst files for ANSYS revision "
                << compiledAnsysRev_ << "\nbut your libbin.so has revision "
                << ansysBinlibRev_ << ".\n"
                << "Please set\n"
                << "'ansys/vXXX/ansys/syslib/" << machineId << "' and\n"
                << "'ansys/vXXX/ansys/customize/misc/" << machineId << "'\n"
                << "in your LD_LIBRARY_PATH variable!");
    }
    
    binlibIface_->SetFNAndOutputNode( fileName_, formatName_,
                                      dirName_, sysPathSep_,
                                      outputNode_ );

    binlibIface_->SetResultFileName( filename );


    binlibIface_->Init(ptGrid, printGridOnly);
  }


  void SimOutputRST::BeginMultiSequenceStep( UInt step,
                                             BasePDE::AnalysisType type,
                                             UInt numSteps)
  {
    msStep_++;
    // Variables for construction of file name
    std::string filename;
    std::stringstream fnstream;
    
    // Generate basename for output file
    // If just the grid is written the file will end with
    // "_grid.rst". If we are in multistep N the file will end
    // with "_msN.rst".
    fnstream << dirName_ << sysPathSep_ << fileName_;
    fnstream << "_ms" << msStep_;
    fnstream << ".rst";
    
    filename = fnstream.str();

    RemoveFile(filename, "Error while removing previous .rst file.");
    
    binlibIface_->SetResultFileName( filename );

    binlibIface_->BeginMultiSequenceStep( step, type, numSteps);
  }

  //! Register result (within one multisequence step)
  void SimOutputRST::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd,
                                     bool isHistory )
  {
    binlibIface_->RegisterResult( sol,
                                  saveBegin,
                                  saveInc,
                                  saveEnd,
                                  isHistory );
  }

  void SimOutputRST::BeginStep( UInt stepNum, Double stepVal ) {
    binlibIface_->BeginStep( stepNum, stepVal );
  }
 
  void SimOutputRST::AddResult( shared_ptr<BaseResult> sol ) {
    binlibIface_->AddResult( sol );
  }
  
  void SimOutputRST::FinishStep( ) {
    binlibIface_->FinishStep( );
  }
  
  void SimOutputRST::FinishMultiSequenceStep( )
  {
    binlibIface_->FinishMultiSequenceStep();
  }
  

  void SimOutputRST::Finalize() 
  {
  }

  void SimOutputRST::GetAnsysRuntimeInfos() 
  {
    std::string TestFN;
    std::string revision;
    std::string compiled_rev;
    std::stringstream sstr;

    sstr.clear();
    sstr.str("");
    sstr << dirName_ << sysPathSep_ << fileName_;
    sstr << "_test";
    sstr << ".rst";
    TestFN = sstr.str();
    
    RemoveFile(TestFN, "Error while writing test .rst file.");

    binlibIface_->GetANSYSRevisionInfos(revision,
                                        ansysMachineId_,
                                        compiled_rev,
                                        TestFN);

    sstr.clear(); sstr.str("");
    sstr << revision;
    sstr >> ansysBinlibRev_;

    sstr.clear(); sstr.str("");
    sstr << compiled_rev;
    sstr >> compiledAnsysRev_;
    
    boost::algorithm::trim(ansysMachineId_);
    
    RemoveFile(TestFN, "Error while removing test .rst file.");

    //    std::cout.precision(1);
    //    std::cout << std::fixed << ansysBinlibRev_ << " "
    //              << ansysMachineId_ << " compiled rev: "
    //              << std::fixed << compiledAnsysRev_
    //              << std::endl;
  }

  void SimOutputRST::RemoveFile(const std::string& fileName,
                                const std::string& exception) 
  {
    try 
    {
      if(fs::exists(fileName))
      {
        fs::remove( fileName );
      }
    } catch (std::exception &ex)
    {
      EXCEPTION(exception);
    }
  }
  
  
} // end of namespace
