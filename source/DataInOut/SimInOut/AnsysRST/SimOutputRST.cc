// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iomanip>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include <boost/algorithm/string/trim.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/ElemMapping/Elem.hh"

#include "SimOutputRST.hh"

#include "Utils/DynamicLibrary.hh"
#include "Utils/DynamicObject.hh"
#include "AnsysBinlibIface.hh"

#include <def_use_ansysrst.hh>
#include <def_cfs_stats.hh>

namespace CoupledField {
  
  SimOutputRST::SimOutputRST( const std::string& fileName,
                              PtrParamNode outputNode,
                              PtrParamNode infoNode)
    : SimOutput( fileName, outputNode, infoNode ),
      dynLibrary_(NULL),
      binlibIface_(NULL)
  {

    // Initialize variables
    formatName_ = "rst";
    fileName_ = fileName;
    dirName_ = "results_" + formatName_;
    outputNode->GetValue("directory", dirName_, ParamNode::INSERT );
    outputNode_ = outputNode;

    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

  }


  SimOutputRST::~SimOutputRST() {
    // Get rid of the binlibIface_ object
    if(binlibIface_) 
    {
      binlibIface_->deleteSelf();
      binlibIface_ = NULL;
    }
    
    // Close the dynamic library
    delete dynLibrary_;
  }

  void SimOutputRST::Init( Grid* ptGrid, bool printGridOnly ) {
    std::stringstream sstr;

    std::map<std::string, std::string> revisionDLLMap;
    std::map<std::string, std::string> archMachIdMap;

    ptGrid_ = ptGrid;
    
    try 
    {
      sysPathSep_ = fs::path("/").string();
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

#ifndef _WIN32
    revisionDLLMap["10.0"] = "libCFSAnsysBinlibIfaceRev100.so";
    revisionDLLMap["11.0"] = "libCFSAnsysBinlibIfaceRev110.so";
    revisionDLLMap["12.0"] = "libCFSAnsysBinlibIfaceRev120.so";
    revisionDLLMap["12.1"] = "libCFSAnsysBinlibIfaceRev121.so";
    revisionDLLMap["13.0"] = "libCFSAnsysBinlibIfaceRev130.so";
    revisionDLLMap["14.0"] = "libCFSAnsysBinlibIfaceRev140.so";

    archMachIdMap["I386"] = "linia32";
    archMachIdMap["X86_64"] = "linop64";
#else
    revisionDLLMap["10.0"] = "libCFSAnsysBinlibIfaceRev100.dll";
    revisionDLLMap["11.0"] = "libCFSAnsysBinlibIfaceRev110.dll";
    revisionDLLMap["12.0"] = "libCFSAnsysBinlibIfaceRev120.dll";
    revisionDLLMap["12.1"] = "libCFSAnsysBinlibIfaceRev121.dll";
    revisionDLLMap["13.0"] = "libCFSAnsysBinlibIfaceRev130.dll";
    revisionDLLMap["14.0"] = "libCFSAnsysBinlibIfaceRev140.dll";
    
    archMachIdMap["I386"] = "winia32";
    archMachIdMap["X86_64"] = "winx64";
#endif
    
    // The default revision for ANSYS is now set to 11.0. Note however, that
    // not all features of the 11.0, 12.1 and even later revisions are supported.
    // This is especially true for animation of harmonic results in ANSYS Classic.
    // Moreover, revisions later than 11.0 may cause problems, when trying to load
    // .rst files into Tecplot or EnSight.
    std::string ansysRev = "11.0";
    outputNode_->GetValue("ansysRevision", ansysRev, ParamNode::INSERT);
    sstr.clear(); sstr.str("");
    sstr << ansysRev;
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
    
    std::string machineId = archMachIdMap[CFS_ARCH];

    // First try to load the binlib interface DLL from the LD_LIBRARY_PATH
    sstr.clear(); sstr.str("");
    sstr << revisionDLLMap[requestedRev];
    int dynLoadFlags = 0;
#ifndef _WIN32
    dynLoadFlags = RTLD_NOW;
#endif

    // Load the dynamic library with the high-pass filter class
    dynLibrary_ = 
      DynamicLoader::loadObjectFile(sstr.str().c_str(), dynLoadFlags);
      
    if(dynLibrary_ == NULL) {

      // Now try to load from builtin library path.
      sstr.clear(); sstr.str("");
      sstr << CFS_ANSYS_LD_PATH << sysPathSep_ << revisionDLLMap[requestedRev];

      dynLibrary_ = 
        DynamicLoader::loadObjectFile(sstr.str().c_str(), dynLoadFlags);
      if(dynLibrary_ == NULL) {
        EXCEPTION("Couldn't load the dynamic library '" << sstr.str() 
                  << "'.\n"
                  << "Please make sure you have the necessary paths:\n"
                  << "'" << CFS_ANSYS_LD_PATH << "'\n and"
                  << "'ansys/vXXX/ansys/syslib/" << machineId << "' and\n"
                  << "'ansys/vXXX/ansys/customize/misc/" << machineId << "'\n"
                  << "in your LD_LIBRARY_PATH variable!");
      }
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

    WriteComponentsFile(printGridOnly);
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
    WriteComponentsFile(false);
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
  
  void SimOutputRST::WriteComponentsFile(bool printGridOnly) 
  {
    // Finally, let's  write the  regions, named nodes  and named  elements as
    // components into an ANSYS components .cm file.
    std::stringstream sstr;
    StdVector<Elem*> elems;
    StdVector<RegionIdType> volRegions, surfRegions;
    // UInt numRegions;
    StdVector< std::string > regionNames;
    std::string rstFileName;
    std::string compFileName;
    
    sstr.clear();
    sstr.str("");
    sstr << dirName_ << sysPathSep_ << fileName_;
    if(printGridOnly) {
        sstr << "_grid";
    } else {
        sstr << "_ms" << msStep_;
    }
    rstFileName = sstr.str() + ".rst";
    compFileName = sstr.str() + ".cm";

    compFile_.open( compFileName.c_str(), std::ios::binary );

    if(compFile_.fail())
      EXCEPTION("Cannot open file '" << sstr.str()
                <<"' for writing components!");

    // Write file header with a little help info for the user.
    compFile_ << "! Components file for '" << rstFileName << "'." << std::endl;
    compFile_ << "! Read it into the ANSYS database using the following sequence of APDL commands:" << std::endl;
    compFile_ << "!" << std::endl << "! /POST1" << std::endl;
    compFile_ << "! INRES,ALL" << std::endl;
    compFile_ << "! FILE,'" << sstr.str() << "','rst','.'" << std::endl;
    compFile_ << "! SET,FIRST" << std::endl;
    compFile_ << "! CDREAD,DB,'" << sstr.str() << "','cm'" << std::endl;
    compFile_ << "! CMLIST" << std::endl << "!" << std::endl;
    compFile_ << "! The component/region information should also be available in ANSYS CFD-Post" << std::endl;
    compFile_ << "! when loading the .rst file." << std::endl;

    // numRegions = ptGrid_->GetNumRegions();
    ptGrid_->GetVolRegionIds( volRegions );
    ptGrid_->GetSurfRegionIds( surfRegions );

    ptGrid_->GetRegionNames(regionNames);
    // loop over regions and write out elements
    for(UInt i = 0, n=volRegions.GetSize(); i < n; i++)
    {
      UInt idx = volRegions[i];
      std::string regionName = regionNames[idx];
      
      ptGrid_->GetElems(elems, idx);
      UInt nElems = elems.GetSize();
      StdVector<UInt> elemNums(nElems);

      UInt j=0;
      for(j=0; j<nElems; j++) {
        elemNums[j] = elems[j]->elemNum;
      }

      WriteComponent(regionName, "ELEM", elemNums);
    }

    for(UInt i = 0, n=surfRegions.GetSize(); i < n; i++)
    {
      UInt idx = surfRegions[i];
      std::string regionName = regionNames[idx];
      
      ptGrid_->GetElems(elems, idx);
      UInt nElems = elems.GetSize();
      StdVector<UInt> elemNums(nElems);

      std::cout << "region: " << regionName << " " << nElems << std::endl;

      UInt j=0;
      for(j=0; j<nElems; j++) {
        elemNums[j] = elems[j]->elemNum;
      }

      WriteComponent(regionName, "ELEM", elemNums);
    }

    StdVector< UInt > entities;
    StdVector<std::string> entityNames;
    UInt numEntityGroups = 0;

    // obtain list with names of nodes
    ptGrid_->GetListNodeNames(entityNames);
    numEntityGroups = entityNames.GetSize();

    for(UInt i = 0; i < numEntityGroups; i++ ) {
      ptGrid_->GetNodesByName(entities, entityNames[i]);

      WriteComponent(entityNames[i], "NODE", entities);
    }

    // obtain list with names of nodes
    ptGrid_->GetListElemNames(entityNames);
    numEntityGroups = entityNames.GetSize();

    for(UInt i = 0; i < numEntityGroups; i++ ) {
      ptGrid_->GetElemsByName(elems, entityNames[i]);

      UInt j=0;
      UInt nElems = elems.GetSize();
      entities.Resize(nElems);
      for(j=0; j<nElems; j++) {
        entities[j] = elems[j]->elemNum;
      }

      WriteComponent(entityNames[i], "ELEM", entities);
    }

    compFile_.close();
  }
  
  void SimOutputRST::WriteComponent(const std::string& compName,
                                    const std::string& compType,
                                    const StdVector<UInt>& entities) 
  {
    char namebuf[33];
    UInt nEntities = entities.GetSize();

    int byteswritten = snprintf(namebuf, 32, "%s", compName.c_str());
    if(byteswritten < 32) 
    {
      std::fill(namebuf+byteswritten, namebuf+sizeof(namebuf), ' ');        
    }    
    namebuf[32] = 0;
    
    compFile_ << "CMBLOCK," << namebuf << "," << compType << ","
              << std::setw(8) << nEntities << std::endl;
    compFile_ << "(8i10)" << std::endl;
    
    UInt j=0;
    UInt m;
    for(j=0, m=nEntities/8; j<m; j++) {
      for(UInt k=0, o=8; k<o; k++) {
        compFile_ << std::setw(10) << entities[j*8+k];
      }
      compFile_ << std::endl;
    }
    for(UInt k=0, o=nEntities%8; k<o; k++) {
      compFile_ << std::setw(10) << entities[j*8+k];
    }
    compFile_ << std::endl;
  }
  
} // end of namespace
