#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <fstream>
#include <filesystem>

#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Exception.hh"
#include "SimOutputXDMF.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "FeBasis/BaseFE.hh"



namespace CoupledField {

  SimOutputXDMF::SimOutputXDMF(std::string fileName, PtrParamNode inputNode,
                               PtrParamNode infoNode, bool isRestart ) :
    SimOutput(fileName, inputNode, infoNode, isRestart ),
    simOutHDF5_(NULL)
  {
    // The restart case is currently not implemented, i.e. resuls from a 
    // partial simulation get overwritten.
    if( isRestart_ ) {
      WARN( "The XDMF-Writer is currently not adapted to write restarted "
            "results correctly, thus the results of the previous run get"
            " overwritten." );
    }

    fileName_ = fileName;
    formatName_ = "xdmf";
    dirName_ = "results_hdf5";
    
    std::string dirString = "results_" + formatName_; 
    inputNode->GetValue("directory", dirString, ParamNode::PASS );
    dirName_ = dirString; 
    
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    gridWritten_ = false;
    printGridOnly_ = false;

    currMS_ = 0;
    currStep_ = 0;

//     myParam_->GetValue("compressionLevel", compressionLevel, ParamNode::PASS );
//     if( compressionLevel > 9) {
//       EXCEPTION( "Value for compressionLevel must be betwen 1 and 9" );
//     }
//     myParam_->GetValue("maxChunkSize", maxChunkSize, ParamNode::PASS );
//     // Change defaults according to XML file
//     myParam_->GetValue("externalFiles", externalFiles_, ParamNode::PASS);
  }


  SimOutputXDMF::~SimOutputXDMF() {
    //    CloseFile();
  }

  void SimOutputXDMF::Init( Grid* ptGrid, bool printGridOnly ) {
    ptGrid_ = ptGrid;
    printGridOnly_ = printGridOnly;
    WriteGrid();
  }

  void SimOutputXDMF::BeginMultiSequenceStep( UInt step,
                                              BasePDE::AnalysisType type,
                                              UInt numSteps  ) {
//     std::stringstream msName;
//     H5::Group resultDescGroup;


//     currMSNumSteps_ = numSteps;

//     // If it does not exist, create Group for Grid / Volume data
//     try {
//       resultsGroup_ = mainGroup_.openGroup("Results");
//     } H5_CATCH( "Could open group for results" );

//     // Map analysistype
//     std::string analysisType = BasePDE::analysisType.ToString(type);
//     currAnalysisType_ = type;

//     // Assemble name of multistep
//     msName << "MultiStep_" << step;

//     // 1) Write group for mesh results
//     if( registeredMeshResults_.size() != 0 ) {
//       try{
//         meshResultsGroup_ = resultsGroup_.openGroup("Mesh");
//       } catch (H5::Exception&) {
//         try {
//           meshResultsGroup_ = resultsGroup_.createGroup("Mesh");
//           // write attribute indicating use of external files for simlation
//           // steps
//           H5IO::WriteAttribute( meshResultsGroup_, "ExternalFiles", externalFiles_ );
//         } H5_CATCH( "Could not create group for mesh results" );
//       }

//       try {
//         // create new multistep group.
//         try{ // first try to use an existing one, so we can call BeginMultiSequenceStep multiple times
//           currMSMeshGroup_ = meshResultsGroup_.openGroup(msName.str());
//         }catch (H5::Exception&){
//           currMSMeshGroup_ = meshResultsGroup_.createGroup(msName.str());

//           // add analysistype and number of steps to group
//           H5IO::WriteAttribute( currMSMeshGroup_, "AnalysisType", analysisType );

//           // add a group for the result description datasets.
//           resultDescGroup = currMSMeshGroup_.createGroup("ResultDescription");

//           // write result meta information to file
//           WriteResultDescriptions( resultDescGroup, numSteps, false);
//           resultDescGroup.close();
//         }
//       } H5_CATCH( "Could not create group for multi sequence step " << step );
//     }

//     // 2) Write group for history results
//     if( registeredHistResults_.size() != 0 ) {
//       try{
//         histResultsGroup_ = resultsGroup_.openGroup("History");
//       } catch (H5::Exception&) {
//         try {
//           histResultsGroup_ = resultsGroup_.createGroup("History");
//         } H5_CATCH( "Could not create group for history results" );
//       }

//       try {
//         try { // first try to use an existing one, so we can call BeginMultiSequenceStep multiple times
//           currMSHistGroup_ = histResultsGroup_.openGroup(msName.str());
//         }catch (H5::Exception&){
//           // create new multistep group.
//           currMSHistGroup_ = histResultsGroup_.createGroup(msName.str());

//           // add analysistype and number of steps to group
//           H5IO::WriteAttribute( currMSHistGroup_, "AnalysisType", analysisType );

//           // add a group for the result description datasets.
//           resultDescGroup = currMSHistGroup_.createGroup("ResultDescription");

//           // write result meta information to file
//           WriteResultDescriptions( resultDescGroup, numSteps, true );
//           resultDescGroup.close();

//           // iterate over all results
//           ResDescType::iterator it;
//           for( it = registeredHistResults_.begin();
//               it != registeredHistResults_.end();
//               it++ ) {

//             // create for each result a group within the ms group
//             H5::Group resultGroup;
//             try {
//               resultGroup = currMSHistGroup_.openGroup(it->first);
//             } catch( H5::Exception& ) {
//               resultGroup = currMSHistGroup_.createGroup(it->first);
//             }

//             // calculate number of "real steps"
//             UInt saveBegin, saveEnd, saveInc, numRealSteps;
//             saveBegin = std::max( (UInt) 1, histResultSaveBegin_[it->first] );
//             saveEnd = std::min( currMSNumSteps_, histResultSaveEnd_[it->first] );
//             saveInc = histResultSaveInc_[it->first];
//             numRealSteps = (UInt) (saveEnd-saveBegin) / saveInc + 1;

//             // create subgroup for entitytype
//             ResultInfo::EntityUnknownType definedOn
//             = it->second[0]->GetResultInfo()->definedOn;
//             std::string entityString = H5IO::MapUnknownTypeAsString(definedOn );
//             H5::Group entityTypeGroup = resultGroup.createGroup( entityString );

//             // iterate over all entitylists of result and create sub-subgroup
//             std::vector<shared_ptr<BaseResult> > const & lists = it->second;
//             for( UInt iList = 0; iList < lists.size(); iList++ )  {

//               // iterate over all entities in this list
//               EntityIterator entIt = lists[iList]->GetEntityList()->GetIterator();
//               UInt numDofs = lists[iList]->GetResultInfo()->dofNames.GetSize();
//               for( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
//                 H5::Group entityGroup; 
//                 try {
//                   entityGroup = entityTypeGroup.openGroup( entIt.GetIdString() );
//                   WARN("You are trying to add history entity '" << entIt.GetIdString()
//                        << "' under group '"
//                        << "History/" << msName.str() << "/" << it->first << "/" << entityString 
//                        << "'\nwhich already exists under a different name! Please check your mesh and XML files.");
//                   entityGroup.close();

//                   continue;
//                 } catch( H5::Exception& h5Ex ) {
//                   entityGroup = entityTypeGroup.createGroup( entIt.GetIdString() );
//                 }

//                 H5IO::Reserve2DArray<Double>(entityGroup, "Real", numRealSteps,
//                     numDofs, dPropList_ );
//                 if( lists[iList]->GetEntryType() == BaseMatrix::COMPLEX){
//                   H5IO::Reserve2DArray<Double>(entityGroup, "Imag", numRealSteps,
//                       numDofs, dPropList_ );
//                 }
//                 entityGroup.close();
//               }
//             }
//             entityTypeGroup.close();
//           }


//         }
//       } H5_CATCH( "Could not create group for multi sequence step " << step );
//     }
//     currMS_ = step;
  }

  void SimOutputXDMF::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd,
                                      bool isHistory ) {

//     std::string resultName = sol->GetResultInfo()->resultName;

//     if( !isHistory ) {
//       registeredMeshResults_[resultName].push_back(sol);
//       meshResultSaveBegin_[resultName] = saveBegin;
//       meshResultSaveEnd_[resultName] =  saveEnd;
//       meshResultSaveInc_[resultName] = saveInc;
//       usedCapabilities_.insert(MESH_RESULTS);
//     } else {
//       registeredHistResults_[resultName].push_back(sol);
//       histResultSaveBegin_[resultName] = saveBegin;
//       histResultSaveEnd_[resultName] = saveEnd;
//       histResultSaveInc_[resultName] = saveInc;
//       usedCapabilities_.insert(HISTORY);
//     }
  }

  void SimOutputXDMF::BeginStep( UInt stepNum, Double stepVal ) {
    
//     // always also call BeginMultiSequenceStep
//     // the MultiSequenceStep is closed at FinishStep
//     // the following was changed to reload the current MultiSequenceStep if it exists and only else create a new one
//     BeginMultiSequenceStep(currMS_, currAnalysisType_, currMSNumSteps_); // this recreates the multi-sequence step environment

//     currStep_ = stepNum;
//     currStepValue_ = stepVal;

  }


  void SimOutputXDMF::AddResult( shared_ptr<BaseResult> sol )
  {
//     std::string resultName = sol->GetResultInfo()->resultName;

//     // try to determine, if current result is a history or mesh result
//     bool isHistory = false;
//     if( registeredHistResults_.find(resultName) !=
//       registeredHistResults_.end() ) {
//       if( std::find( registeredHistResults_[resultName].begin(),
//                      registeredHistResults_[resultName].end(),
//                      sol) !=
//                        registeredHistResults_[resultName].end() ) {
//         isHistory = true;
//       }
//     }

//     if( !isHistory) {
//       AddMeshResult( sol );
//     } else {
//       AddHistResult( sol );
//     }
  }

  void SimOutputXDMF::AddMeshResult( shared_ptr<BaseResult> sol) {

//     H5::Group resultGroup, subGroup, regionGroup;
//     UInt numDOFs;
//     Vector<Double> realVec, imagVec;
//     std::vector<std::string> resultNames;

//     std::string regionName = sol->GetEntityList()->GetName();
//     shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
//     std::string resultName = resInfo->resultName;
//     numDOFs = resInfo->dofNames.GetSize();

//     // check, if step group is already open
//     if( currMeshStepGroup_.getId() <= 0 ) {
//       std::stringstream stepName;
//       stepName << "Step_" << currStep_;

//       // Create new step group.
//       try {
//         // TODO: strieben - Creation of step groups fails when converting from h5 to h5 (see cube2d). Do something about it!
//         currMeshStepGroup_= currMSMeshGroup_.createGroup(stepName.str());
//         H5IO::WriteAttribute( currMeshStepGroup_, "StepValue", currStepValue_ );

//         if(externalFiles_ )
//           CreateExternalFile();
//       } H5_CATCH( "Can not create dataset for step " << currStep_ );
//     }

//     // Add current stepvalue to related group in result description,
//     // if not yet present
//     bool writeStep = false;
//     StdVector<UInt> & myStepNums = meshResultStepNums_[resultName];
//     StdVector<Double> & myStepVals = meshResultStepVal_[resultName];
//     if( myStepNums.GetSize() == 0 ) {
//       writeStep = true;
//     } else {
//       if( myStepNums.Last() != currStep_ )
//         writeStep = true;
//     }


//     if( writeStep ) {
//       myStepNums.Push_back( currStep_ );
//       myStepVals.Push_back( currStepValue_ );
//       try {
//         H5::Group resDescGroup =
//           currMSMeshGroup_.openGroup( "ResultDescription").openGroup(resultName);
//         H5IO::SetEntries1DArray(resDescGroup, "StepValues", myStepNums.GetSize()-1,
//                                 myStepNums.GetSize()-1,
//                                 &myStepVals[myStepNums.GetSize()-1] );
//         H5IO::SetEntries1DArray(resDescGroup, "StepNumbers", myStepNums.GetSize()-1,
//                                         myStepNums.GetSize()-1,
//                                         &myStepNums[myStepNums.GetSize()-1] );
//         resDescGroup.close();
//       } H5_CATCH( "Could not write current step value for result '"
//                   << resultName << "'" );
//     }

//     // check, if result was already written
//     try {
//       resultGroup = currMeshStepGroup_.openGroup( resultName );
//     } catch( H5::Exception& ) {
//       resultGroup = currMeshStepGroup_.createGroup( resultName );
//     }

//     // determine, on whicht type of entity the result is defined
//     std::string entityString;
//     switch( resInfo->definedOn ) {
//     case ResultInfo::NODE:
//     case ResultInfo::PFEM:
//       entityString = "Nodes";
//       break;
//     case ResultInfo::ELEMENT:
//     case ResultInfo::SURF_ELEM:
//       entityString = "Elements";
//       break;
//     default:
//       EXCEPTION( "Result of type '" << resInfo->resultName
//                   << "'can not be written as mesh result" );
//     }

//     // try to create regionGroup
//     try {
//       regionGroup = resultGroup.openGroup( regionName );
//     } catch( H5::Exception& ) {
//       regionGroup = resultGroup.createGroup( regionName );
//     }

//     // try to create subgroup for entity
//     try {
//       subGroup = regionGroup.createGroup( entityString );
//     } H5_CATCH( "Could not create subgroup " << entityString
//                 << " for result " << resultName << " on region "
//                 << regionName );

//     if( sol->GetEntryType() == BaseMatrix::DOUBLE ) {

//       Vector<Double> & resultVec = dynamic_cast<Result<Double>&>
//       (*sol).GetVector();

//       WriteResults(subGroup, resultVec, numDOFs, false);
//     } else {
//       Vector<Complex> & resultVec = dynamic_cast<Result<Complex>&>
//       (*sol).GetVector();


//       realVec.Resize( resultVec.GetSize() );
//       imagVec.Resize( resultVec.GetSize() );

//       for(UInt i = 0; i < resultVec.GetSize(); i++) {
//         realVec[i] = resultVec[i].real();
//         imagVec[i] = resultVec[i].imag();
//       }
//       WriteResults(subGroup, realVec, numDOFs, false);
//       WriteResults(subGroup, imagVec, numDOFs, true);

//     }

//     // close groups
//     subGroup.close();
//     regionGroup.close();
//     resultGroup.close();
  }

  void SimOutputXDMF::AddHistResult( shared_ptr<BaseResult> sol) {

//     shared_ptr<ResultInfo> resInfo = sol->GetResultInfo();
//     std::string resultName = resInfo->resultName;
//     UInt numDofs = resInfo->dofNames.GetSize();
//     std::string entityString = H5IO::MapUnknownTypeAsString(resInfo->definedOn );


//     // Add current stepvalue to related group in result description,
//     // if not yet present
//     bool writeStep = false;
//     StdVector<UInt> & myStepNums = histResultStepNums_[resultName];
//     StdVector<Double> & myStepVals = histResultStepVal_[resultName];
//     if( myStepNums.GetSize() == 0 ) {
//       writeStep = true;
//     } else {
//       if( myStepNums.Last() != currStep_ )
//         writeStep = true;
//     }
//     if( writeStep ) {
//       myStepNums.Push_back( currStep_ );
//       myStepVals.Push_back( currStepValue_ );

//       try {
//         H5::Group resDescGroup =
//           currMSHistGroup_.openGroup( "ResultDescription").openGroup(resultName);
//         H5IO::SetEntries1DArray(resDescGroup, "StepValues", myStepNums.GetSize()-1,
//                                 myStepNums.GetSize()-1,
//                                 &myStepVals[myStepNums.GetSize()-1] );
//         H5IO::SetEntries1DArray(resDescGroup, "StepNumbers", myStepNums.GetSize()-1,
//                                         myStepNums.GetSize()-1,
//                                         &myStepNums[myStepNums.GetSize()-1] );
//         resDescGroup.close();
//       } H5_CATCH( "Could not write current step value for result '"
//                   << resultName << "'" );
//     }
//     // ---------------------
//     //  Write result itself
//     // ---------------------
//     try {
//       H5::Group resultGroup = currMSHistGroup_.openGroup( resultName);
//       H5::Group entityTypeGroup = resultGroup.openGroup( entityString );


//       // iterate over all entities in this list
//       EntityIterator entIt = sol->GetEntityList()->GetIterator();
//       UInt pos = 0;
//       for( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
//         H5::Group entityGroup =
//           entityTypeGroup.openGroup( entIt.GetIdString() );

//         if( sol->GetEntryType() == BaseMatrix::DOUBLE ){
//           Vector<Double> & resultVec = dynamic_cast<Result<Double>&>
//           (*sol).GetVector();
//           H5IO::SetEntries2DArray( entityGroup, "Real", myStepNums.GetSize()-1,
//                                    myStepNums.GetSize()-1, 0, numDofs-1,
//                                    &resultVec[pos] );
//           pos += numDofs;
//         } else {
//           Vector<Complex> & resultVec = dynamic_cast<Result<Complex>&>
//           (*sol).GetVector();
//           Vector<Double> realVec(numDofs), imagVec(numDofs);

//           for( UInt i = 0; i < numDofs; i++ ) {
//             realVec[i] = resultVec[pos+i].real();
//             imagVec[i] = resultVec[pos+i].imag();
//           }

//           H5IO::SetEntries2DArray( entityGroup, "Real", myStepNums.GetSize()-1,
//                                    myStepNums.GetSize()-1, 0, numDofs-1,
//                                    &realVec[0] );
//           H5IO::SetEntries2DArray( entityGroup, "Imag", myStepNums.GetSize()-1,
//                                    myStepNums.GetSize()-1, 0, numDofs-1,
//                                    &imagVec[0] );
//           pos += numDofs;
//         }
//         entityGroup.close();
//       }
//       entityTypeGroup.close();
//       resultGroup.close();
//     } H5_CATCH( "Could not open history result group for result '"
//                 << resultName << "'" );
  }


  void SimOutputXDMF::FinishStep( )
  {
//     if(externalFiles_)
//     {
//       PtrParamNode in = info->Get("analysis/output/externalFile");
//       try {
//         in->Get("name")->SetValue(currStepFile_.getFileName());
//         in->Get("size")->SetValue((int) currStepFile_.getFileSize());
//         info->ToFile();
//       } catch (H5::FileIException &h5ex) {}
//     }

//    // we close everything here, so that the file is in a consistent state
//     currStepFile_.close();

//     if( currMeshStepGroup_.getId() > 0 )
//       currMeshStepGroup_.close();

//     if( currHistStepGroup_.getId() > 0 )
//       currHistStepGroup_.close();
    
//     if( registeredMeshResults_.size() > 0 ) {

//       // write attributes containing stepNumber and stepValue
//       // of last simulation step
//       H5IO::WriteAttribute( currMSMeshGroup_, "LastStepNum", currStep_ );
//       H5IO::WriteAttribute( currMSMeshGroup_, "LastStepValue", currStepValue_ );

//       currMSMeshGroup_.close();
//       meshResultsGroup_.close();
//     }

//     if( registeredHistResults_.size() > 0 ) {
//       // write attributes containing stepNumber and stepValue
//       // of last simulation step
//       H5IO::WriteAttribute( currMSHistGroup_, "LastStepNum", currStep_ );
//       H5IO::WriteAttribute( currMSHistGroup_, "LastStepValue", currStepValue_ );

//       currMSHistGroup_.close();
//       histResultsGroup_.close();
//     }
    
//     resultsGroup_.close();
    
//     CloseFile();
    
//     // everything is closed
//     // now reopen the file and basic groups
//     // the multisequencestep is reload in BeginStep again
    
//     OpenFile(false);
  }

  void SimOutputXDMF::FinishMultiSequenceStep( ) {
    
//     // all closing is already done in FinishStep
    
//     registeredMeshResults_.clear();

//     // reset all data per sequence step
//     meshResultSaveBegin_.clear();
//     meshResultSaveEnd_.clear();
//     meshResultSaveInc_.clear();
//     histResultSaveBegin_.clear();
//     histResultSaveEnd_.clear();
//     histResultSaveInc_.clear();
//     meshResultStepNums_.clear();
//     histResultStepNums_.clear();
//     meshResultStepVal_.clear();
//     histResultStepVal_.clear();
    
  }

  void SimOutputXDMF::Finalize() {

//     // Write file header
//     WriteFileInfoHeader();

//     // return, if only the grid is to be printed
//     if( printGridOnly_ )
//       return;

//     // return, if no commandLine handler or
//     // global root ParaemNode are present
//     if( !progOpts || !param )
//       return;

//     std::vector<std::string> fileNames;
//     std::vector<std::string> dataSetNames;
//     std::ifstream fin;
//     std::ostringstream dumpStr;

//     fileNames.push_back(progOpts->GetParamFileStr());
//     fileNames.push_back(myParam_->GetRoot()->Get("fileFormats")
//                         ->Get("materialData")
//                         ->Get("file")->As<std::string>());

//     dataSetNames.push_back("ParameterFile");
//     dataSetNames.push_back("MaterialFile");

//     for(UInt i=0; i<fileNames.size(); i++)
//     {
//       fin.open( fileNames[i].c_str(), std::ios::binary );

//       if(fin.fail())
//         EXCEPTION("Cannot open file '" << fileNames[i]
//                   <<"' to dump into XDMF!");

//       // seek to the end of the file
//       fin.seekg (0, std::ios::end);
//       UInt numBytes = fin.tellg();
//       fin.seekg (0, std::ios::beg);

//       std::string str;
//       str.resize(numBytes);
//       fin.read(&str[0], numBytes);
//       WriteStringToUserData(dataSetNames[i], str);
//       fin.close();
//     }

//     progOpts->GetVersionString( dumpStr, false );
//     WriteStringToUserData( "ProgramStats", dumpStr.str() );
  }

  void SimOutputXDMF::InitModule() {
//     std::string pathsep;
//     std::string fileName;
//     std::ostringstream strBuffer;

//     // concatenate output file name
//     try {
//       fs::create_directory( dirName_ );
//       pathsep = fs::path("/").native_directory_string();
//     } catch (std::exception &ex) {
//       EXCEPTION(ex.what());
//     }

//     strBuffer << dirName_ << pathsep << fileName_ << ".h5";
//     currFileName_ = strBuffer.str();
    
//     OpenFile(true);
  }
  
  void SimOutputXDMF::WriteGrid() {
    std::string gridName = fileName_+std::string("_grid.xmf");
    fs::path gridFN = dirName_ / gridName;
    //std::string gridFN = dirName_ + "/" + fileName_ + 
    std::ofstream gridFile;

    // ensure that grid gets only written once
    if(!gridWritten_)
      InitModule();
    else
      return;

    WriteDTD();

    gridFile.open(gridFN);
    gridFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    gridFile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl << std::endl;
    gridFile << "<Xdmf>" << std::endl;
    gridFile << "  <Domain>" << std::endl;
    gridFile << "    <!-- Global node coordinate array -->" << std::endl;
    gridFile << "    <DataItem Name=\"NodeCoords\" Dimensions=\"" << ptGrid_->GetNumNodes()
             << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
    gridFile << "      " << fileName_ << ".h5:/Mesh/Nodes/Coordinates" << std::endl;
    gridFile << "    </DataItem>" << std::endl << std::endl;

    WriteRegions(gridFile);
    
    gridFile << "  </Domain>" << std::endl;
    gridFile << "</Xdmf>" << std::endl;

    gridFile.close();
    


//     // write the dimension of the grid.
//     H5IO::WriteAttribute( meshGroup_, "Dimension", ptGrid_->GetDim() );

//     // ================
//     //  Node Locations
//     // ================
//     UInt nNodes = ptGrid_->GetNumNodes();
//     H5::Group nodeGroup;
//     try {
//       nodeGroup = meshGroup_.createGroup( "Nodes" );
//     } H5_CATCH( "Could not create node group" );

//     H5IO::WriteAttribute( nodeGroup, "NumNodes", nNodes );

//     // collect all nodal coordinates
//     std::vector<Double> locs( nNodes * 3 );
//     for (UInt i = 0; i < nNodes; i++) {
//       Point p;
//       ptGrid_->GetNodeCoordinate(p, i+1);
//       locs[i*3+0] = p[0];
//       locs[i*3+1] = p[1];
//       locs[i*3+2] = p[2];
//     }

//     H5IO::Write2DArray( nodeGroup, "Coordinates", nNodes,
//                         3, &locs[0], dPropList_ );

//     nodeGroup.close();


//     // =====================
//     //  Element definitions
//     // =====================
//     UInt nElems = ptGrid_->GetNumElems();
//     H5::Group elemGroup;
//     try{
//       elemGroup = meshGroup_.createGroup("Elements");
//     } H5_CATCH( "Could not create element group" );

//     UInt maxNumNodes = ptGrid_->GetMaxNumNodesPerElem();
//     std::vector<UInt> connect (nElems * maxNumNodes);
//     std::vector<UInt> elConnect;
//     std::vector<Integer> feTypes (nElems);
//     std::vector< UInt > numElemsOfDim ( 3 );

//     // Fill connectivity array
//     std::fill( connect.begin(), connect.end(), 0 );
//     UInt offset;
//     Elem::FEType eType;
//     RegionIdType region;

//     // iterate over all elements
//     for( UInt i = 0; i < nElems; i++ ) {
//       elConnect.resize( maxNumNodes );
//       std::fill(elConnect.begin(), elConnect.end(),
//                 0 );
//       ptGrid_->GetElemData( i+1, eType, region, &elConnect[0] );
//       numElemsOfDim[Elem::GetElemDim(eType)-1]++;
//       feTypes[i] = eType;

//       // insert connectivity into global array
//       offset = i * maxNumNodes;
//       for( UInt j = 0; j < elConnect.size(); j++ ) {
//         connect[offset + j] = elConnect[j];
//       }
//     }

//     // write connectivity
//     H5IO::Write2DArray( elemGroup, "Connectivity", nElems,
//                         maxNumNodes, &connect[0], dPropList_ );

//     // write element types
//     H5IO::Write1DArray( elemGroup, "Types", nElems,
//                         &feTypes[0], dPropList_ );


//     // ==========================
//     //  Grid Meta Information
//     // ==========================

//     H5IO::WriteAttribute( elemGroup, "NumElems", nElems );

//     H5IO::WriteAttribute( elemGroup, "QuadraticElems",
//                           ptGrid_->IsQuadratic() );

//     // number of elements per dimension
//     for(UInt i=0; i<3; i++) {
//       std::stringstream attrName;
//       attrName << "Num" << (i+1) << "DElems";
//       H5IO::WriteAttribute( elemGroup, attrName.str(), numElemsOfDim[i] );
//     }

//     // number of elements per type
//     UInt numElemTypes = Elem::feType.map.size();
//     for(UInt i=0; i<numElemTypes; i++) {
//       std::stringstream attrName;
//       attrName << "Num_" << Elem::feType.ToString((Elem::FEType)i);
//       H5IO::WriteAttribute( elemGroup, attrName.str(),
//                             ptGrid_->GetNumElemOfType((Elem::FEType)i) );
//     }

//     // close element group
//     elemGroup.close();

//     // ============================================
//     //  Write Regions, Node Groups, Element Groups
//     // ============================================

//     WriteRegions( meshGroup_ );

//     // create new group for entity groups
//     H5::Group groupsGroup;
//     try{
//       groupsGroup = meshGroup_.createGroup("Groups");
//     } H5_CATCH( "Could not create mesh regiongroup" );

//     WriteNodeGroups( groupsGroup );
//     WriteElemGroups( groupsGroup );
//     groupsGroup.close();

//     gridWritten_ = true;
//     usedCapabilities_.insert(MESH);

//     // close meshGroup
//     meshGroup_.close();
  }


  void SimOutputXDMF::WriteRegions(std::ofstream& gridFile) {
//     H5::Group regionListGroup;
     StdVector< std::string > regionNames;
     StdVector< UInt > regionDims;
     StdVector< std::vector<UInt> > regionElems;
     StdVector< StdVector<UInt> > regionNodes;
     StdVector<Elem*> elems;
     StdVector<RegionIdType> surfRegionIds, volRegionIds, regionIds;
     UInt dim, numRegions;
     Integer idx;

     numRegions = ptGrid_->GetNumRegions();
     if(!numRegions)
       return;
     
     dim = ptGrid_->GetDim();
     ptGrid_->GetVolRegionIds(volRegionIds);
     ptGrid_->GetSurfRegionIds(surfRegionIds);
     ptGrid_->GetRegionNames(regionNames);
     regionDims.Resize(numRegions);
     regionElems.Resize(numRegions);
     regionNodes.Resize(numRegions);
     regionDims.Resize(numRegions);

     for(UInt i=0, n=volRegionIds.GetSize(); i<n; i++) {
       regionIds.Push_back(volRegionIds[i]);
     }
     for(UInt i=0, n=surfRegionIds.GetSize(); i<n; i++) {
       regionIds.Push_back(surfRegionIds[i]);
     }
     
//     // obtain nodes and elements per surface region
     for(UInt i=0, n=regionIds.GetSize(); i<n; i++) {
       idx = regionIds[i];
       regionDims[idx] = dim-1;

       ptGrid_->GetElems(elems, idx);
       UInt nElems = elems.GetSize();
       UInt numConn = 0;
       regionElems[idx].resize(nElems);

       UInt nElemNodes = 0;
       bool wNumNodes = false;
       Elem::FEType elemType = Elem::ET_UNDEF;
       UInt et=0;
       RegionIdType elemRegion=0;
       std::vector<UInt> connect;
       connect.resize(100);

       for(UInt j=0; j<nElems; j++) {
         regionElems[idx][j] = elems[j]->elemNum;

         ptGrid_->GetElemData(elems[j]->elemNum, elemType, elemRegion, &connect[0]);

         MapElemType(elemType, et, nElemNodes, wNumNodes);

         numConn += nElemNodes;

         if(wNumNodes) numConn++;
       }

       ptGrid_->GetNodesByRegion(regionNodes[idx], idx);

//       // determine dimensionality
//       if( ptGrid_->GetDim() == 3 ) {
//         regionDims[idx] = 2;
//       } else {
//         regionDims[idx] = 1;
//       }

       gridFile << "    <!-- Region: " << regionNames[idx] << " -->" << std::endl;

       UInt numRegionNodes = regionNodes[idx].GetSize();
       gridFile << "    <!-- Global node numbers belonging to region '" << regionNames[idx] << "' -->" << std::endl;
       gridFile << "    <DataItem Name=\"Nodes_" << regionNames[idx] << "\" Dimensions=\"" << numRegionNodes << "\" ItemType=\"Function\" Function=\"$0 - 1\">" << std::endl;
       gridFile << "      <DataItem Dimensions=\"" << numRegionNodes << "\" NumberType=\"UInt\" Precision=\"4\" Format=\"HDF\">" << std::endl;
       gridFile << "        " << fileName_ << ".h5:/Mesh/Regions/" << regionNames[idx] << "/Nodes" << std::endl;
       gridFile << "      </DataItem>" << std::endl;
       gridFile << "    </DataItem>" << std::endl << std::endl;

       gridFile << "    <!-- Node coordinate array for region '" << regionNames[idx] << "' -->" << std::endl;
       gridFile << "    <DataItem Name=\"NodeCoords_" << regionNames[idx] << "\" Dimensions=\"" << numRegionNodes << " 3\" ItemType=\"Function\" Function=\"JOIN($0, $1, $2)\">" << std::endl;
       for(UInt d=0; d<3; d++) 
       {
         gridFile << "      <DataItem Dimensions=\"" << numRegionNodes << "\" ItemType=\"Function\" Function=\"$0[($1 * 3) + " << d << "]\">" << std::endl;
         gridFile << "        <DataItem Reference=\"XML\" Dimensions=\"" << numRegionNodes << "\">" << std::endl;
         gridFile << "          /Xdmf/Domain/DataItem[@Name=\"NodeCoords\"]" << std::endl;
         gridFile << "        </DataItem>" << std::endl;
         gridFile << "        <DataItem Reference=\"XML\" Dimensions=\"" << numRegionNodes << "\">" << std::endl;
         gridFile << "          /Xdmf/Domain/DataItem[@Name=\"Nodes_" << regionNames[idx] << "\"]" << std::endl;
         gridFile << "        </DataItem>" << std::endl;
         gridFile << "      </DataItem>" << std::endl;
       }
       gridFile << "    </DataItem>" << std::endl << std::endl;


       std::map<UInt, UInt> glob2LocNode;       
       for(UInt n=0; n<numRegionNodes; n++) 
       {
         glob2LocNode[regionNodes[idx][n]] = n;         
       }

       gridFile << "    <!-- Element connectivities for region '" << regionNames[idx] << "' -->" << std::endl;
       gridFile << "    <DataItem Name=\"Connectivity_" << regionNames[idx] << "\" Format=\"XML\" NumberType=\"UInt\" Dimensions=\"" << (numConn + nElems) << "\">" << std::endl;
       for(UInt j=0; j<nElems; j++) {
         ptGrid_->GetElemData(elems[j]->elemNum, elemType, elemRegion, &connect[0]);

         MapElemType(elemType, et, nElemNodes, wNumNodes);

         gridFile << "      " << et << "   ";
         if(wNumNodes) gridFile << nElemNodes << "   ";
         for(UInt n=0; n<nElemNodes; n++) 
         {
           gridFile << glob2LocNode[elems[j]->connect[n]] << " ";
         }
         gridFile << std::endl;
       }
       
       gridFile << "    </DataItem>" << std::endl << std::endl;

       gridFile << "    <!-- Grid for region '" << regionNames[idx] << "' -->" << std::endl;
       gridFile << "    <Grid Name=\"" << regionNames[idx] << "\" GridType=\"Uniform\">" << std::endl;
       gridFile << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << nElems << "\">" << std::endl;
       gridFile << "        <DataItem Reference=\"XML\" Dimensions=\"" << nElems << "\">" << std::endl;
       gridFile << "          /Xdmf/Domain/DataItem[@Name=\"Connectivity_" << regionNames[idx] << "\"]" << std::endl;
       gridFile << "        </DataItem>" << std::endl;
       gridFile << "      </Topology>" << std::endl;
       gridFile << "      <Geometry GeometryType=\"XYZ\">" << std::endl;
       gridFile << "        <DataItem Reference=\"XML\" Dimensions=\"" << numRegionNodes << " 3\">" << std::endl;
       gridFile << "          /Xdmf/Domain/DataItem[@Name=\"NodeCoords_" << regionNames[idx] << "\"]" << std::endl;
       gridFile << "        </DataItem>" << std::endl;
       gridFile << "      </Geometry>" << std::endl;
       gridFile << "    </Grid>" << std::endl << std::endl;
       
     }

     
//     // obtain nodes and elements per volume region
//     for(UInt i=0, n=volRegionIds.GetSize(); i<n; i++) {
//       idx = volRegionIds[i];
//       regionDims[idx] = dim;

//       ptGrid_->GetElems(elems, idx);
//       UInt nElems = elems.GetSize();
//       regionElems[idx].resize(nElems);
//       for(UInt j=0; j<nElems; j++) {
//         regionElems[idx][j] = elems[j]->elemNum;
//       }

//       ptGrid_->GetNodesByRegion(regionNodes[idx], idx);

//       // determine dimensionality
//       regionDims[idx] = ptGrid_->GetDim();
//     }


//     // loop over regions and write out nodes and elements
//     for(UInt i = 0; i < numRegions; i++)
//     {
//       // create new region group
//       H5::Group actRegionGroup;
//       try {
//         actRegionGroup = regionListGroup.createGroup(regionNames[i] );
//       } H5_CATCH( "Could not create region group for region '"
//                   << regionNames[i] << "'" );
//       H5IO::WriteAttribute( actRegionGroup, "Dimension",
//                             regionDims[i] );

//       // create new node group
//       H5IO::Write1DArray<UInt>( actRegionGroup, "Nodes",
//                           regionNodes[i].GetSize(),
//                           (const UInt*)&regionNodes[i][0], dPropList_ );

//       // create new element group
//       H5IO::Write1DArray( actRegionGroup,
//                           "Elements",
//                           regionElems[i].size(),
//                           (const Integer*)&regionElems[i][0],
//                           dPropList_);

//       // create new face group
//       // .. to be implemented ..

//       // create new edge group
//       // .. to be implemented ..

//       // close current region group
//       actRegionGroup.close();
//     }

//     // close regionlist group
//     regionListGroup.close();


  }

  void SimOutputXDMF::WriteNodeGroups() {
//     H5::Group myGroup;
//     StdVector< UInt > nodes;
//     StdVector<std::string> nodeNames;
//     UInt numNodeGroups = 0;

//     // obtain list with names of nodes
//     ptGrid_->GetListNodeNames(nodeNames);
//     numNodeGroups = nodeNames.GetSize();

//     for(UInt i = 0; i < numNodeGroups; i++ ) {
//       ptGrid_->GetNodesByName(nodes, nodeNames[i]);

//       // try to open group with given name
//       try {
//         myGroup = meshGroup.openGroup( nodeNames[i] );
//       } catch (H5::Exception& ) {
//         myGroup = meshGroup.createGroup( nodeNames[i] );
//       }
//       H5IO::WriteAttribute( myGroup, "Dimension", (Integer) 0 );
//       H5IO::Write1DArray( myGroup, "Nodes",
//                           nodes.GetSize(), &nodes[0], dPropList_ );

//       // close nodes array of current group
//       myGroup.close();
//     }
  }

  void SimOutputXDMF::WriteElemGroups() {
//     H5::Group myGroup;
//     StdVector< UInt > elemNums, elemNodes;
//     StdVector<Elem*> elems;
//     StdVector<std::string> elemNames;
//     std::set<UInt> nodeSet;

//     // obtain list with names of elements
//     ptGrid_->GetListElemNames(elemNames);
//     UInt numElemGroups = elemNames.GetSize();

//     for(UInt i = 0; i < numElemGroups; i++ ) {
//       ptGrid_->GetElemsByName(elems, elemNames[i]);
//       elemNums.Resize( elems.GetSize() );
//       nodeSet.clear();
      
//       std::set<UInt> dims;

//       for( UInt j = 0; j < elems.GetSize(); j++ ) {
//         elemNums[j] = elems[j]->elemNum;
//         dims.insert( elems[j]->ptElem->GetDim() );
//         nodeSet.insert( elems[j]->connect.Begin(),
//                         elems[j]->connect.End() );
//       }
//       if( dims.size() > 1 ) {
//         EXCEPTION( "Element group '" << elemNames[i]
//                     << "' contains elements of different dimensions" );
//       }
//       try {
//         myGroup = meshGroup.openGroup( elemNames[i] );
//       } catch (H5::Exception& ) {
//         myGroup = meshGroup.createGroup( elemNames[i] );
//       }
//       H5IO::WriteAttribute( myGroup, "Dimension", *dims.begin() );
//       H5IO::Write1DArray( myGroup, "Elements",
//                           elemNums.GetSize(), &elemNums[0],
//                           dPropList_);

//       // Write nodes of element group
//       elemNodes.Resize( nodeSet.size() );
//       std::copy( nodeSet.begin(), nodeSet.end(), elemNodes.Begin() );
//       H5IO::Write1DArray( myGroup, "Nodes",
//                           elemNodes.GetSize(), &elemNodes[0],
//                           dPropList_);

//       // close nodes array of current group
//       myGroup.close();
//     }
  }

//   void SimOutputXDMF::WriteResultDescriptions( const H5::Group& resGroup,
//                                                UInt numSteps,
//                                                bool isHistory ) {
//     std::string resultName, unit;
//     UInt definedOn, numDofs, entryType, saveBegin, saveEnd, saveInc;
//     std::vector<std::string> entityNames;
//     ResDescType::const_iterator it, end;
//     std::vector< shared_ptr<BaseResult> >::const_iterator solIt, solEnd;
//     shared_ptr<ResultInfo> resInfo, actResInfo;

//     if( !isHistory ) {
//       it = registeredMeshResults_.begin();
//       end = registeredMeshResults_.end();
//     } else {
//       it = registeredHistResults_.begin();
//       end = registeredHistResults_.end();
//     }

//     for( ; it != end; it++ ) {
//       resultName = it->first;

//       solIt = it->second.begin();
//       solEnd = it->second.end();

//       resInfo = (*solIt)->GetResultInfo();
//       numDofs = resInfo->dofNames.GetSize();
//       unit = resInfo->unit;
//       definedOn = H5IO::MapUnknownType( resInfo->definedOn );
//       entryType = H5IO::MapEntryType( resInfo->entryType );
//       if( !isHistory ) {
//         saveBegin = std::max( (UInt)1, meshResultSaveBegin_[resultName]);
//         saveEnd = std::min( currMSNumSteps_, meshResultSaveEnd_[resultName] );
//         saveInc = meshResultSaveInc_[resultName];
//       } else {
//         saveBegin = std::max( (UInt) 1, histResultSaveBegin_[resultName] );
//         saveEnd = std::min( currMSNumSteps_, histResultSaveEnd_[resultName] );
//         saveInc = histResultSaveInc_[resultName];
//       }

//       // Generate list of entityNames for the current result.
//       entityNames.clear();
//       for( ; solIt != solEnd; solIt++ ) {
//         actResInfo = (*solIt)->GetResultInfo();

//         entityNames.
//           push_back((*solIt)->GetEntityList()->GetName());
//       }

//       // Reset solIt to beginning of result vector
//       solIt = it->second.begin();

//       try {

//         // Generate compound datatype

//         /* // First version: Compound data type
//            H5IO::CompoundType resInfo;
//            typedef std::pair<std::string, boost::any> CEntryType;
//            resInfo.push_back( CEntryType( "DefinedOn", definedOn ) );
//            resInfo.push_back( CEntryType( "Regions", regions ) );
//            resInfo.push_back( CEntryType( "NumDOFs", numDOFs ) );
//            resInfo.push_back( CEntryType( "DOFNames", dofNames ) );
//            resInfo.push_back( CEntryType( "EntryType", entryType ) );
//            resInfo.push_back( CEntryType( "Unit", unit ) );

//            H5IO::WriteCompound( currAttrDescGroup_, resNames[0], resInfo );
//         */

//         // === Second version: Separate datasets for each entry
//         H5::Group actGroup = resGroup.createGroup(resultName );

//         H5IO::Write1DArray( actGroup, "DefinedOn", 1, &definedOn, dPropList_ );
//         H5IO::Write1DArray( actGroup, "EntityNames", entityNames.size(),
//                             &entityNames[0], dPropList_ );
//         H5IO::Write1DArray( actGroup, "NumDOFs", 1, &numDofs, dPropList_ );
//         H5IO::Write1DArray( actGroup, "DOFNames", resInfo->dofNames.GetSize(),
//                             &(resInfo->dofNames[0]), dPropList_ );
//         H5IO::Write1DArray( actGroup, "EntryType", 1, &entryType, dPropList_ );
//         H5IO::Write1DArray( actGroup, "Unit", 1, &unit, dPropList_ );

//         UInt numRealSteps = (UInt) (saveEnd-saveBegin) / saveInc + 1;

//         // do not reserve 1D array for StepValues and StepNumbers but initialize with 0
//         // otherwise strange memory bugs occur! :(
//         StdVector<Double> tmp_double;
//         tmp_double.Resize(numRealSteps, 0.0);
//         H5IO::Write1DArray( actGroup, "StepValues", numRealSteps, tmp_double.GetPointer(), dPropList_ );

//         StdVector<UInt> tmp_uint;
//         tmp_uint.Resize(numRealSteps, 0);
//         H5IO::Write1DArray( actGroup, "StepNumbers", numRealSteps, tmp_uint.GetPointer(), dPropList_ );

//         actGroup.close();

//       } H5_CATCH( "Could not write result description for result '"
//                   << resultName << "'" );
//     }
//   }

  void SimOutputXDMF::WriteDTD() 
  {
    fs::path dtdFN = dirName_ / "Xdmf.dtd";
    std::ofstream dtdFile(dtdFN);
    
    dtdFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    dtdFile << "<!--Root element of dataset-->" << std::endl;
    dtdFile << "<!ELEMENT Xdmf (Information*, Domain+)>" << std::endl;
    dtdFile << "<!ATTLIST Xdmf" << std::endl;
    dtdFile << "    Version CDATA #IMPLIED" << std::endl;
    dtdFile << ">" << std::endl;
    dtdFile << "<!--Contains domain data information-->" << std::endl;
    dtdFile << "<!ELEMENT Domain (Information*, Grid+)>" << std::endl;
    dtdFile << "<!ATTLIST Domain" << std::endl;
    dtdFile << "	Name CDATA #IMPLIED" << std::endl;
    dtdFile << ">" << std::endl;
    dtdFile << "<!--Contains a collection of homogeneous elements-->" << std::endl;
    dtdFile << "<!ELEMENT Grid (Information*, Time*, Grid*, Topology*, Geometry*, Attribute*)>" << std::endl;
    dtdFile << "<!ATTLIST Grid" << std::endl;
    dtdFile << "	Name CDATA #IMPLIED" << std::endl;
    dtdFile << "	GridType (Uniform | Collection | Tree | Subset) \"Uniform\"" << std::endl;
    dtdFile << "    CollectionType   (Spatial | Temporal) \"Spatial\"" << std::endl;
    dtdFile << "    Section (DataItem | All) \"All\"" << std::endl;
    dtdFile << ">" << std::endl;
    dtdFile << "<!-- Described Temporal Relationship -->" << std::endl;
    dtdFile << "<!ELEMENT Time (Information*, DataItem*)>" << std::endl;
    dtdFile << "<!ATTLIST Time" << std::endl;
    dtdFile << "	TimeType (Single | HyperSlab | List | Range) \"Single\"" << std::endl;
    dtdFile << "	Value CDATA #IMPLIED" << std::endl;
    dtdFile << ">" << std::endl;
    dtdFile << "<!--Describes the general organization of the data-->" << std::endl;
    dtdFile << "<!ELEMENT Topology (Information*, DataItem*)>" << std::endl;
    dtdFile << "<!ATTLIST Topology" << std::endl;
    dtdFile << "	TopologyType (Polyvertex | Polyline | Polygon | Triangle | Quadrilateral | Tetrahedron | Pyramid| Wedge | Hexahedron | Edge_3 | Triagle_6 | Quadrilateral_8 | Tetrahedron_10 | Pyramid_13 | Wedge_15 | Hexahedron_20 | Mixed | 2DSMesh | 2DRectMesh | 2DCoRectMesh | 3DSMesh | 3DRectMesh | 3DCoRectMesh) #REQUIRED" << std::endl;
    dtdFile << "	Dimensions CDATA #IMPLIED" << std::endl;
    dtdFile << "	Order CDATA #IMPLIED" << std::endl;
    dtdFile << "	BaseOffset CDATA \"0\"" << std::endl;
    dtdFile << "	NodesPerElement CDATA #IMPLIED" << std::endl;
    dtdFile << "	NumberOfElements CDATA #IMPLIED" << std::endl;
    dtdFile << ">" << std::endl;
    dtdFile << "<!--Describes the XYZ values of the mesh-->" << std::endl;
    dtdFile << "<!ELEMENT Geometry (Information*, DataItem+)>" << std::endl;
    dtdFile << "<!ATTLIST Geometry" << std::endl;
    dtdFile << "	Name CDATA #IMPLIED" << std::endl;
    dtdFile << "	GeometryType (XYZ | XY | X_Y_Z | VXVYVZ | ORIGIN_DXDYDZ) \"XYZ\"" << std::endl;
    dtdFile << ">" << std::endl;
    dtdFile << "<!--Lowest level element, describes the data that is present in the XDMF dataset-->" << std::endl;
    dtdFile << "<!ELEMENT DataItem (#PCDATA | DataItem)*>" << std::endl;
    dtdFile << "<!ATTLIST DataItem" << std::endl;
    dtdFile << "	Name CDATA #IMPLIED" << std::endl;
    dtdFile << "	ItemType (Uniform | Collection | Tree | HyperSlab | Coordinates | Function) \"Uniform\"" << std::endl;
    dtdFile << "	Dimensions CDATA #REQUIRED" << std::endl;
    dtdFile << "	NumberType (Char | UChar | Float | Int | UInt) \"Float\"" << std::endl;
    dtdFile << "	Precision (1 | 4 | 8) \"4\"" << std::endl;
    dtdFile << "	Reference CDATA #IMPLIED" << std::endl;
    dtdFile << "	Endian (Big | Little | Native) \"Native\"" << std::endl;
    dtdFile << "	Compression (Zlib |BZip2 | Raw) \"Raw\"" << std::endl;
    dtdFile << "	Format (XML | HDF | Binary) \"XML\"" << std::endl;
    dtdFile << ">" << std::endl;
    dtdFile << "<!--Describes the values on the mesh-->" << std::endl;
    dtdFile << "<!ELEMENT Attribute (Information*, DataItem)>" << std::endl;
    dtdFile << "<!ATTLIST Attribute" << std::endl;
    dtdFile << "	Name CDATA #IMPLIED" << std::endl;
    dtdFile << "	Center (Node | Cell | Grid | Face | Edge) \"Node\"" << std::endl;
    dtdFile << "	AttributeType (Scalar | Vector | Tensor | Tensor6 | Matrix) \"Scalar\"" << std::endl;
    dtdFile << ">" << std::endl;
    dtdFile << "<!-- Application Dependent -->" << std::endl;
    dtdFile << "<!ELEMENT Information (#PCDATA | Information | EMPTY)*>" << std::endl;
    dtdFile << "<!ATTLIST Information " << std::endl;
    dtdFile << "	Name CDATA #IMPLIED" << std::endl;
    dtdFile << "	Value CDATA #IMPLIED" << std::endl;
    dtdFile << ">" << std::endl;

    dtdFile.close();
  }

#define XDMF_NOTOPOLOGY     0x0
#define XDMF_POLYVERTEX     0x1
#define XDMF_POLYLINE       0x2
#define XDMF_POLYGON        0x3
#define XDMF_TRI            0x4
#define XDMF_QUAD           0x5
#define XDMF_TET            0x6
#define XDMF_PYRAMID        0x7
#define XDMF_WEDGE          0x8
#define XDMF_HEX            0x9
#define XDMF_EDGE_3         0x0022
#define XDMF_TRI_6          0x0024
#define XDMF_QUAD_8         0x0025
#define XDMF_TET_10         0x0026
#define XDMF_PYRAMID_13     0x0027
#define XDMF_WEDGE_15       0x0028
#define XDMF_WEDGE_18       0x0029
#define XDMF_HEX_20         0x0030
#define XDMF_HEX_24         0x0031
#define XDMF_HEX_27         0x0032

  void SimOutputXDMF::MapElemType(Elem::FEType eType, UInt& et,
                                  UInt& numNodes, bool& wNodeNodes)
  {
    wNodeNodes = false;
    et = XDMF_NOTOPOLOGY;
    numNodes = 0;
    
    switch(eType) 
    {
    case Elem::ET_UNDEF:
      break;
    case Elem::ET_POINT:
      et = XDMF_POLYVERTEX;
      numNodes = 1;
      wNodeNodes = true;
      break;
    case Elem::ET_LINE2:
      et = XDMF_POLYLINE;
      numNodes = 2;
      wNodeNodes = true;
      break;
    case Elem::ET_LINE3:
      et = XDMF_EDGE_3;
      numNodes = 3;           
      break;
    case Elem::ET_TRIA3:
      et = XDMF_TRI;
      numNodes = 3;
      break;
    case Elem::ET_TRIA6:
      et = XDMF_TRI_6;
      numNodes = 6;
      break;
    case Elem::ET_QUAD4:
      et = XDMF_QUAD;
      numNodes = 4;
      break;
    case Elem::ET_QUAD8:
      et = XDMF_QUAD_8;
      numNodes = 8;
      break;
    case Elem::ET_QUAD9:
      et = XDMF_QUAD_8;
      numNodes = 8;      
      break;
    case Elem::ET_TET4:
      et = XDMF_TET;
      numNodes = 4;
      break;
    case Elem::ET_TET10:
      et = XDMF_TET_10;
      numNodes = 10;
      break;
    case Elem::ET_HEXA8:
      et = XDMF_HEX;
      numNodes = 8;
      break;
    case Elem::ET_HEXA20:
      et = XDMF_HEX_20;
      numNodes = 20;
      break;
    case Elem::ET_HEXA27:
      et = XDMF_HEX_20;
      numNodes = 20;
      break;
    case Elem::ET_PYRA5:
      et = XDMF_PYRAMID;
      numNodes = 5;
      break;
    case Elem::ET_PYRA13:
      et = XDMF_PYRAMID_13;
      numNodes = 13;
      break;
    case Elem::ET_PYRA14:
      et = XDMF_PYRAMID_13;
      numNodes = 13;
      break;
    case Elem::ET_WEDGE6:
      et = XDMF_WEDGE;
      numNodes = 6;
      break;
    case Elem::ET_WEDGE15:
      et = XDMF_WEDGE_15;
      numNodes = 15;
      break;
    case Elem::ET_WEDGE18:
      et = XDMF_WEDGE_18;
      numNodes = 18;
      break;
    }
    
  }

} // end of namespace
