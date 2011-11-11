#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "TextSimOutput.hh"
#include <fstream>


#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/BaseFE.hh"


namespace CoupledField {


  // initialize delimiter string
  std::string SimOutputText::delim_ = "  ";

  SimOutputText::SimOutputText( const std::string& fileName,
                                PtrParamNode outputNode )
    : SimOutput( fileName, outputNode ){

    // initialize variables
    formatName_ = "text";
    dirName_ = "history";
    fileName_ = fileName;
    coordSys_ = NULL;
    stepNumOffset_ = 0;
    stepValOffset_ = 0.0;
    globalNumbering_ = true;

    capabilities_.insert( HISTORY );

    // Create subdirectory 'history'
    try {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }
    
    // Get comment char 
    cmChar_ = '#';
    if( outputNode )
      cmChar_ = (outputNode->Get("commentChar")->As<std::string>())[0];

    // Get type of collection
    collecType_ = ENTITY;
    if ( outputNode ) {
      if( outputNode->Get("fileCollect")->As<std::string>() == "timeFreq" )
        collecType_ = TIMEFREQ;
      if( outputNode->Get("fileCollect")->As<std::string>() == "altogether" )
        collecType_ = ALTOGETHER;
    }
    
    // Get type of entity numbering
    if( outputNode ) {
      std::string numType;
      outputNode->GetValue( "entityNumbering", numType, ParamNode::EX );
      if( numType == "global" ) 
        globalNumbering_ = true;
      else
        globalNumbering_ = false;
    }
  }

  SimOutputText::~SimOutputText() {
    
    // Delete all 'outFiles', which might be still open
    std::map<shared_ptr<BaseResult>,
      StdVector<std::ofstream*> >::iterator it = outFiles_.begin();
    for( ; it != outFiles_.end(); it++ ) {
      
      StdVector<std::ofstream*> & actVec  = it->second;
      for( UInt i = 0; i < actVec.GetSize(); i++ ) {
        actVec[i]->close();
        delete actVec[i];
      }
    }
  }
  
  void SimOutputText::Init( Grid * ptGrid, bool printGridOnly ) {
    ptGrid_ = ptGrid;

    // Get reference coordinate system
    if( myParam_ ) {
      std::string coordSysName;
      myParam_->GetValue("coordSysId", coordSysName);
      coordSys_ = domain->GetCoordSystem( coordSysName );
    }
  }
  
  void SimOutputText::BeginMultiSequenceStep( UInt step, 
                                              BasePDE::AnalysisType type,
                                              UInt numSteps ) {
      
      actAnalysis_ = type;
    }
  
  void SimOutputText::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd,
                                      bool isHistory ) {

  }
  

  void SimOutputText::BeginStep( UInt stepNum, Double stepVal ) {

    actStep_ = stepNum;
    actStepVal_ = stepVal;

    // add  offset to step value to account for multisequence steps
    if( actAnalysis_ == BasePDE::TRANSIENT ||
        actAnalysis_ == BasePDE::STATIC  ) { 
      actStep_ += stepNumOffset_;
      actStepVal_ += stepValOffset_;
    }

    resultMap_.clear();
  }


  void SimOutputText::AddResult( shared_ptr<BaseResult> sol ) {

    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );
  }
  
  void SimOutputText::FinishStep( ) {

    // call correct function according to  fileCollectionType
    if( collecType_ == TIMEFREQ ) {
      WriteStepCollectTimeFreq();
    } else if ( collecType_ == ALTOGETHER ) {
      WriteStepCollectAltogether();
    } else {
      WriteStepCollectEntity();
    }

  }

  void SimOutputText::FinishMultiSequenceStep( ) {
    
    // set offset for step value and number to last values
    if( actAnalysis_ == BasePDE::TRANSIENT ||
         actAnalysis_ == BasePDE::STATIC ) {
      stepNumOffset_ = actStep_; 
      stepValOffset_ = actStepVal_;
     }
   }
  
  void SimOutputText::WriteStepCollectTimeFreq() {

    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {

      // get result info object and results for current result type
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;
      
      // Iterate over all 'BaseResults'
      for( UInt iSol = 0; iSol < actResults.GetSize(); iSol++ ) {
        
        // Check, if output stream for entity is already existant
        if( outFiles_.find( actResults[iSol]) == outFiles_.end() ) {
          CreateFiles( actResults[iSol], actStep_, actStepVal_ );
        } 
        
        std::ofstream  & actFile  = *outFiles_[actResults[iSol]][0];
        
        // Create function pointer, which returns for each entity
        // a number (node: nodenumber, element: element number,
        // region: region number...) and eventually the coordinates
        // w.r.t. to a local coordinate system
        StdVector<std::string> (SimOutputText::*pt2Func)(EntityIterator&) const;
        switch (actInfo.definedOn) {
        case ResultInfo::NODE:
          pt2Func = &SimOutputText::GetNodeInfo;
          break; 
        case ResultInfo::ELEMENT:
          pt2Func = &SimOutputText::GetElemInfo;
          break;
        case ResultInfo::SURF_ELEM:
          pt2Func = &SimOutputText::GetElemInfo;
          break;
        case ResultInfo::REGION:
          pt2Func = &SimOutputText::GetRegionInfo;
          break;
        case ResultInfo::SURF_REGION:
          pt2Func = &SimOutputText::GetRegionInfo;
          break;
        default:
          EXCEPTION( "Case not implemented" );
        }

        if( actResults[iSol]->GetEntryType() == BaseMatrix::DOUBLE ) {

          // --- Real part ---
          Result<Double> & actRes = 
            dynamic_cast<Result<Double>& >( *(actResults[iSol]) );
          Vector<Double> & vec = actRes.GetVector();

          // get entity list of current result object
          shared_ptr<EntityList> actList  = actResults[iSol]->GetEntityList();
          
          // Iterate over all 'entities' of particular result
          EntityIterator it = actList->GetIterator();
          UInt numDofs = actInfo.dofNames.GetSize();
          
          // Iterate over entities
          for( it.Begin(); !it.IsEnd(); it++ ) {

            // write node/elem/region number
            StdVector<std::string> loc = (this->*pt2Func)(it);
            actFile << loc[0];
            for( UInt i=1; i < loc.GetSize(); i++) {
              actFile << delim_ << loc[i];
            }
            
            // print value(s)
            for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
              actFile << delim_ << vec[it.GetPos()*numDofs + iDof];
            }
            actFile << std::endl;
            
          }
        } else {

          // --- Complex part ---
          Result<Complex> & actRes = 
            dynamic_cast<Result<Complex>& >( *(actResults[iSol]) );
          Vector<Complex> & vec = actRes.GetVector();

          // get entity list of current result object
          shared_ptr<EntityList> actList  = actResults[iSol]->GetEntityList();

          // create function pointer to function, which extract complex value
          // in correct format
          std::string (SimOutputText::*ptComplexOutput)(const Complex&) const;
          if( actInfo.complexFormat == REAL_IMAG ) {
            ptComplexOutput = &SimOutputText::ComplexAsRealImag;
          } else {
            ptComplexOutput = &SimOutputText::ComplexAsAmplPhase;
          }
          
          // Iterate over all 'entities' of particular result
          EntityIterator it = actList->GetIterator();
          UInt numDofs = actInfo.dofNames.GetSize();
              
          for( it.Begin(); !it.IsEnd(); it++ ) {
            // write node/elem/region info (number + evtl. coordinates )
            StdVector<std::string> loc = (this->*pt2Func)(it);
            actFile << loc[0];
            for( UInt i=1; i < loc.GetSize(); i++) {
              actFile << delim_ << loc[i];
            }
            
            // print value(s)
            for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
              actFile << delim_ << 
                (this->*ptComplexOutput)( vec[it.GetPos()*numDofs + iDof] );
            }
            actFile << std::endl;
          }
        }
        
        // close file and delete it
        actFile.close();
        delete &actFile;
        outFiles_.clear();
      }
      
    }
  }


  void SimOutputText::WriteStepCollectEntity() {

    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {

      // get result info object and results for current result type
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;
      
      // Iterate over all 'BaseResults'
      for( UInt iSol = 0; iSol < actResults.GetSize(); iSol++ ) {
        
        // Check, if output stream for entity is already existant
        if( outFiles_.find( actResults[iSol]) == outFiles_.end() ) {
          CreateFiles( actResults[iSol], actStep_, actStepVal_ );
        } 
        
        StdVector<std::ofstream*> & ptFiles  = outFiles_[actResults[iSol]];
        
        // get entity list of current result object
        shared_ptr<EntityList> actList  = actResults[iSol]->GetEntityList();
        
        // Iterate over all 'entities' of particular result
        EntityIterator it = actList->GetIterator();
        UInt numDofs = actInfo.dofNames.GetSize();
        
        // *** Transient part ***
        if( actResults[iSol]->GetEntryType() == BaseMatrix::DOUBLE ) {
          Result<Double> & actRes = 
            dynamic_cast<Result<Double>& >( *(actResults[iSol]) );
          Vector<Double> & vec = actRes.GetVector();
        
          for( it.Begin(); !it.IsEnd(); it++ ) {
            std::ofstream& actOut = *ptFiles[it.GetPos()];
            actOut << actStepVal_;;
            for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
              actOut << delim_ << vec[it.GetPos()*numDofs + iDof];
            }
            actOut << std::endl;
          
          }
        
          // *** Harmonic part ***
        } else {
          Result<Complex> & actRes = 
            dynamic_cast<Result<Complex>& >( *(actResults[iSol]) );
          Vector<Complex> & vec = actRes.GetVector();
          
          // create function pointer to function, which extract complex value
          // in correct format
          std::string (SimOutputText::*ptComplexOutput)(const Complex&) const;
          if( actInfo.complexFormat == REAL_IMAG ) {
            ptComplexOutput = &SimOutputText::ComplexAsRealImag;
          } else {
            ptComplexOutput = &SimOutputText::ComplexAsAmplPhase;
          }
          for( it.Begin(); !it.IsEnd(); it++ ) {
            std::ofstream& actOut = *ptFiles[it.GetPos()];
            actOut << actStepVal_;;
            for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
              actOut << delim_  << (this->*ptComplexOutput)
                ( vec[it.GetPos()*numDofs + iDof] ); 
            }
            actOut << std::endl;
          } // entities         
        } // harmonic part
      } // loop over results
    } // loop over result types
  }
  
  void SimOutputText::WriteStepCollectAltogether() {

    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {

      // get result info object and results for current result type
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;
      
      // Iterate over all 'BaseResults'
      for( UInt iSol = 0; iSol < actResults.GetSize(); iSol++ ) {
        
        // Check, if output stream for entity is already existant
        if( outFiles_.find( actResults[iSol]) == outFiles_.end() ) {
          CreateFiles( actResults[iSol], actStep_, actStepVal_ );
        } 
        
        std::ofstream  & actFile  = *outFiles_[actResults[iSol]][0];
        
        // get entity list of current result object
        shared_ptr<EntityList> actList  = actResults[iSol]->GetEntityList();
        
        // Iterate over all 'entities' of particular result
        EntityIterator it = actList->GetIterator();
        UInt numDofs = actInfo.dofNames.GetSize();
        
        // *** Transient part ***
        if( actResults[iSol]->GetEntryType() == BaseMatrix::DOUBLE ) {
          Result<Double> & actRes = 
            dynamic_cast<Result<Double>& >( *(actResults[iSol]) );
          Vector<Double> & vec = actRes.GetVector();
        
          actFile << actStepVal_;
          for( it.Begin(); !it.IsEnd(); it++ ) {
            for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
              actFile << delim_ << vec[it.GetPos()*numDofs + iDof];
            }
          }
          actFile << std::endl;
        
        // *** Harmonic part ***
        } else {
          Result<Complex> & actRes = 
            dynamic_cast<Result<Complex>& >( *(actResults[iSol]) );
          Vector<Complex> & vec = actRes.GetVector();
          
          // create function pointer to function, which extract complex value
          // in correct format
          std::string (SimOutputText::*ptComplexOutput)(const Complex&) const;
          if( actInfo.complexFormat == REAL_IMAG ) {
            ptComplexOutput = &SimOutputText::ComplexAsRealImag;
          } else {
            ptComplexOutput = &SimOutputText::ComplexAsAmplPhase;
          }
          
          actFile << actStepVal_;
          for( it.Begin(); !it.IsEnd(); it++ ) {
            for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
              actFile << delim_  << (this->*ptComplexOutput)
                ( vec[it.GetPos()*numDofs + iDof] ); 
            }
          }
          actFile << std::endl;
        } // harmonic part
      } // loop over results
    } // loop over result types
  }


  
  void SimOutputText::CreateFiles( shared_ptr<BaseResult> res,
                                   UInt step,
                                   Double stepVal ) {
    
    std::string namePrefix="history/" + fileName_ + "-";
    std::string totalName;
    
    // determine type of entity the result is defined on
    ResultInfo & actInfo = *(res->GetResultInfo());
    namePrefix += actInfo.resultName;
    
    std::string entityString, entTypeString;
    ResultInfo::Enum2String( actInfo.definedOn, entTypeString );
    
    
    // =========================
    //  ENTITY - CollectionType
    // =========================
    if( collecType_ == ENTITY ) {

      // create new vector with ofstreams
      StdVector<std::ofstream*> outFiles;
    
      shared_ptr<EntityList> list = res->GetEntityList();
      outFiles.Resize( list->GetSize() );
      EntityIterator it = list->GetIterator();

      UInt entityNum = 1;
      for( it.Begin(); !it.IsEnd(); it++, entityNum++ ) {

        std::string idString = it.GetIdString();
        std::string numString = lexical_cast<std::string>(entityNum);

        // now we have to switch the naming scheme
        switch (list->GetType() ) {
          case EntityList::NODE_LIST:
          case EntityList::ELEM_LIST:
          case EntityList::SURF_ELEM_LIST:
            // check naming scheme
            if( !globalNumbering_ ) {
              entityString = numString; // local numbering 
            } else {
              entityString = idString; // global node / element number
            }
            entityString += "-" + list->GetName();
            break;
          default:
            entityString = idString;
        }
        
        //entityString = it.GetIdString();
        //entityString = lexical_cast<std::string>(entityNum);
//        
//        if(list->GetType() == EntityList::NODE_LIST){
//          //ok so it seems we have history nodes.
//          //lets add the name of the nodelist to the filename
//          entityString += "-" + list->GetName();
//        }
//        switch( actInfo.definedOn ) {
//          
//        case ResultInfo::NODE:
//          entityString = lexical_cast<std::string>(it.GetNode() );
//          break;
//        case ResultInfo::PFEM:
//          entityString = lexical_cast<std::string>(it.GetNode() );
//          entTypeString="node";
//          break;
//        case ResultInfo::ELEMENT:
//          entityString = lexical_cast<std::string>(it.GetElem()->elemNum );
//          break;
//        case ResultInfo::SURF_ELEM:
//          entityString = lexical_cast<std::string>(it.GetSurfElem()->elemNum );
//          break;
//        case ResultInfo::REGION:
//          entityString = ptGrid_->RegionIdToName( it.GetRegion() );
//          break;
//        case ResultInfo::SURF_REGION:
//          entityString = ptGrid_->RegionIdToName( it.GetRegion() );
//          break;
//        default:
//          EXCEPTION( "Not implemented!" );
//        }

        // Concatenate 'simName-resultName-NODE/ELEM-#.hist'
        totalName = namePrefix + "-";
        totalName += entTypeString;
        totalName += "-";
    
        totalName += entityString + ".hist";
      
        if( usedFileNames_.find( totalName ) != usedFileNames_.end() ) {
          outFiles[it.GetPos()] = new std::ofstream( totalName.c_str(), 
                                                     std::ios::out | 
                                                     std::ios::app );
        } else {
          outFiles[it.GetPos()] = new std::ofstream( totalName.c_str() );
          usedFileNames_.insert( totalName );
        }

        // write header to file
        std::ofstream & actFile = *outFiles[it.GetPos()];
        actFile << cmChar_ << " Result: '" << actInfo.resultName 
                << "' on " << entTypeString << "(s) '" 
                << list->GetName() << "'";
        
        // now we have to switch the naming scheme
        switch (list->GetType() ) {
          case EntityList::NODE_LIST:
          case EntityList::ELEM_LIST:
          case EntityList::SURF_ELEM_LIST:
            actFile << " number #" << it.GetIdString();
            break;
          default:
            break;
        }
        actFile << "\n" << cmChar_;
        
        if( res->GetEntryType() == BaseMatrix::DOUBLE ) {
          actFile << " t (s)  ";        
        } else {
          actFile << " f (Hz)  ";
        }
        actFile << ResultDofString( res) << "\n" << cmChar_ 
                << " ---------------------------------------------------------------------------\n";
      }      
      // CHECK, if creation of file succeeded
      outFiles_[res] = outFiles;

    } else if( collecType_ == TIMEFREQ ) {

      // =========================
      //  TIMEFREQ - CollectionType
      // =========================

      // create output stream;
      std::ofstream* outFile;
      
      // Compose name as 'simName-result-LISTNAME-STEP.hist'
      shared_ptr<EntityList> list = res->GetEntityList();
      totalName = namePrefix + "-";
      totalName += list->GetName();
      totalName += "-" + lexical_cast<std::string>(step);
//      if( res->GetEntryType() == BaseMatrix::DOUBLE ) {
//        totalName += "s";
//      } else {
//        totalName += "Hz";
//      }
      totalName +=+ ".hist";

      // open stream
      outFile = new std::ofstream( totalName.c_str() );
      
      // Write header
      *outFile << cmChar_ << " Result: '" << actInfo.resultName 
              << "' on " << entTypeString << "(s) '" 
              << list->GetName();
      if( res->GetEntryType() == BaseMatrix::DOUBLE ) {
        *outFile << "' at time " << stepVal << " s";
      } else {
        *outFile << "' at frequency " << stepVal << " Hz";
      }
      *outFile << "\n" << cmChar_ << "\n";

      // write entityTypeString
      *outFile <<  cmChar_ << " " << entTypeString << delim_;
      
      // write coordinate system entries(only for node/elem/surfElem results)
      if( actInfo.definedOn == ResultInfo::NODE 
          || actInfo.definedOn == ResultInfo::ELEMENT
          || actInfo.definedOn == ResultInfo::SURF_ELEM ) {
        for (UInt iDim = 0; iDim < ptGrid_->GetDim(); iDim++ ) {
          *outFile << coordSys_->GetDofName(iDim+1) << delim_;
        }
      }
      
      // write dof result string
      *outFile << ResultDofString( res) << "\n" << cmChar_ 
              << " ---------------------------------------------------------------------------\n";

       // Push back file to outFiles_
      outFiles_[res].Resize(1);
      outFiles_[res][0] = outFile;
    } else {
      
      // ============================
      //  ALTOGETHER - CollectionType
      // ============================

      // create output stream;
      std::ofstream* outFile;
      
      // Compose name as 'simName-result-LISTNAME.hist'
      shared_ptr<EntityList> list = res->GetEntityList();
      totalName = namePrefix + "-";
      totalName += list->GetName();
      totalName +=+ ".hist";

      // open stream
      outFile = new std::ofstream( totalName.c_str() );
      
      // write header to file
      *outFile << cmChar_ << " Result: '" << actInfo.resultName 
               << "' on " << entTypeString << "(s) '" 
               << list->GetName() << "' with\n"
               << cmChar_ << "   line 1: " << entTypeString 
               << " number, line 2-3/2-4: " << entTypeString
               << " coordinates\n" << cmChar_ 
               << "   1st column is placeholder for equal no. of columns in all rows!\n";
      
      // write node/element number and coordinates
      StdVector<StdVector<std::string> > entityMatrix;
      StdVector<std::string> entityVector;
      
      EntityIterator it = list->GetIterator();
      for( it.Begin(); !it.IsEnd(); it++ ) {

        // pointer to function, which retrieves node/element number and coordinates
        StdVector<std::string> (SimOutputText::*pt2Func)(EntityIterator&) const;
        switch (actInfo.definedOn) {
        case ResultInfo::NODE:
          pt2Func = &SimOutputText::GetNodeInfo;
          break;
        case ResultInfo::ELEMENT:
          pt2Func = &SimOutputText::GetElemInfo;
          break;
        case ResultInfo::SURF_ELEM:
          pt2Func = &SimOutputText::GetElemInfo;
          break;
        default:
          EXCEPTION( "Case not implemented" );
        }
        // node/element number and coordinates 
        entityVector = (this->*pt2Func)(it);
        entityMatrix.Push_back(entityVector);
      }
      for( UInt iDim=0; iDim<entityVector.GetSize(); iDim++) {
        *outFile << " 0 ";
        for( UInt numEnt=0; numEnt<entityMatrix.GetSize(); numEnt++) {
          *outFile << delim_ << entityMatrix[numEnt][iDim];
        }
        *outFile << "\n";
      }
      
      // write dof result string
      if( res->GetEntryType() == BaseMatrix::DOUBLE ) {
        *outFile << cmChar_ << " t (s)  ";        
      } else {
        *outFile << cmChar_ << " f (Hz)  ";
      }
      *outFile << ResultDofString( res) << "\n" << cmChar_ 
              << " ---------------------------------------------------------------------------\n";

      // Push back file to outFiles_
      outFiles_[res].Resize(1);
      outFiles_[res][0] = outFile;
    }
  }

  std::string SimOutputText::ResultDofString( shared_ptr<BaseResult> res ) {
    
    // extract result info object
    ResultInfo & actInfo = *(res->GetResultInfo());

    std::stringstream ret;

    if( res->GetEntryType() == BaseMatrix::DOUBLE ) {
      for( UInt i = 0; i < actInfo.dofNames.GetSize(); i++ ) {
        ret << actInfo.dofNames[i]
            << " (" << actInfo.unit << ")  ";
      }
    } else {
      if( actInfo.complexFormat == REAL_IMAG ) {
        for( UInt i = 0; i < actInfo.dofNames.GetSize(); i++ ) {
          if( actInfo.dofNames[i] != "" ) {
            ret << actInfo.dofNames[i] << "-real"
                << " (" << actInfo.unit << ")  "
                << actInfo.dofNames[i] << "-imag"
                << " (" << actInfo.unit << ")  ";
          } else {
            ret << "real"
                << " (" << actInfo.unit << ")  "
                << "imag"
                << " (" << actInfo.unit << ")  ";              
          }
        }
      } else {
        for( UInt i = 0; i < actInfo.dofNames.GetSize(); i++ ) {
          if( actInfo.dofNames[i] != "" ) {
            ret << actInfo.dofNames[i] << "-ampl"
                    << " (" << actInfo.unit << ")  "
                    << actInfo.dofNames[i] << "-phase"
                << " (deg)  ";
          } else {
            ret << "ampl"
                << " (" << actInfo.unit << ")  "
                << "phase"
                << " (deg)  ";              
          }
        }
      }
    }

    return ret.str();
  }

  StdVector<std::string> SimOutputText::GetNodeInfo( EntityIterator & it ) const {
    
    StdVector<std::string> ret;

    ret.Push_back( lexical_cast<std::string>( it.GetNode() ) );

    Vector<Double> locCoord, globCoord;
    ptGrid_->GetNodeCoordinate( globCoord, it.GetNode() );
    if( coordSys_ != NULL ) {
      coordSys_->Global2LocalCoord( locCoord, globCoord );
      for( UInt iDim = 0; iDim < locCoord.GetSize(); iDim++ ) {
        ret.Push_back( lexical_cast<std::string>( locCoord[iDim] ) );
      }
    }
  return ret;

  }

  StdVector<std::string> SimOutputText::GetElemInfo( EntityIterator & it ) const {

    StdVector<std::string> ret;
    
    ret.Push_back( lexical_cast<std::string>( it.GetElem()->elemNum ) );
    Vector<Double> globCoord, locCoord;
    
    ptGrid_->GetElemShapeMap(it.GetElem())->GetGlobMidPoint(globCoord);
    
    if( coordSys_ != NULL ) {
      coordSys_->Global2LocalCoord( locCoord, globCoord );
      for( UInt iDim = 0; iDim < locCoord.GetSize(); iDim++ ) {
        ret.Push_back( lexical_cast<std::string>( locCoord[iDim] ));
      }
    }
    return ret;
  }

  StdVector<std::string> SimOutputText::GetRegionInfo( EntityIterator & it ) const {
    
    StdVector<std::string> ret;
    ret.Resize(1);
    
    ret[0] = lexical_cast<std::string>( it.GetRegion() );
    return ret;
  }
  
  std::string SimOutputText::ComplexAsAmplPhase( const Complex& val ) const {
    
    return lexical_cast<std::string>( std::abs(val) ) + delim_ + lexical_cast<std::string>( CPhase(val) );
    
  }

  std::string SimOutputText::ComplexAsRealImag( const Complex& val ) const {
    return lexical_cast<std::string>( val.real()) + delim_ + lexical_cast<std::string>( val.imag() );
  }
}
