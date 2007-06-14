// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "textSimOutput.hh"


#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include "Elements/basefe.hh"

namespace CoupledField {


  // initialize delimiter string
  std::string SimOutputText::delim_ = "  ";

  SimOutputText::SimOutputText( const std::string& fileName,
                                ParamNode * outputNode )
    : SimOutput( fileName, outputNode ){
    ENTER_FCN( "SimOutputText::SimOutputText" );

    // initialize variables
    formatName_ = "text";
    dirName_ = "history";
    fileName_ = fileName;
    coordSys_ = NULL;

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
      cmChar_ = (outputNode->Get("commentChar")->AsString())[0];

    // Get type of collection
    collecType_ = ENTITY;
    if ( outputNode ) {
      if( outputNode->Get("fileCollect")->AsString() == "timeFreq" )
        collecType_ = TIMEFREQ;
    }
  }

  SimOutputText::~SimOutputText() {
    ENTER_FCN( "SimOutputText::~SimOutputText" );
    
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
  
  void SimOutputText::Init( Grid * ptGrid ) {
    ENTER_FCN( "SimOutputText::Init" );
    ptGrid_ = ptGrid;

    // Get reference coordinate system
    if( myParam_ ) {
      std::string coordSysName;
      myParam_->Get("coordSysId", coordSysName);
      coordSys_ = domain->GetCoordSystem( coordSysName );
    }
  }

  void SimOutputText::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd ) {
    ENTER_FCN( "SimOutputText::RegisterResult" );

  }
  

  void SimOutputText::BeginStep( UInt stepNum, Double stepVal ) {
    ENTER_FCN( "SimOutputText::BeginStep" );

    actStep_ = stepNum;
    actStepVal_ = stepVal;
    resultMap_.clear();
  }


  void SimOutputText::AddResult( shared_ptr<BaseResult> sol ) {
    ENTER_FCN( "SimOutputText::AddResult" );

    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );
  }
  
  void SimOutputText::FinishStep( ) {
    ENTER_FCN( "SimOutputText::FinishStep" );

    // call correct function accordint fileCollectionType
    if( collecType_ == TIMEFREQ ) {
      WriteStepCollectTimeFreq();
    } else {
      WriteStepCollectEntity();
    }

  }

  void SimOutputText::WriteStepCollectTimeFreq() {
    ENTER_FCN( "SimInputText::WriteStepCollectTimeFreq" );
    
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
        std::string (SimOutputText::*pt2Func)(EntityIterator&) const;
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

        if( actResults[iSol]->GetEntryType() == EntryType::DOUBLE ) {

          // --- Real part ---
          Result<Double> & actRes = 
            dynamic_cast<Result<Double>& >( *(actResults[iSol]) );
          Vector<Double> & vec = actRes.GetVector();

          // get entity list of current result object
          shared_ptr<EntityList> actList  = actResults[iSol]->GetEntityList();
          
          // Iterate over all 'entities' of particular result
          ResultInfo::EntryType entryType =  actInfo.entryType;
          EntityIterator it = actList->GetIterator();
          UInt numDofs = actInfo.dofNames.GetSize();
          
          // Iterate over entities
          for( it.Begin(); !it.IsEnd(); it++ ) {

            // write node/elem/region number
            actFile << (this->*pt2Func)(it);
            
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
          ResultInfo::EntryType entryType =  actInfo.entryType;
          EntityIterator it = actList->GetIterator();
          UInt numDofs = actInfo.dofNames.GetSize();
              
          for( it.Begin(); !it.IsEnd(); it++ ) {
            // write node/elem/region info (number + evtl. coordinates )
            actFile << (this->*pt2Func)(it);
            
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
    ENTER_FCN( "SimInputText::WriteStepCollectEntity");

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
        ResultInfo::EntryType entryType =  actInfo.entryType;
        EntityIterator it = actList->GetIterator();
        UInt numDofs = actInfo.dofNames.GetSize();
        
        // *** Transient part ***
        if( actResults[iSol]->GetEntryType() == EntryType::DOUBLE ) {
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
  
  
  void SimOutputText::CreateFiles( shared_ptr<BaseResult> res,
                                   UInt step,
                                   Double stepVal ) {
    ENTER_FCN(" SimOutputText::CreateFiles" );
    
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

      for( it.Begin(); !it.IsEnd(); it++ ) {

        switch( actInfo.definedOn ) {
          
        case ResultInfo::NODE:
          entityString = GenStr(it.GetNode() );
          break;
        case ResultInfo::ELEMENT:
          entityString = GenStr(it.GetElem()->elemNum );
          break;
        case ResultInfo::SURF_ELEM:
          entityString = GenStr(it.GetSurfElem()->elemNum );
          break;
        case ResultInfo::REGION:
          entityString = ptGrid_->RegionIdToName( it.GetRegion() );
          break;
        case ResultInfo::SURF_REGION:
          entityString = ptGrid_->RegionIdToName( it.GetRegion() );
          break;
        default:
          EXCEPTION( "Not implemented!" );
        }

        // Concatenate 'simName-resultName-NODE/ELEM-#.hist'
        totalName = namePrefix + "-";
        totalName += entTypeString;
        totalName += "-";
    
        totalName += entityString + ".hist";
      
        outFiles[it.GetPos()] = new std::ofstream( totalName.c_str() );

        // write header to file
        std::ofstream & actFile = *outFiles[it.GetPos()];
        actFile << cmChar_ << " Result: '" << actInfo.resultName 
                << "' on " << entTypeString << "(s) '" 
                << list->GetName() << "'\n" << cmChar_ << "\n" << cmChar_;
        
        if( res->GetEntryType() == EntryType::DOUBLE ) {
          actFile << " t (s)  ";        
        } else {
          actFile << " f (Hz)  ";
        }
        actFile << ResultDofString( res) << "\n" << cmChar_ 
                << " ---------------------------------------------------------------------------\n";
      }      
      // CHECK, if creation of file succeeded
      outFiles_[res] = outFiles;

    } else {

      // =========================
      //  TIMEFREQ - CollectionType
      // =========================

      // create output stream;
      std::ofstream* outFile;
      
      // Compose name as 'simName-result-LISTNAME-STEPVAL_Hz/s.hist'
      shared_ptr<EntityList> list = res->GetEntityList();
      totalName = namePrefix + "-";
      totalName += list->GetName();
      totalName += "-" + GenStr(stepVal);
      if( res->GetEntryType() == EntryType::DOUBLE ) {
        totalName += "s";
      } else {
        totalName += "Hz";
      }
      totalName +=+ ".hist";

      // open stream
      outFile = new std::ofstream( totalName.c_str() );
      
      // Write header
      *outFile << cmChar_ << " Result: '" << actInfo.resultName 
              << "' on " << entTypeString << "(s) '" 
              << list->GetName();
      if( res->GetEntryType() == EntryType::DOUBLE ) {
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
    }
  }

  std::string SimOutputText::ResultDofString( shared_ptr<BaseResult> res ) {
    
    // extract result info object
    ResultInfo & actInfo = *(res->GetResultInfo());

    std::stringstream ret;

    if( res->GetEntryType() == EntryType::DOUBLE ) {
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

  std::string SimOutputText::GetNodeInfo( EntityIterator & it ) const {
    
    std::stringstream ret;

    ret << it.GetNode()  << delim_;

    Vector<Double> locCoord, globCoord;
    ptGrid_->GetNodeCoordinate( globCoord, it.GetNode() );
    if( coordSys_ != NULL ) {
      coordSys_->Global2LocalCoord( locCoord, globCoord );
      for( UInt iDim = 0; iDim < locCoord.GetSize()-1; iDim++ ) {
        ret <<  locCoord[iDim] << delim_;
      }
      ret << locCoord[locCoord.GetSize()-1];
    }
  return ret.str();

  }

  std::string SimOutputText::GetElemInfo( EntityIterator & it ) const {

    std::stringstream ret;
    ret << it.GetElem()->elemNum  << delim_;
    
    Vector<Double> elemLocCoord, locCoord, globCoord;
    
    Matrix<Double> cornerCoords;
    ptGrid_->GetElemNodesCoord( cornerCoords, it.GetElem()->connect, false);
    BaseFE * ptelem = it.GetElem()->ptElem;
    
    ptelem->GetCoordMidPoint( elemLocCoord );
    ptelem->Local2GlobalCoord(globCoord, elemLocCoord,  
                              cornerCoords, it.GetElem() );
    if( coordSys_ != NULL ) {
      coordSys_->Global2LocalCoord( locCoord, globCoord );
      for( UInt iDim = 0; iDim < locCoord.GetSize()-1; iDim++ ) {
        ret << locCoord[iDim] << delim_;
      }
      
      ret << locCoord[locCoord.GetSize()-1];
    }
    return ret.str();
  }

  std::string SimOutputText::GetRegionInfo( EntityIterator & it ) const {
    
    std::string ret;
    ret += GenStr( it.GetRegion() );
    return ret;
  }
  
  std::string SimOutputText::ComplexAsAmplPhase( const Complex& val ) const {
    
    return GenStr( std::abs(val) ) + delim_ + GenStr( CPhase(val) );
    
  }

  std::string SimOutputText::ComplexAsRealImag( const Complex& val ) const {
    return GenStr( val.real()) + delim_ + GenStr( val.imag() );
  }
}
