// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "textSimOutput.hh"


namespace CoupledField {


  SimOutputText::SimOutputText( const std::string& fileName,
                                ParamNode * outputNode )
    : SimOutput( fileName, outputNode ){
    ENTER_FCN( "SimOutputText::SimOutputText" );

    // initialize variables
    formatName_ = "text";
    dirName_ = "history";
    fileName_ = fileName;

    capabilities_.insert( HISTORY );

    // Create subdirectory 'history'
    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    // Get comment char 
    cmChar_ = '#';

  }


  SimOutputText::~SimOutputText() {
    ENTER_FCN( "SimOutputText::~SimOutputText" );

  }

  void SimOutputText::Init( Grid * ptGrid ) {
    ENTER_FCN( "SimOutputText::Init" );
    ptGrid_ = ptGrid;

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
              actOut << "  " << vec[it.GetPos()*numDofs + iDof];
            }
            actOut << std::endl;
          
          }
        
          // *** Harmonic part ***
        } else {
          Result<Complex> & actRes = 
            dynamic_cast<Result<Complex>& >( *(actResults[iSol]) );
          Vector<Complex> & vec = actRes.GetVector();
          if( actInfo.complexFormat == REAL_IMAG ) {
            for( it.Begin(); !it.IsEnd(); it++ ) {
              std::ofstream& actOut = *ptFiles[it.GetPos()];
              actOut << actStepVal_;;
              for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
                actOut << "  " << vec[it.GetPos()*numDofs + iDof].real() 
                       << " " << vec[it.GetPos()*numDofs + iDof].imag();
              }
              actOut << std::endl;
            
            }
          
          } else {
            for( it.Begin(); !it.IsEnd(); it++ ) {
              std::ofstream& actOut = *ptFiles[it.GetPos()];
              actOut << actStepVal_;;
              for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
                Complex & actVal = vec[it.GetPos()*numDofs + iDof];
                actOut << "  " << std::abs( actVal ) << " ";
                actOut << ((std::abs(actVal.imag()) > 1e-16) ?
                           std::atan2(actVal.imag(),actVal.real() )*180/PI : 0.0) ;
              }
              actOut << std::endl;
            
            } 
          }         
        
        }
      }
    }
  }
  
  
  void SimOutputText::CreateFiles( shared_ptr<BaseResult> res,
                                   UInt step,
                                   Double stepVal ) {
    ENTER_FCN(" SimOutputText::CreateFiles" );
    
    // Concatenate 'simName-resultName-NODE/ELEM-#.hist'
    std::string namePrefix="history/" + fileName_ + "-";
    std::string totalName, entityName;
    
    // determine type of entity the result is defined on
    ResultInfo & actDof = *(res->GetResultInfo());
    namePrefix += actDof.resultName;
    
    if( actDof.definedOn == ResultInfo::NODE
        || actDof.definedOn == ResultInfo::PFEM ) {
      entityName = "node";
    } else if( actDof.definedOn == ResultInfo::ELEMENT
               || actDof.definedOn == ResultInfo::SURF_ELEM ) {
      entityName = "elem";
    } else if( actDof.definedOn == ResultInfo::REGION ) {
      entityName = "region";
    } else if( actDof.definedOn == ResultInfo::SURF_REGION ) {
      entityName = "surfRegion";
    } else {
      EXCEPTION( "Not implemented!" );
    }

    
    // create new vector with ofstreams
    StdVector<std::ofstream*> outFiles;
    
    shared_ptr<EntityList> list = res->GetEntityList();
    outFiles.Resize( list->GetSize() );
    EntityIterator it = list->GetIterator();

    for( it.Begin(); !it.IsEnd(); it++ ) {
      
      totalName = namePrefix + "-";
      totalName += entityName;
      totalName += "-";
      std::string entityString, entTypeString;
      
      switch( actDof.definedOn ) {
        
      case ResultInfo::NODE:
        entityString = GenStr(it.GetNode() );
        entTypeString = "node(s)";
        break;
      case ResultInfo::ELEMENT:
        entityString = GenStr(it.GetElem()->elemNum );
        entTypeString = "element(s)";
        break;
      case ResultInfo::SURF_ELEM:
        entityString = GenStr(it.GetSurfElem()->elemNum );
        entTypeString = "surface element(s)";
        break;
      case ResultInfo::REGION:
        entityString = ptGrid_->RegionIdToName( it.GetRegion() );
        entTypeString = "region(s)";
        break;
      case ResultInfo::SURF_REGION:
        entityString = ptGrid_->RegionIdToName( it.GetRegion() );
        entTypeString = "surface region(s)";
        break;
      default:
        EXCEPTION( "Not implemented!" );
      }
      totalName += entityString + ".hist";
      
      outFiles[it.GetPos()] = new std::ofstream( totalName.c_str() );

      // write header to file
      std::ofstream & actFile = *outFiles[it.GetPos()];
      actFile << cmChar_ << " Result: '" << actDof.resultName 
              << "' on " << entTypeString << " '" 
              << list->GetName() << "'\n" << cmChar_ << "\n" << cmChar_;
      if( res->GetEntryType() == EntryType::DOUBLE ) {
        actFile << " t (s)  ";
        for( UInt i = 0; i < actDof.dofNames.GetSize(); i++ ) {
          actFile << actDof.dofNames[i]
                    << " (" << actDof.unit << ")  ";
        }
        actFile << "\n" 
                << cmChar_ << " ---------------------------------------------------------------------------\n";
      } else {
        actFile << " f (Hz)  ";
        if( actDof.complexFormat == REAL_IMAG ) {
          for( UInt i = 0; i < actDof.dofNames.GetSize(); i++ ) {
            if( actDof.dofNames[i] != "" ) {
              actFile << actDof.dofNames[i] << "-real"
                      << " (" << actDof.unit << ")  "
                      << actDof.dofNames[i] << "-imag"
                      << " (" << actDof.unit << ")  ";
            } else {
              actFile << "real"
                      << " (" << actDof.unit << ")  "
                      << "imag"
                      << " (" << actDof.unit << ")  ";              
            }
          }
        } else {
          for( UInt i = 0; i < actDof.dofNames.GetSize(); i++ ) {
            if( actDof.dofNames[i] != "" ) {
              actFile << actDof.dofNames[i] << "-ampl"
                      << " (" << actDof.unit << ")  "
                      << actDof.dofNames[i] << "-phase"
                      << " (deg)  ";
            } else {
              actFile << "ampl"
                      << " (" << actDof.unit << ")  "
                      << "phase"
                      << " (deg)  ";              
            }
          }
        }
        actFile << "\n"
                << cmChar_ << " ---------------------------------------------------------------------------\n";
      }
    }      
    // CHECK, if creation of file succeeded
    outFiles_[res] = outFiles;

  }
}
