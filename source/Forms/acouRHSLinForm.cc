// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <map>
#include <math.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>


// #include "DataInOut/infiles.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/simInput.hh"
#include "Domain/entityList.hh"
#include "Driver/singleDriver.hh"

#include "Utils/result.hh"
#include "Utils/boost-serialization.hh"
#include "DataInOut/resultHandler.hh"
#include "acouRHSLinForm.hh"

namespace CoupledField {

  // ====================================================================
  // volume source integration
  // ====================================================================

  AcouRHSLinForm::AcouRHSLinForm(std::string id,
                                 std::string regionName,
                                 std::string factor,
                                 bool interpolate,
                                 std::string srcInputId,
                                 std::string srcRegions,
                                 std::string coordSysId,
                                 Double globalEpsilon,
                                 Double localEpsilon,
                                 std::string restartFileMode,
                                 std::string asynchSteps,
                                 bool node_warnings)
    : LinearForm(),
      id_(id),
      regionName_(regionName),
      interpolate_(interpolate),
      srcInputId_(srcInputId),
      restartFileMode_(restartFileMode),
      coordSysId_(coordSysId),
      globalEpsilon_(globalEpsilon),
      localEpsilon_(localEpsilon),
      asynchSteps_(asynchSteps),
      node_warnings_(node_warnings),
      destNodeList_(NULL),
      destElemList_(NULL),
      step_(0),
      lastStep_(0)
  {
    name_ = "AcouRHSLinForm";
    
    // Determine analysis type
    analysistype_ = domain->GetSingleDriver()->GetAnalysisType();
    
    // For asynchronous steps we need the time values
    if (asynchSteps_ != "no") {
      mph_tf_ = mParser_->GetNewHandle();
      if (analysistype_ == BasePDE::TRANSIENT) {
        // Get a handle to math parser's time variable
        mParser_->SetExpr(mph_tf_, "t");
      }
      else {
        // Get a handle to math parser's frequency variable
        mParser_->SetExpr(mph_tf_, "f");
        
        // linear interpolation does not make sense here
        if (asynchSteps_ == "interpolate") {
          *warning << "For a harmonic analysis it makes no sense to use "
            << "linear frequency step interpolation; switching to nearest "
            << "neighbor.";
          Warning(__FILE__, __LINE__);
          asynchSteps_ = "nearest";
        }
      }
        
      // Find the ResultInfo for acouRhsValues
      shared_ptr<SimInput> inputReader = domain->GetResultHandler()
          ->GetInputReader(interpolate_ ? srcInputId_ : id_);
      StdVector< shared_ptr<ResultInfo> > resInfos;
      inputReader->GetResultTypes(1, resInfos);
      UInt numResults = resInfos.GetSize();
      UInt iRes;
      for (iRes = 0; iRes < numResults; ++iRes) {
        if (resInfos[iRes]->resultName == "acouRhsLoad") {
          break;
        }
      }
      if (iRes >= resInfos.GetSize()) {
        EXCEPTION("acouRhsLoad not found in input file (id=" << id_ << ")");
      }
      
      // Retrieve the time values of the stored acouRhsValues
      inputReader->GetStepValues(1, resInfos[iRes], stepValues_);
    }
    
    // Set Expressions for math parser
    
    // Specify that we want to use the current step number.
    mParser_->SetExpr( mHandle_, "step" );

    // Specify the ramping factor function for the RHS to the parser.
    mHandle2_ =  mParser_->GetNewHandle();
    mParser_->SetExpr( mHandle2_, factor );
    //std::cout << "********** factor is: " << factor << std::endl;

    typedef boost::tokenizer<char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");
    Tok t(srcRegions, sep);
    Tok::iterator it, end;
    it = t.begin();
    end = t.end();
    
    for( ; it != end; it++)
      srcRegions_.Push_back(*it);

    UInt numSourceRegions = srcRegions_.GetSize();

    sourceNodeLists_.Resize(numSourceRegions);
    sourceNodeLists_.Init(NULL);

    consInterpWeights_.Resize(numSourceRegions);
    
    /*
    std::cout << "id " <<              id_ << std::endl;
    std::cout << "regionName " <<       regionName_ << std::endl;
    std::cout << "interpolate " <<      interpolate_ << std::endl;
    std::cout << "srcInputId " <<       srcInputId_ << std::endl;
    std::cout << "srcRegion " <<        srcRegion_ << std::endl;
    std::cout << "coordSysId " <<       coordSysId_ << std::endl;
    std::cout << "globalEpsilon " <<    globalEpsilon_ << std::endl;
    std::cout << "localEpsilon " <<     localEpsilon_ << std::endl;
    */
  }

  AcouRHSLinForm::~AcouRHSLinForm()
  {
    mParser_->ReleaseHandle( mHandle2_ );
    if (asynchSteps_ != "no")
      mParser_->ReleaseHandle(mph_tf_);
  }

  void AcouRHSLinForm::CalcElemVector( Vector<Double> & elemVec,
                                       EntityIterator& ent ) 
  {
    ResultHandler* resultHandler = domain->GetResultHandler();
    Grid* source;
    Grid* dest;
    shared_ptr<BaseResult> acouRHSVal, acouRHSVal2;
    bool timeInterp = (asynchSteps_ == "interpolate");
    Double factor;
    Double intFactor = 1.0;
    Double t = 0.0;
    const Double timeTol = 0.001;
    UInt prevStep = 0;
    
    step_ = (UInt) mParser_->Eval(mHandle_);
    factor = mParser_->Eval( mHandle2_ );
    
    // Determine asynchronous time step
    if (asynchSteps_ != "no") {
      // get the current time
      t = mParser_->Eval(mph_tf_);
      
      // we have to use an iterator in order to tackle non-contiguous step
      // numbering (saveInc>1)
      std::map< UInt, Double >::iterator it = stepValues_.lower_bound(lastStep_);
      step_ = it->first;
      
      // find first step with a larger time value than the current time
      while (it != stepValues_.end()) {
        if (it->second >= t)
          break;
        step_ = (++it)->first;
      }
   
      if (it == stepValues_.end()) {
        step_ = (--it)->first;
      }

      if (asynchSteps_ == "nearest") { // nearest neighbor interpolation
        // check if the previous step is closer to the current time
        Double dt = fabs(it->second - t); 
        if (dt > fabs(t - (--it)->second))
          step_ = it->first;
      }
      else { // compute factor for linear interpolation
        Double t2 = it->second;
        if (it == stepValues_.begin()) {
          intFactor = 1.0; // no interpolation before first time step
        }
        else {
          --it;
          prevStep = it->first;
          intFactor = (t - it->second) / (t2 - it->second);
          
          // choose just one step, if we are very close to it
          if (intFactor < timeTol) {
            step_ = prevStep;
            intFactor = 1.0;
          }
          else if (intFactor > 1.0) {
            intFactor = 1.0;
            if (intFactor > 1.0 + timeTol) {
              *warning << "Current time/frequency step is beyond the last "
                       << "step of input file (id='" << srcInputId_ << "')";
              Warning(__FILE__, __LINE__);
            }
          }
        }
      }
    }
    
    
    // Determine if values for a new step have to be read.
    
    if (step_ != lastStep_)
    {
      if (interpolate_)
      {
#ifdef USE_INTERPOLATION
        UInt numSourceRegions = srcRegions_.GetSize();
        StdVector<UInt> regionNodes;

        // do we need data of previous time step for interpolation?
        if (timeInterp) {
          // are the data still in memory?
          if (prevStep == lastStep_) {
            // only if previous step gives a significant contribution
            if (fabs(1.0-intFactor) > timeTol) {
              rhsValues2_.Resize(rhsValues_.GetSize());
              rhsValues2_ = rhsValues_; // just copy it
            }
          }
        }
        
        // initialize vector for current time step
        dest = domain->GetGrid();
        dest->GetNodesByRegion( regionNodes, dest->RegionNameToId(regionName_));
        rhsValues_.Resize(regionNodes.GetSize());
        std::fill(rhsValues_.Begin(), rhsValues_.End(), 0.0);
        
        if (timeInterp) {
          rhsValues2_.Resize(regionNodes.GetSize());
          std::fill(rhsValues2_.Begin(), rhsValues2_.End(), 0.0);
        }
        
        for(UInt i=0; i<numSourceRegions; i++)
        {
          acouRHSVal = resultHandler->GetResult( srcInputId_,
                                                 1,
                                                 step_,
                                                 ACOU_RHS_LOAD,
                                                 srcRegions_[i] );

          if(!destElemList_)
          {
            destElemList_ = new ElemList(dest);
            destElemList_->SetRegion(dest->RegionNameToId(regionName_));
          }
          
          if(!sourceNodeLists_[i])
          {
            source = acouRHSVal->GetEntityList()->GetGrid();
            sourceNodeLists_[i] = new NodeList(source);
            RegionIdType srcRegionId = source->RegionNameToId(srcRegions_[i]);
            sourceNodeLists_[i]->SetNodesOfRegion(srcRegionId);
          }

          // read data from archive
          if(restartFileMode_ == "r" || restartFileMode_ == "rw") {
            // open a character archive for input
            std::ostringstream fNameStr;
            fNameStr << "cons_interpol_weights_" << srcRegions_[i] << ".dat"; 
            std::ifstream ifs(fNameStr.str().c_str(), std::ios::binary);
            if ( ifs.good() ) {
              boost::archive::binary_iarchive ia(ifs);
            
              // read conservative interpolation weights from archive
              ia >> consInterpWeights_[i];
              // archive and stream closed when destructors are called
            }
            else {
              (*warning) << "An error occured while reading the restart file "
                         << "for conservative interpolation weights. All "
                         << "weights will be recalculated.";
              Warning(__FILE__, __LINE__);
            }
          }

          dest->ComputeConservativeInterpolationWeights(*destElemList_,
              *sourceNodeLists_[i],
              coordSysId_,
              globalEpsilon_,
              localEpsilon_,
              consInterpWeights_[i],
              node_warnings_);

          // save data to archive
          if(restartFileMode_ == "w" || restartFileMode_ == "rw") {
            // create and open a character archive for output
            std::ostringstream fNameStr;
            fNameStr << "cons_interpol_weights_" << srcRegions_[i] << ".dat"; 
            std::ofstream ofs(fNameStr.str().c_str(), std::ios::binary);
            boost::archive::binary_oarchive oa(ofs);
            
            // write class instance to archive
            oa << ((const std::vector< std::map<UInt, Double> >&) consInterpWeights_[i]);
            // archive and stream closed when destructors are called
          }
          
          // Get data of previous time step for asynchronous time stepping
          // with interpolation
          if (timeInterp) {
            // do we need to load the data?
            if (prevStep != lastStep_) {
              // only if previous step gives a significant contribution
              if (fabs(1.0-intFactor) > timeTol) {
                // load data from file
                acouRHSVal2 = resultHandler->GetResult(srcInputId_, 1,
                    prevStep, ACOU_RHS_LOAD, srcRegions_[i]);
                
                Result<Double> *result2 =
                  dynamic_cast<Result<Double>*>(&(*acouRHSVal2));
                if (result2 == NULL) {
                  EXCEPTION("Cannot read result 'acouRhsLoad' from input id '"
                      << id_ << "'");
                }
                Vector<Double>& resVec2 = result2->GetVector();
                
                std::map<UInt, Double>::const_iterator it, end;
                
                for (UInt j=0; j < consInterpWeights_[i].size(); ++j) {
                  it = consInterpWeights_[i][j].begin();
                  end = consInterpWeights_[i][j].end();
  
                  if(!std::distance(it, end))
                    continue;
                  
                  for (; it != end; ++it) {
                    rhsValues2_[it->first] += it->second * resVec2[j];
                  }
                }
              } // fabs(1.0-intFactor) > timeTol
            } // prevStep != lastStep_
          } // timeInterp

          std::map<UInt, Double>::const_iterator it, end;
          Result<Double> *result =
            dynamic_cast<Result<Double>*>(&(*acouRHSVal)); 
          if (result == NULL) {
              EXCEPTION("Cannot read result 'acouRhsLoad' from input id '"
                  << id_ << "'");
            }
          Vector<Double>& resVec = result->GetVector();

          //Integer pos;
          for(UInt j=0; j<consInterpWeights_[i].size(); j++)
          {
            it = consInterpWeights_[i][j].begin();
            end = consInterpWeights_[i][j].end();

            if(!std::distance(it, end))
              continue;

            for(; it != end; ++it)
            {
              /*
              std::cout << "Node: " << (i+1) << " -> " << it->first
              << ": " << it->second << std::endl;
               */

              //            pos = regionNodes.Find(it->first);
              rhsValues_[it->first] += it->second * resVec[j];
            }
          }
        }
#else // USE_INTERPOLATION
        EXCEPTION("This executable does not support interpolation!");
#endif // USE_INTERPOLATION
      }
      else // if (interpolate_)
      {
        // Get data of previous time step for asynchronous time stepping
        // with interpolation
        if (timeInterp) {
          // only if previous step gives a significant contribution
          if (fabs(1.0-intFactor) > timeTol) {
            // are the data still in memory?
            if (prevStep == lastStep_) {
              rhsValues2_.Resize(rhsValues_.GetSize());
              rhsValues2_ = rhsValues_; // just copy it
            }
            else { // load data from file
              acouRHSVal2 = resultHandler->GetResult(id_, 1, prevStep,
                  ACOU_RHS_LOAD, regionName_);
              
              Result<Double> *result2 =
                dynamic_cast<Result<Double>*>(&(*acouRHSVal2));
              if (result2 == NULL) {
                EXCEPTION("Cannot read result 'acouRhsLoad' from input id '"
                    << id_ << "'");
              }
              Vector<Double>& resVec2 = result2->GetVector();
              UInt numValues = resVec2.GetSize();
              
              rhsValues2_.Resize(numValues);
              for (UInt i=0; i < numValues; ++i) {
                rhsValues2_[i] = resVec2[i];
              }
            }
          }
        } // if (timeInterp)

        // Get reference result
        acouRHSVal = resultHandler->GetResult( id_,
                                               1,
                                               step_,
                                               ACOU_RHS_LOAD,
                                               regionName_ );

        Result<Double> *result =
          dynamic_cast<Result<Double>*>(&(*acouRHSVal));
        if (result == NULL) {
          EXCEPTION("Cannot read result 'acouRhsLoad' from input id '"
              << id_ << "'");
        }
        Vector<Double>& resVec = result->GetVector();

        rhsValues_.Resize(resVec.GetSize());
        for(UInt i=0, n=resVec.GetSize();
            i < n;
            i++)
          rhsValues_[i] = resVec[i];
      } // if (interpolate_)

      lastStep_ = step_;
    }

    elemVec.Resize(1);
    // do we use asynchronous time stepping with interpolation?
    if (timeInterp && (fabs(1.0-intFactor) > timeTol)) {
      UInt iEnt = ent.GetPos();
      // linear interpolation
      elemVec[0] = factor*(intFactor*rhsValues_[iEnt]
                   + (1.0-intFactor)*rhsValues2_[iEnt]);
    }
    else {
      elemVec[0] = rhsValues_[ent.GetPos()]*factor;
    }
    
  } // if (step_ != lastStep_)

  void AcouRHSLinForm::CalcElemVector( Vector<Complex> & elemVec,
                                       EntityIterator& ent ) 
  {
    ResultHandler* resultHandler = domain->GetResultHandler();
    shared_ptr<BaseResult> acouRHSVal;
    Double factor;
    Double f = 0.0;

    step_ = (UInt) mParser_->Eval(mHandle_);
    factor = mParser_->Eval( mHandle2_ );
    
    if (asynchSteps_ != "no") {
      // get current frequency
      f = mParser_->Eval(mph_tf_);
      
      // we have to use an iterator in order to tackle non-contiguous step
      // numbering (saveInc>1)
      std::map< UInt, Double >::iterator it = stepValues_.lower_bound(lastStep_);
      step_ = it->first;
      
      // find first step with a larger time value than the current time
      while (it != stepValues_.end()) {
        if (it->second >= f)
          break;
        step_ = (++it)->first;
      }
   
      if (it == stepValues_.end()) {
        step_ = (--it)->first;
      }

      // check if the previous step is closer to the current frequency
      // (nearest neighbor interpolation)
      Double df = fabs(it->second - f); 
      if (df > fabs(f - (--it)->second))
        step_ = it->first;
    }
    
    if(step_ != lastStep_)
    {
      acouRHSVal = resultHandler->GetResult( id_,
                                             1,
                                             step_,
                                             ACOU_RHS_LOAD,
                                             regionName_ );
      Result<Complex> *result = dynamic_cast<Result<Complex>*>(&(*acouRHSVal));
      if (result == NULL) {
        EXCEPTION("Cannot read result 'acouRhsLoad' from input id '"
                  << id_ << "'");
      }
      Vector<Complex> &resVec = result->GetVector();
      
      rhsValuesComplex_.Resize(resVec.GetSize());
      for(UInt i=0, n=resVec.GetSize();
          i < n;
          i++)
        rhsValuesComplex_[i] = resVec[i];
      
      lastStep_ = step_;
    }
    
    Complex complexValue( rhsValuesComplex_[ent.GetPos()].real() * factor,
                          rhsValuesComplex_[ent.GetPos()].imag() * factor );

    elemVec.Resize(1);
    elemVec[0] = complexValue;

  }

  /*
    valAmpl = rhsValues_[ent.GetNode()-1] * factor;
    valPhase = rhsValues2_[ent.GetNode()-1];
    Complex complexValue( valAmpl * cos( valPhase / 180 * PI ),
                          valAmpl * sin( valPhase / 180 * PI ) );

    elemVec.Resize(1);
    elemVec[0] = complexValue;
  }
  */

} // end of namespace
