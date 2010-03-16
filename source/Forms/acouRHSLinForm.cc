// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>


// #include "DataInOut/infiles.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
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

  AcouRHSLinForm::AcouRHSLinForm(PtrParamNode rhsValuesNode)
    : LinearForm(),
      id_("default"),
      interpolate_(false),
      restartFileMode_("rw"),
      overwriteOldSrcs_(true),
      coordSysId_("default"),
      globalEpsilon_(1.0e-4, 0.0, 0.0),
      localEpsilon_(1.0e-3, 0.0, 0.0),
      z_(0.0), zEpsilon_(1.0e-10),
      asyncSteps_("no"),
      node_warnings_(CI_WARN_YES),
      ptGrid_(domain->GetGrid()),
      globCoSy_(domain->GetCoordSystem()),
      step_(0), lastStep_(0),
      destNodeList_(NULL),
      destElemList_(NULL),
      isFirstStep_(true)
  {
    name_ = "AcouRHSLinForm";

    std::string srcRegions, factor = "1.0";

    // Read parameters from XML file
    PtrParamNode tmpNode;

    try {
      // destination region
      rhsValuesNode->GetValue("region", regionName_);
      // input ID of destination region
      rhsValuesNode->GetValue("inputId", id_);
      // factor (is multiplied with RHS)
      rhsValuesNode->GetValue("factor", factor, ParamNode::PASS);
      
      // interpolation parameters (if any)
      PtrParamNode intNode;
      if(rhsValuesNode->Has("interpolation"))
        intNode = rhsValuesNode->Get("interpolation", ParamNode::PASS);
      
      if (intNode) {
        interpolate_ = true;
        
        tmpNode = intNode->Get("srcRegions", ParamNode::PASS);
        if (tmpNode->HasChildren()) {
          // names of source regions
          tmpNode->GetValue("names", srcRegions);
          // input ID of source regions
          tmpNode->GetValue("inputId", srcInputId_);
          // coordinate system ID of source regions
          tmpNode->GetValue("coordSysId", coordSysId_);
        }
        
        // if 3D data are to be interpolated to a 2D grid,
        // where should the xy-plane be?
        tmpNode = intNode->Get("xyPlane", ParamNode::PASS);
        if (tmpNode->HasChildren()) {
          tmpNode->GetValue("z", z_, ParamNode::PASS);
          tmpNode->GetValue("tol", zEpsilon_, ParamNode::PASS);
        }
        
        tmpNode = intNode->Get("tolerances", ParamNode::PASS);
        if (tmpNode->HasChildren()) {
          // tolerance in global coordinates
          PtrParamNode tolNode = tmpNode->Get("global", ParamNode::PASS);
          if (tolNode) {
            tolNode->GetValue("start", globalEpsilon_.start, ParamNode::PASS);
            tolNode->GetValue("end", globalEpsilon_.end, ParamNode::PASS);
            tolNode->GetValue("inc", globalEpsilon_.inc, ParamNode::PASS);
          }
          
          // tolerance in local coordinates
          tolNode = tmpNode->Get("local", ParamNode::PASS);
          if (tolNode) {
            tolNode->GetValue("start", localEpsilon_.start, ParamNode::PASS);
            tolNode->GetValue("end", localEpsilon_.end, ParamNode::PASS);
            tolNode->GetValue("inc", localEpsilon_.inc, ParamNode::PASS);
          }
        }
        
        // restart file mode
        intNode->GetValue("restartFileMode", restartFileMode_, ParamNode::PASS);
        
        // verbosity of warnings
        tmpNode = intNode->Get("nodeWARNs", ParamNode::PASS);
        if (tmpNode->HasChildren()) {
          std::string dispStr = tmpNode->Get("display")->As<std::string>();
          if (dispStr == "verbose")
            node_warnings_ = (ciWarnFlags) (CI_WARN_YES | CI_WARN_VERBOSE);
          else if (dispStr == "yes")
            node_warnings_ = CI_WARN_YES;
          else
            node_warnings_ = CI_WARN_NO;
          
          if (tmpNode->Get("writeNodes")->As<bool>()) {
            node_warnings_ = (ciWarnFlags) (node_warnings_ | CI_WARN_LIST);
          }
        }

        // overwrite old source terms present on destination mesh?
        intNode->GetValue("overwriteOldSrcs", overwriteOldSrcs_, ParamNode::PASS);
      }
      
      // do we use asynchronous time/frequency steps?
      rhsValuesNode->GetValue("asyncSteps", asyncSteps_, ParamNode::PASS);
      
    } catch (Exception& ex) 
    {
      RETHROW_EXCEPTION(ex, "Error while trying to read parameters for AcouRHSLinForm.");
    }
    
    // Determine analysis type
    analysistype_ = domain->GetSingleDriver()->GetAnalysisType();

    // For asynchronous steps we need the time values
    if (asyncSteps_ != "no") {
      mph_tf_ = mParser_->GetNewHandle();
      if (analysistype_ == BasePDE::TRANSIENT) {
        // Get a handle to math parser's time variable
        mParser_->SetExpr(mph_tf_, "t");
      }
      else {
        // Get a handle to math parser's frequency variable
        mParser_->SetExpr(mph_tf_, "f");

        // linear interpolation does not make sense here
        if (asyncSteps_ == "interpolate") {
          WARN("For a harmonic analysis it makes no sense to use "
               << "linear frequency step interpolation; switching to nearest "
               << "neighbor.");
          asyncSteps_ = "nearest";
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
    mphFactor_ =  mParser_->GetNewHandle();
    mParser_->SetExpr( mphFactor_, factor );

    // split the string given in srcRegions into region names
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
  }

  AcouRHSLinForm::~AcouRHSLinForm()
  {
    mParser_->ReleaseHandle( mphFactor_ );
    if (asyncSteps_ != "no")
      mParser_->ReleaseHandle(mph_tf_);
  }

  void AcouRHSLinForm::CalcElemVector( Vector<Double> & elemVec,
                                       EntityIterator& ent )
  {
    ResultHandler* resultHandler = domain->GetResultHandler();
#ifdef USE_INTERPOLATION
    Grid* source = NULL;
#endif
    shared_ptr<BaseResult> acouRHSVal, acouRHSVal2;
    bool timeInterp = (asyncSteps_ == "interpolate");
    Double factor;
    Double intFactor = 1.0;
    Double t = 0.0;
    const Double timeTol = 0.001;
    UInt prevStep = 0;
    Vector<Double> nodeCoord;
    
    step_ = (UInt) mParser_->Eval(mHandle_);

    ptGrid_->GetNodeCoordinate(nodeCoord, ent.GetNode());
    mParser_->SetCoordinates(mphFactor_, *globCoSy_, nodeCoord);
    factor = mParser_->Eval( mphFactor_ );
    
    // Determine asynchronous time step
    if (asyncSteps_ != "no") {
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

      if (asyncSteps_ == "nearest") { // nearest neighbor interpolation
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
              WARN("Current time/frequency step is beyond the last "
                   << "step of input file (id='" << srcInputId_ << "')");
            }
          }
        }
      }
    }


    // Determine if values for a new step have to be read.

    if (step_ != lastStep_)
    {
      // Load data before doing interpolation
      StdVector<UInt> regionNodes;

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
              Result<Double> *result2 = NULL;
              try {
                // Try to read in RHS values from main grid file.
                acouRHSVal2 = resultHandler->GetResult(id_, 1, prevStep,
                                                       ACOU_RHS_LOAD, regionName_);

                result2 = dynamic_cast<Result<Double>*>(&(*acouRHSVal2));
              } catch (Exception& ex) {};
              

              if (result2 == NULL) {
                if(!interpolate_)
                {
                  EXCEPTION("Cannot read result 'acouRhsLoad' from input id '"
                            << id_ << "'");
                }
                else
                {
                  ptGrid_->GetNodesByRegion( regionNodes, ptGrid_->GetRegion().Parse(regionName_));
                  rhsValues2_.Resize(regionNodes.GetSize());
                  std::fill(rhsValues2_.Begin(), rhsValues2_.End(), 0.0);
                }
              }
              else
              {
                Vector<Double>& resVec2 = result2->GetVector();
                UInt numValues = resVec2.GetSize();

                rhsValues2_.Resize(numValues);
                for (UInt i=0; i < numValues; ++i) {
                  rhsValues2_[i] = resVec2[i];
                }
              }
            }
          }
        } // if (timeInterp)

        // Get reference result
        Result<Double> *result = NULL;
        try {
          // Try to read in RHS values from main grid file.
          acouRHSVal = resultHandler->GetResult( id_,
                                                 1,
                                                 step_,
                                                 ACOU_RHS_LOAD,
                                                 regionName_ );
          
          result = dynamic_cast<Result<Double>*>(&(*acouRHSVal));
        } catch (Exception& ex) {};
        
        if (result == NULL) {
          if(!interpolate_)
          {
//            EXCEPTION("Cannot read result 'acouRhsLoad' from input id '"
//                      << id_ << "'");
              std::fill(rhsValues_.Begin(), rhsValues_.End(), 0.0);
          }
          else
          {
            ptGrid_->GetNodesByRegion( regionNodes, ptGrid_->GetRegion().Parse(regionName_));
            rhsValues_.Resize(regionNodes.GetSize());
            std::fill(rhsValues_.Begin(), rhsValues_.End(), 0.0);
          }
          
        }
        else 
        {
          Vector<Double>& resVec = result->GetVector();
          
          rhsValues_.Resize(resVec.GetSize());
          for(UInt i=0, n=resVec.GetSize();
              i < n;
              i++)
            rhsValues_[i] = resVec[i];
        }
        
        




      if (interpolate_)
      {
#ifdef USE_INTERPOLATION
        UInt numSourceRegions = srcRegions_.GetSize();
        //        StdVector<UInt> regionNodes;
        StdVector<UInt> unmapped_nodes;

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
        //        ptGrid_->GetNodesByRegion( regionNodes, ptGrid_->GetRegion().Parse(regionName_));
        //        rhsValues_.Resize(regionNodes.GetSize());
        if(overwriteOldSrcs_)
          std::fill(rhsValues_.Begin(), rhsValues_.End(), 0.0);


        if (timeInterp) {
          //          rhsValues2_.Resize(regionNodes.GetSize());
          if(overwriteOldSrcs_)
            std::fill(rhsValues2_.Begin(), rhsValues2_.End(), 0.0);
        }

        if (isFirstStep_) {
          std::cout << "Calculating conservative interpolation weights..."
                    << std::endl;
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
            destElemList_ = new ElemList(ptGrid_);
            destElemList_->SetRegion(ptGrid_->GetRegion().Parse(regionName_));
          }

          if(!sourceNodeLists_[i])
          {
            source = acouRHSVal->GetEntityList()->GetGrid();
            sourceNodeLists_[i] = new NodeList(source);
            RegionIdType srcRegionId = source->GetRegion().Parse(srcRegions_[i]);
            sourceNodeLists_[i]->SetNodesOfRegion(srcRegionId);
          }
          
          // compute interpolation weights only in first step
          if (isFirstStep_) {
            // print name of current region
            std::cout << srcRegions_[i] << " ";
            
            std::ostringstream fNameStr;
            fNameStr << progOpts->GetSimName() << "_ms"
                     << domain->GetSingleDriver()->GetActSequenceStep()
                     << "_" << srcRegions_[i] << "_ciw.dat";
            // read data from archive
            if(restartFileMode_ == "r" || restartFileMode_ == "rw") {
              // open a character archive for input
              std::ifstream ifs(fNameStr.str().c_str(), std::ios::binary);
              if ( ifs.good() ) {
                try {
                  boost::archive::binary_iarchive ia(ifs);
                  // read conservative interpolation weights from archive
                  ia >> consInterpWeights_[i];
                  // archive and stream closed when destructors are called
                } catch (std::exception &ex) {
                  EXCEPTION("The following problem occurred while trying to "
                            << "read conservative interpolation weights from '"
                            << fNameStr << "': " << ex.what()
                            << "\nTry to set restartMode to 'w' in XML file.")
                }
              }
              else {
                WARN("An error occured while reading the restart file "
                     << "for conservative interpolation weights. All "
                     << "weights will be recalculated.");
              }
            }
  
            ptGrid_->ComputeConservativeInterpolationWeights(*destElemList_,
                *sourceNodeLists_[i],
                coordSysId_,
                globalEpsilon_,
                localEpsilon_,
                z_, zEpsilon_,
                consInterpWeights_[i],
                unmapped_nodes);
  
            // save data to archive
            if(restartFileMode_ == "w" || restartFileMode_ == "rw") {
              std::ofstream ofs(fNameStr.str().c_str(), std::ios::binary);
              if ( ofs.good() ) {
                try {
                  // create and open a character archive for output
                  boost::archive::binary_oarchive oa(ofs);
              
                  // write class instance to archive
                  oa << ((const std::vector< std::map<UInt, Double> >&) consInterpWeights_[i]);
                  // archive and stream closed when destructors are called
                } catch (std::exception &ex) {
                  EXCEPTION("The following problem occurred while trying to "
                            << "write conservative interpolation weights to '"
                            << fNameStr << "': " << ex.what())
                }
              }
              else {
                WARN("An error occured while writing the restart file "
                     << "for conservative interpolation weights.");
              }
            }
            
            // print a newline for a proper status display
            std::cout << std::endl;
          } // if (isFirstStep_)

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

          Double sum_orig = 0;            
          Double sum = 0;            
          UInt n = consInterpWeights_[i].size();
          //UInt modval = n / 100;
          
          //Integer pos;
          for(UInt j=0; j<n; j++)
          {
            it = consInterpWeights_[i][j].begin();
            end = consInterpWeights_[i][j].end();

            Double sum_weights = 0;            

            for(; it != end; ++it)
            {
              /*
              std::cout << "Node: " << (i+1) << " -> " << it->first
              << ": " << it->second << std::endl;
               */

              //            pos = regionNodes.Find(it->first);
              rhsValues_[it->first] += it->second * resVec[j];

              if (node_warnings_ & CI_WARN_YES) {
                sum += it->second * resVec[j];
                sum_weights += it->second;
              }
            }
            if (node_warnings_ & CI_WARN_YES) {
              sum_orig += resVec[j];

              // Check for conservitveness every modval source nodes or at least once in the end
              //            if(j % modval == 0 || j == n - 1)
              if(j == n - 1)
              {
                Double ratio = (sum - sum_orig) / sum_orig;
                if( abs(ratio) > 0.01 ) 
                {
                  // If this condition occurs it means that some source nodes
                  // do not have coservative interpolation weights in consInterpWeights_[i]
                  // so that the for loop above does not do anything.
                  WARN("Sum of interpolated source terms (" << sum << ") off by "
                       << (ratio * 100) << "% of sum of original source terms ("
                       << sum_orig << ") from node index " << j
                       << " in source region '" << srcRegions_[i] << "'.");
                }
              }
            }
          }
          
        }

        // interpolation weights are calculated in the first time step only,
        // so we need to check for unmapped nodes only in this case.
        if (isFirstStep_) {
          UInt numBadNodes = unmapped_nodes.GetSize();
          if (numBadNodes > 0)  {
            std::ostringstream fNameStr;
            std::ofstream badNodesFile;

            if (node_warnings_ & CI_WARN_YES) {
              WARN("During conservative interpolation " << numBadNodes
                   << " nodes in total were not mapped.");
            }

            if (node_warnings_ & CI_WARN_LIST) {
              fNameStr << progOpts->GetSimName() << "_ci_unmapped_nodes.txt";
              badNodesFile.open(fNameStr.str().c_str(),
                  std::ios_base::out | std::ios_base::trunc);
              if (!badNodesFile)
                EXCEPTION("Error writing to file " << fNameStr.str());
              badNodesFile.width(1);
              badNodesFile << std::ios_base::left;
            }

            for(UInt i=0; i<numBadNodes; ++i) {
              if (node_warnings_ & CI_WARN_VERBOSE) {
                WARN("Node " << unmapped_nodes[i] << " from source "
                     << "grid could not be mapped to destination grid.");
              }
              if ((node_warnings_ & CI_WARN_LIST) && badNodesFile) {
                badNodesFile << unmapped_nodes[i] << std::endl;
                if (!badNodesFile)
                  EXCEPTION("Error writing to file " << fNameStr.str());
              }
            }

            if (node_warnings_ & CI_WARN_LIST) {
              badNodesFile.close();
            }
          }
        }

        isFirstStep_ = false;

#else // USE_INTERPOLATION
        
        EXCEPTION("This executable does not support interpolation!");
        
#endif // USE_INTERPOLATION
        
      }
      //      else // if (interpolate_)
      //      {
      //      } // if (interpolate_)

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
    Vector<Double> nodeCoord;

    step_ = (UInt) mParser_->Eval(mHandle_);

    ptGrid_->GetNodeCoordinate(nodeCoord, ent.GetNode());
    mParser_->SetCoordinates(mphFactor_, *globCoSy_, nodeCoord);
    factor = mParser_->Eval( mphFactor_ );
    
    if (asyncSteps_ != "no") {
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
      for(UInt i=0, n=resVec.GetSize(); i < n; i++)
        rhsValuesComplex_[i] = resVec[i];

      lastStep_ = step_;
    }

    Complex complexValue( rhsValuesComplex_[ent.GetPos()].real() * factor,
                          rhsValuesComplex_[ent.GetPos()].imag() * factor );

    elemVec.Resize(1);
    elemVec[0] = complexValue;

  }

} // end of namespace
