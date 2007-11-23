// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <math.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>


// #include "DataInOut/infiles.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/entityList.hh"

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
                                 std::string restartFileMode)
    : LinearForm(),
      id_(id),
      regionName_(regionName),
      interpolate_(interpolate),
      srcInputId_(srcInputId),
      coordSysId_(coordSysId),
      globalEpsilon_(globalEpsilon),
      localEpsilon_(localEpsilon),
      restartFileMode_(restartFileMode),
      step_(0),
      lastStep_(0),
      destNodeList_(NULL),
      destElemList_(NULL)
  {
    name_ = "AcouRHSLinForm";

    // Set Expressions for math parser

    // Specify that we want to use the current step number.
    mParser_->SetExpr( mHandle_, "step" );

    // Specify the ramping factor function for the RHS to the parser.
    mHandle2_ =  mParser_->GetNewHandle();
    mParser_->SetExpr( mHandle2_, factor );
    std::cout << "********** factor is: " << factor << std::endl;

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
  }

  void AcouRHSLinForm::CalcElemVector( Vector<Double> & elemVec,
                                       EntityIterator& ent ) 
  {
    ResultHandler* resultHandler = domain->GetResultHandler();
    Grid* source;
    Grid* dest;
    shared_ptr<BaseResult> acouRHSVal;
    Double factor;

    step_ = mParser_->Eval(mHandle_);
    factor = mParser_->Eval( mHandle2_ );
    
    // Determine if values for a new step have to be read.
    if(step_ != lastStep_)
    {
      if(interpolate_)
      {
#ifdef USE_INTERPOLATION
        UInt numSourceRegions = srcRegions_.GetSize();
        StdVector<UInt> regionNodes;
        dest = domain->GetGrid();
        dest->GetNodesByRegion( regionNodes, dest->RegionNameToId(regionName_));
        rhsValues_.resize(regionNodes.GetSize());
        std::fill(rhsValues_.begin(), rhsValues_.end(), 0);

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
            boost::archive::binary_iarchive ia(ifs);
            
            // read conservative interpolation weights from archive
            ia >> consInterpWeights_[i];
            // archive and stream closed when destructors are called
          }

          dest->ComputeConservativeInterpolationWeights(*destElemList_,
              *sourceNodeLists_[i],
              coordSysId_,
              globalEpsilon_,
              localEpsilon_,
              consInterpWeights_[i]);

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
          
          std::map<UInt, double>::const_iterator it, end;
          Vector<double>* resVec = dynamic_cast<Vector<Double>*>(acouRHSVal->GetCFSVector());

          Integer pos;
          for(UInt j=0; j<consInterpWeights_[i].size(); j++)
          {
            it = consInterpWeights_[i][j].begin();
            end = consInterpWeights_[i][j].end();

            if(!std::distance(it, end))
              continue;

            for(; it != end; it++)
            {
              /*
              std::cout << "Node: " << (i+1) << " -> " << it->first
              << ": " << it->second << std::endl;
               */

              //            pos = regionNodes.Find(it->first);
              rhsValues_[it->first] += it->second * (*resVec)[j];
            }
          }
        }
#else // USE_INTERPOLATION
        EXCEPTION("This executable does not support interpolation!");
#endif // USE_INTERPOLATION
      }
      else
      {
        // Get reference result
        acouRHSVal = resultHandler->GetResult( id_,
                                               1,
                                               step_,
                                               ACOU_RHS_LOAD,
                                               regionName_ );
        Vector<double>& resVec = dynamic_cast<Result<Double>&>(*acouRHSVal).GetVector();
//        Vector<double>* resVec = dynamic_cast<Vector<Double>*>(acouRHSVal->GetCFSVector());

        rhsValues_.resize(resVec.GetSize());
        for(UInt i=0, n=resVec.GetSize();
            i < n;
            i++)
          rhsValues_[i] = resVec[i];
      }

      lastStep_ = step_;
    }

    elemVec.Resize(1);
    elemVec[0] = rhsValues_[ent.GetNode()-1]*factor;

  }

  void AcouRHSLinForm::CalcElemVector( Vector<Complex> & elemVec,
                                       EntityIterator& ent ) 
  {
    ResultHandler* resultHandler = domain->GetResultHandler();
    shared_ptr<BaseResult> acouRHSVal;
    Double factor;

    step_ = (UInt) mParser_->Eval(mHandle_);
    factor = mParser_->Eval( mHandle2_ );
    
    
    if(step_ != lastStep_)
    {
      acouRHSVal = resultHandler->GetResult( id_,
                                             1,
                                             step_,
                                             ACOU_RHS_LOAD,
                                             regionName_ );
      
      Vector<Complex> &resVec = dynamic_cast<Result<Complex>&>(*acouRHSVal).GetVector();
      
      rhsValuesComplex_.resize(resVec.GetSize());
      for(UInt i=0, n=resVec.GetSize();
          i < n;
          i++)
        rhsValuesComplex_[i] = resVec[i];
      
      lastStep_ = step_;
    }
    
    Complex complexValue( rhsValuesComplex_[ent.GetNode()-1].real() * factor,
                          rhsValuesComplex_[ent.GetNode()-1].imag() * factor );

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
