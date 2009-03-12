// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include "Domain/domain.hh"
#include "Domain/grid.hh"

#include "dampLayer.hh"

namespace CoupledField
{

  DampLayer::DampLayer(std::string type) {

    //check correct damping type
    if ( type == "linear" ) {
      dampFncType_ = type;
    }
    else if ( type == "quadDist" ) {
      dampFncType_ = type;
    }
    else if ( type == "inverseDist" ) {
      dampFncType_ = type;
    }
    else if ( type == "exponential" ) {
      dampFncType_ = type;
    }
    else {
      EXCEPTION("Damping type for Sponge Layer not known");
    }
  }


 
  DampLayer::~DampLayer()
  {
  }


  void DampLayer::CalcDampFactor( Double& factor, EntityIterator& ent)
  {

    domain->GetGrid()->GetElemNodesCoord( ptCoord_, 
                                          ent.GetElem()->connect,
                                          false );

    UInt numVals = ptCoord_.GetNumRows();
    Vector<Double> pos;
    pos.Resize(numVals);
    pos.Init(0);
    for (UInt i=0; i<ptCoord_.GetNumRows(); i++) {
      for (UInt j=0; j<ptCoord_.GetNumCols(); j++) {
	pos[i] += ptCoord_[i][j];
      }
      pos[i] /= (Double) ptCoord_.GetNumCols();
    }

    //    std::cout << "pos:\n" << pos << std::endl;

    factor = EvalDampFnc(pos);
    //    std::cout << "DampFactor: " << factor << "\n" << std::endl;
  }


  Double DampLayer::EvalDampFnc(Vector<Double>& pos)
  {
		// initialize or receive compiler warning
    Double factor = 0.0;

    if ( dampFncType_ == "linear" ) {
      factor = dampingFactor_;
    }

    else if ( dampFncType_ == "quadDist" ) {
      factor = dampingFactor_;
    }

    else if ( dampFncType_  == "inverseDist" ) {
      EXCEPTION("PML damping inverseDist divides by factor smaller 1E-12");
    }
    else if ( dampFncType_ == "exponential" ) {
      UInt dim = pos.GetSize(); 
      Vector<Double> diff(dim);

      diff =  pos - midPoint_;
      Double dist = diff.NormL2() - startRadius_;

      //      std::cout << "dist= " << dist << std::endl;

      if ( dist > 0 ) {
	factor = std::exp( dampingFactor_*dist / layerThickness_ );
	if ( factor > dampingFactorMax_ ) 
	  factor = dampingFactorMax_;
      }
      else {
	factor = 1.0;
      }
    }

    return factor;
  }


  void DampLayer::SetDampingParams(Vector<Double> point, 
				   Double& dampFactor, 
				   Double& dampFactorMax, 
				   Double& startRadius, 
				   Double& endRadius) {

    // center point
    midPoint_ = point;

    startRadius_ = startRadius;
    endRadius_   = endRadius;

    dampingFactor_   = dampFactor;
    dampingFactorMax_ = dampFactorMax;
    layerThickness_ =  endRadius - startRadius;
    if ( layerThickness_ < 0 ) {
      EXCEPTION("Start-Radius has to be larger then End-Radius");
    }
  }


} // end namespace CoupledField
