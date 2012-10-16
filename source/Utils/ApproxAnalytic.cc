// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <vector>

#include "Domain/domain.hh"
#include "ApproxAnalytic.hh"

namespace CoupledField
{ 

  ApproxAnalytic::ApproxAnalytic(std::string fncStr, MaterialType matType )
    : ApproxData( "", matType )
  {

  }

  //! set the analytic expressions
  void ApproxAnalytic::SetAnalyticExpressions( std::string fncStr, 
                                               std::string fncDerivStr ) {

    mParser_ = domain->GetMathParser();

    // Specify the function 
    mphFactor_ =  mParser_->GetNewHandle();
    mParser_->SetExpr( mphFactor_,  fncStr );

    mphFactorDeriv_ =  mParser_->GetNewHandle();
    mParser_->SetExpr( mphFactorDeriv_, fncDerivStr );


  }

  //! computes the magnetic reluctivity
  Double ApproxAnalytic::EvaluateFuncNu( Double t ) {     
    //limit the maximal value
    if (t > yMax_ ) {
      t = yMax_;
    }

    //set the value
    mParser_->SetValue( mphFactor_, "b", t );
    return mParser_->Eval( mphFactor_ );

  }


  //! computes the magnetic reluctivity
  Double ApproxAnalytic::EvaluatePrimeNu( Double t ) {     
    //limit the maximal value
    if (t > yMax_ ) {
      t = yMax_;
    }

    //set the value
    mParser_->SetValue( mphFactorDeriv_, "b", t );
    return mParser_->Eval( mphFactorDeriv_ );
  }

}

