// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "singleEntryInt.hh"
#include "Domain/domain.hh"


namespace CoupledField {


 SingleEntryInt::SingleEntryInt(  const std::string& real,
                                  UInt dof, UInt numDofs )
    : BaseForm( NULL, FULL, false ) {

    name_ = "SingleEntryInt";
    imag = "0.0";
    dof_ = dof;
    numDofs_ = numDofs;
    isComplex_ = false;

    // obtain handles
    rHandle_ =  domain->GetMathParser()->GetNewHandle();
    iHandle_ =  domain->GetMathParser()->GetNewHandle();
    mParser_->SetExpr( rHandle_, real );
    mParser_->SetExpr( iHandle_, imag );

  }


  SingleEntryInt::SingleEntryInt(  const std::string& real,
                                   const std::string& imag,
                                  UInt dof, UInt numDofs )
    : BaseForm( NULL, FULL, false ) {

    name_ = "SingleEntryInt";
    dof_ = dof;
    numDofs_ = numDofs;
    isComplex_ = true;

    // obtain handles
    rHandle_ =  domain->GetMathParser()->GetNewHandle();
    iHandle_ =  domain->GetMathParser()->GetNewHandle();
    mParser_->SetExpr( rHandle_, real );
    mParser_->SetExpr( iHandle_, imag );

  }
  
  SingleEntryInt::~SingleEntryInt() {
  }
  
  void SingleEntryInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                          EntityIterator& ent1, 
                                          EntityIterator& ent2 ) {

    // Resize and init elemMat
    elemMat.Resize( numDofs_ );
    elemMat.Init();
    elemMat[dof_-1][dof_-1] = mParser_->Eval( rHandle_ );
    
  }

  void SingleEntryInt::CalcElementMatrix( Matrix<Complex>& elemMat,
                                          EntityIterator& ent1, 
                                          EntityIterator& ent2 ) {

    // Resize and init elemMat
    elemMat.Resize( numDofs_ );
    elemMat.Init();
    
    Complex entry = Complex(mParser_->Eval(rHandle_ ),
                            mParser_->Eval( iHandle_ ) );
    elemMat[dof_-1][dof_-1] = entry;
  }
  
}
