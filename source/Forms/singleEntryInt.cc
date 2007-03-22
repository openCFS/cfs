// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "singleEntryInt.hh"


namespace CoupledField {



  SingleEntryInt::SingleEntryInt( Double value, UInt dof, UInt numDofs )
    : BaseForm( NULL, FULL, false ) {
    ENTER_FCN( "SingleEntryInt::SingleEntryInt" );

    name_ = "SingleEntryInt";
    entry_ = value;
    dof_ = dof;
    numDofs_ = numDofs;
  }
  
  SingleEntryInt::~SingleEntryInt() {
    ENTER_FCN( "SingleEntryInt::~SingleEntryInt" );
  }
  
  void SingleEntryInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                          EntityIterator& ent1, 
                                          EntityIterator& ent2 ) {
    ENTER_FCN( "SingleEntryInt::CalcElementMatrix" );

    // Resize and init elemMat
    elemMat.Resize( numDofs_ );
    elemMat.Init();
    elemMat[dof_-1][dof_-1] = entry_;
    
  }

}
