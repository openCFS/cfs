// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "basedriver.hh"
#include "General/exception.hh"

namespace CoupledField
{

  BaseDriver :: BaseDriver( )
  {
    ENTER_FCN( "BaseDriver::BaseDriver" );

    actSequenceStep_ = 1;
    nummeshes_=0;
  }

  BaseDriver :: ~BaseDriver()
  {
    ENTER_FCN( "BaseDriver::~BaseDriver" );
    //delete ptdomain_;
  }

  UInt BaseDriver::GetActSequenceStep() {
    return actSequenceStep_;
  }

  // for computation with adaptivity
  bool BaseDriver::printMeshesOrNot() {

    ENTER_FCN( "BaseDriver::DefinePrintMeshesOrNot" );
  
    EXCEPTION( "Currently not working, need change to XML-Standard" );
    bool meshesInfo=false;
  
    return meshesInfo;
  }

  void BaseDriver :: PrintSeqMeshes()
  {
    ENTER_FCN( "BaseDriver::PrintSeqMeshes" );

    Warning( "Not implemented anymore", __FILE__, __LINE__ );
  }

}
