// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     cfsdat.cc
 *       \brief    <Description>
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================

//System includes
#include <iostream>
#include <stdio.h>


//CFS includes
#include <def_cfs_stats.hh>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"

#include "General/defs.hh"
#include "General/Environment.hh"

//CFSDatIncludes
#include "cfsdat/Filters/BaseFilter.hh"
#include "cfsdat/Filters/Input/InputFilterGeneric.hh"

namespace CFSDat {

  namespace po = boost::program_options;

  void GenerateParamNode(int argc, char** argv, PtrParamNode param){

  }

}

int main(int argc, char** argv)
{

  //Create dummy program options
  const char * myOpt[] = { "CFSDatComputation" };
  CoupledField::progOpts = new CoupledField::ProgramOptions(1,myOpt);

  shared_ptr<CFSDat::BaseFilter> newFilter(new CFSDat::InputFilterGeneric() );
return 0;
}
