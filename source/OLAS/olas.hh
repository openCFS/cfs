// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_OLAS_HEADERS_HH
#define OLAS_OLAS_HEADERS_HH


//! \file This file represents the interface from OLAS to CFS++.
//! It makes the important classes avilable in the 'CoupledField' namespace
//! and imports some enum-types for communicating with OLAS_Params and
//! OLAS_Report.

#include "algsys/basesystem.hh"
#include "algsys/standardsys.hh"
#include "algsys/sbmsystem.hh"

// Note: OLAS uses a different namespace for writing out 
// the trace-, debug- and cla-stream. Therefore This
// declaration here has to be done
namespace OutInfo{

  // ------------------- Files for debug, trace and information -------------
  extern std::ostream * trace;     // basename.trace
  extern std::ostream * debug;     // basename.deb
  extern std::ostream * cla;       // basename.las
  extern std::ostream * memtrace;  // basename.mem 
  extern std::ostream * data;      // basename.data
}


namespace CoupledField {

  // import main classes
  using OLAS::BaseSystem;
  using OLAS::StandardSystem;
  using OLAS::SBM_System;
  using OLAS::OLAS_Params;
  using OLAS::OLAS_Report;


  // import some important enum typs for communication

  //! Type for identifying PDEs in OLAS
  //! \note This type is guaranteed to be of an Integer or
  //! enumeration typem, so it can be used as array index
  using OLAS::PdeIdType;


  //! Type of finite element matrix 
  using OLAS::FEMatrixType;

  using OLAS::NOTYPE;
  using OLAS::SYSTEM;
  using OLAS::STIFFNESS;
  using OLAS::DAMPING;
  using OLAS::CONVECTION;
  using OLAS::MASS;



}

#endif
