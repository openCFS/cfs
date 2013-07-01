// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

// SuperLU_MT header
// #include "pcsp_defs.h"

#include "MatVec/BaseMatrix.hh"
#include "MatVec/StdMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"

#include "def_use_superlu.hh"
#include "SuperLUSolver.hh"

namespace CoupledField {

  // Declare logging stream and make sure that it is also available in
  // release mode by using BOOST_DECLARE_LOG() instead of DECLARE_LOG()
  BOOST_DECLARE_LOG(superluSolver)
  DEFINE_LOG(superluSolver, "olas.solvers.superlu")

  // ***********************
  //   Default Constructor
  // ***********************
  template<typename T>
  SuperLUSolver<T>::SuperLUSolver() {
    EXCEPTION( "Default constructor of SuperLUSolver is forbidden!" );
  }


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  SuperLUSolver<T>::SuperLUSolver( PtrParamNode solverNode,
                                   PtrParamNode olasInfo ) {


    // Set pointers to communication objects
    xml_ = solverNode;
    infoNode_ = olasInfo->Get("superlu",ParamNode::APPEND);

    // Set default solver type to direct sparse solver
    PtrParamNode sNode = xml_->Get("superlu", ParamNode::INSERT);
  }


  // **************
  //   Destructor
  // **************
  template<typename T>
  SuperLUSolver<T>::~SuperLUSolver() {
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void SuperLUSolver<T>::Setup( BaseMatrix &sysMat, PtrParamNode analysis_step ) {

    PtrParamNode sNode;
    sNode = xml_->Get("superlu", ParamNode::INSERT);
    
    // Determine, whether we are expected to be verbose
    LOG_TRACE(superluSolver) << " -----------------------------------------"
                             << "-------------------------------------";

    EXCEPTION("SuperLU::Setup not implemented!");
  }



  // *************************
  //   Solve linear system
  // *************************
  template<typename T>
  void SuperLUSolver<T>::Solve( const BaseMatrix &sysmat,
                                const BaseVector &rhs,
                                BaseVector &sol,
                                PtrParamNode analysis_id ) {

    LOG_TRACE(superluSolver) << " -----------------------------------------"
                             << "-------------------------------------";
    LOG_TRACE(superluSolver) << " Solving linear system ...";

    EXCEPTION("SuperLU::Solve not implemented!");

  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class SuperLUSolver<Double>;
  template class SuperLUSolver<Complex>;
#endif

}
