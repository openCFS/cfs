#include <stdio.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "ParamIdent/piezoParamIdent.hh"
#include "newmark.hh"
#include "Elements/basefe.hh"
#include "blocknodeEQN.hh"
#include "scalarblockEQN.hh"
#include "Forms/forms_header.hh"
#include "Forms/linPressureInt.hh"
#include "Driver/solveStepPiezo.hh"

#include "piezoPDE.hh" 

namespace CoupledField {

  PiezoPDE::PiezoPDE( Grid *aptgrid, TimeFunc *aptTimeFunc,
                      WriteResults *aptOut ) :
    SinglePDE( aptgrid, aptOut, aptTimeFunc ) {

    ENTER_FCN( "PiezoPDE::PiezoPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================


    Error("We no longer support piezpPDE, use direct coupling!!",
	  __FILE__,__LINE__);

  }

  void PiezoPDE::ReadDampingInformation( )
  {
    ENTER_FCN( "PiezoPDE::ReadDampingInformation" );
    
  }

  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoPDE::DefineIntegrators() {
    ENTER_FCN( "PiezoPDE::DefineIntegerators" );

  }


  // **********************
  //   GetStiffIntegrator
  // **********************
  BaseForm *PiezoPDE::GetStiffIntegrator( BaseMaterial& actSDMat,
                                          Boolean reducedInt,
                                          Boolean isdamping ) {

    ENTER_FCN( "PiezoPDE::GetStiffIntegrator" );

  }

  void PiezoPDE::DefineSolveStep()
  {
    ENTER_FCN( "PiezoPDE::DefineSolveStep" );
  }


  // ======================================================
  // TRANSIENT SOLVING SECTION
  // ======================================================


  void PiezoPDE::InitTimeStepping() {
    ENTER_FCN( "PiezoPDE::InitTimeStepping" );
  }


  void PiezoPDE::WriteResultsInFile( const UInt kstep,
                                     const Double asteptime,
                                     UInt stepOffset,
                                     Double timeOffset ) {

    ENTER_FCN( "PiezoPDE::WriteResultsInFile" );

  }


  void PiezoPDE::WriteHistoryInFile( const UInt kstep,
                                     const Double asteptime,
                                     UInt stepOffset,
                                     Double timeOffset ) {

    ENTER_FCN( "PiezoPDE::WriteHIstoryInFile" );

  }



  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void PiezoPDE::ReadStoreResults() {

    ENTER_FCN( "PiezoPDE::ReadStoreResults" );

  }

  // ************************************************************
  //   PostProcess
  // ************************************************************

  void PiezoPDE::PostProcess() {

    ENTER_FCN( "PiezoPDE::PostProcess" );

  }


  void PiezoPDE::CalcEfield(){
    ENTER_FCN( "PiezoPDE::CalcEfield" );
  }


  void PiezoPDE::CalcComplexValuedEfield(){
    ENTER_FCN( "PiezoPDE::CalcComplexValuedEfield" );
  }


  void PiezoPDE::CalcStress(){
    ENTER_FCN( "PiezoPDE::CalcSress" );
  
  }
  
  void PiezoPDE::CalcComplexValuedStress(){
    ENTER_FCN( "PiezoPDE::CalcComplexValuedStress" );
  }
  
  
  void PiezoPDE::CalcCharges(){
  //   ENTER_FCN( "PiezoPDE::CalcCharges" );
  }



  void PiezoPDE::CalcComplexValuedCharges(){
    ENTER_FCN( "PiezoPDE::CalcComplexValuedCharges" );

  }

} // end of namespace
