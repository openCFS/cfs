#ifdef OLAS

#include <string.h>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "CFSOLASParams.hh"

namespace CoupledField {

  // *******************
  //   SetSolverParams
  // *******************
  void CFSOLASParams::SetSolverParams( std::string pdename,
				       BaseParamHandler *cfs,
				       OLAS_Params *olas ) {

    ENTER_FCN( "CFSOLASParams::SetSolverParams" );

    // First determine the type of solver for this PDE
    std::string sTypeString;
    SolverType sType;
    cfs->Get( "type", sTypeString, pdename, "solver" );
    if ( sTypeString == "PCG" ) {
      sType = CG;
    }
    if ( sTypeString == "hyprePCG" ) {
      sType = HYPRE_PCG;
    }
    if ( sTypeString == "hypreGMRES" ) {
      sType = HYPRE_GMRES;
    }
    if ( sTypeString == "hypreBICGSTAB" ) {
      sType = HYPRE_BICGSTAB;
    }
    else {
      std::string errmsg = "Solver " + sTypeString;
      errmsg += " not supported yet.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Now determine which parameters have been set by the user
    // and insert them into the olasParams object.
    std::vector<std::string> list;

    switch( sType ) {

    case CG:
      cfs->GetList( "tol", list, pdename, "PCG" );
      if( list.size() == 1 ) {
	olas->SetValue( "CG_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "hyprePCG" );
      if( list.size() == 1 ) {
	olas->SetValue( "CG_maxIter", atoi(list[0].c_str()) );
      }
      break;

    case HYPRE_PCG:
      cfs->GetList( "tol", list, pdename, "hyprePCG" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREPCG_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "hyprePCG" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREPCG_maxIter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "hyprePCG" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREPCG_logging", atoi(list[0].c_str()) );
      }
      cfs->GetList( "printLevel", list, pdename, "hyprePCG" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREPCG_printLevel", atoi(list[0].c_str()) );
      }
      break;

    case HYPRE_GMRES:
      cfs->GetList( "tol", list, pdename, "hypreGMRES" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREGMRES_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "hypreGMRES" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREGMRES_maxIter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "maxKrylovDim", list, pdename, "hypreGMRES" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREGMRES_maxKrylovDim", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "hypreGMRES" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREGMRES_logging", atoi(list[0].c_str()) );
      }
      cfs->GetList( "printLevel", list, pdename, "hypreGMRES" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREGMRES_printLevel", atoi(list[0].c_str()) );
      }
      break;

    case HYPRE_BICGSTAB:
      cfs->GetList( "tol", list, pdename, "hypreBICGSTAB" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREBICGSTAB_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "hypreBICGSTAB" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREBICGSTAB_maxIter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "hypreBICGSTAB" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREBICGSTAB_logging", atoi(list[0].c_str()) );
      }
      cfs->GetList( "printLevel", list, pdename, "hypreBICGSTAB" );
      if( list.size() == 1 ) {
	olas->SetValue( "HYPREBICGSTAB_printLevel", atoi(list[0].c_str()) );
      }
      break;

    }

  }


  // ********************
  //   SetPrecondParams
  // ********************
  void CFSOLASParams::SetPrecondParams( std::string pdename,
					BaseParamHandler *cfs,
					OLAS_Params *olas ) {

    ENTER_FCN( "CFSOLASParams::SetPrecondParams" );

    // First determine the type of preconditioner for this PDE
    std::string pTypeString;
    PrecondType pType;
    cfs->Get( "precond", pTypeString, pdename, "solver" );
    if ( pTypeString == "hypreSPAI" ) {
      pType = HYPRE_SPAI;
    }
    else if ( pTypeString == "hypreILU" ) {
      pType = HYPRE_ILU;
    }
    else if ( pTypeString == "hypreAMG" ) {
      pType = HYPRE_AMG;
    }
    else {
      std::string errmsg = "Preconditioner " + pTypeString;
      errmsg += " not supported yet.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Now determine which parameters have been set by the user
    // and insert them into the olasParams object.
    std::vector<std::string> list;

    switch( pType ) {

    case HYPRE_ILU:
      cfs->GetList( "level", list, pdename, "hypreILU" );
      if( list.size() == 1 ) {
	olas->SetValue( "EUCLID_level",  atoi(list[0].c_str()) );
      }
      cfs->GetList( "stats", list, pdename, "hypreILU" );
      if( list.size() == 1 ) {
	bool stats = (list[0] == "yes");
	olas->SetValue( "EUCLID_stats", stats );
      }
      cfs->GetList( "memory", list, pdename, "hypreILU" );
      if( list.size() == 1 ) {
	Integer memory = (list[0] == "yes");
	olas->SetValue( "EUCLID_memory", memory );
      }
      break;

    case HYPRE_SPAI:
      cfs->GetList( "thresh", list, pdename, "hypreSPAI" );
      if( list.size() == 1 ) {
	olas->SetValue( "PARASAILS_thresh", atof(list[0].c_str()) );
      }
      cfs->GetList( "levels", list, pdename, "hypreSPAI" );
      if( list.size() == 1 ) {
	olas->SetValue( "PARASAILS_levels", atoi(list[0].c_str()) );
      }
      cfs->GetList( "filter", list, pdename, "hypreSPAI" );
      if( list.size() == 1 ) {
	olas->SetValue( "PARASAILS_filter", atof(list[0].c_str()) );
      }
      cfs->GetList( "symmetry", list, pdename, "hypreSPAI" );
      if( list.size() == 1 ) {
	olas->SetValue( "PARASAILS_symmetry", atoi(list[0].c_str()) );
      }
      cfs->GetList( "loadBalance", list, pdename, "hypreSPAI" );
      if( list.size() == 1 ) {
	Integer balance = (list[0] == "yes");
	olas->SetValue( "PARASAILS_loadBalance", balance );
      }
      cfs->GetList( "reuse", list, pdename, "hypreSPAI" );
      if( list.size() == 1 ) {
	Integer balance = (list[0] == "yes");
	olas->SetValue( "PARASAILS_reuse", balance );
      }
      break;

    case HYPRE_AMG:
      cfs->GetList( "maxLevels", list, pdename, "hypreAMG" );
      if( list.size() == 1 ) {
	olas->SetValue( "BOOMERAMG_maxLevels", atoi(list[0].c_str()) );
      }
      cfs->GetList( "alpha", list, pdename, "hypreAMG" );
      if( list.size() == 1 ) {
	olas->SetValue( "BOOMERAMG_alpha", atof(list[0].c_str()) );
      }
      cfs->GetList( "numSweeps", list, pdename, "hypreAMG" );
      if( list.size() == 1 ) {
	olas->SetValue( "BOOMERAMG_numSweeps", atoi(list[0].c_str()) );
      }
      cfs->GetList( "cycleType", list, pdename, "hypreAMG" );
      if( list.size() == 1 ) {
	Integer cycleType = list[0] == "W-cycle" ? 2 : 1;
	olas->SetValue( "BOOMERAMG_cycleType", cycleType );
      }
      break;
    }

    olas->ShowPool( OLAS_Params::INT_POOL     , std::cerr );
    olas->ShowPool( OLAS_Params::DOUBLE_POOL  , std::cerr );
    olas->ShowPool( OLAS_Params::BOOLEAN_POOL , std::cerr );
    olas->ShowPool( OLAS_Params::STRING_POOL  , std::cerr );
    olas->ShowPool( OLAS_Params::ENUM_POOL    , std::cerr );

  }

}
#endif
