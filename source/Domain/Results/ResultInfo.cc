#include "ResultInfo.hh"
#include <iostream>

#include "Domain/ElemMapping/EntityLists.hh"
#include "General/Environment.hh"

namespace CoupledField {

  ResultInfo::ResultInfo() {

    resultType = NO_SOLUTION_TYPE;
    dofNames.Clear();
    unit = "";
    complexFormat = AMPLITUDE_PHASE;
    entryType = UNKNOWN;
    definedOn = FREE;
    fromOptimization = false;
    isStatic = false;
  }


  UInt ResultInfo::GetDofIndex( const std::string & dof ) const {

    Integer pos = dofNames.Find( dof );
    if( pos < 0  ) {
      EXCEPTION( "Dof with name '" << dof << "' not found!" );
    }
    return (UInt) (pos);
  }


  void ResultInfo::SetVectorDOFs(UInt dim, bool is_axi)
  {
    if(dim == 3)
      dofNames = "x", "y", "z";
    if(dim == 2 && !is_axi)
      dofNames = "x", "y";
    if(dim == 2 && is_axi)
      dofNames = "r", "z";
  }

  void ResultInfo::SetVectorDOFs(UInt dim, bool is_axi, bool is2p5)
  {
    if(dim == 3 || is2p5)
      dofNames = "x", "y", "z";
    else
    {
      if(dim == 2 && !is_axi)
        dofNames = "x", "y";
      if(dim == 2 && is_axi)
        dofNames = "r", "z";
    }
  }

  std::string ResultInfo::GetDofName( const UInt dof ) const {
    if( dof <= 0 || dof > dofNames.GetSize()+1 ) {
      EXCEPTION( "'dof' must be in the range of [1.." 
                 << dofNames.GetSize()+1 << "]!" );
    }
    
    return dofNames[dof];
  }

  std::string ResultInfo::ToString() const 
  {
    std::ostringstream os;
    os << "name = " << resultName << " dofs " << dofNames.GetSize();
    return os.str();
  }



  ResultInfo& ResultInfo::operator=( const ResultInfo& data ) {

    resultType = data.resultType;
    resultName = data.resultName;
    dofNames = data.dofNames;
    unit = data.unit;
    isStatic = data.isStatic;
    complexFormat = data.complexFormat;
    entryType = data.entryType;
    definedOn = data.definedOn;
    return *this;
  }


  void ResultInfo::Enum2String(EntityUnknownType in, 
                               std::string& out ) {
    
    switch( in ) {
      
    case NODE: 
      out = "node";
      break;
    case EDGE:
      out = "edge";
      break;
    case FACE:  
      out = "face";
      break;
    case ELEMENT:  
      out = "element";
      break;
    case SURF_ELEM:
      out = "surfElement";
      break;
    case REGION:
      out = "region";
      break;
    case REGION_AVERAGE:
      out = "regionAverage";
      break;
    case SURF_REGION:  
      out = "surfRegion";
      break;
    case NODELIST:
      out = "nodeList";
      break;
    case COIL:
      out = "coil";
      break;
    case FREE:
      out = "free";
      break;
    default:
      EXCEPTION( "Conversion of " << in
                 << " to EntityUnkwownType not possible" );
    }
    
  }
  
  void  ResultInfo::String2Enum( const std::string& in,
                                  EntityUnknownType& out ) {

    if( in == "node")
      out =  NODE;
      else if( in == "edge")
        out = EDGE;
      else if( in == "face")
        out = FACE;
      else if( in == "element")
        out = ELEMENT; 
      else if( in == "surfElement")
        out = SURF_ELEM;
      else if( in == "region")
        out = REGION;
      else if( in == "regionAverage")
        out = REGION_AVERAGE;
      else if( in == "surfRegion")
        out = SURF_REGION;
      else if( in == "nodeList")
        out = NODELIST;
      else if( in == "coil")
        out = COIL;
      else if( in == "free")
        out = FREE;
      else {
        EXCEPTION( "Can not convert '" << in << "' to EntityUnknownType");
      }
  }

  bool operator==( const ResultInfo& a, const ResultInfo& b ) {
    return a.resultType == b.resultType && a.dofNames == b.dofNames &&
           a.definedOn == b.definedOn && a.entryType == b.entryType;
  }

  bool operator!=( const ResultInfo& a, const ResultInfo& b ) {
    return !(a == b);
  }

  bool operator<( const ResultInfo& a, const ResultInfo& b ) {
    //return !(a==b);
    return ( a.resultType  < b.resultType );
  }
  // ************************************************************************
  // ENUM INITIALIZATION
  // ************************************************************************

  // Definition of unknown entity types
  static EnumTuple unknownTypeTuples[] = 
  {
   EnumTuple(ResultInfo::NODE,  "NODE"), 
   EnumTuple(ResultInfo::EDGE,  "EDGE"),
   EnumTuple(ResultInfo::FACE,  "FACE"),
   EnumTuple(ResultInfo::ELEMENT,  "ELEMENT"),
   EnumTuple(ResultInfo::SURF_ELEM,  "SURF_ELEM"),
   EnumTuple(ResultInfo::REGION,  "REGION"),
   EnumTuple(ResultInfo::REGION_AVERAGE,  "REGION_AVERAGE"),
   EnumTuple(ResultInfo::SURF_REGION,  "SURF_REGION"),
   EnumTuple(ResultInfo::NODELIST,  "NODELIST"),
   EnumTuple(ResultInfo::COIL,  "COIL"),
   EnumTuple(ResultInfo::FREE,  "FREE")
  };

  Enum<ResultInfo::EntityUnknownType> ResultInfo::EntityUnknownTypeEnum_ = \
      Enum<ResultInfo::EntityUnknownType>("Entry type of Results",
                                          sizeof(unknownTypeTuples) / sizeof(EnumTuple),
                                          unknownTypeTuples);

  // Definition of entry types of resultinfo
  static EnumTuple entryTypeTuples[] = 
  {
   EnumTuple(ResultInfo::UNKNOWN,  "NO_DIM"), 
   EnumTuple(ResultInfo::SCALAR,  "SCALAR"),
   EnumTuple(ResultInfo::VECTOR,  "VECTOR"),
   EnumTuple(ResultInfo::TENSOR,  "TENSOR"),
   EnumTuple(ResultInfo::STRING,  "STRING")
  };

  Enum<ResultInfo::EntryType> ResultInfo::EntryTypeEnum_ = \
      Enum<ResultInfo::EntryType>("Entry type of Results",
                                  sizeof(entryTypeTuples) / sizeof(EnumTuple),
                                  entryTypeTuples);

  ResultInfo::EntityUnknownType ResultInfo::MapSolTypeToDefinedOn(SolutionType solType) // static
  {
    // Map each SolutionType to its canonical EntityUnknownType (definedOn)
    // Organized by physical domain and entity type
    
    switch(solType)
    {
      // ============================================
      // PRIMARY SOLUTIONS - NODAL
      // ============================================
      // MECHANICAL
      case MECH_DISPLACEMENT:
      case MECH_VELOCITY:
      case MECH_ACCELERATION:
      case MECH_RHS_LOAD:
        return ResultInfo::NODE;
      
      // ACOUSTIC
      case ACOU_POTENTIAL:
      case ACOU_PRESSURE:
      case ACOU_PRESSURE_DERIV_1:
      case ACOU_PRESSURE_DERIV_2:
      case ACOU_VELOCITY:
      case ACOU_ACCELERATION:
      case ACOU_FORCE:
      case ACOU_RHS_LOAD:
      case ACOU_RHSVAL:
        return ResultInfo::NODE;
      
      // ACOUSTIC - TDEF auxiliary (nodal)
      case ACOU_TDEF_PHI_C_1: case ACOU_TDEF_PHI_C_2: case ACOU_TDEF_PHI_C_3:
      case ACOU_TDEF_PHI_C_4: case ACOU_TDEF_PHI_C_5: case ACOU_TDEF_PHI_C_6:
      case ACOU_TDEF_PHI_C_7: case ACOU_TDEF_PHI_C_8: case ACOU_TDEF_PHI_C_9:
      case ACOU_TDEF_PHI_C_10: case ACOU_TDEF_PHI_C_11: case ACOU_TDEF_PHI_C_12:
      case ACOU_TDEF_PHI_C_13: case ACOU_TDEF_PHI_C_14: case ACOU_TDEF_PHI_C_15:
      case ACOU_TDEF_PSI_C_1: case ACOU_TDEF_PSI_C_2: case ACOU_TDEF_PSI_C_3:
      case ACOU_TDEF_PSI_C_4: case ACOU_TDEF_PSI_C_5: case ACOU_TDEF_PSI_C_6:
      case ACOU_TDEF_PSI_C_7: case ACOU_TDEF_PSI_C_8: case ACOU_TDEF_PSI_C_9:
      case ACOU_TDEF_PSI_C_10: case ACOU_TDEF_PSI_C_11: case ACOU_TDEF_PSI_C_12:
      case ACOU_TDEF_PSI_C_13: case ACOU_TDEF_PSI_C_14: case ACOU_TDEF_PSI_C_15:
      case ACOU_TDEF_PHI_V_1: case ACOU_TDEF_PHI_V_2: case ACOU_TDEF_PHI_V_3:
      case ACOU_TDEF_PHI_V_4: case ACOU_TDEF_PHI_V_5: case ACOU_TDEF_PHI_V_6:
      case ACOU_TDEF_PHI_V_7: case ACOU_TDEF_PHI_V_8: case ACOU_TDEF_PHI_V_9:
      case ACOU_TDEF_PHI_V_10: case ACOU_TDEF_PHI_V_11: case ACOU_TDEF_PHI_V_12:
      case ACOU_TDEF_PHI_V_13: case ACOU_TDEF_PHI_V_14: case ACOU_TDEF_PHI_V_15:
      case ACOU_TDEF_PSI_V_1: case ACOU_TDEF_PSI_V_2: case ACOU_TDEF_PSI_V_3:
      case ACOU_TDEF_PSI_V_4: case ACOU_TDEF_PSI_V_5: case ACOU_TDEF_PSI_V_6:
      case ACOU_TDEF_PSI_V_7: case ACOU_TDEF_PSI_V_8: case ACOU_TDEF_PSI_V_9:
      case ACOU_TDEF_PSI_V_10: case ACOU_TDEF_PSI_V_11: case ACOU_TDEF_PSI_V_12:
      case ACOU_TDEF_PSI_V_13: case ACOU_TDEF_PSI_V_14: case ACOU_TDEF_PSI_V_15:
      case ACOU_PMLAUXSCALAR:
      case ACOU_PMLAUXVEC:
        return ResultInfo::NODE;
      
      // HEAT
      case HEAT_TEMPERATURE:
      case HEAT_TEMPERATURE_D1:
      case HEAT_RHS_LOAD:
        return ResultInfo::NODE;
      
      // ELECTRIC
      case ELEC_POTENTIAL:
      case ELEC_POTENTIAL_DERIV1:
      case ELEC_POTENTIAL_DERIV_1:
        return ResultInfo::NODE;
      
      // ELECTRIC QUASI-STATIC (nodal)
      case ELEC_RHS_LOAD:
        return ResultInfo::NODE;
      
      // MAGNETIC
      case MAG_POTENTIAL:
      case MAG_POTENTIAL_DERIV1:
      case MAG_TOTAL_POTENTIAL:
      case MAG_REDUCED_POTENTIAL:
      case MAG_RHS_LOAD:
        return ResultInfo::NODE;
      
      // COIL
      case COIL_INDUCTANCE:
      case COIL_INDUCED_VOLTAGE:
      case COIL_CURRENT:
      case COIL_VOLTAGE:
      case COIL_VOLTAGE_INTEGRAL:
      case COIL_LINKED_FLUX:
      case COIL_CURRENT_DERIV1:
        return ResultInfo::COIL;
      
      // FLUID MECHANICS
      case FLUIDMECH_PRESSURE:
      case FLUIDMECH_VELOCITY:
      case FLUIDMECH_MESH_VELOCITY:
      case FLUIDMECH_MESH_VELOCITY_NODE:
      case FLUIDMECH_ZERO_PRESSURE:
      case MEAN_FLUIDMECH_VELOCITY:
      case MEAN_FLUIDMECH_VELOCITY_NORMAL:
      case LAMBDA_K:
        return ResultInfo::NODE;
      
      // SMOOTH
      case SMOOTH_DISPLACEMENT:
      case SMOOTH_VELOCITY:
      case SMOOTH_ZERO_PRESSURE:
        return ResultInfo::NODE;
      
      // WATER WAVES
      case WATER_PRESSURE:
      case WATER_PMLAUXSCALAR:
      case WATER_PMLAUXVEC:
      case WATER_RHS_LOAD:
        return ResultInfo::NODE;
      
      // SPLITTING (nodal primary)
      case SPLIT_SCALAR:
      case SPLIT_VECTOR:
      case SPLIT_RHS_LOAD:
      case SPLIT_SCALAR_VELOCITY:
      case SPLIT_VECTOR_VELOCITY:
        return ResultInfo::NODE;
      
      // TEST
      case TEST_DOF:
      case TEST_RHS_LOAD:
        return ResultInfo::NODE;
      
      // PERTURBATION (nodal primary)
      // Most pressure/velocity primary fields are NODE
      
      // ============================================
      // ELEMENT-BASED RESULTS
      // ============================================
      // MECHANICAL - VELOCITY, STRESSES AND STRAINS
      case MECH_VELOCITY_ELEM: 
      case MECH_STRESS:
      case MECH_PRINCIPAL_STRESS:
      case MECH_PRINCIPAL_STRESS_MIN:
      case MECH_PRINCIPAL_STRESS_MAX:
      case MECH_PRINCIPAL_STRESS_MED:
      case MECH_PRINCIPAL_STRESS_MIN_SCAL:
      case MECH_PRINCIPAL_STRESS_MAX_SCAL:
      case MECH_PRINCIPAL_STRESS_MED_SCAL:
      case MECH_STRAIN:
      case MECH_PRINCIPAL_STRAIN:
      case MECH_PRINCIPAL_STRAIN_MIN:
      case MECH_PRINCIPAL_STRAIN_MAX:
      case MECH_PRINCIPAL_STRAIN_MED:
      case MECH_PRINCIPAL_STRAIN_MIN_SCAL:
      case MECH_PRINCIPAL_STRAIN_MAX_SCAL:
      case MECH_PRINCIPAL_STRAIN_MED_SCAL:
      case MECH_THERMAL_STRAIN:
      case VON_MISES_STRESS:
      case VON_MISES_STRAIN:
      case MECH_DYADIC_STRAIN:
      case MECH_QUAD_DISP:
      case MECH_IRR_STRESS:
      case MECH_IRR_STRAIN:
      case MECH_THERMAL_STRESS:
        return ResultInfo::ELEMENT;
      
      // MECHANICAL - ENERGY DENSITIES (element-based)
      case MECH_KIN_ENERGY_DENS:
      case MECH_DEFORM_ENERGY_DENS:
      case MECH_TOTAL_ENERGY_DENS:
      case MECH_STRUCT_INTENSTIY:
        return ResultInfo::ELEMENT;
      
      // ACOUSTIC - ELEMENT BASED
      //case ACOU_VELOCITY:  // element variant
      case ACOU_POSITION:
      case ACOU_ELEM_SPEED_OF_SOUND:
      case ACOU_INTENSITY:
      case ACOU_POWERDENSITY:
      case ACOU_DIV_LH_TENSOR:
      case ACOU_RHS_LOAD_DENSITY:
      case ACOU_RHS_LOAD_DENSITY_VECTOR:
      case ACOU_RHS_LOAD_DENSITY_VECTOR_IN_NORMAL:
      case ACOU_RHS_LOADP:
      case ACOU_APE_RHS_LOAD:
      case ACOU_VORTEX_RHS_LOAD:
        return ResultInfo::ELEMENT;
      
      // HEAT - ELEMENT BASED
      case HEAT_FLUX_DENSITY:
      case HEAT_SOURCE_DENSITY:
      case HEAT_FLUX_INTENSITY:
        return ResultInfo::ELEMENT;
      
      // ELECTRIC - ELEMENT BASED
      case ELEC_FIELD_INTENSITY:
      case ELEC_ENERGY_DENSITY:
      case ELEC_FLUX_DENSITY:
      case ELEC_POLARIZATION:
      case ELEC_POWER_DENSITY:
      case ELEC_CURRENT_DENSITY:
      case ELEC_NORMAL_CURRENT_DENSITY:
      case ELEC_GRAD_V_INT:
      case ELEC_ELEM_PERMITTIVITY:
      case ELEC_COND_TENSOR:
      //case DISPLACEMENT_CURRENT_DENSITY:
      case DISPLACEMENT_CURRENT_FIELD_INTENSITY:
      case DISPLACEMENT_NORMAL_CURRENT_DENSITY:
      case ELECTRIC_AND_DISPLACEMENT_CURRENT_DENSITY:
      case ELECTRIC_AND_DISPLACEMENT_NORMAL_CURRENT_DENSITY:
      case ELEC_FIELD_INTENSITY_TRANSVERSAL:
      case ELEC_FIELD_INTENSITY_LONGITUDINAL:
        return ResultInfo::ELEMENT;
      
      // ELECTRIC QUASI-STATIC
      case ELEC_PSEUDO_POLARIZATION:
      case ELEC_PHYSICAL_PSEUDO_DENSITY:
        return ResultInfo::ELEMENT;
      
      // MAGNETIC - ELEMENT BASED
      case MAG_FLUX_DENSITY:
      case MAG_FIELD_INTENSITY:
      case MAG_EDDY_CURRENT_DENSITY:
      case MAG_TOTAL_CURRENT_DENSITY:
      case MAG_COIL_CURRENT_DENSITY:
      case MAG_POTENTIAL_DIV:
      case MAG_FORCE_LORENTZ_DENSITY:
      case MAG_FORCE_LORENTZ_DENSITY_STATIC:
      case MAG_FORCE_LORENTZ_DENSITY_HARMONIC:
      case MAG_FORCE_MAXWELL_DENSITY:
      case MAG_NORMALFORCE_MAXWELL_DENSITY:
      case MAG_TANGENTIALFORCE_MAXWELL_DENSITY:
      case MAG_EDDY_POWER_DENSITY:
      case MAG_ENERGY_DENSITY:
      case MAG_CORE_LOSS_DENSITY:
      case MAG_JOULE_LOSS_POWER_DENSITY:
      case MAG_ELEM_PERMEABILITY:
      case MAG_ELEM_RELUCTIVITY:
      case MAG_MAGNETIZATION:
      case MAG_POLARIZATION:
      case FLUX_INDUCED_STRAIN:
        return ResultInfo::ELEMENT;
      
      // FLUID MECHANICS - ELEMENT BASED
      case FLUIDMECH_MESH_VELOCITY_ELEM:
      case FLUIDMECH_TOTAL_VELOCITY_ELEM:
      case FLUIDMECH_DENSITY:
      case FLUIDMECH_STRESS:
      case FLUIDMECH_COMP_STRESS:
      case FLUIDMECH_VISC_STRESS:
      case FLUIDMECH_PRES_TENS:
      case FLUIDMECH_TOTAL_STRESS:
      case FLUIDMECH_STRAINRATE:
      case FLUIDMECH_VELOCITY_DERIV_1:
      case FLUIDMECH_VELOCITY_DERIV_2:
      case DIV_MEAN_FLUIDMECH_VELOCITY:
      case FLUIDMECH_TKE:
      case FLUIDMECH_TDR:
      case FLUIDMECH_TEF:
      case FLUIDMECH_WVT:
      case FLUIDMECH_WVT_DENSITY:
      case FLUIDMECH_VORTICITY:
      case FLUIDMECH_WVT_PHI:
      case FLUIDMECH_WVT_DENSITY_PHI:
      case FLUIDMECH_WVT_U1:
      case FLUIDMECH_WVT_U2:
      case FLUIDMECH_WVT_F:
      case FLUIDMECH_VISCOUS_DISS_POWER_DENS:
      case FLUIDMECH_VISCOUS_DISS_POWER_DENS_DIV:
      case FLUIDMECH_VISCOUS_DISS_POWER_DENS_STRAIN:
      case FLUIDMECH_STABILPARAM:
        return ResultInfo::ELEMENT;
      
      // SMOOTH - ELEMENT BASED
      case SMOOTH_STRAIN:
      case SMOOTH_DEFORM_ENERGY_DENS:
        return ResultInfo::ELEMENT;
      
      // WATER WAVES - ELEMENT BASED
      case WATER_POSITION:
      case WATER_PRES_TENS:
      case WATER_SURFACE_TORQUE_DENSITY:
      case WATER_TDT:
        return ResultInfo::ELEMENT;
      
      // TEST
      case TEST_FIELD:
        return ResultInfo::ELEMENT;
      
      // OPTIMIZATION (element-based pseudo density)
      case MECH_PSEUDO_DENSITY:
      case PSEUDO_DENSITY:
      case PHYSICAL_PSEUDO_DENSITY:
      case MECH_TENSOR_TRACE:
      case MECH_TENSOR:
      case MECH_TENSOR_HILL_MANDEL:
      case MECH_ELEM_VOL:
      case MECH_ELEM_POROSITY:
      case ELEC_TENSOR_TRACE:
      case ELEC_TENSOR:
      case RHS_PSEUDO_DENSITY:
      case PHYSICAL_RHS_PSEUDO_DENSITY:
      case LBM_PHYSICAL_PSEUDO_DENSITY:
      case ELEM_DENSITY:
        return ResultInfo::ELEMENT;
      
      // PML
      case PML_DAMP_FACTOR:
      case PML_TENSOR:
      case PML_DETERMINANT:
      case PML_DISTANCE:
        return ResultInfo::ELEMENT;
      
      // GEOMETRY
      case ELEM_LOC_DIR:
      case JACOBIAN:
      case ASPECT_RATIO:
      case VOLUME:
        return ResultInfo::ELEMENT;
      
      // LBM
      case LBM_PROBABILITY_DISTRIBUTION:
      case LBM_VELOCITY:
      case LBM_DENSITY:
      case LBM_PRESSURE:
        return ResultInfo::ELEMENT;
      
      // ============================================
      // SURFACE ELEMENT - BOUNDARY QUANTITIES
      // ============================================
      case MECH_NORMAL_DISPLACEMENT:
      case MECH_NORMAL_STRESS:
      case MECH_NORMAL_VELOCITY:
      case MECH_NORMAL_STRUCT_INTENSITY:
      case MECH_FORCE:
        return ResultInfo::SURF_ELEM;
      
      // ACOUSTIC - SURFACE
      case ACOU_SURFPRESSURE:
      case ACOU_NORMAL_VELOCITY:
      case ACOU_NORMAL_INTENSITY:
      case ACOU_SURFINTENSITY:
      case ACOU_MIXED_MASS_LOAD:
      case ACOU_MIXED_MOMENTUM_LOAD:
      //case ACOUSSURF_RHSVAL:
        return ResultInfo::SURF_ELEM;
      
      // HEAT - SURFACE
      //case HEAT_FLUX_NORMAL:
      //  return ResultInfo::SURF_ELEM;
      
      // ELECTRIC - SURFACE
      case ELEC_FORCE_DENSITY:
      case ELEC_SURFACE_CHARGE_DENSITY:
      case ELEC_FORCE:
      case DISPLACEMENT_CURRENT_SURF:
        return ResultInfo::SURF_ELEM;
      
      // MAGNETIC - SURFACE
      case MAG_FLUX_DENSITY_SURF:
      case MAG_NORMAL_FLUX_DENSITY:
        return ResultInfo::SURF_ELEM;
      
      // FLUID MECHANICS - SURFACE
      case FLUIDMECH_NORMAL_VELOCITY:
      case FLUIDMECH_SURFACE_TRACTION:
      case FLUIDMECH_NORMAL_INTENSITY:
      case FLUIDMECH_NORMAL_INTENSITY_PRESSURE_ONLY:
      case FLUIDMECH_SURFIMPEDANCE:
      case FLUIDMECH_SURFINTENSITY:
      case FLUIDMECH_SURFINTENSITY_PRESSURE_ONLY:
        return ResultInfo::SURF_ELEM;
      
      // WATER WAVES - SURFACE
      case WATER_SURFACE_TRACTION:
      case WATER_SURFACE_TORQUE:
      case WATER_SURFACE_FORCE:
        return ResultInfo::SURF_ELEM;
      
      // ============================================
      // REGION - INTEGRATED QUANTITIES
      // ============================================
      case MECH_KIN_ENERGY:
      case MECH_DEFORM_ENERGY:
      case MECH_TOTAL_ENERGY:
      case MECH_QUAD_DISP_SUM:
      case MECH_DYADIC_STRAIN_SUM:
        return ResultInfo::REGION;
      
      // ACOUSTIC - REGION
      case ACOU_ENERGY:
      case ACOU_POT_ENERGY:
      case ACOU_KIN_ENERGY:
      case ACOU_LAMB_RHS:
        return ResultInfo::REGION;
      
      // HEAT - REGION
      case HEAT_MEAN_TEMPERATURE:
      case HEAT_CONDUCTIVITY_TENSOR_HOM:
        return ResultInfo::REGION;
      
      // ELECTRIC - REGION
      case ELEC_ENERGY:
        return ResultInfo::REGION;
      
      // MAGNETIC - REGION
      case MAG_ENERGY:
      case MAG_EDDY_POWER:
      case MAG_CORE_LOSS:
        return ResultInfo::REGION;
      
      // FLUID MECHANICS - REGION
      case FLUIDMECH_ENERGY:
      case FLUIDMECH_VISCOUS_DISS_POWER:
        return ResultInfo::REGION;
      
      // SMOOTH - REGION
      case SMOOTH_DEFORM_ENERGY:
        return ResultInfo::REGION;
      
      // ============================================
      // SURFACE REGION - INTEGRATED SURFACE QUANTITIES
      // ============================================
      case MECH_POWER:
      case MECH_DEF_SURF_VOLUME:
      case MECH_WEIGHT:
        return ResultInfo::SURF_REGION;
      
      // ACOUSTIC - SURFACE REGION
      case ACOU_POWER:
      case ACOU_POWER_PLANEWAVE:
        return ResultInfo::SURF_REGION;
      
      // HEAT - SURFACE REGION
      case HEAT_FLUX:
        return ResultInfo::SURF_REGION;
      
      // ELECTRIC - SURFACE REGION
      case ELEC_CHARGE:
      case ELEC_CURRENT:
      case ELEC_POWER:
        return ResultInfo::SURF_REGION;
      
      // ELECTRIC QUASI-STATIC
      case ELEC_CURRENT_SURF:
        return ResultInfo::SURF_REGION;
      
      // MAGNETIC - SURFACE REGION
      case MAG_FLUX:
      case MAG_FORCE_VWP:
      case MAG_FORCE_LORENTZ:
      case MAG_FORCE_MAXWELL:
      case MAG_EDDY_CURRENT:
      case MAG_EDDY_CURRENT1:
      case MAG_EDDY_CURRENT2:
      case MAG_TOTAL_CURRENT:
        return ResultInfo::SURF_REGION;
      
      // FLUID MECHANICS - SURFACE REGION
      case FLUIDMECH_FORCE:
      case FLUIDMECH_PRESSURE_DERIV_1:
      case FLUIDMECH_PRESSURE_DERIV_2:
      case FLUIDMECH_PRESSURE_TIME_DERIV_1:
      case FLUIDMECH_PRESSURE_TIME_DERIV_2:
      case FLUIDMECH_INTENSITY:
      case FLUIDMECH_INTENSITY_PRESSURE_ONLY:
      case FLUIDMECH_POWER:
      case FLUIDMECH_POWER_PRESSURE_ONLY:
      case FLUIDMECH_IMPEDANCE:
        return ResultInfo::SURF_REGION;
      
      // WATER WAVES - SURFACE REGION
      //case WATER_SURFACE_FORCE:
      //  return ResultInfo::SURF_REGION;
      
      // SMOOTH - SURFACE REGION
      case SMOOTH_CONTACT_FORCE:
      case SMOOTH_CONTACT_FORCE_DENSITY:
        return ResultInfo::SURF_REGION;
      
      // NODE_LIST for nodal list results
      //case NODELIST:
      //  return ResultInfo::NODELIST;
      
      // GENERIC and misc results default to FREE
      case HOMOGENIZED_TENSOR:
      case OPT_RESULT_1: case OPT_RESULT_2: case OPT_RESULT_3: case OPT_RESULT_4: case OPT_RESULT_5:
      case OPT_RESULT_6: case OPT_RESULT_7: case OPT_RESULT_8: case OPT_RESULT_9: case OPT_RESULT_10:
      case OPT_RESULT_11: case OPT_RESULT_12: case OPT_RESULT_13: case OPT_RESULT_14: case OPT_RESULT_15:
      case OPT_RESULT_16: case OPT_RESULT_17: case OPT_RESULT_18: case OPT_RESULT_19: case OPT_RESULT_20:
      case OPT_RESULT_21: case OPT_RESULT_22: case OPT_RESULT_23: case OPT_RESULT_24: case OPT_RESULT_25:
      case OPT_RESULT_26: case OPT_RESULT_27: case OPT_RESULT_28: case OPT_RESULT_29: case OPT_RESULT_30:
      case OPT_RESULT_31: case OPT_RESULT_32: case OPT_RESULT_33: case OPT_RESULT_34: case OPT_RESULT_35:
      case OPT_RESULT_36: case OPT_RESULT_37: case OPT_RESULT_38: case OPT_RESULT_39: case OPT_RESULT_40:
      case OPT_RESULT_41: case OPT_RESULT_42: case OPT_RESULT_43: case OPT_RESULT_44: case OPT_RESULT_45:
      case OPT_RESULT_46: case OPT_RESULT_47: case OPT_RESULT_48: case OPT_RESULT_49: case OPT_RESULT_50:
      case OPT_RESULT_51: case OPT_RESULT_52: case OPT_RESULT_53: case OPT_RESULT_54: case OPT_RESULT_55:
      case OPT_RESULT_56: case OPT_RESULT_57: case OPT_RESULT_58: case OPT_RESULT_59: case OPT_RESULT_60:
      case OPT_RESULT_61: case OPT_RESULT_62: case OPT_RESULT_63: case OPT_RESULT_64: case OPT_RESULT_65:
      case OPT_RESULT_66:
      case LAGRANGE_MULT:
      case LAGRANGE_MULT_DERIV_1:
      case LAGRANGE_MULT_DERIV_2:
      case LAGRANGE_MULT_1:
      case THERMOMECH_FORCE:
      case THERMOELEC_FORCE:
      case GRAD_ACOU_SOLUTION:
      case GRAD_ELEC_POTENTIAL:
      case GRAD_ELEC_POTENTIAL_DERIV1:
      case GENERIC_RESULT_0: case GENERIC_RESULT_1: case GENERIC_RESULT_2: case GENERIC_RESULT_3: case GENERIC_RESULT_4:
      case GENERIC_RESULT_5: case GENERIC_RESULT_6: case GENERIC_RESULT_7: case GENERIC_RESULT_8: case GENERIC_RESULT_9:
      case NO_SOLUTION_TYPE:
      case INVALID_SOLUTION_TYPE:
      case SPLIT_POT_ENERGY:
      case SPLIT_LAMB:
      case SPLIT_DIVLAMB:
      case COOSY_X: case COOSY_Y: case SURFACE_NORMAL: case AREA: case ETA: case XI: case NODE_NORMAL:
      default:
        // Fallback for unmapped types
        return ResultInfo::FREE;
    }
  }

}
