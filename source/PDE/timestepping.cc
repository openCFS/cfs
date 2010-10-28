// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "timestepping.hh"

#include "OLAS/algsys/basesystem.hh"

namespace CoupledField {

  TimeStepping::TimeStepping(BaseSystem * algebraicsystem )
  {
    algsys_  = algebraicsystem;
    rhsSize_ = 0 ;
  }

  const std::map<FEMatrixType,Double>  &
  TimeStepping::GetEffSysMatFactors( ) const {
    return matrix_factors_;
  }
                                              

  TimeStepping::~TimeStepping()
  {
  }

  bool TimeStepping::FeMatrixPresent( FEMatrixType type) {
    std::set<FEMatrixType> matTypes;    
    algsys_->GetFEMatrixTypes(matTypes);
    
    if ( matTypes.find(type) != matTypes.end() ) {
      return true;
    } else {
      return false;
    }
  }

  SolutionType
  TimeStepping::mapDerivToSolutionType(const SolutionType solType, const DERIVType derivType) const
  {
    SolutionType retSolType = NO_SOLUTION_TYPE;
    switch (derivType)
    {
      case FIRST_DERIV:
        switch (solType)
        {
          case FLUIDMECH_VELOCITY:
            retSolType = FLUIDMECH_VELOCITY_DERIV_1;
            break;
          case FLUIDMECH_PRESSURE:
            retSolType = FLUIDMECH_PRESSURE_DERIV_1;
            break;
          case MECH_DISPLACEMENT:
            retSolType = MECH_DISPLACEMENT_DERIV_1;
            break;
          case ACOU_PRESSURE:
            retSolType = ACOU_PRESSURE_DERIV_1;
            break;
          case ACOU_POTENTIAL:
            retSolType = ACOU_POTENTIAL_DERIV_1;
            break;
          case BUBBLE_RADIUS:
            retSolType = BUBBLE_RADIUS_DERIV_1;
            break;
          default:
            EXCEPTION("not implemented for that deriv: " \
                << SolutionTypeEnum.ToString(solType));
        }
        break;
      case SECOND_DERIV:
        switch (solType)
        {
          case FLUIDMECH_VELOCITY:
            retSolType = FLUIDMECH_VELOCITY_DERIV_2;
            break;
          case FLUIDMECH_PRESSURE:
            retSolType = FLUIDMECH_PRESSURE_DERIV_2;
            break;
          case MECH_DISPLACEMENT:
            retSolType = MECH_DISPLACEMENT_DERIV_2;
            break;
          case ACOU_PRESSURE:
            retSolType = ACOU_PRESSURE_DERIV_2;
            break;
          case ACOU_POTENTIAL:
            retSolType = ACOU_POTENTIAL_DERIV_2;
            break;
          default:
            EXCEPTION("not implemented for that deriv: " \
                << SolutionTypeEnum.ToString(solType));
        }
        break;
      default:
        EXCEPTION("not implemented for that deriv: " \
            << SolutionTypeEnum.ToString(solType));
    }
    return retSolType;
  }
     


} // end of namespace
