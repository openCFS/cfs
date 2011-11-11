#include <fstream>
#include <iostream>
#include <string>

#include "TimeStepping.hh"

#include "OLAS/algsys/AlgebraicSys.hh"

namespace CoupledField {

  TimeStepping::TimeStepping(AlgebraicSys * algebraicsystem )
  {
    algsys_  = algebraicsystem;
    rhsSize_ = 0 ;
    omitFirstPredictor_ = false;
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
            retSolType = MECH_VELOCITY;
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
          case SMOOTH_DISPLACEMENT:
            retSolType = SMOOTH_VELOCITY;
            break;
          default:
            EXCEPTION("not implemented for first derivative of: " \
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
            retSolType = MECH_ACCELERATION;
            break;
          case ACOU_PRESSURE:
            retSolType = ACOU_PRESSURE_DERIV_2;
            break;
          case ACOU_POTENTIAL:
            retSolType = ACOU_POTENTIAL_DERIV_2;
            break;
          default:
            EXCEPTION("not implemented for second derivative of: " \
                << SolutionTypeEnum.ToString(solType));
        }
        break;
      default:
        EXCEPTION("not implemented for that deriv: " \
            << SolutionTypeEnum.ToString(solType));
    }
    return retSolType;
  }

  bool
  TimeStepping::isDeriv(const SolutionType solType) const
  {
    switch (solType)
    {
      case MECH_VELOCITY:
        return true;
      case MECH_ACCELERATION:
        return true;
      case SMOOTH_VELOCITY:
        return true;
      case BUBBLE_RADIUS_DERIV_1:
        return true;
      case ACOU_PRESSURE_DERIV_1:
        return true;
      case ACOU_PRESSURE_DERIV_2:
        return true;
      case FLUIDMECH_VELOCITY_DERIV_1:
        return true;
      case FLUIDMECH_VELOCITY_DERIV_2:
        return true;
      case FLUIDMECH_PRESSURE_DERIV_1:
        return true;
      case FLUIDMECH_PRESSURE_DERIV_2:
        return true;
      case ACOU_POTENTIAL_DERIV_1:
        return true;
      default:
        break;
    }
    switch (solType)
    {
      case MECH_DISPLACEMENT:
        return false;
      case SMOOTH_DISPLACEMENT:
        return false;
      case BUBBLE_RADIUS:
        return false;
      case ACOU_PRESSURE:
        return false;
      case FLUIDMECH_VELOCITY:
        return false;
      case FLUIDMECH_PRESSURE:
        return false;
      case ACOU_POTENTIAL:
        return false;
      default:
        EXCEPTION("TimeStepping::isDeriv() is not implemented for this solution Type:\n" \
            << SolutionTypeEnum.ToString(solType))
    }
  }



} // end of namespace
