//==============================================================================
/*!
 *       \file     CoefFunctionStabParams.cc
 *       \brief    Coefficient function which computes the stabilization
 *                 parameters for SUPG
 *
 *       \date     07/04/2013
 *       \author   Manfred Kaltenbacher
 */
//==============================================================================

#include "CoefFunctionStabParams.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField{

  // Definition of stabilization type mappings
  static EnumTuple stabTypeTuples[] = {
    EnumTuple(CoefFunctionStabParams::NO_TYPE, "no stabilization"),
    // Streamline Upwind Petrov Galerkin
    EnumTuple(CoefFunctionStabParams::SUPG,    "SUPG"),
    // Pressure Stabilized Petrov Galerkin
    EnumTuple(CoefFunctionStabParams::PSPG,    "PSPG"),
    // Least Squared Incompressibility Condition
    EnumTuple(CoefFunctionStabParams::LSIC,    "LSIC")
  };
  Enum<CoefFunctionStabParams::StabType> CoefFunctionStabParams::stabType =
    Enum<CoefFunctionStabParams::StabType>(
      "Stabilization types for perturbed flow",
      sizeof(stabTypeTuples) / sizeof(EnumTuple),
      stabTypeTuples);


  CoefFunctionStabParams::CoefFunctionStabParams(PtrCoefFct density,
                                                 PtrCoefFct viscosity,
                                                 StabType type) 
    : CoefFunction() 
  {
    PerformInitialization( density, viscosity, type );
  }
  
  CoefFunctionStabParams::CoefFunctionStabParams(PtrCoefFct density,
                                                 PtrCoefFct viscosity,
                                                 PtrCoefFct meanFlow,
                                                 BaseBOperator* opt,
                                                 shared_ptr<BaseFeFunction> feFnc,
                                                 StabType type) 
    : CoefFunction()
  {    
    PerformInitialization(density, viscosity, type );
    meanFlowCoef_ = meanFlow;
    IsSetMeanFlow_ = true;
    bOperator_ = opt;
    feFct_ = feFnc;    
  }

  void CoefFunctionStabParams::GetScalar(Double& scal,
                                         const LocPointMapped& lpm ) 
  {
    Double dens, vis;
    density_->GetScalar(dens,lpm);
    viscosity_->GetScalar(vis,lpm);
    
    Double vol = lpm.shapeMap->CalcVolume();
    
    
    UInt dim = lpm.jac.GetNumRows();
    Double edgeLength, edgeLengthVel;
    if ( dim == 2 ) 
    {
      edgeLength =  std::pow(vol, 1.0/2.0);
    }
    else
    {
      edgeLength =  std::pow(vol, 1.0/3.0);
    }
    
    Double velL2, stabValVec;
    
    edgeLengthVel = 0.0;
    velL2 = 0.0;
    if ( IsSetMeanFlow_ ) {
      Vector<Double> myVec;
      meanFlowCoef_->GetVector(myVec,lpm);
      velL2 = myVec.NormL2();
      
      if ( velL2 > 1e-10 ) {
        Matrix<Double> xiDx;
        BaseFE * ptFe = feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
        bOperator_->CalcOpMat(xiDx, lpm, ptFe);
        Double help;
        edgeLengthVel = 0.0;
        for ( UInt i=0; i<xiDx.GetNumCols(); i++ ) {
          help = 0.0;
          for ( UInt j=0; j<xiDx.GetNumRows(); j++ )
            help += myVec[j]*xiDx[j][i];
          edgeLengthVel += std::abs(help);
        }
        edgeLengthVel = 2.0 * myVec.NormL2() / edgeLengthVel ;
      }
    }
    
    if ( stabType_ == LSIC ) {
      Double stabVal = edgeLengthVel * velL2 * dens / (6.0 * vis);
      if ( stabVal < 1.0 )
        scal = velL2 * edgeLengthVel * stabVal / 2.0;
      else
        scal = velL2 * edgeLengthVel / 2.0;
    }
    else if ( stabType_ == SUPG ){
      Double stabValDiff = ( dens / (12.0*vis) ) * edgeLengthVel * edgeLengthVel;
      if ( velL2 < 1e-10 ) {
        scal = stabValDiff;
      }
      else {
        stabValVec    = edgeLengthVel / (2.0 * velL2);
        if ( stabValVec < stabValDiff )
          scal = stabValVec;
        else
          scal = stabValDiff;
      }
    }
    else if ( stabType_ == PSPG ){
      Double stabValDiff = ( dens / (12.0*vis) ) * edgeLength * edgeLength;
      if ( IsSetMeanFlow_ ) {
        if ( velL2 < 1e-10 )
          scal = stabValDiff;
        else {
          stabValVec = edgeLength / (2.0 * velL2);
          if ( stabValVec < stabValDiff )
            scal = stabValVec;
          else
            scal = stabValDiff;
        }
      }
      else
        scal = stabValDiff;
    }
    //	std::cout << "\n StabVal:" << scal << std::endl;
  }
  
  void CoefFunctionStabParams::GetTensor( Matrix<Double>& tensor, 
                                          const LocPointMapped& lpm ) {
    EXCEPTION( "Sabaradi nomoi nei!" );

  }
  
  void CoefFunctionStabParams::PerformInitialization(PtrCoefFct density,
                                                     PtrCoefFct viscosity,
                                                     StabType type) {
    
    dimType_ = SCALAR;
    isAnalytic_ = false;
    dependType_ = GENERAL;
    
    density_ = density;
    viscosity_ = viscosity;
    
    IsSetMeanFlow_ = false;
    stabType_ = type;
  }
  
  
}// end of namespace
