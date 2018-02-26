#ifndef FILE_SOLVESTEPHYST
#define FILE_SOLVESTEPHYST

#include <map>
#include <fstream>

#include "StdSolveStep.hh"

#include "Driver/Assemble.hh"
#include "MatVec/SBM_Vector.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/BaseMaterial.hh"
#include "DataInOut/ResultHandler.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"

namespace CoupledField
{
  // forward class declarations
  class StdPDE;
  class WriteResults;
  struct ResultInfo;
  class SingleDriver;
  class IDBC_Handler;
  class BaseIDBC_Handler;
  class FeSpace;
  class SolStrategy;

  //  class Domain;
  
  //! Derived class for step-wise solving of StdPDEs
  class SolveStepHyst : public StdSolveStep {

  public:

    //! Constructor
    SolveStepHyst(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepHyst();

    void SetupRESorRHS(Integer timeLevelMaterial, Integer timeLevelDeltaMat, 
          Integer timeLevelRHSHyst, SBM_Vector& solutionForUpdate, SBM_Vector& returnVector);
    
    void SetupRHS();
    void SetupRES();
    void AssembleSystem(UInt evalCase);
    void AssembleSystem(Integer timeLevelMaterial, Integer timeLevelDeltaMat);
    void InitTimeStepping();
    void SolveStepTrans();
    void ReadNonLinDataHyst();
    
    virtual void StepTransNonLinHysteresisNew();
    
    //----------------------- HYSTERESIS -------------------------------------
    //! solves for one nonlinear transient step
    //! consideres hystreresis nonlinearities in direct coupled PDEs
    virtual void StepTransNonLinHysteresis();
    virtual void StepTransNonLinHysteresisTotal();

    void SetLastItOrLastTSSBMVectors(bool lastTS);
    /*!
     * Helper funciton for setting up the equation system during
     *            StepTransNonLinHysteresis()
     * Background: During the solve step, the matrices and the rhs have to be
     *  assembled multiple times during linesearch, for the calculation of the
     *  residual error and of course to get a system to be solved;
     *  for simplification, encapsulate that sequence of function calls
     *  in a separate function
     */
    /*
     * for residual computation we need a slightly different version -> see .cc file
     */
    virtual void CalcResidualAndConfigSystemForHysteresis(SBM_Vector& oldSolution,SBM_Vector& solIncrement, SBM_Vector& stageSol, Double usedEta,
                                                          UInt stage, UInt callingCnt, UInt evalVersion, bool trans, bool forceReevaluation,
                                                          bool skipReassembly, bool debugOutput, bool reset);

    //! does a line search and returns the optimal residual norm
    Double LineSearchHyst(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double& etaLineSearch, UInt evalVersion, UInt callingCnt,
                      bool trans=false, bool performLineSearch=true, bool forceReevaluation=false, bool debugOutput=false, bool reset=false,UInt allowedSteps=5);
    

    //! returns the hysteresis operator
    Hysteresis * GetHystOperator(UInt idx) {
      return hyst_[idx];
    };
    
    
  private:
    bool flagsInitialized_;
    bool trans_;
    
    bool systemRequiresReassembly_;
    bool LinRHSRequiresReassembly_;
    bool NonLinRHSRequiresReassembly_;
    bool PredictorCorrectorRequiresReassembly_;
    
    Integer timeLevelRHSHystCurrent_;
    Integer timeLevelMaterialCurrent_;
    Integer timeLevelDeltaMatCurrent_;
    
    SBM_Vector LinRhsVec_;
    SBM_Vector NonLinRhsVec_;
    SBM_Vector predictorCorrectorUpdate_;

    SBM_Vector RhsVec_;
    SBM_Vector ResVec_;
    
    SBM_Vector stageRHSUpdate_;
    SBM_Vector solVecLastTS_;
            
    UInt evalCase_;
    
  };


    
  
} // end of namespace

#endif

