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
          Integer timeLevelRHSHyst, SBM_Vector& solutionForUpdate, SBM_Vector& solutionForMatrixAssembly, 
          SBM_Vector& solutionForRHSAssembly, SBM_Vector& returnVector);
    
    void SetupRHS(SBM_Vector& solutionForUpdate, SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForRHSAssembly, UInt iterationCounter);
    void SetupRES(SBM_Vector& solutionForUpdate, SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForRHSAssembly, UInt iterationCounter);
    void ComputeKeffTimesSolution(SBM_Vector& solution, SBM_Vector& returnVector);
    
    void AssembleSystem(SBM_Vector& solutionForAssembly, UInt iterationCounter);
    void AssembleSystem(Integer timeLevelMaterial, Integer timeLevelDeltaMat, SBM_Vector& solutionForAssembly);
    
    void InitTimeStepping();
    void SolveStepTrans();
    void ReadNonLinDataHyst();
    
    virtual void StepTransNonLinHysteresis();
    
    void WriteHystIterToInfoXML(const std::string& pdeName,
          const std::string& nonLinSolveMethod,
          const UInt coupledIterStep,
          const UInt solStep,
          const UInt iterationCounter,
          const UInt iterationCounterTotal,
          const Double residualErr, 
          const Double incrementalErr, double etaLineSearch);
    
    void WriteHystIterToInfoXML(const std::string& pdeName,
          const std::string& nonLinSolveMethod,
          const UInt solStep,
          const UInt iterationCounter,
          const UInt iterationCounterTotal,
          const Double residualErr, 
          const Double incrementalErr, double etaLineSearch);
    
    //! does a line search and returns the optimal residual norm
    Double LineSearchHyst(SBM_Vector& currentSol, SBM_Vector& solIncrement, UInt iterationCounter);
        
  private:
    // new flag indicating if material tensors (mu,eps,e,d,...) depend on hysteresis or not
    // if not and deltaMat == false > system matrix needs no reassembly!
    bool materialTensorsHystDepenedent_;
    bool flagsInitialized_;
    bool trans_;
    
    bool LinRHSRequiresAssembly_;
    bool PredictorCorrectorRequiresAssembly_;
    
    Integer timeLevelRHSHystCurrent_;
    Integer timeLevelMaterialCurrent_;
    Integer timeLevelDeltaMatCurrent_;
    
    SBM_Vector LinRhsVec_;
    SBM_Vector NonLinRhsVec_;
    SBM_Vector predictorCorrectorUpdate_;

    SBM_Vector resVec_;
    
    SBM_Vector stageRHSUpdate_;
    SBM_Vector solVecLastTS_;
    SBM_Vector solVecLastIT_;
    SBM_Vector solVecForUpdate_; 
    
    SBM_Vector lastUpdateVecForNonLinRHS_; 
    SBM_Vector lastUpdateVecForMatrices_; 
    
    /*
     * NEW: describe behavior via new nonlin parameter from xml file
     */
    UInt evalDepth_;
    UInt measurePerformance_;
    UInt testInversion_;
    PtrParamNode testNode_;

    // + solVec which is defined above

    // oldTS > values after last iteration of previous TS
    // to be stored after iteration suceeded
    //SBM_Vector oldTSLinRhsVec_;
    //SBM_Vector oldTSNonLinRhsVec_;
    //SBM_Vector oldTSSolVec_;

    // oldIt > values of the current TS but from previous It
    // during first iterartion of a new TS, these vectors contain the values
    // after the last iteration of the previous TS (similar as oldTS...)
    //SBM_Vector oldItNonLinRhsVec_;
    //SBM_Vector oldItResVec_;

    //! Additional flags and parameter for hyst
    bool forceReevaluation_;
    bool debugOutput_;
    
    // new failback tolerance that might also be ok if real criterion fails
    Double failBackCrit_;
    
    bool allowReset_;
    UInt minNumberItTillReset_;
    Double minResDecreaseTillReset_;
    std::string nonLinMethodAfterReset_;
    std::string lineSearchAfterReset_;
    
  };


    
  
} // end of namespace

#endif

