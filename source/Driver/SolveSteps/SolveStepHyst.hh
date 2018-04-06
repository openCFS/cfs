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
          Integer timeLevelRHSHyst, SBM_Vector& solutionForUpdate, SBM_Vector& returnVector);
    
    void SetupRHS(SBM_Vector& solutionForUpdate);
    void SetupRES(SBM_Vector& solutionForUpdate);
    
    void AssembleSystem();
    void AssembleSystem(Integer timeLevelMaterial, Integer timeLevelDeltaMat);
    
    void InitTimeStepping();
    void SolveStepTrans();
    void ReadNonLinDataHyst();
    
    virtual void StepTransNonLinHysteresis();
    
    //! does a line search and returns the optimal residual norm
    Double LineSearchHyst(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double& etaLineSearch, UInt evalVersion, UInt callingCnt,
                      bool trans=false, bool performLineSearch=true, bool forceReevaluation=false, bool debugOutput=false, bool reset=false,UInt allowedSteps=5);
        
  private:
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

    SBM_Vector RhsVec_;
    SBM_Vector ResVec_;
    
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
    UInt testInvesion_;

   
    //! Vectors used for NonLinHysteresis
    // current > current timestep and iteration
    SBM_Vector currentLinRhsVec_;
    //SBM_Vector currentRhsVec_;
    SBM_Vector currentResVec_;


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
    
  };


    
  
} // end of namespace

#endif

