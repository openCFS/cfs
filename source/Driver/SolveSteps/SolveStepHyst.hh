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
#include "Domain/Domain.hh"
#include <deque>

namespace CoupledField
{
  // forward class declarations
  class StdPDE;
  struct ResultInfo;
  struct StoppingCriterion;
  struct LinesearchParameter;
  class SingleDriver;
  class IDBC_Handler;
  class BaseIDBC_Handler;
  class FeSpace;
  class SolStrategy;
  class MathParser;

  /*
   * Local enums
   */

  /*
   * Local structs > TODO: make use of them or delete
   */
  struct flagsForHystOperator {
    /*
     * Struct collecting flags that get passed to coefFunction hyst
     */
    /*
     * timeLevelMaterial
     * -2: material tensors (eps, e, h, ...) are independent of hysteresis
     * -1: material tensors are to be evaluated with hyst output from previous ts (n)
     *  0: material tensors are to be evaluated with current hyst output (k)
     *  1: material tensors are to be evaluated with hyst output from previous it (k-1)
     */
    int timeLevelMaterial;

    /*
     * scalingOfMaterialTensor (NEW, only for Fixpoint cases)
     *  0: no scaling; use standard material tensors (eps, mu, ...)
     *  1: apply scaling for B-version of fixpoint iteration (TODO: has to be implemented in CoefFunctionHyst->getTensor)
     *  2: apply scaling for H-version of fixpoint iteration
     */
    int scalingOfMaterialTensor;

    /*
     * timeLevelDeltaMat
     * -2: no delta matrix formulation (deactivates ALL deltaMatrices
     * -1: delta matrix formulation with delta towards last ts (n)
     *			e.g. deltaP = P_k - P_n
     *  0: deactivate specific deltaMat
     *  1: delta matrix formulation with delta towards previous it (k-1)
     *			e.g. deltaP = P_k - P_k-1
     *  2: return last computed value of delta matrix
     *  3: compute Finite-Difference approximation of Jacobian instead of delta matrix
     */
    int timeLevelDeltaMatPol;
    int timeLevelDeltaMatStrain;
    int timeLevelDeltaMatCoupling;

    /*
     * timeLevelRHSHyst
     * -2: hyst operator will not be added to rhs
     * -1: output of hyst operator from last ts (n) will be added to rhs
     *  0: current output of hyst operator (k) will be added to rhs
     *  1: output of hyst operator from last it (k-1) will be added to rhs
     */
    int timeLevelRHSPol;
    int timeLevelRHSStrain;
    int timeLevelRHSCoupling;

    /*
     * timeLevelRESIDUALHyst
     * -2: hyst operator will not be added to residual
     * -1: output of hyst operator from last ts (n) will be added to residual
     *  0: current output of hyst operator (k) will be added to residual
     *  1: output of hyst operator from last it (k-1) will be added to residual
     */
    int timeLevelRESIDUALPol;
    int timeLevelRESIDUALStrain;
    int timeLevelRESIDUALCoupling;

    /*
     * inversionApproach (NEW, not required in electrostatics)
     *  0: invert hyst operator using a local inversion scheme (as specified in mat.xml)
     *      > for B-iteration, i.e. polarization/hystvalue = function(B)
     *  1: evaluate hyst operator with previously obtained value of H; no local inversion
     *      > for H-iteration, i.e. polarization/hystvalue = fucntion(H)
     */
    int inversionApproach;

  };

  // simple struct for ls-results
  struct LSResults {
    std::string name;
    linesearchFlag lsFlag;
    bool success;
    Double alphaFound;
    Double resultingRes;
    Double totalEvalTime;
    UInt numResidualEvals;
  };

  // simple struct for solution approach
  struct SolApproachParameter {
    solveFlag usedSolveFlag;
    linesearchFlag usedLSFlag;
    linesearchFlag failbackLSFlag;
    Double minEtaLS;

    // if true: compute update towards previous time step result
    // if false: compute update towards zero
    bool incremental;
  };

  //! Derived class for step-wise solving of StdPDEs
  class SolveStepHyst : public StdSolveStep {

  public:

    SolveStepHyst(StdPDE& apde);
    virtual ~SolveStepHyst();

    void InitTimeStepping();
    void SolveStepTrans();

  private:

    void Init_Step1_ReadNonLinDataHyst();
    void Init_Step1a_ReadSolutionMethod(PtrParamNode& methodNode);
    void Init_Step1b_ReadLinesearches(PtrParamNode& linesearchNode);
    LinesearchParameter Init_Step1b_ReadSingleLinesearch(PtrParamNode& singleLinesearchNode,std::string nameTag);
    void Init_Step1c_ReadStoppingCriteria(PtrParamNode& stoppingNode,bool isFailback);
    StoppingCriterion Init_Step1c_ReadSingleStoppingCriterion(PtrParamNode& singleStoppingNode,std::string nameTag);
    void Init_Step2_InitGlobalArrays();

    void Solve_Step1_PrepareSolveStep();
    bool Solve_Step2_NonLinIteration(UInt resetCnt);
    void Solve_Step3_FinishSolveStep();

    SBM_Vector ComputeSolutionUpdate_BGM(SBM_Vector& z0);
    SBM_Vector ComputeSolutionUpdate_ALGSYS(solveFlag flagForSolving, bool forceStepTowardsPrevTS,
        bool forceFullStep, bool forceJacobian, bool forceReuseOfJacobian);//, int iterationCounterLog);

    bool AssembleSystem(matrixFlag requestedAssemblyFlag, SBM_Vector& solutionForMatrixAssembly, bool forceReuseOfJacobian);

    void BuildSystem(bool prohibitStorage);

    bool LineSearchHyst_GetInitialSearchInterval(matrixFlag flag_for_matrix_assembly, vectorFlag flag_for_vector_assembly,
      SBM_Vector& currentSol, SBM_Vector& searchDirection, Double alpha, Double h, Double t, UInt kMax, Double& a, Double& b, Double& m,
      UInt iterationCounter);

    Double LineSearchHyst(matrixFlag flag_for_matrix_assembly, vectorFlag flag_for_vector_assembly,
      SBM_Vector& currentSol, SBM_Vector& solIncrement_n_k_local,UInt iterationCounter);

    Double SetupRESVector_Linesearch(matrixFlag flag_for_matrix_assembly, vectorFlag flag_for_vector_assembly, SBM_Vector& targetVector,
      SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForNonLinRHSAssembly, SBM_Vector& oldSolution, SBM_Vector& solutionUpdate, Double currentEta);

    Double SetupRESorRHS(matrixFlag flag_for_matrix_assembly, vectorFlag flag_for_vector_assembly,
      SBM_Vector& solutionForUpdate, SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForRHSAssembly,
      SBM_Vector& returnVector, bool comingFromLinesearch);

    /*
     * Returnvalue:
     * -1 = no success
     * 0 = success
     * 1 = failback
     *
     * loggingLevel:
     * 0: just name of criteria and if it was passed or not
     * 1: name of criteraion, target value, obtained value, passed or not
     */
    int checkStoppingCriteria(UInt iterationCounter, std::ostringstream& logger, UInt loggingLevel);

    void Helper_GetParamsFromSolveFlag(solveFlag flag_for_systemSolution, matrixFlag& flag_for_systemMatrix, matrixFlag& flag_for_resEvalMatrix,
      vectorFlag& flag_for_rhsVector, vectorFlag& flag_for_resVector, bool& estimatorUseful);

    void Helper_GetHystOperatorParamsFromMatrixFlag(matrixFlag flag_for_assembly, int& timeLevelMaterial, int& FPMaterialTensor,
          int& FPRHSCorrectionTensor, int& timeLevelDeltaMatPol, int& timeLevelDeltaMatStrain, int& timeLevelDeltaMatCoupling, int& inversionApproach);

    void Helper_GetHystOperatorParamsFromVectorFlag(vectorFlag flag_for_res_or_rhs, int& timeLevelPol, int& timeLevelSrain,
      int& timeLevelCoupling);

    void Helper_ComputeKeffTimesSolution(SBM_Vector& solutionForUpdate, SBM_Vector& returnVector,bool useBackup, matrixFlag flag_for_matrix_assembly);
    void Helper_CreateSystemBackup(bool backupSOL, bool backupRHS, bool backupHystValues, bool backupSYSTEM, FEMatrixType backupSLOT);
    void Helper_LoadSystemBackup(bool loadSOL, bool loadRHS, bool loadHystValues, bool loadSYSTEM, FEMatrixType loadSLOT);
    Vector<UInt> Helper_GetStorageIndices();

    void Helper_MatrixTimesSBMVector(Matrix<Double>& matIn, SBM_Vector& vecIn, SBM_Vector& vecOut);
    void Helper_GetSizeInfoFromSBMVector(SBM_Vector& vecIn, UInt& totalSize,
    UInt& numSubvectors, Vector<UInt>& sizeSubvectors, Vector<UInt>& posOffsets);
    void Helper_DyadicProductFromSBMVectors(Matrix<Double>& dyadicProduct, SBM_Vector& vec1,
    SBM_Vector& vec2,Double scaling = 1.0,bool addIdentity = false);

    void Output_WriteHystIterToInfoXML(const std::string& pdeName,
          const std::string& nonLinSolveMethod,
          const std::string& linesearchMethod,
          const UInt coupledIterStep,
          const UInt solStep,
          const UInt iterationCounter,
          const std::string& convergenceLogMSG);

    void Output_WriteHystIterToInfoXML(const std::string& pdeName,
    const std::string& nonLinSolveMethod,
    const std::string& linesearchMethod,
    const UInt coupledIterStep,
    const UInt solStep,
    const UInt iterationCounter,
    const UInt iterationCounterTotal,
    const Double residualErr, const Double residualErrRel, bool useRelativeRelErr,
    const Double incrementalErr, const Double incrementalErrRel, bool useRelativeIncErr,
    double etaLineSearch, bool reset, bool failback);

    /*
     * Enum to string
     */
    std::string linesearchFlagToString(linesearchFlag flag);

    std::string solveFlagToString(solveFlag flag);

    std::string matrixFlagToString(matrixFlag flag, int verbose = 0);

    std::string vectorFlagToString(vectorFlag flag, int verbose = 0);

    std::string FEMatrixTypeToString(FEMatrixType matType);

    /*
     * NEW class variables (previously local variables)
     */
    SBM_Vector solInc;
    SBM_Vector bestSol;
    SBM_Vector resInc;
    SBM_Vector resVecInitial_;

    std::map<FEMatrixType,Integer> matrices;
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;

    bool flagsInitialized_;
    bool trans_;

    bool LinRHSRequiresAssembly_;
    bool PredictorCorrectorRequiresAssembly_;

    Integer timeLevelMaterialCurrent_;
    Integer timeLevelRHSPolCurrent_;
    Integer timeLevelRHSStrainCurrent_;
    Integer timeLevelRHSCouplingCurrent_;

    Integer timeLevelDeltaMatPolCurrent_;
    Integer timeLevelDeltaMatStrainCurrent_;
    Integer timeLevelDeltaMatCouplingCurrent_;

    //    Integer timeLevelRHSHystCurrent_;
    //    Integer timeLevelDeltaMatCurrent_;

    SBM_Vector resVec_;

    SBM_Vector stageRHSUpdate_;
    SBM_Vector solVecLastTS_;
    SBM_Vector solVecLastIT_;
    SBM_Vector solVecForUpdate_;

    SBM_Vector lastUpdateVecForNonLinRHS_;


    /*
     * for non-monotonic-Armijo
     */
    std::deque<Double> resNormHistoryArmijo_;
    UInt nonmonotonicArmijo_numEntries_;

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

    // new failback tolerance that might also be ok if real criterion fails
//    Double failBackCrit_;

//    UInt minNumberItTillReset_;
//    Double minResDecreaseTillReset_;
//    std::string nonLinMethodAfterReset_;
//    std::string lineSearchAfterReset_;
//    bool towardsLastIterationReset_;
//    bool deltaMatCouplingReset_;
//    bool deltaMatStrainReset_;

//    bool towardsLastIteration_;
    bool deltaMatCoupling_;
    bool deltaMatStrain_;

    bool useRelativeNormForRes_;
    bool useRelativeNormForInc_;
    bool reuseDeltaMat_;
    bool disableDeltaMat_;

    Double residualErrAbsPrev_;
    bool assemblyApproximatesJacobian;
    bool systemApproximatesJacobian;

    /*
     * -2: linear system
     * -1: DeltaMat Last TS
     * 1: DeltaMat Last Iteration
     * 3: FixPoint system
     */
    int currentAssemblyFlag;
    int currentBuildFlag;

    bool linearSystemBackuped_;

    /*
     * =====================
     * ==== NEW SECTION ====
     * =====================
     */
    // dummy vector that will be used throughout the solve step to
    // setup local (and global) SBM_Vectors
    SBM_Vector zeroVec_;

    /*
     * NEW nomenclature
     *
     * n = current timestep
     * k = current iteration
     * kPrev = previou iteration
     *
     */
    bool resVec_n_0_Norm_;
    SBM_Vector resVec_n_0_;
    SBM_Vector resVec_n_kPrev_;
    SBM_Vector resVec_n_k_;
    SBM_Vector resIncrement_n_k_;

    SBM_Vector solVec_n_0_;
    SBM_Vector solVec_n_kPrev_;
    SBM_Vector solVec_n_k_;
    SBM_Vector solIncrement_n_k_;
    SBM_Vector bestSol_;

    /*
     * Backup vectors for algsys
     * > system itself is backuped in algsys class
     */
    SBM_Vector solVecBAK_;
    SBM_Vector rhsVecBAK_;

    /*
     * linesearch
     */
    std::map< linesearchFlag, LinesearchParameter > specifiedLinesearchesMap_;
    std::set< linesearchFlag > specifiedLinesearches_;
//    std::set< LinesearchParameter > specifiedLinesearchParameter_;
    UInt LSSelectionCriterion_;
    Double etaLinsearch_n_k_;

    /*
     * normal and failback stopping criterion
     * 1. stoppingCriteria:
     *    iteration successful if ALL stopping criteria are satisfied
     * 2. failbackCriterai:
     *    iteration is stopped (with failback) if ANY failback criterion is satisfied
     */
    std::set< StoppingCriterion > stoppingCriteria_;
    std::set< StoppingCriterion > failbackCriteria_;
    Double systemSize_; // currently length of solution vector; should be restricted to hysteretic / nonlinear parts

    /*
     * for BGM
     */
    UInt historyLength_; // current number of stored residual/solution increments
    UInt hisotryLengthMax_; // max numboer of stored residual/solution incrmenets
    Vector<Double> etaLinsearch_history_;
    SBM_Vector* solIncrements_history_;
    SBM_Vector* resIncrements_history_;

    /*
     * flag indicating if there are further non-linearities
     * to be considered in the system matrix
     * if false: residual computation can be done with the initial linear system matrix
     * if true: a new nonlinear but hysteresis free system matrix has to be setup for residual computation
     */
    bool nonHystereticNonLinearitiesPresent_;
    /*
     * flag indicating if material tensors (mu,eps,e,d,...) depend on hysteresis
     * this is e.g. the case in magnetostriction and piezos where the coupling tensors are rotated and scaled
     * according to the current value of polarization
     * if true: a new nonlinear but hysteresis free system matrix has to be setup for residual computation
     */
    bool materialTensorsHystDepenedent_;

    SBM_Vector linRhsVec_n_0_;
    SBM_Vector nonLinRhsVec_n_k_;
    SBM_Vector predictorCorrectorUpdate_n_0_;
    bool linRhsVec_n_0_precomputed_;
    bool predictorCorrectorUpdate_n_0__precomputed_;

    /*
     * in pure hysteresis case, the residual is computed
     * via the linear system matrix and hysteresis is treated
     * via rhs vector; during linesearch, we can thus precompute
     * some values
     * --> CAN BE MADE LOCAL
     */
    SBM_Vector linSysMatrix_n_0_TIMES_solVec_n_kPrev_;
    bool linSysMatrix_n_0_TIMES_solVec_n_kPrev_precomputed_;

    SBM_Vector linSysMatrix_n_0_TIMES_solIncrement_n_k_;
    bool linSysMatrix_n_0_TIMES_solIncrement_n_k_precomputed_;

    /*
     * for each different matrix assembly type, keep track of
     * the last build counter and the corresponding solution vector;
     * idea: check if re-assembly is required or not
     */
    std::map< matrixFlag, SBM_Vector > latestUpdateVectors_;
    std::map< matrixFlag, int > latestBuildCounter_;


    /*
     * In case of additional non-linearities (e.g. due to B-H-saturation curves)
     * we need to build the nonlinear system matrix for residual computation,
     * but without approximating the Jacobian (during residual computation, hysteresis
     * is treated via non-linear rhs vector)
     */
    //    SBM_Matrix& hystFreeSysMatrix_n_k_;
    SBM_Vector lastUpdateVecForHystFreeSystem_;
    int hystFreeSysMatrix_n_k_lastBuildCounter_;

    /*
     * lastBuildCounter_
     *  -1:  no build
     * >=0:  iteration at which matrix was setup for the last time
     */
    //    SBM_Matrix& FDJacobianSysMatrix_n_k_;
    SBM_Vector lastUpdateVecForFDJacobianSystem_;
    int FDJacobianSysMatrix_n_k_lastBuildCounter_;

    //    SBM_Matrix& DeltaMatJacobianSysMatrix_n_k_;
    SBM_Vector lastUpdateVecForDeltaMatJacobianSystem_;
    int DeltaMatJacobianSysMatrix_n_k_lastBuildCounter_;

    // local enums
    //    enum nonLinMethodFlag { NOMETHODSET = -1, FIXPOINT_B_STEPPING, FIXPOINT_H_STEPPING,
    //    QUASI_NEWTON_DELTAMAT, QUASI_NEWTON_FDJACOBIAN, QUASI_NEWTON_BGM };




    /*
     * -1: system has not been assembled/built yet (during this timestep)
     *  0: hystFree/Fixpoint
     *  1: FDJacobian
     *  2: DeltaMatJacobian
     *  3: ReusedOldApproximation
     */
    matrixFlag currentAssemblyFlag_;
    matrixFlag requestedBuildFlag_;
    matrixFlag currentBuildFlag_;

    vectorFlag currentRESandRHSFlag_;

    solveFlag solutionApproach_; // solution approach as defined in input.xml
    solveFlag solutionApproachFailback_; // solution approach for failback as defined in input.xml
    solveFlag currentSolutionApproach_; // currently used solution approach
    UInt initialNumberFPSteps_; // for Quasi-Newton methods it might be reasonable to have some initial fixpoint steps
                                // (B-stepping) to come into the contractive region; usually the line-search methods
                                // will bring Newton into a contractive region BUT the initial iteration in each time
                                // time step has to be a full step (unfortunately!) because of the IDBC!
    bool ReEvalFPMaterialTensors_; // if we switch from globalFP to localFP, we have to reevaluate the localFP slope!

    UInt iterationsTillReset_; // readded option to reset solution vector to 0; might help here and there to get a better
                               // starting point for locally converging methods; does not work so great however;
                               // per default, this value is maxIter + 1 such that it is not triggered

    /*
     * solver parameter
     */
    Double FPLocalC_;
    Double FPGlobalC_;
    UInt numIterationsTillUpdateOfJacobian_;
    bool towardsPreviousTimestep_;

    linesearchFlag currentLinesearchApproach_;

    fixpointFlag FPApproach_;
    fixpointFlag currentFPApproach_;
    std::string currentNonLinMethod_;

    bool forceRHSevaluation_;
    bool forceReassembly_;
    bool systemAssemblyRequired_;
    bool systemBuildRequired_;
    bool systemLoadRequired_;

    UInt timestepCnt_;
    UInt totalIterationCnt_;
    UInt maxIterationPerTSCnt_;
    UInt iterationCntCurrentTS_;
    UInt iterationOfLatestReset_;


    std::string nonLinBackup_;
    std::string lineSearchBackup_;
//    bool towardsLastIterationBackup_;
    bool deltaMatCouplingBackup_;
    bool deltaMatStrainBackup_;

    /*
     * actual solution related flags from input.xml
     * > new ones
     */
    bool initStep1Done_;
    UInt checkDirection_;
    bool performEstimatorStep_;
    Double minEta_;
    UInt allowedNumResets_;

    /*
     * Debug flags from input.xml
     */
    bool DEBUG_testJacobianApproximations_;
    UInt DEBUG_compareLinesearches_;
    bool DEBUG_measureLSPerformance_;
  };

} // end of namespace

#endif
