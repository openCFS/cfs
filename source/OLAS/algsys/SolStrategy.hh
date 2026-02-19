#ifndef OLAS_SOL_STRATEGY_HH
#define OLAS_SOL_STRATEGY_HH

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "MatVec/BaseMatrix.hh"
#include "General/Enum.hh"


namespace CoupledField {


//! Base class for defining a solution strategy

//! This class encapsulates complex solution strategies, which
//! involves passing around non-local parameters within CFS and OLAS.
//! Initially this class gets generated within the AlgebraicSystem,
//! getting the pointer to a specific <solutionStrategy> sub-element.
//! 
//! The AlgebraicSystem uses this class for querying parameters of
//! SBM-Matrix blocks and also about the compound preconditioners.
//! In complex setups with several SBM-blocks and a compound
//! preconditioner, which acts differently on each block, the
//! SolStrategy class encapsulates the specific layout of the 
//! XML-tree and provides abstract information like number of SBM-blocks,
//! layout of each sub-block and preconditioner type for each SBMBlock.
//! 
//! The FeSpace uses this class to determine, how the equation numbering
//! scheme has to be set up (e.g. splitting of the equations in equations
//! of lower and higher order or for non-matching grids or spatial
//! splitting for domain-decomposition methods) and how the 
//! FeSpace registers its SBM- and submatrix blocks.
//!
//! The SolveStep class can use this object to query the number of
//! solution steps required (i.e. we start on a system of lower order,
//! solve it until convergence and continue on the full system).
class SolStrategy {
  
public:

  //! Enum for different strategy types
  typedef enum {NO_STRATEGY, STD_STRATEGY, TWO_LEVEL_STRATEGY} StrategyType;
  static Enum<StrategyType> strategyType;
  
  //! Constructor
  SolStrategy( PtrParamNode param );
  
  //! Destructor
  virtual ~SolStrategy();
  
  //! Factory methods
  static shared_ptr<SolStrategy> Generate( PtrParamNode );
  
  //! Return type of solution strategy
  StrategyType GetType() { 
    return type_;
  }
  
  // ========================================================================
  //  SOLUTION STEPS / SPLITTING LEVELS
  // ========================================================================
  //! Get number of solution steps
  virtual UInt GetNumSolSteps() = 0;
  
  //! Set current solution step (1-based)
  //! \param stepNum current solution step number (1-based)
  virtual void SetActSolStep(UInt stepNum) = 0;
  
  //! Get current solution step (1-based)
  virtual UInt GetActSolStep() = 0;
  
  //! Get number of SBM-blocks
  virtual UInt GetNumSBMBlocks() = 0;
  
  // ========================================================================
  //  OLAS-PARAMETER HANDLING
  // ========================================================================
  
  //! Return use of static condensation
  virtual bool UseStaticCondensation() = 0;
  
  //! Return use of penalty approach for Dirichlet handling
  virtual bool UseDirichletPenalty() = 0;

  //! Return calculation of conditionNumber
  virtual bool CalcConditionNumber() = 0;
  
  //! Return pointer to <setup> element
  virtual PtrParamNode GetSetupNode() = 0;
  
  //! Return storage type for matrix block
  virtual BaseMatrix::StorageType GetStorageType(UInt blockNum) = 0;
  
  //! Return <matrix> nodes
  virtual PtrParamNode GetMatrixNode(UInt level) = 0;
  
  //! Return pointer to <solver
  virtual std::string GetSolverId() = 0;
  
  //! Return pointer to <solver
  virtual std::string GetEigenSolverId() = 0;
  
  //! Return pointer to <precond>
  virtual std::string GetPrecondId() = 0;
  
  //! Return pointer to <nonLinear> element
  virtual PtrParamNode GetNonLinNode() = 0;
  
  //! Return pointer to <hysteresis> element
  virtual PtrParamNode GetHystNode() = 0;

  //! Return pointer to <timeStepping> element
  virtual PtrParamNode GetTimeSteppingNode() = 0;
  
  PtrParamNode GetParamNode() { return param_; }

  unsigned int GetEstimatedRowSize() { return estimated_row_size_; }  

  // ========================================================================
  //  MULTIHARMONIC SECTION
  // ========================================================================
  //! Set flag is multiharmonic analysis is used
  void SetMultHarm(bool isMultHarm){ isMultHarm_ = isMultHarm;}

  bool IsMultHarm(){ return isMultHarm_; }

  void SetMultHarm(const UInt& bF, const UInt& nN, const UInt& nM, const UInt& numFFT, bool fullSystem){
    baseFreq_ = bF;
    numHarmN_ = nN;
    numHarmM_ = nM;
    numFFT_ = numFFT;
    fullSystem_ = fullSystem;
  }

  UInt GetBaseFreq(){ return baseFreq_; }
  UInt GetNumHarmN(){ return numHarmN_; }
  UInt GetNumHarmM(){ return numHarmM_; }
  UInt GetNumFFT(){ return numFFT_; }
  bool IsFullSystem(){ return fullSystem_; }

  // ========================================================================
  //  2.5D HARMONIC SECTION
  // ========================================================================
  //! Set flag if 2.5d harmonic analysis is used
  void Set25dHarm(bool is25dHarm){ is25dHarm_ = is25dHarm;}
  bool Is25dHarm(){ return is25dHarm_; }

protected:
  
  //! Paramnode for strategy element
  PtrParamNode param_;
  
  //! Strategy type of this class
  StrategyType type_;
  
  //! Number of solution steps
  UInt numSolSteps_;
  
  //! Current solution step (0-based)
  UInt curSolStep_;
  
  //! Special matrix element for static condensation
  PtrParamNode statCondMatNode_;

  //! Flag if multiharmonic analysis is used
  bool isMultHarm_;

  //! Flag if 2.5d harmonic analysis is used
  bool is25dHarm_;

  //! Base frequency for multiharmonic excitation
  //! Also used for 2.5d harmonic analysis
  UInt baseFreq_;
  //! Number of harmonics for solution
  UInt numHarmN_;
  //! Number of harmonics for nonlinearity
  UInt numHarmM_;

  //! Number of considered time evaluation points for FFT and iFFT
  UInt numFFT_;

  //! Boolean, which tells us if we need to incorporate the zero harmonic
  bool fullSystem_;
 
  /** for BaseGraph, 100 in .xsd. If too small we have re-hasing. Check .info.xml for "graph" */
  unsigned int estimated_row_size_ = 100;

};


//! Standard solution strategy

//! This class represents the "standard" solution strategy, i.e. we have a
//! SBM-system with only one block and we use a standard solver / 
//! preconditioner combination with optional static condensation.
//! The FeSpace does not have to perform a specific numbering strategy 
//! (exception: static condensation for interior block) and the SolveStep
//! class just has to trigger once the assembly of the system and 
//! the solution of the system.
class SolStrategyStd : public SolStrategy {
public:

  //! Constructor
  SolStrategyStd( PtrParamNode node );

  //! Destructor
  ~SolStrategyStd();

  // ========================================================================
  //  SOLUTION STEPS
  // ========================================================================
  //! Get number of solution steps
  virtual UInt GetNumSolSteps();

  //! Set current solution step
  virtual void SetActSolStep(UInt stepNum);
  
  //! Get current solution step
  virtual UInt GetActSolStep();

  //! Return use of static condensation
  virtual bool UseStaticCondensation();

  //! Get number of SBM-blocks
  virtual UInt GetNumSBMBlocks();
   
  // ========================================================================
  //  OLAS-PARAMETER HANDLING
  // ========================================================================

  //! Return use of penalty approach for Dirichlet handling
  bool UseDirichletPenalty();
  
  //! Return calculation of conditionNumber
  bool CalcConditionNumber();
  
  //! Return pointer to <setup> element
  virtual PtrParamNode GetSetupNode();
  
  //! Return storage type for matrix block
  virtual BaseMatrix::StorageType GetStorageType(UInt blockNum);
  
  //! Return <matrix> node
  virtual PtrParamNode GetMatrixNode(UInt level);

  //! Return eigenSolverId
  virtual std::string GetEigenSolverId();
  
  //! Return solver id
  virtual std::string GetSolverId();

  //! Return pointer to <precond>
  virtual std::string GetPrecondId();

  //! Return pointer to <nonLinear> element
  virtual PtrParamNode GetNonLinNode();

  //! Return pointer to <hysteresis> element
  virtual PtrParamNode GetHystNode();

  //! Return pointer to <timeStepping> element
  virtual PtrParamNode GetTimeSteppingNode();

protected:
  
  //! <setup> element
  PtrParamNode setupNode_;
  
  //! <matrix> element
  PtrParamNode matrixNode_;

  //! <eigenSolver> element
  PtrParamNode esNode_;
  
  //! <solver> element
  PtrParamNode solverNode_;
  
  //! <precond> element
  PtrParamNode precondNode_;
  
  //! <nonLinear> element
  PtrParamNode nonlinNode_;
  
  //! <hysteresis> element
  PtrParamNode hystNode_;

  //! <timeStepping> element
  PtrParamNode tsNode_;
  
};

//! Two-Level Solution Strategy

//! This class represents the "standard" solution strategy, i.e. we have a
//! SBM-system with only one block and we use a standard solver / 
//! preconditioner combination with optional static condensation.
//! The FeSpace does not have to perform a specific numbering strategy 
//! (exception: static condensation for interior block) and the SolveStep
//! class just has to trigger once the assembly of the system and 
//! the solution of the system.
class SolStrategyTwoLevel : public SolStrategy {
public:

  //! Constructor
  SolStrategyTwoLevel( PtrParamNode node );

  //! Destructor
  ~SolStrategyTwoLevel();

  // ========================================================================
  //  SOLUTION STEPS
  // ========================================================================
  //! Get number of solution steps
  virtual UInt GetNumSolSteps();

  //! Set current solution step
  virtual void SetActSolStep(UInt stepNum );
  
  //! Get current solution step
  virtual UInt GetActSolStep();

  //! Return use of static condensation
  virtual bool UseStaticCondensation();

    //! Get number of SBM-blocks
  virtual UInt GetNumSBMBlocks();

  // ========================================================================
  //  OLAS-PARAMETER HANDLING
  // ========================================================================

  //! Return use of penalty approach for Dirichlet handling
  bool UseDirichletPenalty();
  
  //! Return calculation of conditionNumber
  bool CalcConditionNumber();
  
  //! Return pointer to <setup> element
  virtual PtrParamNode GetSetupNode();
  
  //! Return storage type for matrix block
  virtual BaseMatrix::StorageType GetStorageType(UInt blockNum);
  
  //! Return <matrix> node
  virtual PtrParamNode GetMatrixNode(UInt level);

  //! Return eigenSolverId
  virtual std::string GetEigenSolverId();
  
  //! Return solver id
  virtual std::string GetSolverId();

  //! Return pointer to <precond>
  virtual std::string GetPrecondId();

  //! Return pointer to <nonLinear> element
  virtual PtrParamNode GetNonLinNode();

  //! Return pointer to <hysteresis> element
  virtual PtrParamNode GetHystNode();

  //! Return pointer to <timeStepping> element
  virtual PtrParamNode GetTimeSteppingNode();

protected:
  
  //! <setup> node
  PtrParamNode setupNode_;

  //! Matrix nodes s per level
  ParamNodeList matrixNodes_;
  
  //! EigenSolver nodes per solution step
  ParamNodeList eigenSolverNodes_;
    
  //! Solver nodes per solution step
  ParamNodeList solverNodes_;
  
  //! Preconditioner nodes per solution step
  ParamNodeList precondNodes_;
  
  //! Nonlinear nodes per solution step
  ParamNodeList nonLinNodes_;
  
    //! Nonlinear nodes per solution step
  ParamNodeList hystNodes_;

  //! Timestepping nodes per solution step
  ParamNodeList timeStepNodes_;

};

} // end of namespace
#endif
