#ifndef OLAS_SOL_STRATEGY_HH
#define OLAS_SOL_STRATEGY_HH

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "MatVec/BaseMatrix.hh"


namespace CoupledField {


//! Base class for defining a solution strategy

//! This class encapsulates complex solution strategies, which
//! involves passing around non-local parameters within CFS and OLAS.
//! Initially this class gets generated within the AlgebraicSystem,
//! getting the pointer to a specific <solutionStrategy> sub-element.
//! 
//! The AlgebraicSystem uses this class for querying parameters of
//! SBM-Matrix blocks and also about the compound preconiditioners.
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

  //! Typedef for splitting strategy
  typedef enum {NO_SPLIT_STRAT, ORDER_TWO_LEVEL} SplittingType;
  Enum<SplittingType> SplittingTypeEnum;
  
  //! Constructor
  SolStrategy( PtrParamNode param );
  
  //! Destructor
  virtual ~SolStrategy();
  
  //! Factory methods
  static shared_ptr<SolStrategy> Generate( PtrParamNode );
  
  // ========================================================================
  //  SOLUTION STEPS
  // ========================================================================
  //! Get number of solution steps
  virtual UInt GetNumSolSteps() = 0;
  
  //! Set current solution step
  virtual void SetActSolStep() = 0;
  
  //! Return splitting strategy
  virtual SplittingType GetSplitStrategy() = 0;
  
  // ========================================================================
  //  OLAS-PARAMEER HANDLING
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
  virtual std::string GetPrecondId(UInt level) = 0;
  
  //! Return pointer to <exportLinSys> element
  virtual PtrParamNode GetExportLinSysNode() = 0;
  
  //! Return pointer to <nonLinear> element
  virtual PtrParamNode GetNonLinNode() = 0;
  
  //! Return pointer to <timeStepping> element
  virtual PtrParamNode GetTimeSteppingNode() = 0;
  
protected:
  
  //! Paramnode for strategy element
  PtrParamNode param_;
  
  //! Number of solution steps
  UInt numSolSteps_;
};


//! Standard solution strategy

//! This class represents the "stanard" solution strategy, i.e. we have a
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
  virtual void SetActSolStep();

  //! Return use of static condensation
  virtual bool UseStaticCondensation();

  //! Return splitting strategy
  virtual SplittingType GetSplitStrategy();

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
  virtual std::string GetPrecondId(UInt level);

  //! Return pointer to <exportLinSys> element
  virtual PtrParamNode GetExportLinSysNode();

  //! Return pointer to <nonLinear> element
  virtual PtrParamNode GetNonLinNode();

  //! Return pointer to <timeStepping> element
  virtual PtrParamNode GetTimeSteppingNode();

protected:
  
  // <setup> element
  PtrParamNode setupNode_;
  
  // <exportLinSys> element
  PtrParamNode exportNode_;
  
  // <matrix> element
  PtrParamNode matrixNode_;

  // <eigenSolver> element
  PtrParamNode esNode_;
  
  // <solver> element
  PtrParamNode solverNode_;
  
  // <precond> element
  PtrParamNode precondNode_;
  
  // <nonLinear> element
  PtrParamNode nonlinNode_;
  
  // <timeStepping> element
  PtrParamNode tsNode_;
  
};

////! Solution strategy for two level solver
//class SolStrategyTwoLevel : public SolStrategy {
//  
//public:
//  
//protected:
//  
//};

} // end of namespace
#endif
