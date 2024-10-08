// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_ASSEMBLE_HH
#define CFS_ASSEMBLE_HH

#include <set>
#include <map>

#include "FormsContexts.hh"
#include "Domain/BCs.hh"
#include "PDE/BasePDE.hh"

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace CoupledField {


  // Forward class declarations
  class TimeFunc;
  class Timer;
  class StdPDE;
  class AlgebraicSys;
  class MathParser;

  //! Class for assembling element/entities matrices and RHS vectors
  class Assemble {

  public:

    //! Constructor
    Assemble( AlgebraicSys* algsys, 
              BasePDE::AnalysisType analysis,
              MathParser* mp,
              PtrParamNode infoNode); 

    //! Destructor
    virtual ~Assemble();
    
    //! Explicitly set the algebraic system
    void SetAlgSys(AlgebraicSys * algsys);

    // ======================================================
    //  REGISTRATION METHODS
    // ======================================================

    //! Add a bilinearform, wrapped in a BilinFormContext
    void AddBiLinearForm( BiLinFormContext* biLinContext );

    //! Add a linearform, wrapped in a LinearFormContext
    void AddLinearForm( LinearFormContext* linContext );

    //! Re-Set status of matrix re-assembly
    void ResetMatrixReassembly();
    
    //! Set an additional custom mapping for equation numbers
    
    //! This method registers an additional map for equations numbers.
    //! This can be used, e.g. to handle bilinearforms, which are defined
    //! on just a sub-set of the original domain.
    void SetEqnCustomMap( const std::map<Integer, Integer>& eqnMap,
                          const std::map<FeFctIdType, FeFctIdType>& fctIdMap );

    // ======================================================
    //  ASSEMBLING METHODS
    // ======================================================

    //! Assemble matrix graph of given pair of pdes
    void SetupMatrixGraph( FeFctIdType fctId1, FeFctIdType fctId2 );

    //! Trigger assembly of the matrices
    void AssembleMatrices(bool isNewtonPart=false);
    
    //! Assemble matrices with static condensation for transient simulations
    void AssembleMatrices_CondTrans(bool isNewtonPart,UInt currentStage,
                                    std::map<FeFctIdType,
                                    std::map<FEMatrixType,Double> > timeStepFactors);

    //! Assemble matrices for multiharmonic analysis
    void AssembleMatrices_MultHarm(Integer harmonic, UInt N, UInt M,
     const std::map<RegionIdType, StdVector<NonLinType> >& regionNonLinTypes,
     const StdVector<Double>& multHarmFreqVec = StdVector<Double>());

    //! Initialize matrices for multiharmonic analysis
    //! Usually this is done in the Assemble method but
    //! for the multiharmonic analysis, we need to do it externally
    void InitMultHarm();

    //! Trigger assembly of all linear right hand side terms
    void AssembleLinRHS();

    //! Trigger assenbly of all non-linear right hand side terms
    void AssembleNonLinRHS();

    void PostAssemble();

    // ======================================================
    //  MISCELLANEOUS METHODS
    // ======================================================
    
    //! Tell Assemble to Reassemble all the matrices the next time (done on init and from Optimization)
    void SetAllReassemble(){ CheckNonLinearities(true); }

    void DisableMatrixAssembly();
    //! Query if matrix related to BiLinearForm is symmetric
    bool IsFEMatSymmetric(  BiLinFormContext* ctx );
    
    //! Query if matrix related to BiLinearForm is complex
    bool IsFEMatComplex(  BiLinFormContext* ctx );

    //! Returns true, if matrices have changed since last call of
    //! AssembleMatrices
    bool IsMatrixUpdated(){ return matrixUpdated_;}

    //! Return if any RHS integrator depends on the solution
    bool IsRhsSolDependent();

    /** Append info about registered (bi)linearforms */
    void ToInfo(PtrParamNode in);

    //Sets a flag to skip element assembly for std matrix case
    void SkipElemAssembly();

    //! Stops the timer in a multiharmonic analysis. It was started in InitMultHarm
    void TimerStop();

    /** search for an integrator.
     * @param integrator name where we do only a startswith, hence "curlCurlIntegrator" returns also a "curlCurlIntegrator-NL"
     * @param pde2 the second pde, note the order -> see debug file.
     * @param pde1/pde2 this is the first and second pde. If NULL not compared.
     * @param silent if false no NULL can be returned
     * @return the form is GetIntegrator() of the context. NULL only for silent true
     * @exception if not silent and nothing found */
    BiLinFormContext* GetBiLinForm(const std::string& integrator, RegionIdType regionId, SinglePDE* pde1 = NULL, SinglePDE* pde2 = NULL, bool silent = false);

    LinearFormContext* GetLinForm(const std::string& integrator, RegionIdType regionId, StdPDE* pde, bool silent);

    bool HasBiLinForm(const std::string& integrator, RegionIdType regionId, SinglePDE* pde1 = NULL, SinglePDE* pde2 = NULL) {
      return GetBiLinForm(integrator, regionId, pde1, pde2, true) != NULL;
    }

    /** Returns the algebraic system
     * TODO check if really used */
    AlgebraicSys* GetAlgSys() { return algsys_; }

    /** Returns the bilinear forms list for Shape Optimization does need to loop these as assemble does */
    std::set<BiLinFormContext*>& GetBiLinForms() { return allBiLinForms_; }

    /** Returns the linear forms list for external modification */
    StdVector<LinearFormContext*>& GetLinForms(bool take_ownership = false)
    {
      if(take_ownership)
        lin_forms_given_ = true; // we won't delete it
      return linForms_;
    }

    /** Do we use the region? */
    bool UseRegion(RegionIdType reg);

    //! Perform re-mapping of equation numbers in case of eqn permutation
    void ReMapEquations( StdVector<Integer>&  eqns, FeFctIdType& fctId );

  protected:

    //! Assemble matrices without static condensation
    void AssembleMatrices_Std(bool isNewtonPart=false);
    
    //! Assemble matrices with static condensation
    void AssembleMatrices_Cond(bool isNewtonPart=false);

    //! Assemble linearForms of right hand side
    void AssembleRHSLinForms(bool nonLin );

    //! Transform real-valued matrix to complex-valued
    void Matrix2Complex(Matrix<Complex>& harmMat, Matrix<Double>& origMat);

    //! Transform real-valued element matrix to harmonic representation
    void Matrix2Harmonic( Matrix<Complex>& harmMat,
                          Matrix<Double>& origMat,
                          FEMatrixType matrixType,
                          Global::ComplexPart matDataType,
                          Double omega );

    /**  Transform complex-valued element matrix to harmonic representation */
    void Matrix2Harmonic( Matrix<Complex>& harmMat,
                          Matrix<Complex>& origMat,
                          FEMatrixType matrixType,
                          Global::ComplexPart matDataType,
                          Double omega );

    //! Create map for mapping general FEMatrixtype to analysis-specific ones
    void CreateMatrixMap();
    
    //! Perform re-mapping of functionId
    void ReMapFctId( FeFctIdType& fctId );

    //! Insert matrix into algebraic system and adapt harmonic matrices
    void InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                       Matrix<Double>& elemMat, StdVector<Integer>& eqnVec1,
                       StdVector<Integer>& eqnVec2,
                       FeFctIdType fctId1, FeFctIdType fctId2,
                       bool preventStaticCondensation = false,
                       const StdVector<UInt>& sbmIndices = StdVector<UInt>(),
                       const Double& f = 0,
                       bool isMultHarmDiag = false);

    //! Insert complex matrix into algebraic system and adapt harmonic matrices
    void InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                       Matrix<Complex>& elemMat, StdVector<Integer>& eqnVec1,
                       StdVector<Integer>& eqnVec2,
                       FeFctIdType fctId1, FeFctIdType fctId2,
                       bool preventStaticCondensation = false,
                       const StdVector<UInt>& sbmIndices = StdVector<UInt>(),
                       const Double& f = 0,
                       bool isMultHarmDiag = false);

    //! Check which integrator is non-linear due to solution-dependent
    //! non-linearities or updated lagrangian formulation
    //! if setall is set, set all to be reassembled, called from SetAllReassemble (on Init or from Optimization)
    void CheckNonLinearities(bool setall = false);

    // ======================================================
    //  DATA MEMBERS
    // ======================================================

    //! Pointer to algebraic system
    AlgebraicSys* algsys_;

    //! Analysistype
    BasePDE::AnalysisType analysisType_;
    
    //! Pointer to math parser instance
    MathParser* mp_;
    
    //! Handle for expression
    unsigned int mHandle_;

    //! Flag indicating if system was already assembled
    bool isFirstTime_;

    //! Map from general FEMatrixType to analysis specific one
    std::map<FEMatrixType,FEMatrixType> matrixMap_;

    //! List of bilinear integrator contexts
    
    //! Associate pair of entity lists with bilinearform
    typedef std::map<
        std::pair<shared_ptr<EntityList>,shared_ptr<EntityList> >, 
        StdVector<BiLinFormContext*> > BiLinContextListType;  
     BiLinContextListType biLinForms_;

    //! Set containing all bilinear integrator contexts
    std::set<BiLinFormContext*> allBiLinForms_;
    
    /** List of linear integrator contexts. They are generated in the PDEs deletet here.
     * When we do muliload optimization Excitations gains ownership and linForms_ is manipulated.
     * @see Excitytion::form */
    StdVector<LinearFormContext*> linForms_;

    //! Set containing all linear integrator contexts
    std::set<LinearFormContext*> allLinForms_;

    /** when set, the destructor won't delete linForms_ (but Excitation will do it) */
    bool lin_forms_given_;

    //! Map with flags if FE matrix has to be reassembled
    std::map<FEMatrixType, bool> matReassemble_;

    //! Map with additional permutation for equations
    std::map<Integer, Integer> customEqnMap_;
    
    //! Map with additional permutation for function Ids
    std::map<FeFctIdType, FeFctIdType> customFctIdMap_;
    
    // ======================================================
    //  MISCELLANEOUS DATA
    // ======================================================
    
    //! Info node for output of list of (bi)linearforms
    PtrParamNode info_;
    
    //! flag indicating if matrices have changed since
    //! last call of AssembleMatrices
    bool matrixUpdated_;

    /** The object is within a ParamNode and deleted there! */
    boost::shared_ptr<Timer> timer_;
    
    //! Flag, if progress bar should be printed
    bool printProgressBar_;

    //flag for skipping element assembly when using petsc geometric multigrid
    bool skipElemAssembly_;

  };
}
#endif
