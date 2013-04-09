// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_ASSEMBLE_HH
#define CFS_ASSEMBLE_HH

#include <set>
#include <map>

#include "FormsContexts.hh"
#include "Domain/BCs.hh"
#include "Utils/mathParser/mathParser.hh"
#include "PDE/BasePDE.hh"

namespace CoupledField {


  // Forward class declarations
  class TimeFunc;
  class Timer;
  class StdPDE;
  class AdjointParameters;
  class AlgebraicSys;

  //! Class for assembling element/entities matrices and RHS vectors
  class Assemble {

  public:

    //! Constructor
    Assemble( AlgebraicSys* algsys, BasePDE::AnalysisType analysis,
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
    void AssembleMatrices();

    //! Trigger assembly of all linear right hand side terms
    void AssembleLinRHS(AdjointParameters* adjointParams = NULL);

    //! Trigger assenbly of all non-linear right hand side terms
    void AssembleNonLinRHS(AdjointParameters* adjointParams = NULL);

    // ======================================================
    //  MISCELLANEOUS METHODS
    // ======================================================
    
    //! Tell Assemble to Reassemble all the matrices the next time (done on init and from Optimization)
    void SetAllReassemble(){ CheckNonLinearities(true); }

    //! Query if matrix related to BiLinearForm is symmetric
    bool IsFEMatSymmetric(  BiLinFormContext* ctx );
    
    //! Query if matrix related to BiLinearForm is complex
    bool IsFEMatComplex(  BiLinFormContext* ctx );

    //! Returns true, if matrices have changed since last call of
    //! AssembleMatrices
    bool IsMatrixUpdated(){ return matrixUpdated_;}

    /** Append info about registered (bi)linearforms */
    void ToInfo(PtrParamNode in);

    /** <p>The PDEs don't know their own Integrators (the Element matrices K_{uu},
     *  ...) but when one wants to use it, we have to get it back from the
     * assemble class.</p>
     * <p>The query needs to define a unique form.</p>
     * @param regionId guess what!
     * @param pde1 this is the first pde
     * @param pde2 the second pde, note the order -> see debug file.
     * @param integrator: linElastInt, MassInt, linElecInt, linPiezoCoupling
     * @param silent exception or NULL if nothing found
     * @return the defined context, never NULL
     * @exception error when nothing found or not unique specification */
    BiLinFormContext* GetBiLinForm(RegionIdType regionId, StdPDE* pde1, 
                                   StdPDE* pde2, const std::string& integrator, bool silent = false);

    /** @see GetBiLinForm() */
    LinearForm* GetLinearForm(RegionIdType regionId, StdPDE* pde,  const std::string& integrator, bool silent = false);

    /** Overwrites the linearForms to implement the multi-load optimization */
    void SetLinForms(StdVector<LinearFormContext*>* linForms) { linForms_ = linForms; }

    /** Returns the algebraic system
     * TODO check if really used */
    AlgebraicSys* GetAlgSys() { return algsys_; }

//    /** Returns the bilinear forms list for Shape Optimization does need to loop these as assemble does */
//    StdVector<BiLinFormContext*>& GetBiLinForms() { return *biLinForms_; }
    

    /** Returns the linear forms list for external modification */
    StdVector<LinearFormContext*>& GetLinForms() { return *linForms_; }

  protected:

    //! Assemble matrices without static condensation
    void AssembleMatrices_Std();
    
    //! Assemble matrices with satic condensation
    void AssembleMatrices_Cond();

    
    //! Assemble linearForms of right hand side
    void AssembleRHSLinForms(bool nonLin );

    //! Transform real-valued element matrix to harmonic representation
    void Matrix2Harmonic( Matrix<Complex>& harmMat,
                          Matrix<Double>& origMat,
                          FEMatrixType matrixType,
                          Global::ComplexPart matDataType,
                          Double omega );

    //! Transform complex-valued element matrix to harmonic representation
    void Matrix2Harmonic( Matrix<Complex>& harmMat,
                          Matrix<Complex>& origMat,
                          FEMatrixType matrixType,
                          Global::ComplexPart matDataType,
                          Double omega );

    //! Create map for mapping general FEMatrixtype to analysis-specific ones
    void CreateMatrixMap();
    
    //! Perform re-mapping of equation numbers in case of eqn permutation
    void ReMapEquations( StdVector<Integer>&  eqns, FeFctIdType& fctId );
    
    //! Perform re-mapping of functionId
    void ReMapFctId( FeFctIdType& fctId );

    //! Insert matrix into algebraic system and adapt harmonic matrices
    void InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                       Matrix<Double>& elemMat, StdVector<Integer>& eqnVec1,
                       StdVector<Integer>& eqnVec2,
                       FeFctIdType fctId1, FeFctIdType fctId2 );

    //! Insert complex matrix into algebraic system and adapt harmonic matrices
    void InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                       Matrix<Complex>& elemMat, StdVector<Integer>& eqnVec1,
                       StdVector<Integer>& eqnVec2,
                       FeFctIdType fctId1, FeFctIdType fctId2 );

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
    
    //! List of linear integrator contexts
    StdVector<LinearFormContext*>* linForms_;

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
  };
}
#endif
