// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_ASSEMBLE_HH
#define CFS_ASSEMBLE_HH

#include <set>
#include <map>

#include "formsContext.hh"
#include "Domain/bcs.hh"
#include "Utils/mathParser/mathParser.hh"
#include "PDE/basePDE.hh"
#include "DataInOut/Scripting/scriptable.hh"

namespace CoupledField {


  // Forward class declarations
  class TimeFunc;
  class StdPDE;
  class InfoNode;

  //! Class for assembling element/entities matrices and RHS vectors
  class Assemble : public Scriptable {

  public:

    //! Constructor
    Assemble( BaseSystem* algsys, BasePDE::AnalysisType analysis, UInt maxTimeDerivOrder );

    //! Destructor
    ~Assemble();
    
    //! Explicitly set the algebraic system
    void SetAlgSys(BaseSystem * algsys);

    // ======================================================
    //  REGISTRATION METHODS
    // ======================================================

    //! Add a bilinearform, wrapped in a BilinFormContext
    void AddBiLinearForm( BiLinFormContext* biLinContext );

    //! Add a linearform, wrapped in a LinearFormContext
    void AddLinearForm( LinearFormContext* linContext );

    //! Add right hand side load definitions
    void AddLoads( LoadList& loads );
    
    //! Re-Set status of matrix re-assembly
    void ResetMatrixReassembly();

    // ======================================================
    //  ASSEMBLING METHODS
    // ======================================================

    //! Assemble matrix graph of given pair of pdes
    void SetupMatrixGraph( FeFctIdType pdeId1, FeFctIdType pdeId2 );

    //! Trigger assembly of the matrices
    void AssembleMatrices();

    void CalcMinMaxStrain();
    
    //! Trigger assembly of all linear right hand side terms
    void AssembleLinRHS();

    //! Trigger assenbly of all non-linear right hand side terms
    void AssembleNonLinRHS();

    //! Assemble nodal load values of right hand side
    void AssembleRHSLoads();

    // ======================================================
    //  MISCELLANEOUS METHODS
    // ======================================================

    //! Query if resulting matrix will be symmetric
    bool IsFEMatSymmetric( FEMatrixType matType = SYSTEM );

    //! Returns true, if matrices have changed since last call of
    //! AssembleMatrices
    bool IsMatrixUpdated(){ return matrixUpdated_;}

    /** Append info about registered (bi)linearforms */
    void ToInfo(InfoNode* in);

    /** <p>The PDEs don't know their own Integrators (the Element matrices K_{uu},
     *  ...) but when one wants to use it, we have to get it back from the
     * assemble class.</p>
     * <p>The query needs to define a unique form.</p>
     * @param regionId guess what!
     * @param pde1 this is the first pde
     * @param pde2 the second pde, note the order -> see debug file.
     * @param integ the integrator: linElastInt, MassInt, linElecInt, linPiezoCoupling
     * @return the defined context, never NULL
     * @exception error when nothing found or not unique specification */
    BiLinFormContext* GetBiLinForm(RegionIdType regionId, StdPDE* pde1, StdPDE* pde2, const std::string& integrator);

    /** Returns the load list for external modification */
    LoadList& GetLoads() { return loads_; }

    /** Overwrites the loads to implmented the adjoint solution for SIMP
     * mechanism optimization */
    void SetLoads(LoadList& new_loads) { loads_ = new_loads; }
    
    /** Overwrites the linearForms to implement the multi-load optimization */
    void SetLinForms(std::set<LinearFormContext*>& linForms) { linForms_ = linForms; }

    /** Returns the algebraic system
     * TODO check if really used */
    BaseSystem* GetAlgSys() { return algsys_; }

    /** Returns the bilinear forms list for Shape Optimization does need to loop these as assemble does */
    std::set<BiLinFormContext*>* GetBiLinForms() { return &biLinForms_; }
    

    /** Returns the linear forms list for external modification */
    std::set<LinearFormContext*>* GetLinForms() { return &linForms_; }

  protected:

    //! Assemb2le linearForms of right hand side
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

    //! Insert matrix into algebraic system and adapt harmonic matrices
    void InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                       Matrix<Double>& elemMat, StdVector<Integer>& eqnVec1,
                       StdVector<Integer>& eqnVec2,
                       FeFctIdType pdeId1, FeFctIdType pdeId2 );

    //! Insert complex matrix into algebraic system and adapt harmonic matrices
    void InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                       Matrix<Complex>& elemMat, StdVector<Integer>& eqnVec1,
                       StdVector<Integer>& eqnVec2,
                       FeFctIdType pdeId1, FeFctIdType pdeId2 );

    //! Check which integrator is non-linear due to solution-dependent
    //! non-linearities or updated lagrangian formulation
    void CheckNonLinearities();

    // ======================================================
    // SCRIPTING SECTION
    // ======================================================
    //@{
    //! \name Scripting Methods

    //! Register scriptable functions
    virtual void RegisterFunctions();
    
    //! Map with element numbers for which element matrix should be written
    std::set<UInt> printElemNums_;
    
    //! Provide list of elements for which matrix / vector is printed
    void Wrap_AddPrintElemNum( );
    
    //! Set all element matrices / vectors to be printed
    void Wrap_PrintAllElems();
    //@}
    
    // ======================================================
    //  DATA MEMBERS
    // ======================================================

    //! Pointer to algebraic system
    BaseSystem* algsys_;

    //! Analysistype
    BasePDE::AnalysisType analysisType_;

    //! Flag indicating if system was already assembled
    bool isFirstTime_;

    //! Map from general FEMatrixType to analysis specific one
    std::map<FEMatrixType,FEMatrixType> matrixMap_;

    //! List of bilinear integrator contexts
    std::set<BiLinFormContext*> biLinForms_;

    //! List of linear integrator contexts
    std::set<LinearFormContext*> linForms_;

    //! Map with flags if FE matrix has to be reassembled
    std::map<FEMatrixType, bool> matReassemble_;

    // ======================================================
    //  BOUNDARIES AND LOADS
    // ======================================================

    //! List of right hand side nodal values
    LoadList loads_;

    // ======================================================
    //  MISCELLANEOUS DATA
    // ======================================================

    //! flag indicating if matrices have changed since
    //! last call of AssembleMatrices
    bool matrixUpdated_;

    //! Maximum order of partial derivatives w.r.t. time
    UInt maxTimeDerivOrder_;

    //! Handle for MathParser object
    MathParser::HandleType mHandle_;
  };
}
#endif
