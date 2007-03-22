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

namespace CoupledField {


  //! Class for assembling element/entities matrices and RHS vectors
  class Assemble {

  public:

    //! Constructor
    Assemble( BaseSystem* algsys, AnalysisType analysis, UInt maxTimeDerivOrder );

    //! Destructor
    ~Assemble();

    // ======================================================
    //  REGISTRATION METHODS
    // ======================================================    
    
    //! Add a bilinearform, wrapped in a BilinFormContext
    void AddBiLinearForm( BiLinFormContext* biLinContext );

    //! Add a linearform, wrapped in a LinearFormContext
    void AddLinearForm( LinearFormContext* linContext );

    //! Add right hand side load definitions
    void AddLoads( LoadList& loads );

    // ======================================================
    //  ASSEMBLING METHODS
    // ======================================================

    //! Assemble matrix graph of given pair of pdes
    void SetupMatrixGraph( PdeIdType pdeId1, PdeIdType pdeId2 );

    //! Trigger assembly of the matrices
    void AssembleMatrices();

    //! Trigger assembly of all linear right hand side terms
    void AssembleLinRHS( Double actTimeFreq );

    //! Trigger assenbly of all non-linear right hand side terms
    void AssembleNonLinRHS( Double actTimeFreq );


    // ======================================================
    //  MISCELLANEOUS METHODS
    // ======================================================

    //! Query if resulting matrix will be symmetric
    bool IsFEMatSymmetric( FEMatrixType matType = OLAS::SYSTEM );

    //! Returns true, if matrices have changed since last call of
    //! AssembleMatrices
    bool IsMatrixUpdated(){ return matrixUpdated_;}
    
    //! Print information about registered (bi)linearforms and general data
    void PrintInfo( std::ostream& out );

  protected:

    //! Assemb2le linearForms of right hand side
    void AssembleRHSLinForms( Double actTimeFreq, bool nonLin );

    //! Assemble noda load values of right hand side
    void AssembleRHSLoads( Double actTimeFreq);

    //! Transform real-valued element matrix to harmonic representation
    void Matrix2Harmonic( Vector<Double>& harmMat,
                          Matrix<Double>& origMat,
                          FEMatrixType matrixType,
                          DataType matDataType,
                          Double omega );
      
    //! Create map for mapping general FEMatrixtype to analysis-specific ones
    void CreateMatrixMap();

    //! Modify damping factor in harmonic case
    void AdjustDamping( BiLinFormContext& context );

    //! Insert matrix into algebraic system and adapt harmonic matrices
    void InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                       Matrix<Double>& elemMat, StdVector<Integer>& eqnVec1,
                       StdVector<Integer>& eqnVec2,
                       PdeIdType pdeId1, PdeIdType pdeId2 );

    //! Check which integrator is non-linear due to solution-dependent
    //! non-linearities or updated lagrangian formulation
    void CheckNonLinearities();

    // ======================================================
    //  DATA MEMBERS
    // ======================================================

    //! Pointer to algebraic system
    BaseSystem* algsys_;

    //! Analysistype
    AnalysisType analysisType_;

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
    //  DAMPING SPECIFIC DATA
    // ======================================================

    //! Current (adjusted) Rayleigh damping factor
    Double raylDampFactor_;
    
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
