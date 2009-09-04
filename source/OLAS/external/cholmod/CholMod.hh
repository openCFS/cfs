// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CHOLMOD_HH
#define CHOLMOD_HH

#include <def_expl_templ_inst.hh>

#include "General/environment.hh"
#include "OLAS/solver/basesolver.hh"

#include "General/Enum.hh"

// include the original cholmod header
#include "cholmod.h" 

namespace CoupledField 
{
  class BaseMatrix;  
  class BaseVector;
  class BasePrecond;
  class Flags;
  
  template<typename T>
  class CholMod : public BaseIterativeSolver 
  {
  public:
    CholMod(ParamNode* param, InfoNode* olasInfo, BaseMatrix::EntryType type);

    ~CholMod();

    /** Every call does the complete factorization
     * @param analysis_id shall be the current info/analysis/progress/step entry and contain an "analysis_id" element */
    void Setup(BaseMatrix &sysmat, InfoNode* analysis_id);
     
    /** solve using a pre computed factorization 
     * @param sysmat shall be the one Setup() is called with
     * @param precond ignored
     * @param analysis_id @see Setup() */
    void Solve( const BaseMatrix &sysmat, const BasePrecond &precond,
                const BaseVector &rhs, BaseVector &sol, InfoNode* analysis_id);

    /** Return what solver is used (here CholMod) */
    SolverType GetSolverType() { return CHOLMOD; }

    typedef enum { ANALYSIS_SIMPLICIAL = CHOLMOD_SIMPLICIAL, ANALYSIS_AUTO = CHOLMOD_AUTO, ANALYSIS_SUPERNODAL = CHOLMOD_SUPERNODAL } Analysis;

  private:

    /** Set up the matrix mat 
     * @param base_mat input */
    void SetMatrix(const BaseMatrix &base_mat);

    /** Initializes the parameter with default settings by CholMod, conditionally
     * overwrite them with the xml settings and print out the complete stuff. */ 
    void InitParameters(); 
  
    /** Calls Cholmod analyze */
    void CholModAnalyze();
    
    /** Calls Factorize */
    void CholModFactorize();
    
    /** Calls Default */
    void CholModDefaultParams();
    
    /** Calls CholMod free on factorization */
    void CholModFreeFactorization();
    
    /** Calls Solve */
    void CholModSolve();
    
    /** Test the return status (in common_) and throw an exception if != CHOLMOD_OK
     * @param where text occuring in exception to identify where the error occured */
    void TestException(std::string where);
    
    /** This method sets up the "enums", it fills them with the string representations. */
    void SetEnums();

    /** are we complex */
    bool isComplex_;
    
    /** CholMod session with all parameters */
    cholmod_common common_;
    
    /** CholMod matrix */
    cholmod_sparse* mat_;
    
    /** Stored CholMod Factorization */
    cholmod_factor* fact_;
    
    Enum<Analysis> analysis;
    
  };

} // end of namespace

#endif
