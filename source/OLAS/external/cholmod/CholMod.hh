#ifndef CHOLMOD_HH
#define CHOLMOD_HH

#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"

#include "General/Enum.hh"

// include the original cholmod header
#include "cholmod.h" 

namespace CoupledField 
{
  class BaseMatrix;  
  class BaseVector;
  class Flags;
  
  template<typename T>
  class CholMod : public BaseIterativeSolver 
  {
  public:
    CholMod(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

    ~CholMod();

    /** Every call does the complete factorization */
    void Setup(BaseMatrix &sysmat);

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented for Cholmod solver. GetRidOfZeros for NCIs will not work.");};
     
    /** solve using a pre computed factorization 
     * @param sysmat shall be the one Setup() is called with */
    void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol);

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
