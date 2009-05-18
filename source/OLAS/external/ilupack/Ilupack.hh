// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ILUPACK_HH
#define ILUPACK_HH

#include <def_expl_templ_inst.hh>

#include "General/environment.hh"
#include "OLAS/solver/basesolver.hh"

#include "General/Enum.hh"

// include the original ilupack header
extern "C"
{
   #include <ilupack.h> 
}

namespace CoupledField 
{
  class BaseMatrix;  
  class BaseVector;
  class BasePrecond;
  class Flags;

  template<typename T>
  class Ilupack : public BaseIterativeSolver 
  {
  public:
    Ilupack(ParamNode* param, InfoNode* olasInfo, BaseMatrix::EntryType type);

    ~Ilupack();

    void Setup( BaseMatrix &sysmat);

     
    /** To satisfy the compiler */
    void Solve( const BaseMatrix &sysmat, const BasePrecond &precond,
                const BaseVector &rhs, BaseVector &sol);


    /** This is not Ilupack intern but Olas stuff. It means that this is
     * a Ilupack, the actual solver called by the lib is in Ilupack::Solver */
    SolverType GetSolverType() { return ILUPACK_SOLVER; }

    /** This are the names of the matrix types in ILUPACK style */
    typedef enum { GNL, SYM, PD, HER } Matrix;

    /** There are combinations with matching from Pardiso and MC64, most symmetric
     * is based on matching but a matching is not mandatory. */
    typedef enum { MUMPS, MC64, MWM } Matching;

        
  private:

    /** Set up the ilupack matrix mat_ 
     * @param base_mat input */
    void SetMatrix(const BaseMatrix &base_mat);

    /** Determine the matrix type for ilupack. This can be set optionally in the xml file.
     * If it is not set, then it is SYM or GNL but neither SPD, HPD nor HER */
    void DetermineMatrixType(BaseMatrix &sysMat);

    /** Initializes the parameter with default settings by ilupack, conditionally
     * overwrite them with the xml settings and print out the complete stuff. 
     * Note that the complex and real versions can be savely
     * casted against each other as its only the pointers. output to param
     * and the ilupack matrix mat_ has to be set in advance by SetMatrix() */ 
    void InitParameters(); 
  
    /** Calls the AMGinit() variant of for the system from the lib */
    void IlupackAMGInit();
    
    /** Calls the AMGFactor() variant of for the system from the lib */
    int IlupackAMGFactor();

    int IlupackAMGSolver(T* sol_ptr, T* rhs_ptr);

    /** Relases the memory allocated for ilupack by calling its delete method. 
     * Does nothing if mat_.a is NULL. */ 
    void IlupackAMGDelete();
    
    /** This method sets up the "enums", it fills them with the string representations. */
    void SetEnums();

     
    /** If the factorization didn't result with 0, this method encrypts the error code. */
    std::string PrecondError(int result, const DILUPACKparam& param, int level);

    /** If the ilupack solver didn't result with 0, this method encrypts the error code. */
    std::string SolverError(int result);     

    /** Helper to Set ilupack char* from ParamNode */
    char* GetString(ParamNode* node);
    
    
    
    /** Helper method for InitParameters. Takes the ilupack value, prints it to the InfoNode (via
     * param_name).
     * If there is a parameter with param_name in the xml file, it is printed.
     * @param param_name simple xpath chain via slash */
    void CheckParameter(InfoNode* out, double* ilupack_val, const char* param_name);
    void CheckParameter(InfoNode* out, char** ilupack_val, const char* param_name);
    void CheckParameter(InfoNode* out, int* ilupack_val, const char* param_name);  
      
    /** The ilupack matching */
    Enum<Matching> matching;
    
    /** The ilupack matrix types */
    Enum<Matrix> matrix; 
    
    
    /** This is the ilupack signature of a real permutation function.
     * ILUPACK needs one for the inital, regular and final step */
    typedef int (*RealPermFunc) (Dmat, double*, double*, int*, int*, int*, DILUPACKparam*);
    
    typedef int (*ComplexPermFunc) (Zmat, doublecomplex *,doublecomplex *, int*, int*, int*, ZILUPACKparam*);    
    
                            
    /** Dump the solver stus, extracts nice parameters and prints them */
    void DumpSolverStatus(std::ostream& out);

    /** Dump the parameter setting */
    void DumpParameters(std::ostream& out);

    /** Dump the ilupack flags */ 
    void DumpFlags(int flag, std::ostream& out);

    /** Read the flag settings as provided by the XML file and
     * apply it on the given flag.
     * @param flag the in/out parameter */ 
    void ApplyFlagSettings(int& flag);
    
    /** Calculates the fill in factor of the preconditioner */
    void CalcFillIn(InfoNode* out);
    
    /** killme */
    void Ilupack_symmessages();
    void Ilupack_symprintperformance();
    void Ilupack_symfinalres(T* sol_ptr, T* rhs_ptr);
    
    /** This is the internal matrix type as optinally read from XML */
    Matrix matrix_; 
    
    /** are we complex */
    bool isComplex_;
    
    /** This is the ilupack matrix structure, it is set in SetMatrix(). The complex cases is only a cast of
     * the double version */    
    Dmat mat_;
    
    /** A pointer to mat_ */
    Dmat* dmat_ptr_;
    /** Also a pointer to (Dmat) dmat_ as complex is same-same here */
    Zmat* zmat_ptr_;
    
    /** Here ilupack stores all the paramters, initialized in InitParameters(). The complex case is just a cast */    
    DILUPACKparam param;

    /** A pointer to param */
    DILUPACKparam* dparam_ptr_;
    /** Also a pointer to (DILUPACKparam) param as complex is same-same here */
    ZILUPACKparam* zparam_ptr_;
    
    
    /** Here ilupack stores the preconditioner. */
    DAMGlevelmat precond_;

    /** A pointer to precond_ */
    DAMGlevelmat* dprecond_ptr_;
    /** Also a pointer to (DILUPACKparam) param as complex is same-same here */
    ZAMGlevelmat* zprecond_ptr_;
    
    /** This counts how many rhs are solved while the last precond call */
    int precondAge_;
    
    /** We can do setup multiple times, so we cout the sutup runs */
    int setupRuns_;
    
    /** Ilupack has ILUPACK_N[I/F]PAR (40) parmeters in *ILUPACKparam.[i/f]par. The entries are not defined
     * but from magic indices in the example the content can be extracted. Check this for ilupack updates!!<br>
     * MAX_S_ROW_ENTRIES = limit maximum number of entries in each row of S<br>
     * MAX_L_COL_U_ROW_ENTRIES = limit maximum number of entries in each column of L, row of U
     * MAX_ITERATIONS = maximum number of iteration steps
     * STOP_CRITERIA = stopping criteria (StoppingCriteria)
     * RESTART_LENGTH = number of restart length for GMRES,FOM,...
     * SOLVER_PARAM = init solver -> param.ipar[20]=0; values < 0 seem to be the error code 
     * PRECOND_PARAM = sides of PrecondSettings
     * ARNDOLDI_SOLVER_PARAM -> Arnoldi-type solvers: if (SOLVER>=7) -> switch (param.ipar[31]) 
     * */
    enum { MAX_L_COL_U_ROW_ENTRIES = 3, SOLVER_TYPE = 5,  MAX_S_ROW_ENTRIES = 9, SOLVER_PARAM = 20, PRECOND_PARAM = 21, STOP_CRITERIA = 22, 
           RESTART_LENGTH = 24, MAX_ITERATIONS = 25, ITERATION_STEPS = 26, ARNDOLDI_SOLVER_PARAM = 31} IntegerIndex; 
    
    /** See ipar_Index 
     * SCHUR_DROP_TOL = drop tolerance for the approximate Schur complement<br>
     * REL_ERROR_TOL = relative error tolerance -> param.fpar[20]=restol; - to set <br>
     * ABS_ERROR_TOL = absolute error tolerance -> param.fpar[21]=0;      - to set <br>
     *   */
    enum { SCHUR_DROP_TOL = 8, REL_ERROR_TOL = 20, ABS_ERROR_TOL = 21, INITIAL_RESIDUAL_NORM = 22,  
           TARGET_RESIDUAL_NORM = 23, CURRENT_RESIDUAL_NORM = 25, CONVERGENCE_RATE = 26} RealIndex;

    /** ZERO_INPUT_VECTOR = due to zero input vector <br>
     *  VEC_CONTAINS_ABNORMAL_NUMBERS = since input vector contains abnormal numbers<br>
     *  LINEAR_COMBINATION =  since input vector is a linear combination of others<br>
     *  NULL_RANK_TRINANGUAL = since triangular system in GMRES/FOM/etc. has null rank<br> */ 
    enum { ZERO_INPUT_VECTOR = -1, VEC_CONTAINS_ABNORMAL_NUMBERS = -2, LINEAR_COMBINATION = -3, 
           NULL_RANK_TRINANGUAL = -4 } ArnoldiSolverParam;
    
    /** stopping criteria, energy norm ipar[2]=3;
     * stopping criteria, backward error : param.ipar[22]=3;<br> */
    enum {ENERGY_NORM = 3 } StoppingCriteria;
    
    enum { NO_PRECOND = 0, LEFT = 1, RIGHT = 2, BOTH_SIDES = 3} PrecondParam;
    
    /** 0 -> everything is fine; -1 -> too many iterations; -2 -> not enough work space provided; -3 -> not enough work space, algorithm breaks down */ 
    enum { TOO_MANY_ITERATIONS = -1, NOT_ENOUGH_WORK_SPACE = -2,  ALGORITHM_BREAKS_DOWN = -3, INIT = 0 } SolverParam;


    /** The flags for parm.flag. The parameter "flag" is bitwise
     * modified by flag|=DROP_INVERSE to set and flag&=~DROP_INVERSE to
     * turn off inverse-based dropping. There is a FL_ before the constants to
     * differentiate them from the original ilupack defines. <br>
     * FL_DROP_INVERSE = switch to indicate inverse-based dropping. It is used in AMGINIT
     *               AMGGETPARAMS and AMGSETPARAMS. In AMGINIT, DROP_INVERSE is set by default <br>
     * FL_NO_SHIFT = switch for not shifting away zero pivots. This switch is used in ILUC
     *            which does not have pivoting prevent small diagonal entries.<br>
     * FL_TISMENETSKY_SC = switch for using Tismenetsky update.<br>
     * FL_REPEAT_FACT = switch for repeated ILU <br>
     * FL_IMPROVED_ESTIMATE =  switch for enhanced estimate for the norm of the inverses<br>
     * FL_DIAGONAL_COMPENSATION = switch for using diagonal compensation<br>
     * FL_COARSE_REDUCE = switch for reducing the partial factorization to the non-coarse part<br>
     * FL_FINAL_PIVOTING = switch for using a different pivoting strategy, if the regular reordering
     *                  fails and before we switch to ILUTP<br>
     * FL_ENSURE_SPD = enforce the positve definite property<br>
     * FL_SIMPLE_SC = switch for the most simple Schur complement update*/
    typedef enum { FL_DROP_INVERSE = 1, FL_NO_SHIFT = 2, FL_TISMENETSKY_SC = 4, FL_REPEAT_FACT = 8, FL_IMPROVED_ESTIMATE = 16, 
          FL_DIAGONAL_COMPENSATION = 32, FL_COARSE_REDUCE = 64, FL_FINAL_PIVOTING = 128, FL_ENSURE_SPD = 256, FL_SIMPLE_SC = 512, 
          FL_PREPROCESS_INITIAL_SYSTEM = 1024, FL_PREPROCESS_SUBSYSTEMS = 2048, FL_MULTI_PILUC = 4096, FL_RE_FACTOR = 8192, 
          FL_AGGRESSIVE_DROPPING = 16384, FL_DISCARD_MATRIX = 32768, FL_SYMMETRIC_STRUCTURE = 65536} Flags;    

    /** The Ilupack flags */
    Enum<Flags> flags;
          
          
   //typedef enum { PREPROCESS_INITIAL_SYSTEM = 1024, PREPROCESS_SUBSYSTEMS = 2048, MULTI_PILUC = 4096, RE_FACTOR = 8192} Flags;
   //typedef enum {PREPROCESS_INITIAL_SYSTE = 1} Flags;
    
//    typedef enum { DL_PCG = 1, FL_SBCG = 2, DL_BCG = 3} Flags;     
    
    // turn of column scaling : param.ipar[7]&=~(2+16+128);
    // transposed system is stored : param.ipar[4]=1; 
    // interchange order of scalings for all three cases: param.ipar[7]|=4+32+256;
     
  };

} // end of namespace

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "Ilupack.cc"
#endif

#endif
