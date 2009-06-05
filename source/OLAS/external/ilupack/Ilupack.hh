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

    /** Every call sets up a new preconditionier.
     * @param analysis_id shall be the current info/analysis/progress/step entry and contain an "analysis_id" element */
    void Setup(BaseMatrix &sysmat, InfoNode* analysis_id);
     
    /** To satisfy the compiler 
     * @param sysmat shall be the one Setup() is called with
     * @param precond ignored
     * @param analysis_id @see Setup() */
    void Solve( const BaseMatrix &sysmat, const BasePrecond &precond,
                const BaseVector &rhs, BaseVector &sol, InfoNode* analysis_id);

    /** This is not Ilupack intern but Olas stuff. It means that this is
     * a Ilupack, the actual solver called by the lib is in Ilupack::Solver */
    SolverType GetSolverType() { return ILUPACK_SOLVER; }

    /** This are the names of the matrix types in ILUPACK style */
    typedef enum { GNL, SYM, PD, HER } Matrix;

  private:

    /** Set up the ilupack matrix mat_ 
     * @param base_mat input */
    void SetMatrix(const BaseMatrix &base_mat);

    /** Determine the matrix type for ilupack. This can be set optionally in the xml file.
     * If it is not set, then it is SYM or GNL but neither SPD, HPD nor HER */
    void DetermineMatrixType(BaseMatrix &sysMat, InfoNode* out);

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

    /** Helper method for InitParameters. Takes the ilupack value, prints it to the InfoNode (via
     * param_name).
     * If there is a parameter with param_name in the xml file, it is printed.
     * @param param_name simple xpath chain via slash */
    void CheckParameter(InfoNode* out, double* ilupack_val, const char* param_name);
    void CheckParameter(InfoNode* out, char** ilupack_val, const char* param_name);
    void CheckParameter(InfoNode* out, int* ilupack_val, const char* param_name);  
    void CheckParameter(InfoNode* out, bool* ilupack_val, const char* param_name);
    
    /** The ilupack matrix types */
    Enum<Matrix> matrix; 
                            
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

  };

} // end of namespace

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "Ilupack.cc"
#endif

#endif
