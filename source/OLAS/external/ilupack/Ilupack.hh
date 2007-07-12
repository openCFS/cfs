// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ILUPACK_HH
#define ILUPACK_HH

#include "utils/environment.hh"
#include "solver/basesolver.hh"

#include "General/Enum.hh"

using CoupledField::Enum;

// include the original ilupack header
extern "C"
{
   #include "ilupack.h" 
}

namespace OLAS 
{
  class BaseMatrix;  
  class BaseVector;
  class BasePrecond;

  template<typename T>
  class Ilupack : public BaseIterativeSolver 
  {
  public:
    Ilupack(ParamNode* param, OLAS_Report *myReport, MatrixEntryType type);

    ~Ilupack();

    void Setup( BaseMatrix &sysmat);

     
    /** To satisfy the compiler */
    void Solve( const BaseMatrix &sysmat, const BasePrecond &precond,
                const BaseVector &rhs, BaseVector &sol)
    {
        Solve(sysmat, precond, rhs, sol, false);
    };                 


    /** @param recursive_call when the solver has a problem we Setup() again and
     *         call Solve() once recursively */
    void Solve( const BaseMatrix &sysmat, const BasePrecond &precond,
                const BaseVector &rhs, BaseVector &sol, bool recursive_call);

    /** This is not Ilupack intern but Olas stuff. It means that this is
     * a Ilupack, the actual solver called by the lib is in Ilupack::Solver */
    SolverType GetSolverType() { return ILUPACK_SOLVER; }

    /** The built-in solvers of ILUPACK, some might need additional libs.
     * Note, that this list mixes symmetric, unsymmetric, spd-solvers, ... */
    typedef enum { PCG = 1, SBCG = 2, BCG = 3, SQMR = 4, BCGSTAB = 5, TFQMR = 6, 
                   FOM = 7, GMRES = 8, FGMRES = 9, DQGMRES = 10} Solver;                     

    /** This are the names of the matrix types in ILUPACK style */
    typedef enum { GNL, SYM, PD, HER } Matrix;

    /** This are general reordering/permuation methods. Note, that the actual
     * Ilupak functions are combinations of Matrix, Reodering, Matching.
     * The combinations are quite special and driven by the available functions
     * in the lib */
    typedef enum { NONE, ND, RCM, AMF, AMD, MMD, METIS_E, METIS_N, INDSET, PP, FC, PQ } Ordering; 

    /** There are combinations with matching from Pardiso and MC64, most symmetric
     * is based on matching but a matching is not mandatory. */
    typedef enum { PURE, MC64, MWM } Matching;
    
    /** This are the permutation types */
    typedef enum { INITIAL, REGULAR, FINAL } PermutationRole;  

        
  private:
    /** Relases the memory allocated for ilupack by calling its delete method. 
     * Does nothing if mat_.a is NULL. */ 
    void ReleaseIlupackMemory();

    /** Set up the ilupack matrix mat_ 
     * @param base_mat input */
    void SetMatrix(const BaseMatrix &base_mat);

    /** Determine the matrix type for ilupack. This can be set optionally in the xml file.
     * If it is not set, then it is SYM or GNL but neither SPD, HPD nor HER */
    void DetermineMatrixType(BaseMatrix &sysMat);

    /** Initializes the parameter with default settings by ilupack, conditionally
     * overwrite them with the xml settings and print out the complete stuff. 
     * No permutation here! Note that the complex and real versions can be savely
     * casted against each other as its only the pointers. output to param_
     * and the ilupack matrix mat_ has to be set in advance by SetMatrix() */ 
    void InitParameters(); 
  
  
    /** This method sets up the "enums", it fills them with the string representations. */
    void SetEnums();

    /** Set the permuations - the list and not yet the actual inital, regualar and final. @see GetPermutation() */
    void SetPermutations();
     
    /** If the factorization didn't result with 0, this method encrypts the error code. */
    std::string PrecondError(int result, const DILUPACKparam& param, int level);

    /** If the ilupack solver didn't result with 0, this method encrypts the error code. */
    std::string SolverError(int result);     

    /** This enum holds the string representation of the ILUPACK solver types */ 
    Enum solver; 
 
    /** The ilupack reordering functions */
    Enum ordering;
    
    /** The ilupack matching */
    Enum matching;
    
    /** The ilupack matrix types */
    Enum matrix; 
    
    /** The Ilupack flags */
    Enum flags;
    
    /** The permutation roles (initial, ...) */
    Enum permutationRole;
    
    /** This is the ilupack signature of a real permutation function.
     * ILUPACK needs one for the inital, regular and final step */
    typedef int (*RealPermFunc) (Dmat, double*, double*, int*, int*, int*, DILUPACKparam*);
    
    typedef int (*ComplexPermFunc) (Zmat, doublecomplex *,doublecomplex *, int*, int*, int*, ZILUPACKparam*);    
    
    /** This structure holds the definition of a ilupack permutation function in the lib. The key fields are matrix, reordering and matching */
    struct Permutation
    {
        /** The matrix from the capital letter code of the ilupack function, doesn't need to
         * be of the real problem type. Part of the key. */
        Matrix     matrix;
        
        /** The reordering version, with or without matching. Part of the key. */
        Ordering ordering;
        
        /** In the symmetric case most reoderings is done by one of the 2 matching libs
         * (MC64 or MWM/Pardiso), PURE for without. Part of the key. */
        Matching   matching;
        
        /** The function address as in ilupack.h */
        RealPermFunc realFunc;
        ComplexPermFunc cplxFunc;
        
        /** The name corresponding to the function omitting D and Z for the type */
        std::string ilupack_name;
        
        /** The code as used in ilupack to help orientation */
        std::string code;
        
        /** The description as taken from the ilupack samples or own stuff */
        std::string description;
    };

    /** A convenience method to fill the permuations list*/
    void AddPermutation(Matrix matrix, Ordering ordering, Matching matching,
                        RealPermFunc realFunc, ComplexPermFunc cplxFunc, 
                        const std::string& ilupack_name, const std::string& code,
                        const std::string& description);
                        

    /** Dump the solver stus, extracts nice parameters and prints them */
    void DumpSolverStatus(std::ostream& out);

    /** Dump the parameter setting */
    void DumpParameters(std::ostream& out);

    /** This method dumps the permutation list */
    void DumpPermutations(std::ostream& out);
     
    /** Dump the ilupack flags */ 
    void DumpFlags(int flag, std::ostream& out);

    /** Read the flag settings as provided by the XML file and
     * apply it on the given flag.
     * @param flag the in/out parameter */ 
     void ApplyFlagSettings(int& flag);

    /** Gets the permutation. Either from XML, or from default if not present there.
     * @param role is INITIAL, ...
     * @param out the result if true is returned
     * @throws exception if parameter set in XML or default is invalid */ 
    void GetPermutation(PermutationRole role, Permutation& out);

    /** Simple helper function to find a specigic permutation.
     * @param out the result if true is returned
     * @throws exception if nothing found */  
    void GetPermutation(Matrix matrix, Ordering ordering, Matching matching, Permutation& out);  

    /** This is a simple expert mode for default permutation. Used implicitly by GetPermutation() 
     * @param role if INITIAL, ...
     * @param out the result - if not an exception is thrown*/
    void GetDefaultPermutation(PermutationRole role, Permutation& out);

    /** This vector holds all interesting permutation functions from the available ones in the lib.
     * E.g. only double and double comples */
    std::vector<Permutation> permutations;
    
    /** The inital permutation. All permutations are set via GetPermutation() */
    Permutation initial_;
    
    /** @see intial_ */
    Permutation regular_;
    
    /** @see initial_ */
    Permutation final_;

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
    DILUPACKparam param_;

    /** A pointer to param_ */
    DILUPACKparam* dparam_ptr_;
    /** Also a pointer to (DILUPACKparam) param_ as complex is same-same here */
    ZILUPACKparam* zparam_ptr_;
    
    
    /** Here ilupack stores the preconditioner. */
    DAMGlevelmat precond_;

    /** A pointer to precond_ */
    DAMGlevelmat* dprecond_ptr_;
    /** Also a pointer to (DILUPACKparam) param_ as complex is same-same here */
    ZAMGlevelmat* zprecond_ptr_;
    
    /** This paramerters is automatically determined by ilupack and required for the solver */
    int nlev_; 

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
     * REL_ERROR_TOL = relative error tolerance -> param.fpar[20]=restol; <br>
     * ABS_ERROR_TOL = absolute error tolerance -> param.fpar[21]=0;<br>
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

   //typedef enum { PREPROCESS_INITIAL_SYSTEM = 1024, PREPROCESS_SUBSYSTEMS = 2048, MULTI_PILUC = 4096, RE_FACTOR = 8192} Flags;
   //typedef enum {PREPROCESS_INITIAL_SYSTE = 1} Flags;
    
//    typedef enum { DL_PCG = 1, FL_SBCG = 2, DL_BCG = 3} Flags;     
    
    // turn of column scaling : param.ipar[7]&=~(2+16+128);
    // transposed system is stored : param.ipar[4]=1; 
    // interchange order of scalings for all three cases: param.ipar[7]|=4+32+256;
     
  };

} // end of namespace

#endif
