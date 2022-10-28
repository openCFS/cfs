#ifndef ILUPACK_HH
#define ILUPACK_HH


#include<def_use_ilupack.hh>
#include "General/Environment.hh"
#include "OLAS/solver/BaseSolver.hh"

#include "General/Enum.hh"

// include the original ilupack header
extern "C"
{
#include <ilupack.h>
#ifdef USE_ILUPACK_PARALLEL
#include <SparseSymmetricNew.h>
#endif

}

// Due to some issues in building currently we are using the serial version of the Ilupack 2.2.1 , but the parallel version
// can also be built along with CFS (ie OPENMP version with mtmetis and ilupack own metis).
// The speedup is not significant yet but as it is being currently developed , it might be worth a try in future
// (In other words it means I spent a few weeks and I don't want to remove the code entirely :) )

// Ilupack parallel is a superset of the original ilupack

namespace CoupledField 
{
  class BaseMatrix;  
  class BaseVector;
  class Flags;

  
  /** Ilupack has variants DILUPACKparam/ZILUPACKParam, Dmat/ZMat and DAMGlevelmat/ZAMGlevelmat.
   * Only DILUPACKparam and ZILUPACKParam differ in size as in two cases there are instances but
   * not pointers of the type used. 
   * The following struct is modification of the original Ilupack stuff. */

#ifdef USE_ILUPACK_PARALLEL
  template<typename TYPE>
  struct TYPE_ILUPACKparam
  {
    integer          ipar[ILUPACK_NIPAR];
    TYPE             fpar[ILUPACK_NFPAR];
    integer          type;
    integer          *ibuff;
    integer          *iaux;
    TYPE             *dbuff;
    TYPE             *daux;
    integer          *ju;
    integer          *jlu;
    TYPE             *alu;
    TYPE             *testvector;
    size_t           nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
    integer          rcomflag, returnlabel;
    TYPE             *tv;
    integer          *ind;
    integer          nindicator;
    integer          *indicator;
    Dmat             A;
    integer          istack[30], *pistack[20];
    double           rstack[30], *prstack[10];
    TYPE             fstack[30], *pfstack[10];
    size_t           ststack[5], *pststack[5];
    Dmat             mstack[1];
    DAMGlevelmat     *amglmstack[1];
    integer          (*intfctstack[3])();
    integer          matching;
    char             *ordering;
    double  droptol;
    double  droptolS;
    double  droptolc;
    double  condest;
    double restol;
    integer          maxit;
    double  elbow;
    integer          lfil;
    integer          lfilS;
    char             *typetv;
    char             *amg;
    integer          npresmoothing;
    integer          npostsmoothing;
    integer          ncoarse;
    char             *presmoother;
    char             *postsmoother;
    char             *FCpart;
    char             *typecoarse;
    integer          nrestart;
    integer          flags;
    char             *solver;
    TYPE  damping;
    TYPE  contraction;
    integer          mixedprecision;
    integer          (*perm0)();
    integer          (*perm)();
    integer          (*permf)();
    TYPE  shift0;
    TYPE  shiftmax;
    integer          nshifts;
    TYPE  *shifts;
    Dmat             *shiftmatrix;
    integer          niter;
    integer          nthreads;
    integer          *indpartial;
    integer          *indexpartial;
    TYPE  *valpartial;
    TYPE  *valuepartial;

  } ;


#else
  template<typename TYPE>
  struct TYPE_ILUPACKparam
    {
       integer              ipar[ILUPACK_NIPAR];
       doubleprecision      fpar[ILUPACK_NFPAR];
       integer              type;
       integer             *ibuff;
       integer             *iaux;
       TYPE                *dbuff;
       TYPE                *daux;
       integer             *ju;
       integer             *jlu;
       TYPE                *alu;
       TYPE                *testvector;
       size_t               nibuff, ndbuff, nju,njlu,nalu, ndaux,niaux,ntestvector;
       integer              rcomflag, returnlabel;
       TYPE                *tv;
       integer             *ind;
       integer              nindicator;
       integer             *indicator;
       Dmat                 A;
       integer              istack[30], *pistack[20];
       doubleprecision      rstack[30], *prstack[10];
       TYPE                 fstack[30], *pfstack[10];
       size_t               ststack[5], *pststack[5];
       Dmat                 mstack[1];
       DAMGlevelmat        *amglmstack[1];
       integer            (*intfctstack[3])();
       integer              matching;
       char                *ordering;
       doubleprecision      droptol;
       doubleprecision      droptolS;
       doubleprecision      condest;
       doubleprecision      restol;
       integer              maxit;
       doubleprecision      elbow;
       integer              lfil;
       integer              lfilS;
       char                *typetv;
       char                *amg;
       integer              npresmoothing;
       integer              npostsmoothing;
       integer              ncoarse;
       char                *presmoother;
       char                *postsmoother;
       char                *FCpart;
       char                *typecoarse;
       integer              nrestart;
       integer              flags;
       char                *solver;
       TYPE                 damping;
       integer            (*perm0)();
       integer            (*perm)();
       integer            (*permf)();
       integer            isreal;
       integer            issingle;
    } ;

#endif

  template<typename T>
  class Ilupack : public BaseIterativeSolver 
  {
  public:
    Ilupack(PtrParamNode param, PtrParamNode olasInfo, BaseMatrix::EntryType type);

    ~Ilupack();

    /** Every call sets up a new preconditionier. */
    void Setup(BaseMatrix &sysmat);
     
    /** To satisfy the compiler 
     * @param sysmat shall be the one Setup() is called with */
    void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol);

    /** This is not Ilupack intern but Olas stuff. It means that this is
     * a Ilupack, the actual solver called by the lib is in Ilupack::Solver */
    BaseSolver::SolverType GetSolverType() { return BaseSolver::ILUPACK; }

    /** This are the names of the matrix types in ILUPACK style */
    typedef enum { GNL, SYM, PD, HER } Matrix;

  private:

    /** Set up the ilupack matrix mat 
     * @param base_mat input */
    void SetMatrix(const BaseMatrix &base_mat);

    /** Determine the matrix type for ilupack. This can be set optionally in the xml file.
     * If it is not set, then it is SYM or GNL but neither SPD, HPD nor HER */
    void DetermineMatrixType(BaseMatrix &sysMat, PtrParamNode out);

    /** Initializes the parameter with default settings by ilupack, conditionally
     * overwrite them with the xml settings and print out the complete stuff. 
     * Note that the complex and real versions can be savely
     * casted against each other as its only the pointers. output to param
     * and the ilupack matrix mat has to be set in advance by SetMatrix() */ 
    void InitParameters(); 
  
    /** Calls the AMGinit() variant of for the system from the lib */
    void IlupackAMGInit();
    
    /** Calls the AMGFactor() variant of for the system from the lib */
    int IlupackAMGFactor();

    int IlupackAMGSolver(T* sol_ptr, T* rhs_ptr);

    /** Relases the memory allocated for ilupack by calling its delete method. 
     * Does nothing if mat.a is NULL. */ 
    void IlupackAMGDelete();
    
    /** This method sets up the "enums", it fills them with the string representations. */
    void SetEnums();

    /** The ilupack matrix types */
    Enum<Matrix> matrix; 
                            
    /** Calculates the fill in factor of the preconditioner */
    void CalcFillIn(PtrParamNode out);
    
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
    Dmat mat;
    
    /** Here ilupack stores all the paramters, initialized in InitParameters(). 
     * See TYPE_ILUPACKparam above. */    
    TYPE_ILUPACKparam<T> param;
    
    /** Here ilupack stores the preconditioner. */
    DAMGlevelmat precond;

#ifdef USE_ILUPACK_PARALLEL
    SparseMatrix spr;
#endif

    //The OpenMP version of ILUPACK creates a task dependency tree to parallelize the computation and nleaves describes number
    // of tasks to create.
    int nleaves=1;
    int index = 0;// This indicates the index used for the sparse matrix which is required in parallel ilupack
    int mtmetis=0; //the number of threads for the multithreaded metis

    bool isParallel=false;
    bool firstSetup=true;
  };

} // end of namespace

#endif
