// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASE_EIGENSOLVER_HH
#define OLAS_BASE_EIGENSOLVER_HH

#include "General/Enum.hh"
#include "MatVec/BaseMatrix.hh"
#include "MatVec/Vector.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Environment.hh"

namespace CoupledField {
  
  class OLAS_BaseMatrix;
  class OLAS_BaseVector;
  template<typename> class Vector;
  
  // forward class declaration
  class SolStrategy;
  class StdMatrix;

  class Timer;

  // =========================================================================
  // BASE EIGENVALUE SOLVER
  // =========================================================================
  
  //! Base class for algebraic system eigenvalue solver
  class BaseEigenSolver {
  
  public:
    //! Type of EigenSolver

    typedef enum {NO_EIGENSOLVER, ARPACK, PHIST, FEAST, PALM, QUADRATIC, EXTERNAL} EigenSolverType;
    //! This enumeration data type describes the type of eigensolver which is
    //! applied to solve a generalized eigenvalue problem. The enumeration
    //! contains the following:
    //! - NOEIGENSOLVER
    //! - ARPACK
    //! - PHIST
    //! - FEAST
    //! - PALM
    //! - QUADRATIC
    //! - EXTERNAL
    static Enum<EigenSolverType> eigenSolverType;

    EigenSolverType eigenSolverType_;

    typedef enum {NONE, MAX, NORM} ModeNormalization;
    //! Defines how to normalized the modes
    //! - NONE : do not change what is returned by the solver
    //! - MAX : normalize such that the maximum absolute entry has a value of 1
    //! - NORM : normalize such that the L2 norm of each mode is 1

    //! Default Constructor
    BaseEigenSolver( shared_ptr<SolStrategy> strat, 
                     PtrParamNode eSolverXML,
                     PtrParamNode solverList,
                     PtrParamNode precondList,
                     PtrParamNode eigenInfo )
      : eigenSolverType_(NO_EIGENSOLVER),
        solStrat_(strat),
        xml_(eSolverXML),
        solverList_(solverList),
        precondList_(precondList),
        info_(eigenInfo),
        numFreq_(0),
        freqShift_(0.0),
        isQuadratic_(false),
        isBloch_(false),
//        isBuckling_(false),
        sort_(false),
//        isReal_(true),
//        isSymmetric_(true),
//        isHermitian_(false),
        a_(NULL),
        eigenProblemType_(NO_TYPE),
        modeNormalization_(NONE),
        eigenSolverName_(NO_EIGENSOLVER)
    {
    }
    
    //! Default Destructor
    virtual ~BaseEigenSolver() {
    }

    /** Does constructor stuff only possible after child constructors are called */
    void PostInit();

    //! Get type of eigenvalue problem to be solved
    bool IsQuadratic() const {return isQuadratic_;}// ToDo: remove due to new structure

    //TODO:For testing, remove
    // BaseEigenSolver* Getirgendwas(){return NULL;}

    /** do we calculate in complex generalized EV mode */
    bool IsBloch() const { return isBloch_; }// ToDo: remove due to new structure

    // bool IsComplex() {return (!isSymmetric_ || !isReal_);} // ToDo: remove due to new structure

    //! Setup routine for standard eigenvalue problem

    //! Setup for a standard EVP
    virtual void Setup(const BaseMatrix & A, bool isHermitian=false) =0;

    //! Setup for a generalised EVP
    virtual void Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian=false) =0;

    //! Setup for a quadratic EVP
    virtual void Setup(const BaseMatrix & K, const BaseMatrix & C, const BaseMatrix & M) =0;

    //! returns if eigenvalues will be complex
    bool HasComplexEigenvalues(){
        switch (eigenProblemType_){
            case REAL_SYMMETRIC: return false;
            case REAL_GENERAL: return true;
            case COMPLEX_SYMMETRIC: return true;
            case COMPLEX_HERMITIAN: return false;
            case COMPLEX_GENERAL: return true;
            default: EXCEPTION("set the problem type"); return true;
        }
    }
    //! returns if Modes are complex
    bool HasComplexModes(){
        switch (eigenProblemType_){
            case REAL_SYMMETRIC: return false;
            case REAL_GENERAL: return true;
            case COMPLEX_SYMMETRIC: return true;
            case COMPLEX_HERMITIAN: return true;
            case COMPLEX_GENERAL: return true;
            default: EXCEPTION("set the problem type"); return true;
        }
    }

    //! Setup routine for various initialization tasks of a standard
    //! eigenvalue problem.
    //! \param mat Reference to matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param sort
    virtual void Setup(const BaseMatrix & mat,  UInt numFreq, double freqShift, bool sort) = 0;// ToDo: remove due to new structure
    
    //! Setup routine for a generalized eigenvalue problem
    
    //! Setup routine for various initialization tasks of a generalized
    //! eigenvalue problem.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    virtual void Setup( const BaseMatrix & stiffMat, const BaseMatrix & massMat,
                        UInt numFreq, double freqShift, bool sort, bool bloch) = 0;// ToDo: remove due to new structure
    
    //! Setup routine for a quadratic eigenvalue problem
    
    //! Setup routine for various initialization tasks of a quadratic 
    //! eigenvalue problem.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param dampMat Reference to damping matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    virtual void Setup( const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
                        UInt numFreq, double freqShift, bool sort) = 0;// ToDo: remove due to new structure

    /** Solve the linear generalized eigenvalue problem.
     *  This method triggers the calculation of the eigenvalue problem.
     * @param sol Vector with converged eigenvalues. The size of sol is the number of converged eigenvalues!
     * @param err Vector with error bound of eigenvalues */
    virtual void CalcEigenFrequencies( BaseVector &sol, BaseVector &err ) = 0;// ToDo: remove due to new structure

    //! compute the eigenvalues based on solver defined settings
    //! the return type StdVector<T> depends on the problem type
    virtual void CalcEigenValues( BaseVector &sol, BaseVector &err) = 0;

    //! compute the eigenvalues in an interval [minVal,maxVal]
    //! the return type StdVector<T> depends on the problem type
    virtual void CalcEigenValues( BaseVector &sol, BaseVector &err, Double minVal, Double maxVal ) = 0;

    //! compute N eigenvalues closest to the shift point
    //! the return type StdVector<T> depends on the problem type
    virtual void CalcEigenValues( BaseVector &sol, BaseVector &err, UInt N, Complex shiftPoint ) = 0;
    virtual void CalcEigenValues( BaseVector &sol, BaseVector &err, UInt N, Double shiftPoint ) = 0;
    //! Calculate a particular eigenmode as a postprocessing solution

    //! This method may be called after the CalcEigenFrequencies() method.
    //! It returns a given eigenmode and stores in a use supplied vector.
    //! \param modeNr Number of the (converged) eigenmode to be calculated
    //! \param mode Vector with the eignmode
    virtual void GetEigenMode( UInt modeNr, Vector<Complex> & mode, bool right=true ) = 0;
    virtual void GetComplexEigenMode( UInt modeNr, Vector<Complex> & mode ) = 0;

    
    //! Calculate condition number

    //! This method calculates the condition number of the given matrix,
    //! as well as the 5 smallest and 5 largest eigenvalues
    virtual void CalcConditionNumber( const BaseMatrix& mat,  
                                      Double& condNumber, 
                                      Vector<Double>& evs,
                                      Vector<Double>& err ) = 0;
    
    void CheckMatrix(bool & isReal, bool & isStoredSymmetric,const BaseMatrix & A ) {

        switch (A.GetStorageType()) {
        case BaseMatrix::SPARSE_SYM: isStoredSymmetric=true; break;
        case BaseMatrix::SPARSE_NONSYM: isStoredSymmetric=false; break;
        default: EXCEPTION("storage type" << A.GetStorageType() << " not handeled");
        }
        switch(A.GetEntryType()){
        case BaseMatrix::DOUBLE: isReal=true; break;
        case BaseMatrix::COMPLEX: isReal=false; break;
        default: EXCEPTION("Matrix should be DOUBLE or COMPLEX");
        }

    }

    //! set the problem type depending on the matrix properties A (real|complex, symmetric|non-symmetric)
    //! optional: specify if the problem is Hermitian
    void SetProblemType(const BaseMatrix & A, bool isHermitian=false) {
        bool isReal, isStoredSymmetric;
        CheckMatrix(isReal, isStoredSymmetric, A);
        EigenValueProblemType newType;
        if (isReal) {
            if (isStoredSymmetric) {
                newType = REAL_SYMMETRIC;
            }
            else {
                newType = REAL_GENERAL;
            }
            if (isHermitian) {
                WARN("isHermitian=true does not make sense for a real matrix")
            }
        }
        else {
            if (isHermitian) {
                newType = COMPLEX_HERMITIAN;
                if (!isStoredSymmetric) {
                    EXCEPTION("non-symmetric matrix storage used for hermitian EVP -> Use symmetric matrix storage")
                }
            }
            else {
                if (isStoredSymmetric) {
                    newType = COMPLEX_SYMMETRIC;
                }
                else {
                    newType = COMPLEX_GENERAL;
                }
            }
        }
        if (eigenProblemType_==NO_TYPE) {// was not set yet so just set it
            eigenProblemType_=newType;
        }
        else { // check if it is consistent
            if (eigenProblemType_!=newType) {
                EXCEPTION("cannot change eigen problem type")
            }
        }
    }

    //! normalize a mode
    virtual void NormalizeMode(Vector<Complex> & mode, ModeNormalization type ) {
      double factor;
      switch (type) {
        case ModeNormalization::NONE :
          factor = 1.0;
          break;
        case ModeNormalization::MAX :
          int loc;
          factor = std::abs(mode.MaxAbs(loc));
          break;
        case ModeNormalization::NORM:
          factor = mode.NormL2();
          break;
        default: EXCEPTION("ModeNormalization not known");
      }
      mode.ScalarDiv(factor);
    }

    void GetNormalizedEigenMode( UInt modeNr, Vector<Complex> & mode, bool right=true){
      GetEigenMode(modeNr, mode, right);
      this->NormalizeMode(mode, modeNormalization_);
    }

    void SetModeNormalization(ModeNormalization normType) {
      modeNormalization_ = normType;
    }

    ModeNormalization GetModeNormalization() { return modeNormalization_; }

    virtual EigenSolverType GetEigenSolverName(){ return eigenSolverName_; }

    /** Gives the timer located within PtrParamNode */
    shared_ptr<Timer> GetSetupTimer() { return setupTimer_; }
    shared_ptr<Timer> GetSolveTimer() { return solveTimer_; }

  protected: 

    //! Pointer to solution strategy object
    shared_ptr<SolStrategy> solStrat_;
    
    //! Pointer to parameter object

    //! This is a pointer to a parameter object containing the steering
    //! parameters for this solver.
    PtrParamNode xml_;
    
    //! Pointer to solverList element
    
    //! Stores the pointer to the <solverList> xml. As the arpackSolver
    //! internally also utilizes a solver, we need to pass this
    //! information to the factory method for solver creation.
    PtrParamNode solverList_;
    
    //! Pointer to precondList element

    //! Stores the pointer to the <precondList> xml. As the arpackSolver
    //! internally also utilizes a solver, we need to pass this
    //! information to the factory method for solver creation.
    PtrParamNode precondList_;
    
    //! Pointer to report object
    
    //! This is a pointer to a report object in which the solver will store
    //! general information about the solution of a linear system.
    PtrParamNode info_;
    
    /** shall we scale the B matrix (mass) as suggested by Jonas? */
    bool scale_B_ = false;

    /** in case we have scale_B_, this is the value */
    double scale_B_val_ = 1.0;

    //! Number of frequencies to be calculated
    UInt numFreq_;// ToDo: remove due to new structure
    
    //! Value of frequency shift
    Double freqShift_;// ToDo: remove due to new structure

    //! Flag indicating if a quadratic eigenvalue problem is solved
    bool isQuadratic_;// ToDo: remove due to new structure

    //! Flag indication if a complex generalized bloch mode EV problem is solved
    bool isBloch_;// ToDo: remove due to new structure

    //! Flag indication if a buckling EV problem is solved
    // bool isBuckling_;// ToDo: remove due to new structure

    /** shall we sort the evs` */
    bool sort_;// ToDo: remove due to new structure

    //! flag to specify if the solution is real
    // bool isReal_;// ToDo: remove due to new structure

    //! flag to specify if all matrices are symmetric
    // bool isSymmetric_;// ToDo: remove due to new structure

    // bool isHermitian_;// ToDo: remove due to new structure

    //! the matrix to solve
    const StdMatrix* a_;

    EigenValueProblemType eigenProblemType_;

    //! defines the mode normalization
    ModeNormalization modeNormalization_;

    //! Solver Type Name
    EigenSolverType eigenSolverName_;

    /** This is a pointer to the setup timer. Located within PtrParamNode*/
    shared_ptr<Timer> setupTimer_;
    shared_ptr<Timer> solveTimer_;

  };
  
}
  
#endif
