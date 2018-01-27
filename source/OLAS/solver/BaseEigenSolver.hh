// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASE_EIGENSOLVER_HH
#define OLAS_BASE_EIGENSOLVER_HH

#include "General/Enum.hh"
#include "MatVec/BaseMatrix.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Environment.hh"

namespace CoupledField {
  
  class OLAS_BaseMatrix;
  class OLAS_BaseVector;
  template<typename> class Vector;
  
  // forward class declaration
  class SolStrategy;
  
  // =========================================================================
  // BASE EIGENVALUE SOLVER
  // =========================================================================
  
  //! Base class for algebraic system eigenvalue solver
  class BaseEigenSolver {
  
  public:
    //! Type of EigenSolver

    //! This enumeration data type describes the type of eigensolver which is
    //! applied to solve a generalized eigenvalue problem. The enumeration
    //! contains the following:
    //! - NOEIGENSOLVER
    //! - ARPACK
    //! - FEAST
    typedef enum {NO_EIGENSOLVER, ARPACK, FEAST} EigenSolverType;
    static Enum<EigenSolverType> eigenSolverType;    
    
  public:
    
    //! Default Constructor
    BaseEigenSolver( shared_ptr<SolStrategy> strat, 
                     PtrParamNode eSolverXML,
                     PtrParamNode solverList,
                     PtrParamNode precondList,
                     PtrParamNode eigenInfo )
      : solStrat_(strat),
        xml_(eSolverXML),
        solverList_(solverList),
        precondList_(precondList),
        info_(eigenInfo),
        numFreq_(0),
        freqShift_(0.0),
        isQuadratic_(false),
        isBloch_(false),
        sort_(false)
    {
    }
    
    //! Default Destructor
    virtual ~BaseEigenSolver() {
    }

    //! Get type of eigenvalue problem to be solved
    bool IsQuadratic() const {return isQuadratic_;}

    /** do we calculate in complex generalized EV mode */
    bool IsBloch() const { return isBloch_; }

    //! Setup routine for standard eigenvalue problem

    //! Setup routine for various initialization tasks of a standard
    //! eigenvalue problem.
    //! \param mat Reference to matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver is used
    //! \param sort
    virtual void Setup(const BaseMatrix & mat,  UInt numFreq, double freqShift, bool sort) = 0;
    
    //! Setup routine for a generalized eigenvalue problem
    
    //! Setup routine for various initialization tasks of a generalized
    //! eigenvalue problem.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver
    //!        is used
    virtual void Setup( const BaseMatrix & stiffMat, const BaseMatrix & massMat,
                        UInt numFreq, double freqShift, bool sort, bool bloch) = 0;
    
    //! Setup routine for a quadratic eigenvalue problem
    
    //! Setup routine for various initialization tasks of a quadratic 
    //! eigenvalue problem.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param dampMat Reference to damping matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver
    //!        is used
    virtual void Setup( const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
                        UInt numFreq, double freqShift, bool sort) = 0;

    /** Solve the linear generalized eigenvalue problem.
     *  This method triggers the calculation of the eigenvalue problem.
     * @param sol Vector with converged eigenvalues. The size of sol is the number of converged eigenvalues!
     * @param err Vector with error bound of eigenvalues */
    virtual void CalcEigenFrequencies( BaseVector &sol, BaseVector &err ) = 0;

    //! Calculate a particular eigenmode as a postprocessing solution

    //! This method may be called after the CalcEigenFrequencies() method.
    //! It returns a given eigenmode and stores in a use supplied vector.
    //! \param modeNr Number of the (converged) eigenmode to be calculated
    //! \param mode Vector with the eignmode
    virtual void GetEigenMode( UInt modeNr, Vector<Complex> & mode ) = 0;
    virtual void GetComplexEigenMode( UInt modeNr, Vector<Complex> & mode ) = 0;
    
    
    //! Calculate condition number

    //! This method calculates the condition number of the given matrix,
    //! as well as the 5 smallest and 5 largest eigenvalues
    virtual void CalcConditionNumber( const BaseMatrix& mat,  
                                      Double& condNumber, 
                                      Vector<Double>& evs,
                                      Vector<Double>& err ) = 0;
    
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
    
    //! Number of frequencies to be calculated
    UInt numFreq_;
    
    //! Value of frequency shift
    Double freqShift_;

    //! Flag indicating if a quadratic eigenvalue problem is solved
    bool isQuadratic_;

    //! Flag indication if a complex generalized bloch mode EV problem is solved
    bool isBloch_;

    /** shall we sort the evs` */
    bool sort_;
  };
  
}
  
#endif
