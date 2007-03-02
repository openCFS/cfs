#ifndef OLAS_BASE_EIGENSOLVER_HH
#define OLAS_BASE_EIGENSOLVER_HH

#include "matvec/matvec.hh"
#include "algsys/olasparams.hh"

namespace OLAS {
  
  
  // forward class declaration
  //class BaseSolver;
  
  // =========================================================================
  // BASE EIGENVALUE SOLVER
  // =========================================================================
  
  //! Base class for algebraic system eigenvalue solver
  class BaseEigenSolver {
    
  public:
    
    //! Default Constructor
    BaseEigenSolver( OLAS_Params *myParams, OLAS_Report *myReport )
      : myParams_(myParams),
        myReport_(myReport),
        numFreq_(0),
        freqShift_(0.0),
        isQuadratic_(false)
    {
      ENTER_FCN( "BaseEigenSolver::BaseEigenSolver" );
    }
    
    //! Default Destructor
    virtual ~BaseEigenSolver() {
      ENTER_FCN( "BaseEigenSolver::~BaseEigenSolver" );
    }

    //! Get type of eigenvalue problem to be solved
    bool IsQuadratic() {return isQuadratic_;}
    
    //! Setup routine for a generalized eigenvalue problem
    
    //! Setup routine for various initialization tasks of a generalized
    //! eigenvalue problem.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver
    //!        is used
    virtual void Setup( const BaseMatrix & stiffMat,
                        const BaseMatrix & massMat,
                        UInt numFreq, Double freqShift ) = 0;
    
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
    virtual void Setup( const BaseMatrix & stiffMat,
                        const BaseMatrix & massMat,
                        const BaseMatrix & dampMat,
                        UInt numFreq, Double freqShift ) = 0;

    //! Solve the linear generalized eigenvalue problem
    
    //! This method triggers the calculation of the eigenvalue problem.
    //! Its return value is the number of converged eigenvalues and the
    //! related error.
    //! \param sol Vector with converged eigenvalues
    //! \param err Vector with error bound of eigenvalues
    //! \return Number of converged eigenvalues
    virtual UInt CalcEigenFrequencies( BaseVector &sol,
                                       BaseVector &err ) = 0;

    //! Calculate a particular eigenmode as a postprocessing solution

    //! This method may be called after the CalcEigenFrequencies() method.
    //! It calculates a given eigenmode and stores in a use supplied vector.
    //! \param modeNr Number of the (converged) eigenmode to be calculated
    //! \param mode Vector with the eignmode
    virtual void CalcEigenMode( UInt modeNr, Vector<Double> & mode ) = 0;
    
    
    //! Method to force instantiation of all public member functions
    
    //! This auxillary method is used in our factory concept for solver
    //! generation. The factory function GenerateEigenSolverObject() will
    //! make a pseudo call to this method in order to force the compiler
    //! to instantiate all public methods of a templated solver class.
    //! \note
    //! - The method must never be actually called, since it does not perform
    //!   any sensible operations.
    //! - If a templated solver offers additional public methods besides
    //!   the ones defined in the BaseSolver class, then it should over-write
    //!   the InstantiateAdditionalPublicMethods() member function which is
    //!   called by this method.
    void InstantiatePublicMethods( BaseMatrix &sysMat ) {
      // To be implemented ...
    }

  protected: 

    //! Pointer to parameter object

    //! This is a pointer to a parameter object containing the steering
    //! parameters for this solver.
    OLAS_Params *myParams_;
    
    //! Pointer to report object
    
    //! This is a pointer to a report object in which the solver will store
    //! general information about the solution of a linear system.
    OLAS_Report *myReport_;
    
    //! Number of frequencies to be calculated
    UInt numFreq_;
    
    //! Value of frequency shift
    Double freqShift_;

    //! Flag indicating if a quadratic eigenvalue problem is solved
    bool isQuadratic_;
  };
  
}
  
#endif
