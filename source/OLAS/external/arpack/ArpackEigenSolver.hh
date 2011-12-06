#ifndef OLAS_ARPACK_EIGENSOLVER_HH
#define OLAS_ARPACK_EIGENSOLVER_HH

#include "OLAS/solver/BaseEigenSolver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "ArpackSolver.hh"
#include "ArpackMatInterface.hh"

namespace CoupledField {
  
  class StdMatrix;
  class ParamNode;
  class ParamNode;
  
  // =========================================================================
  //   ARPACK SOLVER
  // =========================================================================
  
  
  //! Class for interfacing to ARPACK solvers
  class ArpackEigenSolver : public BaseEigenSolver {
    
  public:
    
    //! Default Constructor
    ArpackEigenSolver( shared_ptr<SolStrategy> strat,
                       PtrParamNode xml, 
                       PtrParamNode solverList,
                       PtrParamNode precondList,
                       PtrParamNode eigenInfo );
    
    //! Default Destructor
    virtual ~ArpackEigenSolver();

    //! Setup routine for standard eigenvalue problem

    //! Setup routine for various initialization tasks of a standard
    //! eigenvalue problem.
    //! \param mat Reference to matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver
    //!        is used
    void Setup( const BaseMatrix & mat,
                UInt numFreq, Double freqShift );

    //! Setup routine for a generalized eigenvalue problem
    
    //! Setup routine for various initialization tasks.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver
    //!        is used
    void Setup( const BaseMatrix & stiffMat,
                const BaseMatrix & massMat,
                UInt numFreq, Double freqShift );
    
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
    void Setup( const BaseMatrix & stiffMat,
                const BaseMatrix & massMat,
                const BaseMatrix & dampMat,
                UInt numFreq, Double freqShift );


    //! Solve the linear generalized eigenvalue problem
    
    //! This method triggers the calculation of the eigenvalue problem.
    //! Its return value is the number of converged eigenvalues and the
    //! related error.
    //! \param sol Vector with converged eigenvalues
    //! \param err Vector with error bound of eigenvalues
    //! \return Number of converged eigenvalues
    UInt CalcEigenFrequencies( BaseVector &sol,
                               BaseVector &err );
    
    //! Calculate a particular eigenmode as a postprocessing solution

    //! This method may be called after the CalcEigenFrequencies() method.
    //! It calculates a given eigenmode and stores in a use supplied vector.
    //! \param modeNr Number of the (converged) eigenmode to be calculated
    //! \param mode Vector with the eignmode
    void CalcEigenMode( UInt modeNr, Vector<Double> & mode );


    //! Calculate condition number

    //! This method calculates the condition number of the given matrix,
    //! as well as the 5 smallest and 5 largest eigenvalues
    void CalcConditionNumber( const BaseMatrix& mat,  
                              Double& condNumber, 
                              Vector<Double>& evs,
                              Vector<Double>& err );

  private:

    //! Method for writing logging information into .las file
    void PrintInfo();

    //! Pointer to matrix interface
    ArpackSolver * arpackSolver_;

    //! Pointer to matrix interface
    ArpackMatInterface * interface_;
    
    //! Pointer to matrices
    const StdMatrix * matrixA_, * matrixB_, *matrixD_;

    //! Pointer to solver ovbject
    BaseSolver * solver_;

    //! Pointer to precond object
    BasePrecond * precond_;

    /** Attribute for xml paramnode of <solver> section */
    PtrParamNode xml_;
    
    //! Flag indicating if problem is generalized problem
    bool isGeneralized_;

    //! Flag for shift-and-invert mode
    bool shiftAndInvert_;
    
    //! Flag for use of logging
    bool logging_;

    //! Character string for 'which' setting of  arpack
    char * which_;
  };
}

#endif
