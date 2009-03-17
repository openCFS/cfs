// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_ARPACK_SOLVER_HH
#define OLAS_ARPACK_SOLVER_HH

#include "arpackMatInterface.hh"
#include "arpackFortranInterface.hh"

namespace CoupledField {
  // =========================================================================
  //   ARPACK SOLVER FOR SYMMETRIC REAL SYSTEMS
  // =========================================================================
  
  
  //! Class for interfacing with ARPACK FORTRAN library
  //class ArpackSolver : public BaseEigenSolver {
  class ArpackSolver {
    
  public:
    
    //! Default Constructor
    ArpackSolver();
    
    //! Default Destructor
    virtual ~ArpackSolver();

    //! General setup routine
    
    //! Setup routine for various initialization tasks.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver
    //!        is used
    void Setup( ArpackMatInterface *apInterface, 
                UInt size, UInt numFreq, Double freqShift, char* which, 
                char* type, bool shiftMode );

    //! Solve the linear generalized eigenvalue problem
    
    //! This method triggers the calculation of the eigenvalue problem.
    //! Its return value is the number of converged eigenvalues.
    //! \param sol Vector with converged eigenvalues
    //! \return Number of converged eigenvalues
    UInt FindEigenvalues();

    //! This method returns the n-th converged eigenvalue
    //! \param n Number of requested converged eigenvalue
    //! \return Calculated eigenvalue
    Double Eigenvalue(UInt n);

    //! This method returns tolerance obtained for
    //! the n-th converged eigenvalue
    //! \param n Number of requested converged eigenvalue
    //! \return Calculated tolerance for eigenvalue
    Double Tolerance(UInt n);

    //! Calculate a particular eigenmode as a postprocessing solution

    //! This method may be called after the CalcEigenFrequencies() method.
    //! It calculates a given eigenmode and stores in a use supplied vector.
    //! \param modeNr Number of the (converged) eigenmode to be calculated
    //! \param mode Vector with the eignmode
    Double* GetEigenvector( UInt modeNr );

    //! Method to switch arpack debug information off
    void DebugOff();

    //! Method to switch arpack debug information on
    void DebugOn();

    //! Method to set user specified tolerance
    void SetTolerance(Double tol);
    //! Method to set user specified number of iterations
    void SetIterations(Integer maxIt);
    //! Method to set user specified number of Arnoldi vectors
    void SetNumVectors(Integer numVec);

    //! Report access functions to parameters 
    //! number of wanted eigenvalues
    UInt GetNev();
    //! arpack 'which' setting for eigenvalue calculation
    char* GetWhich();
    //! tolerance setting
    Double GetShift();
    //! tolerance setting
    Double GetTol();
    //! number of iterations allowed
    UInt GetMaxit();
    //! number of Arnoldi vectors
    UInt GetNcv();

  private:

    void InitTempSpace(Double*, Double*, Double*, Double*, Double*);
    //! Pointer to parameter object

    //! Pointer to matrix interface
    ArpackMatInterface * interface_;
    
    //! Flag for shift-and-invert mode
    bool shiftAndInvert_;

    //! Character string for 'which' setting of  arpack
    char * which_;
    
    //! Character string denoting generalized or standard eigenvalue problem
    char* type_;

    //! Tolerance to be achieved in solution
    Double tolerance_;

    //! Shift to be applied
    Double freqShift_;

    //! Maximum number of iterations
    UInt numFreq_;

    //! Maximum number of iterations
    UInt maxIterations_;

    //! Number of Arnoldi vectors
    UInt numArnoldiVec_;

    //! Size of equation system
    UInt size_;

    //! stores the calculated eigenvalues
    Double *eigenValues_;

    //! stores the tolerances for the calculated eigenvalues
    Double *eigenTolerances_;

    //! stores the calculated eigenvectors
    Double *eigenVectors_;

  };
}

#endif
