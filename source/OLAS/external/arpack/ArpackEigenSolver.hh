#ifndef OLAS_ARPACK_EIGENSOLVER_HH
#define OLAS_ARPACK_EIGENSOLVER_HH

#include "OLAS/solver/BaseEigenSolver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "ArpackSolver.hh"
#include "ArpackMatInterface.hh"

namespace CoupledField {
  
  class StdMatrix;
  

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

    //! Setup for a standard EVP
    void Setup(const BaseMatrix & A, bool isHermitian=false);

    //! Setup for a generalised EVP
    void Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian=false);

    //! Setup for a quadratic EVP
    void Setup(const BaseMatrix & A, const BaseMatrix & B, const BaseMatrix & C);

    //! Setup routine for standard eigenvalue problem

    //! Setup routine for various initialization tasks of a standard
    //! eigenvalue problem.
    //! \param mat Reference to matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    void Setup(const BaseMatrix & mat, UInt numFreq, double freqShift, bool sort);

    //! Setup routine for a generalized eigenvalue problem
    
    //! Setup routine for various initialization tasks.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    void Setup(const BaseMatrix & stiffMat, const BaseMatrix & massMat,
               UInt numFreq, double freqShift, bool sort, bool bloch);
    
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
    void Setup(const BaseMatrix & stiffMat, const BaseMatrix & massMat, const BaseMatrix & dampMat,
               UInt numFreq, double freqShift, bool sort );

    //! Sets the computational mode for solving the problem.
    //! Must be called before Setup() to have an effect!
    //! \param mode Computational mode
    void SetComputeMode(ArpackMatInterface::ComputeMode mode) { this->computeMode_ = mode; }

    //! Solve the linear generalized eigenvalue problem
    
    //! This method triggers the calculation of the eigenvalue problem.
    //! \param sol Vector with converged eigenvalues. The size is the number of converged evs
    //! \param err Vector with error bound of eigenvalues
    void CalcEigenFrequencies(BaseVector &sol, BaseVector &err);

    //! This method triggers the calculation of the eigenvalue problem.
    //! \param sol Vector with converged eigenvalues. The size is the number of converged evs
    //! \param err Vector with error bound of eigenvalues
    //! \param minVal unused
    //! \param maxVal unused
    void CalcEigenValues(BaseVector& sol, BaseVector& err, Double minVal, Double maxVal) {
      EXCEPTION( "Not implemented yet!" );
    };

    //! This method triggers the calculation of the eigenvalue problem.
    //! \param sol Vector with converged eigenvalues. The size is the number of converged evs
    //! \param err Vector with error bound of eigenvalues
    //! \param N number of requested eigenvalues
    //! \param shiftPoint shiftpoint
    void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Double shiftPoint);
    void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Complex shiftPoint);
    
    //! Calculate a particular eigenmode as a postprocessing solution

    //! This method may be called after the CalcEigenFrequencies() method.
    //! It calculates a given eigenmode and stores it in a supplied vector.
    //! \param modeNr Number of the (converged) eigenmode to be calculated
    //! \param mode Vector with the eigenmode
    void GetEigenMode( UInt modeNr, Vector<Complex> & mode, bool right=true );
    void GetComplexEigenMode( UInt modeNr, Vector<Complex> & mode );

    //! Calculate condition number

    //! This method calculates the condition number of the given matrix,
    //! as well as the 5 smallest and 5 largest eigenvalues
    void CalcConditionNumber( const BaseMatrix& mat,  
                              Double& condNumber, 
                              Vector<Double>& evs,
                              Vector<Double>& err );

  private:
    //! Method for generation of complex matrices from real ones
    void SetupComplexMatrices();

    /** Setup idx_. If not sort_ it is set to match 1:1. But it needs to be called! */
    void SetupIndex(unsigned int numev);

    /** for SetupIndex */
    typedef std::pair<double, unsigned int> ev_idx;

    /** for SetupIndex */
    static bool comperator(const ev_idx& one, const ev_idx& two)
    {
      return one.first < two.first;
    }

    /** print setup information */
    void ToInfo();

    //! Pointer to matrix interface
    ArpackSolver* arpackSolver_;

    //! Pointer to matrix interface
    ArpackMatInterface* interface_;
    
    //! Pointer to matrices
    const StdMatrix* matrixA_;
    const StdMatrix* matrixB_;
    const StdMatrix* matrixD_;

    //! Pointer to complex matrices
    StdMatrix *zStiff_, *zMass_, *zDamp_;

    //! Pointer to solver ovbject
    BaseSolver* solver_;

    //! Pointer to precond object
    BasePrecond* precond_;

    /** Attribute for xml paramnode of <solver> section */
    PtrParamNode xml_;
    
    //! Flag indicating if problem is generalized problem
    bool isGeneralized_;

    //! Computational mode / problem transformation used for solving
    ArpackMatInterface::ComputeMode computeMode_;
    
    //! Character string for 'which' setting of  arpack
    char* which_;

    /** this is the permutation matrix which allows sorting. Always used
     * and in the non-sorting case set to 0,1,2, ... */
    StdVector<unsigned int> idx_;
  };
}

#endif
