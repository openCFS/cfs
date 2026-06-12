#ifndef OLAS_FEAST_EIGENSOLVER_HH
#define OLAS_FEAST_EIGENSOLVER_HH

#include "OLAS/solver/BaseEigenSolver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  class StdMatrix;

  /** This is the BaseEigenSolver implementation based on the FEAST Eigensolver.
   * We are here based on the MKL implementation, but it should be easy to use also the original FEAST academic code */
  class FeastEigenSolver : public BaseEigenSolver {

  public:

    //! Default Constructor
    FeastEigenSolver( shared_ptr<SolStrategy> strat,
                       PtrParamNode xml,
                       PtrParamNode solverList,
                       PtrParamNode precondList,
                       PtrParamNode eigenInfo );

    //! Default Destructor
    virtual ~FeastEigenSolver();

    //! Setup for a standard EVP
    void Setup(const BaseMatrix & A, bool isHermitian=false);

    //! Setup for a generalised EVP
    void Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian=false);

    //! Setup for a quadratic EVP
    void Setup(const BaseMatrix & K, const BaseMatrix & C, const BaseMatrix & M);        

    /* Setup routine for standard eigenvalue problem
     * @see BaseEigenSolver::Setup() */
    void Setup(const BaseMatrix& mat, unsigned int numFreq, double freqShift, bool sort);

    /** Setup routine for a generalized eigenvalue problem
     * @see BaseEigenSolver::Setup() */
    void Setup( const BaseMatrix& stiffMat, const BaseMatrix& massMat,
                unsigned int numFreq, double freqShift, bool sort, bool bloch);

    /** Setup routine for a quadratic eigenvalue problem
     * @see BaseEigenSolver::Setup() */
    void Setup( const BaseMatrix& stiffMat,
                const BaseMatrix& massMat,
                const BaseMatrix& dampMat,
                unsigned int numFreq, double freqShift, bool sort );

    void CalcEigenValues(BaseVector& sol, BaseVector& err, Double minVal, Double maxVal);
    void CalcEigenValues(BaseVector& sol, BaseVector& err);
    void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Double shiftPoint){
        EXCEPTION("not implemented yet");
    }
    void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Complex shiftPoint){
        EXCEPTION("not implemented yet");
    }

    /** Solve the linear generalized eigenvalue problem
     * @see BaseEigenSolver::CalcEigenFrequencies() */
    void CalcEigenFrequencies(BaseVector& sol, BaseVector& err){
        EXCEPTION("obsolete - should be removed from interface");
    }

    /**Calculate a particular eigenmode as a postprocessing solution
     * @see BaseEigenSolver::GetEigenMode() */
    void GetEigenMode(unsigned int modeNr, Vector<Complex>& mode, bool right=true);
    void GetComplexEigenMode(unsigned int modeNr, Vector<Complex>& mode);


    /** @see BaseEigenSolver::CalcConditionNumber() */
    void CalcConditionNumber( const BaseMatrix& mat, double& condNumber, Vector<double>& evs, Vector<double>& err);

    //! Translate the info integer to a message
    std::string FeastInfo(Integer info);

  private:
    /** solve the problem, retrying CalculationAttempt on intermittent FEAST info=-3 */
    void DoCalculation(BaseVector& sol, BaseVector& err);

    /** one FEAST solve attempt as we have a random seed */
    void CalculationAttempt(BaseVector& sol, BaseVector& err);

    /** print setup information */
    void ToInfo();

    /** pointer for info xml */
    PtrParamNode eigenInfo_;

    /** Attribute for xml paramnode of <solver> section */
    PtrParamNode xml_;

    int numFreq_;

    double freqShift_;

    /** the stiffness matrix we setup for all cases */
    const StdMatrix* a_;

    /** the mass matrix for generalized and quadratic problems */
    const StdMatrix* b_;

    /** the damping matrix for quadratic problems */
    const StdMatrix* c_;

    /** feast wants a 1-based row pointer of size n+1 with nnz as last value. As OLAS is 0-based we copy in Setup() */
    StdVector<int> ia_;
    StdVector<int> ib_;
    StdVector<int> ic_;
    StdVector<int> isa_; // Stacked vector for quadratic EVP

    /** the copied column indices converted to 1-based */
    StdVector<int> ja_;
    StdVector<int> jb_;
    StdVector<int> jc_;
    StdVector<int> jsa_; // Stacked vector for quadratic EVP


    /** this is the feastinit parameter space of 128 MKL_INT */
    StdVector<int> fpm_;

    /** size of the symmetric matrix */
    int n_;

    /** number of non-zero elements for quadratic EVP 
     *  according to FEAST documentation it must be equal 
     *  to size of column index vectors in CSR format
    */
    int nnza_ = -1; //stiffness matrix
    int nnzb_ = -1; //mass matrix
    int nnzc_ = -1; //damping matrix

    /** maximum number of non-zero elements for quadratic EVP */
    int nnzmax_ = -1;

    /** number of actually computed eigenvalues */
    int m_ = -1;

    /** subspace dimension n >= m0_ >= m_ */
    int m0_ = -1;

    /** precision of the inner linear solver. 0 = double, 1 = single (legacy default) */
    int precision_ = 1;

    /** the last result value -> MKL manual page 1635 */
    int info_ = -1;

    /** the right and left eigenvectors */
    StdVector<Complex> vr_;
    StdVector<Complex> vl_;

    /** order of the polynomial EVP */
    int p_ = -1;

    /** set by Setup() */
    bool bloch_;

    /**compute only the stochastic estimate for the number of eigenvalues*/
    bool stochasticEstimate_=false;

    /**variables for using custom search contour*/
    int numP_ = -1; // number of pieces that make up contour
    StdVector<Double> z_edge_; // endpoints of each contour piece
    StdVector<Integer> t_edge_; // type of contour piece
    StdVector<Integer> n_edge_; // integration intervals per contour piece
    int numNc_ = -1; // total number of integration nodes
    StdVector<Double> integrationNodes_; // custom integration nodes
    StdVector<Double> integrationWeights_; // custom integration weights

    /**search contour parameters*/
    // set List-I parameters depending on problem type
    // I1, I2 = { Emid, r } : complex-general, complex-symmetric, real-general
    // I1, I2 = { Emin, Emax } : complex-hermitian, real-symmetric
    double I1_ [2];
    double I2_;

    /**no idea what this does*/
    char uplo_ = 'U'; // according to CFS doku U

    /** # of subspace iteration*/
    int loop_ = -1; // # of feast subspace iteration

  }; // end of class
} // end of namespace

#endif // OLAS_FEAST_EIGENSOLVER_HH
