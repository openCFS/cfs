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
    virtual void Setup(const BaseMatrix & K, const BaseMatrix & C, const BaseMatrix & M){
        EXCEPTION("not yet implemented")
    };

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
    /** print setup information */
    void ToInfo();

    /** Attribute for xml paramnode of <solver> section */
    PtrParamNode xml_;

    int numFreq_;

    double freqShift_;

    /** the stiffness matrix we setup for all cases */
    //const StdMatrix* a_;

    /** the mass matrix for generalized problems */
    const StdMatrix* b_;

    /** feast wants a 1-based row pointer of size n+1 with nnz as last value. As OLAS is 0-based we copy in Setup() */
    StdVector<int> ia_;
    StdVector<int> ib_;

    /** the copied column indices converted to 1-based */
    StdVector<int> ja_;
    StdVector<int> jb_;

    /** this is the feastinit parameter space of 128 MKL_INT */
    StdVector<int> fpm_;

    /** size of the symmetric matrix */
    int n_;

    /** number of actually computed eigenvalues */
    int m_;

    /** subspace dimension n >= m0_ >= m_ */
    int m0_;

    /** the last result value -> MKL manual page 1635 */
    int info_;

    /** the right and left eigenvectors */
    StdVector<Complex> vr_;
    StdVector<Complex> vl_;

    /** set by Setup() */
    bool generalized_;

    /** set by Setup() */
    bool bloch_;

  }; // end of class
} // end of namespace

#endif // OLAS_FEAST_EIGENSOLVER_HH
