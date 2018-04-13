#ifndef OLAS_Phist_EIGENSOLVER_HH
#define OLAS_Phist_EIGENSOLVER_HH

#include "PhistCore.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "OLAS/solver/BaseEigenSolver.hh"

namespace CoupledField {
  
  class StdMatrix;
  class sparseMat_t; // phist matrix type

  class PhistEigenSolver : public BaseEigenSolver , PhistCore
  {
  public:
    
    //! Default Constructor
    PhistEigenSolver( shared_ptr<SolStrategy> strat,
                       PtrParamNode xml, 
                       PtrParamNode solverList,
                       PtrParamNode precondList,
                       PtrParamNode eigenInfo );
    
    //! Default Destructor
    virtual ~PhistEigenSolver();


    //! Setup routine for standard eigenvalue problem

    //! Setup routine for various initialization tasks of a standard
    //! eigenvalue problem.
    //! \param mat Reference to matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver
    //!        is used
    void Setup(const BaseMatrix & mat, UInt numFreq, double freqShift, bool sort);

    //! Setup routine for a generalized eigenvalue problem
    
    //! Setup routine for various initialization tasks.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param numFreq Number of eigenvalues/frequencies to be calculated
    //! \param freqShift Frequency shift applied to the system
    //! \param shiftMode Flag indicating if shift-and-invert mode of solver
    //!        is used
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


    //! Solve the linear generalized eigenvalue problem
    
    //! This method triggers the calculation of the eigenvalue problem.
    //! \param sol Vector with converged eigenvalues. The size is the number of converged evs
    //! \param err Vector with error bound of eigenvalues
    void CalcEigenFrequencies(BaseVector &sol, BaseVector &err);
    
    //! Calculate a particular eigenmode as a postprocessing solution

    //! This method may be called after the CalcEigenFrequencies() method.
    //! It calculates a given eigenmode and stores in a use supplied vector.
    //! \param modeNr Number of the (converged) eigenmode to be calculated
    //! \param mode Vector with the eigenmode
    void GetEigenMode(unsigned int modeNr, Vector<Complex> & mode);
    void GetComplexEigenMode(unsigned int modeNr, Vector<Complex> & mode) {
      GetEigenMode(modeNr, mode);
    }


    //! Calculate condition number

    //! This method calculates the condition number of the given matrix,
    //! as well as the 5 smallest and 5 largest eigenvalues
    void CalcConditionNumber( const BaseMatrix& mat,  
                              Double& condNumber, 
                              Vector<Double>& evs,
                              Vector<Double>& err );


    /** this structure is forwarded to (Non)SymSparseMatRowFunc as service value */
    typedef struct {
      /** either stiff, mass, or damping */
      const StdMatrix* mat = NULL;
      /** scale value, e.g. to scale the B-Mat by 1/B[0,0]. Controlled by scale_mass*/
      double scale = 1.0;
    } SparseMatRowFuncService;

  private:

    /** templated instance of the overwritten Setup() */
    template<class TYPE>
    void Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian=false);

    /** templated instance of the overwritten CalcEigenValues() */
    template<class TYPE>
    void CalcEigenValues(BaseVector &sol, BaseVector &err, unsigned int N, double shiftPoint);

    /** print setup information */
    void ToInfo();

    void SetupCommon(unsigned int numFreq, double freqShift);

    /** little helper */
    bool IsSymmetric(const BaseMatrix& cfs) const;

    template<class TYPE>
    void SaveModes(typename phist::types<TYPE>::mvec_ptr X, int nEig);

    phist_jadaOpts opts_;

    phist_comm_ptr comm_ = NULL;

    Enum<phist_EeigSort> which;
    Enum<phist_ElinSolv> linSolv;

    //static int Diag(ghost_gidx row, ghost_lidx *rowlen, ghost_gidx *col, void *val, __attribute__((unused)) void *arg);

    /** phist copy of stiffmess matrix */
    sparseMat_t* A_ = NULL;

    /** phist copy of mass matrix */
    sparseMat_t* B_ = NULL;

    /** Attribute for xml paramnode of <solver> section */
    PtrParamNode xml_;

    /** eigenvalues */
    Vector<std::complex<double> > ev_; // always complex,

    /** norms associated to ev */
    Vector<double> resNorm_;

    /** for some strange reason CFS only expects as complex ?! */
    Matrix<std::complex<double> > mode_;

    /** shall we scale the mass matrix as suggested by Jonas? */
    bool scale_mass_;

    /** in case we have scale_mass_, this is the value */
    double scale_mass_val_ = 1.0;

    /** is this a Hermitian system */
    bool hermitian_ = false;

    /** to know which type we use. Not all complex need to be hermitian! */
    bool complex_ = false;

    /** remove when switching to new interface, we then have eigenProblemType_ */
    bool sym_ = false;

    /** sorts the eigenfrequencies and sets the sort_idx_ permutation. */
    void SetupSortIdx(const StdVector<double>& freq); // TODO: move to BaseEigenSolver

    /** for SetupSortIdx */
    typedef std::pair<double, unsigned int> ev_idx;  // TODO: move to BaseEigenSolver

    /** for SetupIndex */
    static bool comperator(const ev_idx& one, const ev_idx& two)  // TODO: move to BaseEigenSolver
    {
      return one.first < two.first;
    }

    /** this is the permutation matrix which allows sorting. Always used
     * and in the non-sorting case set to 0,1,2, ... */
    StdVector<unsigned int> sort_idx_; // TODO: move to BaseEigenSolver

  };
}

#endif
