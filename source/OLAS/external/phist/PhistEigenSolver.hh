#ifndef OLAS_Phist_EIGENSOLVER_HH
#define OLAS_Phist_EIGENSOLVER_HH

#include "OLAS/solver/BaseEigenSolver.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/Matrix.hh"
#include "phist_enums.h"
#include "phist_jadaOpts.h"
#include "phist_void_aliases.h"

namespace CoupledField {
  
  class StdMatrix;
  class sparseMat_t; // phist matrix type

  class PhistEigenSolver : public BaseEigenSolver
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
    //! \param mode Vector with the eignmode
    void GetEigenMode( UInt modeNr, Vector<Complex> & mode );
    void GetComplexEigenMode( UInt modeNr, Vector<Complex> & mode );


    //! Calculate condition number

    //! This method calculates the condition number of the given matrix,
    //! as well as the 5 smallest and 5 largest eigenvalues
    void CalcConditionNumber( const BaseMatrix& mat,  
                              Double& condNumber, 
                              Vector<Double>& evs,
                              Vector<Double>& err );

  private:
    /** print setup information */
    void ToInfo();

    void SetupCommon(bool sym, unsigned int numFreq, double freqShift, bool sort, bool bloch);

    /** provide A_ or B_ not as pointer value but as pointer itself as the pointer value is NULL prior init
     @return the value we set phist to (redundant to phist** */
    sparseMat_t* InitMatrix(const BaseMatrix& cfs, sparseMat_t** phist);

    /** little helper */
    bool IsSymmetric(const BaseMatrix& cfs) const;

    phist_jadaOpts opts_;

    phist_comm_ptr comm_ = NULL;

    Enum<phist_EeigSort> which;
    Enum<phist_ElinSolv> linSolv;

    //static int Diag(ghost_gidx row, ghost_lidx *rowlen, ghost_gidx *col, void *val, __attribute__((unused)) void *arg);

    /** phist copy of stiffmess matrix */
    sparseMat_t* A_;

    /** phist copy of mass matrix */
    sparseMat_t* B_;

    /** Attribute for xml paramnode of <solver> section */
    PtrParamNode xml_;

    /** eigenvalues */
    StdVector<std::complex<double> > ev_; // always complex,

    /** norms associated to ev */
    StdVector<double> resNorm_;

    Matrix<double> mode_;


    // we do not use solver and preconditioners from CFS for Phist
  };
}

#endif
