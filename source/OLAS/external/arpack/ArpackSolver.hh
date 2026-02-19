#ifndef OLAS_ARPACK_SOLVER_HH
#define OLAS_ARPACK_SOLVER_HH

#include <boost/shared_ptr.hpp>

#include "ArpackMatInterface.hh"
#include "MatVec/Vector.hh"
#include "arpackFortranInterface.hh"

namespace CoupledField {
  // =========================================================================
  //   ARPACK SOLVER FOR SYMMETRIC REAL SYSTEMS
  // =========================================================================
  
  //! Class for interfacing with ARPACK FORTRAN library
  // Note class ArpackEigenSolver : public BaseEigenSolver
  // ArpackSolver is the Arpack frontent
  class ArpackSolver {
    
  public:

    //! Default Constructor
    
    ArpackSolver(PtrParamNode xml);

    //! Default Destructor
    virtual ~ArpackSolver();

    //! General setup routine
    
    //! Setup routine for various initialization tasks.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param computeMode Computational mode of solver
    void Setup(ArpackMatInterface *apInterface, UInt size, char* which, char* type,
               ArpackMatInterface::ComputeMode computeMode, bool complex);

    //! Setup routine for various initialization tasks.
    //! \param stiffMat Reference to stiffness matrix
    //! \param massMat Reference to mass matrix
    //! \param computeMode Computational mode of solver
    void QuadSetup( ArpackMatInterface *apInterface, UInt size, char* which, ArpackMatInterface::ComputeMode computeMode );

    //! Solve the linear generalized eigenvalue problem
    
    //! This method triggers the calculation of the eigenvalue problem.
    //! Its return value is the number of converged eigenvalues.
    //! It is real valued for non quadratic normal EV problems and complex for bloch mode.
    //! Quadratic problems are also complex but solved by FindQuadEigenvalues!
    //! \param sol Vector with converged eigenvalues
    //! \return Number of converged eigenvalues
    template <class TYPE>
    UInt FindEigenvalues(UInt numEV, TYPE valueShift);

    //! Method triggers the calculation of the quadratic eigenvalue problem.
    //! Its return value is the number of converged complex eigenvalues.
    //! \param sol Vector with converged eigenvalues (complex)
    //! \return Number of converged eigenvalues (complex)
    template <class TYPE>
    UInt FindQuadEigenvalues(UInt numEV, TYPE valueShift);

    //! This method returns the n-th converged eigenvalue
    //! This is always the original sorting from arpack. It becomes sorted in ArpackEigenSolver.
    //! \param n Number of requested converged eigenvalue
    //! \return Calculated eigenvalue
    Double Eigenvalue(UInt n);
    Complex CmplxEigenvalue(UInt n);

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
    Complex* GetComplexEigenvector( UInt modeNr );

    //! Method to switch arpack debug information off
    void DebugOff();

    //! Method to switch arpack debug information on
    void DebugOn();

    void ToInfo(PtrParamNode info);

    /** counters for statistics */
    int counter_calll_aupd;
    int counter_solve_OP_x;
    int counter_solve_OP_B_x;
    int counter_B_x;

  private:

    /** for the template FindEigenvalues(). The implementations are concrete! */
    template <class TYPE>
    void CallAUPD(Integer* ido, char* bmat, Integer* n, char* which, Integer* nev, Double* tol,
        TYPE *resid, Integer *ncv, TYPE *V, Integer *ldv, Integer *iparam, Integer *ipntr,
        TYPE *workd, TYPE *workl, Integer *lworkl, Double *workDbleD, Integer *info);

    template <class TYPE>
    void CallEUPD(bool *rvec, char *howMny, Double *select, TYPE *d, TYPE *z, Integer *ldz,
        TYPE *shift, TYPE *zwork, char *bmat, Integer* size, char *which, Integer *nev, Double *tol,
        TYPE *resid, Integer *ncv, TYPE *V, Integer *ldv, Integer *iparam, Integer *ipntr,
        TYPE *workd, TYPE *workl, Integer *lworkl, Double  *workDbleD,Integer *info);

    void InitQuadTempSpace(Complex*, Complex*, Complex*, Complex*, Complex*, Double*);
    //! Pointer to parameter object

    //! Translate error number in meaningfull text describtion
    std::string ArpackError( Integer errNo );


    //! Pointer to matrix interface
    ArpackMatInterface* interface_;
    
    //! Computational mode / problem transformation used for solving
    ArpackMatInterface::ComputeMode computeMode_;

    //! Character string for 'which' setting of arpack
    char* which_;
    
    /** here we store the std::string with which from xml if given */
    std::string xml_which_;

    //! Character string denoting generalized or standard eigenvalue problem
    char* type_;

    //! Tolerance to be achieved in solution
    Double tolerance_;

    //! Shift to be applied. Only there as member variable for ToInfo().
    Complex valueShift_;
    //! Number of eigenvalues. Only there as member variable for ToInfo().
    int numEV_;

    //! Maximum number of iterations
    UInt maxIterations_;

    //! Number of Arnoldi vectors
    int numArnoldiVec_;

    /** Facor of Arnolid vectors w.r.t number of eigenvalues (usually 1.1 ... 2) */
    double arnoldiFactor_;

    //! Size of equation system
    int size_;

    bool logging_;

    //! stores the calculated eigenvalues (real or complex)
    SingleVector* eigenValues_;

    //! stores the calculated eigenvectors (real or complex)
    SingleVector* eigenVectors_;

    //! stores the tolerances for the calculated eigenvalues
    Vector<double> eigenTolerances_;
  };
}

#endif
