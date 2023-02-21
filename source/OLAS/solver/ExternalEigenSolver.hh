// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;// ================================================================================================

#ifndef OLAS_EXTERNAL_EIGENSOLVER_HH
#define OLAS_EXTERNAL_EIGENSOLVER_HH

#include "OLAS/solver/BaseEigenSolver.hh"

#include "MatVec/StdMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include <math.h>


namespace CoupledField {
  class StdMatrix;

  /** This class implements an external Eigensolver, called through the terminal. 
   * The system matrices are exported and fed into an external program, e.g. Python script.
   * The computed Eigenvalues and Eigenvectors are imported. 
   */

  class ExternalEigenSolver : public BaseEigenSolver {
  public:
    // Constructor
    ExternalEigenSolver( shared_ptr<SolStrategy> strat,
                        PtrParamNode xml,
                        PtrParamNode solverList,
                        PtrParamNode precondList,
                        PtrParamNode eigenInfo );

    // Destructor
    ~ExternalEigenSolver();

    // Setup for a standard EVP
    void Setup(const BaseMatrix & A, bool isHermitian=false);

    // Setup for a generalised EVP
    void Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian=false);

    // Setup for a quadratic EVP
    void Setup(const BaseMatrix & K, const BaseMatrix & C, const BaseMatrix & M){
        EXCEPTION("not implemented for external eigensolver")
    };

    /** Setup routine for standard eigenvalue problem
     * @see BaseEigenSolver::Setup() */
    void Setup(const BaseMatrix& mat, unsigned int numFreq, double freqShift, bool sort){
      EXCEPTION("obsolete - should be removed from interface")
    };

    /** Setup routine for a generalized eigenvalue problem
     * @see BaseEigenSolver::Setup() */
    void Setup( const BaseMatrix& stiffMat, const BaseMatrix& massMat,
                unsigned int numFreq, double freqShift, bool sort, bool bloch){
      EXCEPTION("obsolete - should be removed from interface")
    };

    /** Setup routine for a quadratic eigenvalue problem
     * @see BaseEigenSolver::Setup() */
    void Setup( const BaseMatrix& stiffMat,
                const BaseMatrix& massMat,
                const BaseMatrix& dampMat,
                unsigned int numFreq, double freqShift, bool sort ){
        EXCEPTION("obsolete - should be removed from interface")
    };

    // Compute eigenvalues within interval
    void CalcEigenValues(BaseVector& sol, BaseVector& err, Double minVal, Double maxVal);

    // Compute N eigenvalues around real shift point
    void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Complex shiftPoint);

    // Compute N eigenvalues around complex shift point
    void CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Double shiftPoint);

    void CalcEigenValues(BaseVector& sol, BaseVector& err);

    /** Solve the linear generalized eigenvalue problem
     * @see BaseEigenSolver::CalcEigenFrequencies() */
    void CalcEigenFrequencies(BaseVector& sol, BaseVector& err){
      EXCEPTION("obsolete - should be removed from interface")
    };

    /**Calculate a particular eigenmode as a postprocessing solution
     * @see BaseEigenSolver::GetEigenMode() */
    void GetEigenMode(unsigned int modeNr, Vector<Complex>& mode, bool right=true);

    void GetComplexEigenMode(unsigned int modeNr, Vector<Complex>& mode){
      EXCEPTION("not implemented")
    };

    /** @see BaseEigenSolver::CalcConditionNumber() */
    void CalcConditionNumber( const BaseMatrix& mat, double& condNumber, Vector<double>& evs, Vector<double>& err){
        EXCEPTION("not yet implemented")
    };


  private:
    // Constructs the default filenames
    std::string ConstructFileName(std::string fileType);

    // Formats complex numbers in a given format
    std::string FormatedComplex(const std::string& formatString, Complex number);

    // string holding the name of the external script
    std::string cmd_;

    // if true, the terminal output of the external eigensolver is shown
    bool logging_ = false;

    // if true all input and output files of external eigensolver are deleted
    bool deleteFiles_ = true;

    // determines type of EVP
    bool generalized_ = false;

    // strings holding the filenames for the system matrix exports
    std::string AFileName_;
    std::string BFileName_;

    // if true, default filenames are used for the system matrix export files
    bool UseDefaultAFileName_ = true;
    bool UseDefaultBFileName_ = true;
    
    // strings holding the names for the imported result files
    std::string EigenValuesFileName_;
    std::string EigenVectorsFileName_;

    // if true, default filenames are used for the result files
    bool UseDefaultEigenValuesFileName_ = true;
    bool UseDefaultEigenVectorsFileName_ = true;

    // Vector holding the user defined additional arguments in the xml-file
    StdVector<std::string> args_;

    // Strings for storing the input information from the sequence step section
    std::string shiftPointStr_;
    std::string NStr_;
    std::string minValStr_;
    std::string maxValStr_;

    // String holding the output format for complex numbers
    std::string formatComplex_ = "%.6f%+.6fj";

    // Stores the user-defined tolerance for the solution of the EVP
    double tolerance_ = 0;
    std::string toleranceStr_;

    /** Vector of pointers to the argument variables specified in the xml-file 
     * for element "arguments". The vector assures that the order of arguments
     * in the terminal command is the same as that provided by the user in the xml-file.*/ 
    StdVector<std::string*> argumentPointers_;

    // Size of the system matrix
    UInt n_;

    // Number of eigenvalues
    UInt N_;

    // the right eigenvectors
    Vector<Complex> vr_;
  
  };
}

#endif // OLAS_EXTERNAL_EIGENSOLVER_HH
