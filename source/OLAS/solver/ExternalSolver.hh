#ifndef EXTERNALSOLVER_HH
#define EXTERNALSOLVER_HH

#include "OLAS/solver/BaseSolver.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

  /** This class implements a external Solver which is called through the terminal.
   * In this solver the linear system is exported and solved in another program, like python, which is callable through the terminal.
   */

  template<typename T>
  class ExternalSolver : public BaseDirectSolver {

  public:

    /** Constructor*/
    ExternalSolver( PtrParamNode solverNode, PtrParamNode olasInfo );

    /** Destructor*/
    ~ExternalSolver();

    void Setup( BaseMatrix& sysmat);

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    void Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol);

    /** Query type of this solver.
     * This method can be used to query the type of this solver. The answer
     * is encoded as a value of the enumeration data type SolverType.
     * @return External_Solver */
    SolverType GetSolverType() {
      return EXTERNAL_SOLVER;
    }

  private:

    /** Contructs the default filenames.
     * The default filenames are having the following structure:
     * JOB_fileType_MultiSequenceStepNumber_StepNumber */
    std::string ConstructFileName(std::string fileType);

    /** This function returns the time or frequency of the current step in a specified format, based on the
     * @param formatString using formaters based on the one used in @see sprintf() */
    std::string FormatedTimeFreq(const std::string& formatString);

    /** Deletes the file @param filename */
    void DeleteFile(std::string filename);

    /** If it is true the terminal output of the external solver is shown*/
    bool logging_ = false;

    /** String which executes the program of the external solver */
    std::string cmd_;

    /** Vector of the simple arguments of the xml-file as strings
     * They can contain numbers and whitespaces*/
    StdVector<std::string> args_;

    /** Vector of pointers to the argument variables specified in the xml-file under <arguments>.
     * They can be simple arguments like <arg> or special argument like <matrixFileName> or <timeFreq> and they can appear in any order.
     * This vector assures that the order of the arguments in the command is equal to the order in the xml_file*/
    StdVector<std::string*> argumentPointers_;

    /** Sting containing the matrixFileName*/
    std::string matrixFileName_;
    /** If is true the matrixFileName is updated with the default filename */
    bool useDefaultMatrixFileName_ = false;
    /** If is true the matrix will be exported. */
    bool exportMatrix_ = false;

    /** Sting containing the rhsFileName*/
    std::string rhsFileName_;
    /** If is true the rhsFileName is updated with the default filename */
    bool useDefaultRhsFileName_ = false;
    /** If is true the rhs will be exported. */
    bool exportRhs_ = false;

    /** Sting containing the solutionFileName*/
    std::string solutionFileName_;
    /** If is true the solutionFileName is updated with the default filename */
    bool useDefaultSolutionFileName_ = false;
    /** If is true the matrix will be imported. */
    bool importSolution_ = false;

    /** Sting containing the current time/frequency.
     * The formatString is read out of the xml-file if the attribute is specified and
     * specifies the format in which the time/frequency is written into the command.*/
    std::string timeFreq_;
    std::string formatString_ ="%.6f";

    /** specifies the step*/
    std::string step_;

    /** If it is true the matrix and vector files will be deleted after the problem is solved */
    bool deleteFiles_ = false;

    /** Export Timer
     * This is timer measures the time it takes to export the system matrix and rhs vector.
     */
    shared_ptr<Timer> exportTimer_;

    /** Execution Timer
     * This timer measures the time it takes for the external Solver to solve the system.
     */
    shared_ptr<Timer> executionTimer_;

    /** Import Timer
     * This timer measures the time it takes to import the solution vector of the external Solver.
     */
    shared_ptr<Timer> importTimer_;
  };
}

#endif
