#include "MatVec/CRS_Matrix.hh"
#include "OLAS/solver/ExternalSolver.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/AnalysisID.hh"
#include "Utils/Timer.hh"
#include "Utils/StdVector.hh"

#include <stdlib.h>

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  ExternalSolver<T>::ExternalSolver( PtrParamNode solverNode, PtrParamNode olasInfo ) {

    // Set pointers to communication objects
    xml_ = solverNode;
    infoNode_ = olasInfo->Get("externalSolver");

    // Reading information from the xml File

    // get command (required)
    xml_->GetValue("cmd", cmd_, ParamNode::EX);

    // Reading in all the arguments and saving them and adding pointers to @var argumentPointers_ in order to save the order of the arguments,
    // as they can appear in any order in the xml-file.
    // The simple <arg> are saved in a vector of strings, @var args_.
    // Special arguments like filenames, step and timeFreq are class variables.
    if(xml_->Has("arguments"))
    {
      ParamNodeList argumentsNode = xml_->Get("arguments")->GetChildren();
      // The size of @var args_ has to be adapted in order to avoid reallocation, which makes the pointers to
      // earlier arguments in the vector invalid, when another element is added with .push_back().
      args_.Reserve(xml_->Get("arguments")->Count("arg"));
      argumentPointers_.Reserve(argumentsNode.GetSize());
      for(auto pn : argumentsNode)
      {
        if((pn->GetName() == "timeFreq"))
        {
          argumentPointers_.push_back(&timeFreq_);
          pn->GetValue("formatString",formatString_,ParamNode::PASS);
        }
        else if((pn->GetName() == "step"))
        {
          argumentPointers_.push_back(&step_);
        }
        else if(pn->GetName() == "matrixFileName")
        {
          matrixFileName_ = pn->As<std::string>();
          argumentPointers_.push_back(&matrixFileName_);
          if(matrixFileName_ == "")
          {
            useDefaultMatrixFileName_ = true;
          }
        }
        else if(pn->GetName() == "rhsFileName")
        {
          rhsFileName_ = pn->As<std::string>();
          argumentPointers_.push_back(&rhsFileName_);
          if(rhsFileName_ == "")
          {
            useDefaultRhsFileName_ = true;
          }
        }
        else if(pn->GetName() == "solutionFileName")
        {
          solutionFileName_ = pn->As<std::string>();
          argumentPointers_.push_back(&solutionFileName_);
          if(solutionFileName_ == "")
          {
            useDefaultSolutionFileName_ = true;
          }
        }
        else if(pn->GetName() == "arg")
        {
          args_.push_back(pn->As<std::string>());
          argumentPointers_.push_back(&args_[args_.size()-1]);
        }
        else
        {
          EXCEPTION("argument tag '"<< pn->GetName() <<"' not implemented!")
        }
      }
    }

    // Check if the filename tags appear in the xml-File, if not, they will not be exported/imported
    exportMatrix_ = xml_->Has("arguments/matrixFileName");
    exportRhs_ = xml_->Has("arguments/rhsFileName");
    importSolution_ = xml_->Has("arguments/solutionFileName");

    // check if the terminal output of the external solver should be shown, the default value is false
    xml_->GetValue("logging",logging_,ParamNode::PASS);

    // check if the files should be deleted after solving the problem, the default value is true
    xml_->GetValue("deleteFiles",deleteFiles_,ParamNode::PASS);

    // Setup the Timers
    // The timers are subtimers because the "solve_externalSolver" timer is
    // already measuring the total time the solver needs.
    exportTimer_ =  infoNode_->Get(ParamNode::SUMMARY)->Get("export_files/timer")->AsTimer();
    exportTimer_->SetSub();

    executionTimer_ =  infoNode_->Get(ParamNode::SUMMARY)->Get("execution/timer")->AsTimer();
    executionTimer_->SetSub();

    importTimer_ = infoNode_->Get(ParamNode::SUMMARY)->Get("import_file/timer")->AsTimer();
    importTimer_->SetSub();
  }

  // **************
  //   Destructor
  // **************
  template<typename T>
  ExternalSolver<T>::~ExternalSolver() {
    // Deleting the matrix file of the last step
    if(deleteFiles_ && exportMatrix_)
      DeleteFile(matrixFileName_);
  }

  // *********
  //   Setup
  // *********
  template<typename T>
  void ExternalSolver<T>::Setup( BaseMatrix& sysmat) {
    // Getting the analysis id because it contains information about the current step.
    AnalysisID& id = domain->GetDriver()->GetAnalysisId();

    // Getting the current time/step
    if(id.step != -1)
      step_ = std::to_string(id.step);

    if(exportMatrix_)
    {
      // Before creating the default filenames and export the system matrix, the matrix-file from the previous step has to be deleted,
      // except it is the first step. The matrix-file of the last step will be deleted in the destructor. This has to be done here
      // because if default filenames are used, the matrixFileName_ will change in ConstructFileName.
      if(deleteFiles_ && id.step > 1)
      {
        // Create the matrix default name and export the matrix, for the naming schema @see ConstructFileName.
        // The default names of the rhs and solution are created, if necessary, in the @see Solve() method
        // because, if the same system is solved for different right hand sides the matrix only need to be named and exported once.
        DeleteFile(matrixFileName_);
      }
      if(useDefaultMatrixFileName_)
        matrixFileName_ = ConstructFileName("mat");
      exportTimer_->Start();
      sysmat.Export(matrixFileName_,BaseMatrix::OutputFormat::MATRIX_MARKET);
      exportTimer_->Stop();
    }
  }

  // *********
  //   Solve
  // *********
  template<typename T>
  void ExternalSolver<T>::Solve( const BaseMatrix &sysmat, const BaseVector &rhs, BaseVector &sol) {
    if(exportRhs_)
    {
      // Create the rhs and solution default names and export the rhs.
      if(useDefaultRhsFileName_)
        rhsFileName_ = ConstructFileName("rhs");
      if(useDefaultSolutionFileName_)
        solutionFileName_ = ConstructFileName("sol");
      exportTimer_->Start();
      rhs.Export(rhsFileName_,BaseMatrix::OutputFormat::MATRIX_MARKET);
      exportTimer_->Stop();
    }
    timeFreq_ = FormatedTimeFreq(formatString_);

    // Assembling the command string
    std::string cmd = cmd_;
    for (auto pn : argumentPointers_)
    {
      cmd += " " + *(pn);
      // If a filename is appended to cmd, ".mtx" has to be added after the filename,
      // because the @see BaseMatrix::Export() method append ".mtx" to the filenames.
      if((pn == &matrixFileName_) || (pn == &rhsFileName_) || (pn == &solutionFileName_))
      {
        cmd += ".mtx";
      }
    }
    // If logging is false the output is written to nul, so it doesn't appear in the terminal
    if(!logging_)
      cmd+= " > nul:";
    // The command-string has to be cast to const char*, because std::system() can only handle this data type.
    const char * command = cmd.c_str();
    std::cout << "++ External Solver - Executing the command: '" << cmd << "'" << std::endl;
    executionTimer_->Start();
    int retValue = std::system(command);
    executionTimer_->Stop();
    if( retValue != 0)
      EXCEPTION("ERROR at the external Program, called as External Solver, return Value: " << retValue);

    if(importSolution_)
    {
      // Reading in the Solution
      importTimer_->Start();
      sol.Import(solutionFileName_, true);
      importTimer_->Stop();
    }

    // Deleting the rhs- and solution-files
    if(deleteFiles_ & exportRhs_)
      DeleteFile(rhsFileName_);
    if(deleteFiles_ & importSolution_)
      DeleteFile(solutionFileName_);
  }

  template<typename T>
  std::string ExternalSolver<T>::ConstructFileName(std::string fileType)
  {
    AnalysisID& id = domain->GetDriver()->GetAnalysisId();
    std::string sequenceStep = std::to_string(domain->GetDriver()->GetActSequenceStep());
    std::string filename = progOpts->GetSimName() + "_" + fileType + "_" + sequenceStep;
    if(id.step >= 0)
    {
      filename += "_" + std::to_string(id.step);
    }
    return filename;
  }

  template<typename T>
  void ExternalSolver<T>::DeleteFile(std::string filename)
  {
    if(remove((filename + ".mtx").c_str()))
      EXCEPTION("File '" << filename + ".mtx" +"' could not been removed");
  }

  template<typename T>
  std::string ExternalSolver<T>::FormatedTimeFreq(const std::string& formatString){
    // Getting the analysis id because it contains information about the current time, frequency,
    AnalysisID& id = domain->GetDriver()->GetAnalysisId();
    char buffer[30];
    if(id.time != -1)
    {
      std::snprintf(buffer, 20, formatString.c_str(),id.time);
    }
    else if(id.freq != -1)
    {
      std::snprintf(buffer, 20, formatString.c_str(),id.freq);
    }
    else
    {
      std::snprintf(buffer, 20, " ");
    }
    std::string formatedTimeFreq = buffer;
    return formatedTimeFreq;
  }

  // Explicit template instantiation
  template class ExternalSolver<Double>;
  template class ExternalSolver<Complex>;

}
