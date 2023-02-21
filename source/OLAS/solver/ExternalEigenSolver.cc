#include <boost/smart_ptr/shared_ptr.hpp>
#include <DataInOut/Logging/LogConfigurator.hh>
#include <DataInOut/ParamHandling/ParamNode.hh>
#include <General/defs.hh>
#include <General/Exception.hh>
#include <MatVec/SparseOLASMatrix.hh>
#include <OLAS/algsys/SolStrategy.hh>
#include <OLAS/solver/ExternalEigenSolver.hh>
#include <OLAS/solver/generateEigensolver.hh>
#include <Utils/StdVector.hh>
#include <cassert>
#include <string>
#include <stdlib.h>
#include "generateEigensolver.hh"

namespace CoupledField {

  // Constructor
  ExternalEigenSolver::ExternalEigenSolver(shared_ptr<SolStrategy> strat,
                                            PtrParamNode xml,
                                            PtrParamNode solverList,
                                            PtrParamNode precondList,
                                            PtrParamNode eigenInfo)
    : BaseEigenSolver(strat, xml, solverList, precondList, eigenInfo){
    
    xml_ = xml; // Attribute pointer to xml object
    eigenSolverName_ = EigenSolverType::EXTERNAL;

    // get the name of the script called through the terminal
    xml_->GetValue("cmd", cmd_, ParamNode::EX);

    // determine if terminal output should be shown, default is false 
    xml_->GetValue("logging",logging_, ParamNode::PASS);

    // delete the input and output files of the external solvers, default is true
    xml_->GetValue("deleteFiles", deleteFiles_, ParamNode::PASS);

    /** Reading in all supplementary arguments from the "arguments" element in the xml-file.
     *  To guarantee that the order of appereance of the arguments in the terminal command is the same
     *  as given in the xml-file, pointers to the arguments are added and those are stored in the vector 
     *  @var argumentPointers_.*/
    if(xml_->Has("arguments"))
    {
      ParamNodeList argumentsNode = xml_->Get("arguments")->GetChildren();

      /**The size of @var args_ has to be adapted in order to avoid reallocation, when another element 
       * is added with .push_back().
       * This would render the pointers to earlier appearances of the arg element entries invalid*/
      args_.Reserve(xml->Get("arguments")->Count("arg"));

      for(auto pn : argumentsNode)
      {
        if(pn->GetName() == "min")
        {
          minValStr_ = "xmin";
          argumentPointers_.push_back(&minValStr_);
        }
        else if(pn->GetName() == "max")
        {
          maxValStr_ = "xmax";
          argumentPointers_.push_back(&maxValStr_);
        }
        else if(pn->GetName() == "shiftPoint")
        {
          shiftPointStr_ = "xshiftPoint";
          argumentPointers_.push_back(&shiftPointStr_);
          pn->GetValue("formatString", formatComplex_,ParamNode::PASS);
        }
        else if(pn->GetName() == "number")
        {
          NStr_ = "xnumber";
          argumentPointers_.push_back(&NStr_);
        }
        else if(pn->GetName() == "tolerance")
        {
          tolerance_ = pn->As<double>();

          // Save tolerance in string in scientific notation
          std::stringstream ssTol;
          ssTol << std::scientific << tolerance_;
          toleranceStr_ = ssTol.str();
          argumentPointers_.push_back(&toleranceStr_);
        }
        else if(pn->GetName() == "AFileName")
        {
          AFileName_ = pn->As<std::string>();
          argumentPointers_.push_back(&AFileName_);
          if(AFileName_ != "")
          {
            UseDefaultAFileName_ = false;
          }  
        }
        else if(pn->GetName() == "BFileName")
        {
          BFileName_ = pn->As<std::string>();
          argumentPointers_.push_back(&BFileName_);
          if(BFileName_ != "")
          {
            UseDefaultBFileName_ = false;
          }
        }
        else if(pn->GetName() == "EigenValuesFileName")
        {
          EigenValuesFileName_ = pn->As<std::string>();
          argumentPointers_.push_back(&EigenValuesFileName_);
          if(EigenValuesFileName_ != "")
          {
            UseDefaultEigenValuesFileName_ = false;
          }
        }
        else if(pn->GetName() == "EigenVectorsFileName")
        {
          EigenVectorsFileName_ = pn->As<std::string>();
          argumentPointers_.push_back(&EigenVectorsFileName_);
          if(EigenVectorsFileName_ != "")
          {
            UseDefaultEigenVectorsFileName_ = false;
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
      };
    }
    BaseEigenSolver::PostInit();
  }

  // Destructor
  ExternalEigenSolver::~ExternalEigenSolver(){
  }

  // Setup for standard EVP
  void ExternalEigenSolver::Setup(const BaseMatrix & A, bool isHermitian){
    
    this->SetProblemType(A,isHermitian); // check for symmetry
    n_ = A.GetNumRows(); // Determine size of system matrices

    // Define default filename for system matrix A
    if (UseDefaultAFileName_) 
    {
      AFileName_ = ConstructFileName("AMatrix");
    }

    // Export the system matrix A as input for external script
    A.Export(AFileName_,BaseMatrix::OutputFormat::MATRIX_MARKET);
  }

  // Setup for generalized EVP
  void ExternalEigenSolver::Setup(const BaseMatrix & A, const BaseMatrix & B, bool isHermitian){
    generalized_ = true; // set the problem type to generalized
    Setup(A,isHermitian); // call setup routine for standard EVP
    this->SetProblemType(B,isHermitian); // check for symmetry

    // Define default filename for system matrix B
    if (UseDefaultBFileName_) 
    {
      BFileName_ = ConstructFileName("BMatrix");
    }

    // Export the system matrix B as input for external script
    B.Export(BFileName_,BaseMatrix::OutputFormat::MATRIX_MARKET);
  }

  // Compute eigenvalues within given interval
  void ExternalEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err, Double minVal, Double maxVal){
    minValStr_ = std::to_string(minVal);
    maxValStr_ = std::to_string(maxVal);
    CalcEigenValues(sol, err);
  }
  // Compute eigenvalues for complex shift point
  void ExternalEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Complex shiftPoint){
    NStr_ = std::to_string(N);
    shiftPointStr_ = FormatedComplex(formatComplex_,shiftPoint); // format the complex number as string
    CalcEigenValues(sol, err);
  }

  // Compute eigenvalues for real shift point
  void ExternalEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err, UInt N, Double shiftPoint){
    NStr_ = std::to_string(N);
    shiftPointStr_ = std::to_string(shiftPoint);
    CalcEigenValues(sol, err);
  }
  
  /**Handle the call of the external eigensolver. The arguments are passed to the external script in 
   * the order specified in the "arguments" section of the XML file. Values defined in the analysis 
   * section have been passed to their corresponding storage strings in the previous overloading of the 
   * CalcEigenValues() method.   */ 
  void ExternalEigenSolver::CalcEigenValues(BaseVector& sol, BaseVector& err){
    if (UseDefaultEigenValuesFileName_)
    {
      EigenValuesFileName_ = ConstructFileName("EigenValues"); 
    }
    if (UseDefaultEigenVectorsFileName_)
    {
      EigenVectorsFileName_ = ConstructFileName("EigenVectors"); 
    }

    // string for storing the complete terminal command, including all input parameters
    std::string cmd = cmd_;

    /* assemble the command string with additional arguments by iterating 
    * over the vector of all argument pointers
    * in the order defined by the user in the xml-file. */
    for (StdVector<std::string*>::iterator it = argumentPointers_.begin(); it != argumentPointers_.end(); it++)
    { 
      std::string argumentStr = **it;
      // Check if the numerical input arguments have been defined in the analysis section
      if( (*it==&minValStr_) || (*it==&maxValStr_) || (*it==&shiftPointStr_) || (*it==&NStr_) )
      { 
        if( argumentStr.at(0) == 'x')
        {
          std::string missing = argumentStr.substr(1);
          EXCEPTION("The numerical value of the " + missing + " argument has not been defined in the analysis section");
        }
      }

      cmd += " " + argumentStr;

      // For filenames ".mtx" has to be added at the end of the argument,
      // because the @see BaseMatrix::Export() method appends ".mtx" to the filenames.
      if((*it == &AFileName_) || (*it == &BFileName_) || (*it == &EigenValuesFileName_) || (*it == &EigenVectorsFileName_))
      {
        cmd += ".mtx";
      }
    }

    // if logging is not set, the output is written to nul, so it is not shown in the terminal
    if(!logging_)
    {
      cmd+= " > nul:";
    }

    // The command string has to be cast to const char *, because std::system() supports only this data type
    const char * command = cmd.c_str();

    std::cout << "++ External Eigensolver - Executing the command '" << cmd << "'\n";
    int retValue = std::system(command); // execute command in terminal
    if( retValue != 0)
    {
      EXCEPTION("Error while calling external eigensolver, return value: " << retValue);
    }
    
    // import eigenvalues and eigenvectors
    sol.Import(EigenValuesFileName_, false);
    N_ = sol.GetSize();
    vr_.Import(EigenVectorsFileName_, false);

    // if no tolerance is specified by user, machine precision is set
    if (tolerance_ == 0)
    {
      tolerance_ = std::numeric_limits<double>::epsilon();
    }  
    err.Resize(N_);
    if(err.GetEntryType() == BaseMatrix::COMPLEX)
    {
      Complex toleranceComplex;
      toleranceComplex.real(tolerance_);
      toleranceComplex.imag(0); 
      for (UInt i=0; i < N_; i++)
      {
        err.SetEntry(i,toleranceComplex);  
      } 
    }
    else 
    {
      for (UInt i=0; i < N_; i++)
      {
        err.SetEntry(i,tolerance_);
      }
    }

    // remove intermediate files depending on deleteFiles_ value
    if(deleteFiles_)
    {
      if( remove((AFileName_ + ".mtx").c_str()) )
      {
        EXCEPTION("File '" << AFileName_ << ".mtx' could not be removed");
      }
      if( generalized_ && remove((BFileName_ + ".mtx").c_str()) )
      {
        EXCEPTION("File '" << BFileName_ << ".mtx' could not be removed");
      }  
      if( remove((EigenValuesFileName_ + ".mtx").c_str()))
      {
        EXCEPTION("File '" << EigenValuesFileName_ << ".mtx' could not be removed");
      }
      if( remove((EigenVectorsFileName_ + ".mtx").c_str()))
      {
        EXCEPTION("File '" << EigenVectorsFileName_ << ".mtx' could not be removed");
      }
    }
  }
    
  //  Construct a default filename containing the job name
  std::string ExternalEigenSolver::ConstructFileName(std::string fileType){
    std::string filename = progOpts->GetSimName() + "_" + fileType;
    return filename;
  }

  /* Format a complex number as a string in a format which can be specified by the user.
  * Default is the Python complex format. */ 
  std::string ExternalEigenSolver::FormatedComplex(const std::string& formatString, Complex number){
    char buffer [50];
    int count;
    count = snprintf(buffer, 50, formatString.c_str(),real(number), imag(number));
    if ((count >= 50)||(count < 0))
    {
      EXCEPTION("The format for complex numbers is unsuitable for the used values");
    }
    std::string complexOutput = buffer;
    return complexOutput;  
  }

  // Get the eigenvector corresponding to @var modeNr
  void ExternalEigenSolver::GetEigenMode(unsigned int modeNr, Vector<Complex>& mode, bool right){
    mode.Resize(n_); // Set the mode size corresponding to system matrix size
    assert(N_*n_ >= modeNr*n_);

    // Fill the mode vector with the corresponding value range from the stacked vector of all eigenvectors  
    for(UInt i = 0; i < n_; i++) 
    {
      mode[i] = vr_[modeNr * n_ + i];
    }
  }

}