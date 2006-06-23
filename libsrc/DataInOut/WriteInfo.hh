#ifndef FILE_WRITEINFO
#define FILE_WRITEINFO

/**************************************************************************/
/* File:   WriteInfo.hh                                                   */
/* Author: Fred Hofer                                                     */
/* Date:   04. Nov. 2003                                                  */
/*                                                                        */
/* Writes formatted output to the info-file.                              */
/**************************************************************************/


#include <sstream>
#include <ostream>
#include <fstream>
#include "Domain/bcs.hh"
#include "General/environment.hh"
#include "Matrix/matrix.hh"

namespace CoupledField {


  // Forward declarations
  class BaseMaterial;
  template <class TYPE> class Vector;
  template <class TYPE> class StdVector;

  // Forward declaration of classes
  class Coil;

  //! Class for writing formatted output to the info-file.
  class WriteInfo {

  private:

    //! file for informational output
    std::ofstream *cfsInfo;

    //! indicates, if a warning occured already
    bool warningOccured_;

    //! indicates, if there`s a progress running
    bool progressRunning_;

    //! flag for acknowledge message for progress
    //! printouts
    bool needAck_;

  public:

    /// constructor
    WriteInfo();
    
    /// destructor
    virtual ~WriteInfo();

    //! Obtain a pointer to the logfile stream
    std::ostream* GetInfoStreamPointer() {
      std::ostream *ofs = NULL;
      ofs = dynamic_cast<std::ostream*>( cfsInfo );
      if ( ofs == NULL ) {
        (*error) << "GetInfoStreamPointer: Some weird shit happened!";
        CoupledField::Error( __FILE__, __LINE__ );
      }
      return ofs;
    }

    //! Open output file

    //! Calling this method triggers the generation/opening of a file for
    //! logging status messages into. The name of the file is compiled by
    //! adding an <em>.info</em> postfix to the name of the current simulation
    //! run. The latter is obtained from the global pointer to the
    //! BaseCommandLineHandler object.
    //! \note Calling this method twice results in an error, since it will
    //!       refuse to open the file another time.
    void CreateFile();

    /// Prints header for both info-file and standard out
    void PrintHeader();

    /// prints full data of a material
    void PrintMaterial(BaseMaterial* material);
    
    /// prints all data of a coil (e.g. current, area, magnetization, ...)
    void PrintCoil( Coil &coil, AnalysisType &analysistype );

    /// prints the process of a nonlinear iteration
    void WriteNonLinIter(const std::string& pdeName, const UInt iterationCounter,
                         const Double residualErr, const Double incrementalErr,
                         double etaLineSearch=0);

    /// prints the process of a mulitSequence Analysis
    void WriteMultiSequenceStep(const UInt sequenceStep, 
                                const AnalysisType analysis);

    /// prints the process of a transient analysis
    void WriteTimeStep(const std::string& pdeName, const UInt timeStep,    
                       const Double time);

    /// prints the process of a harmonic analysis
    void WriteHarmonicStep(const std::string& pdeName, const UInt freqStep,    
                           const Double frequency);

    /// writes definition of homogeneous dirichet boundary conditions
    void WriteHomDirBC( const std::string& pdeName, HdBcList& list );

    /// writes definition of inhomogeneous dirichet boundary conditions
    void WriteInhomDirBC( const std::string& pdeName, IdBcList& list );

    /// writes definition of inhomogeneous neumann boundary conditions
    void WriteInhomNeuBC( const std::string& pdeName, InBcList& list ); 

    /// writes definition of constraints
    void WriteConstraints( const std::string& pdeName, ConstraintList& list );

    /// writes definition of loads
    void WriteLoad( const std::string& pdeName, LoadList& list );
    

    /// write Result values
    void WriteResult(std::string pdename, std::string resulttype,
                     StdVector<std::string> & subdoms,
                     Vector<Double> & results, std::string unit, 
                     std::string analysis, Double analysisVal);

    void WriteResult(std::string pdename, std::string resulttype, 
                     StdVector<std::string> & subdoms,
                     Vector<Complex> & results, std::string unit, 
                     std::string analysis, Double analysisVal);

    /// just prints a vector
    void PrintVec(Vector<Complex>& vec);

    /// just prints a vector
    void PrintVec(Vector<Double>& vec);

    /// prints a standard vector
    void PrintVec(StdVector<Integer>& vec);

    /// prints a standard vector with a leading comment
    void PrintVec(const char * comment, StdVector<Integer>& vec);

    /// prints a standard vector of strings, each string in a new line
    void PrintVec(const char * comment, StdVector<std::string>& vec);

    /// prints a matrix with a leading commnet
    void PrintMatrix(std::string & comment, const Matrix<Double> & mat);

    /// does a formatted print leaded by the PDE name, equal to std::printf(...)
    void PrintF(const std::string& pdeName, const char * formatStr ...);

    /// prints warning to info-file and std::cerr
    void Warning(const std::string & text, const Char * const filename=NULL,
                 const UInt numline=0);
    
    /// prints error to both std::out and info-file
    void Error( const std::string &text, const Char *const filename,
                const UInt numline );

    /// generates a message, that a certain action has started
    void StartProgress(const std::string &name, bool needAck = true);

    
    /// generates a message, that the last action has finished
    void FinishProgress(const bool success = true);

  };
} // end namespace CoupledField
 
#endif
