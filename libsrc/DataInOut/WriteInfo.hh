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
#include <General/environment.hh>
#include "Matrix/matrix.hh"

namespace CoupledField
{
  // Forward declarations
  class MaterialData;
  template <class TYPE> class Vector;
  template <class TYPE> class StdVector;

  // Forward declaration of classes
  class Coil;

  //! Class for writing formatted output to the info-file.
  class WriteInfo
  {
  private:
    // file for informational output
    std::ofstream * cfsInfo;

    // indicates, if a warning occured already
    Boolean warningOccured_;

    // indicates, if there`s a progress running
    Boolean progressRunning_;

    //! flag for acknowledge message for progress
    //! printouts
    Boolean needAck_;

  public:
    /// constructor
    WriteInfo();

    
    /// destructor
    virtual ~WriteInfo();


    /// create output file
    void CreateFile(const Char *name);

    /// Prints header for both info-file and standard out
    void PrintHeader();

    
    /// prints full data of a piezo material
    void PrintPiezoMat(MaterialData& material);

    
    /// prints all fluid data
    void PrintFluidMat(MaterialData& material);


    /// prints all magnetic data
    void PrintMagMat(MaterialData& material);

    /// prints all data of a coil (e.g. current, area, magnetization, ...)
    void PrintCoil( Coil &coil, AnalysisType &analysistype );

    /// prints the process of a nonlinear iteration
    void WriteNonLinIter(const std::string& pdeName, const Integer iterationCounter,    
			 const Double residualErr, const Double incrementalErr, double etaLineSearch=0);

    /// prints the process of a mulitSequence Analysis
    void WriteMultiSequenceStep(const Integer sequenceStep, 
				const AnalysisType analysis);

    /// prints the process of a transient analysis
    void WriteTimeStep(const std::string& pdeName, const Integer timeStep,    
		       const Double time);

    /// prints the process of a harmonic analysis
    void WriteHarmonicStep(const std::string& pdeName, const Integer freqStep,    
			   const Double frequency);

    /// writes domain and dof of homogenous boundary conditions
    void WriteHomBC(const std::string& pdeName, 
		    const std::string& subDom, Integer dof=0);

    /// writes domain and dof of inhomogenous boundary conditions
    void WriteInHomBC(const std::string& pdeName,const std::string& subDom, 
		      const Double& val, const std::string & fnc, const Integer& dof);

    /// writes domain, value and dof of a load conditon
    void WriteLoad(const std::string& pdeName, const std::string& subDom, 
		   Double value, const std::string & fnc, Integer dof=0);
    

    /// write Result values
    void WriteResult(std::string pdename, std::string resulttype, StdVector<std::string> subdoms,
		     Vector<Double> results);

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
		 const Integer numline=0);
    
    /// prints error to both std::out and info-file
    void Error(const std::string & text, const Char * const filename=NULL,
               const Integer numline=0);


    /// generates a message, that a certain action has started
    void StartProgress(const std::string &name, Boolean needAck = TRUE);

    
    /// generates a message, that the last action has finished
    void FinishProgress(const Boolean success = TRUE);

    //! Convert anything to a standard string
    //! This auxilliary method can convert any data type to a standard string,
    //! if the << operator has been overloaded for the respective type.
    template<class T> std::string GenStr( const T &value ) {
      std::ostringstream mystream;
      mystream << value;
      return mystream.str();
    }

  };
} // end namespace CoupledField
 
#endif
