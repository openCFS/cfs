// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
#include "PDE/basePDE.hh"
#include "Domain/bcs.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"

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
        EXCEPTION( "GetInfoStreamPointer: Some weird shit happened!" );
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

    /// prints full data of a material
    void PrintMaterial(BaseMaterial* material);
    
    /// prints all data of a coil (e.g. current, area, magnetization, ...)
    void PrintCoil( Coil &coil, BasePDE::AnalysisType &analysistype );

    /// for combustion noise
    void WriteCombustionNoiseInfo(std::string filename, std::string cplRegion,
				  UInt sos, UInt src1, UInt src2, UInt src3, 
				  UInt src4, UInt src5, UInt src6, UInt src7);
 
    /// prints the process of a nonlinear iteration
    void WriteNonLinIter(const std::string& pdeName, const UInt iterationCounter,
                         const Double residualErr, const Double incrementalErr,
                         double etaLineSearch=0);

    /// prints the process of a mulitSequence Analysis
    void WriteMultiSequenceStep(const UInt sequenceStep, 
                                const BasePDE::AnalysisType analysis);

    // RHS-Src 
    void PrintSrcRhs( UInt node, UInt eqn, Double val);

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
  };
} // end namespace CoupledField
 
#endif
