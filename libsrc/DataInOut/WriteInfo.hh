#ifndef FILE_WRITEINFO
#define FILE_WRITEINFO

/**************************************************************************/
/* File:   WriteInfo.hh                                                   */
/* Author: Fred Hofer                                                     */
/* Date:   04. Nov. 2003                                                  */
/*                                                                        */
/* Writes formatted output to the info-file.                              */
/**************************************************************************/

#include "MaterialData.hh"
#include "PDE/magEdgePDE.hh"

namespace CoupledField
{

  //! Class for writing formatted output to the info-file.                              */

  
  class WriteInfo
  {
  private:
    // file for informational output
    std::ofstream * cfsInfo;

  public:
    /// constructor
    WriteInfo(const Char * name);
    
    /// destructor
    virtual ~WriteInfo();

    /// Prints header for both info-file and standard out
    void PrintHeader();
    
    /// prints full data of a piezo material
    void PrintPiezoMat(MaterialData& material);
    
    /// prints all fluid data
    void PrintFluidMat(MaterialData& material);

    /// prints all magnetic data
    void PrintMagMat(MaterialData& material);

    /// prints all data of a coil (e.g. current, area, magnetization, ...)
    void PrintCoil(std::string& coilDomain, struct MagEdgePDE::coilDefStruct& coilDef,  AnalysisType& analysistype_);

    /// prints the process of a nonlinear iteration
    void WriteNonLinIter(Integer iterationCounter, 
			 Double residualErr, Double incrementalErr);

    /// just prints a vector
    void PrintVec(Vector<Double>& vec);
    
    /// does a formatted print, equal to std::printf(...)
    void PrintF(char * formatStr ...);
    
  };
} // end namespace CoupledField
 
#endif
