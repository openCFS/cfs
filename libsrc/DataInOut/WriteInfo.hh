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
#include "MaterialData.hh"
#include <PDE/pdes_header.hh>
#include <General/environment.hh>

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
    void PrintCoil(std::string& coilDomain, struct coilDefStruct& coilDef,  AnalysisType& analysistype_);


    /// prints the process of a nonlinear iteration
    void WriteNonLinIter(const std::string& pdeName, const Integer iterationCounter,    
			 const Double residualErr, const Double incrementalErr, double etaLineSearch=0);

    /// prints the process of a nonlinear iteration
    void WriteTimeStep(const std::string& pdeName, const Integer timeStep,    
		       const Double time);

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
    void WriteResult(std::string pdename, std::string resulttype, std::vector<std::string> subdoms,
		     std::vector<Double> reults);

    /// just prints a vector
    void PrintVec(Vector<Double>& vec);

    /// prints a standard vector
    void PrintVec(std::vector<Integer>& vec);

    /// prints a standard vector with a leading comment
    void PrintVec(const char * comment, std::vector<Integer>& vec);

    /// prints a standard vector of strings, each string in a new line
    void PrintVec(const char * comment, std::vector<std::string>& vec);

    /// does a formatted print leaded by the PDE name, equal to std::printf(...)
    void PrintF(const std::string& pdeName, char * formatStr ...);

    /// prints warning to info-file
    void Warning(const std::string & text);
    
    /// prints error to both std::out and info-file
    void Error(const std::string & text, const Char * const filename=NULL,
               const Integer numline=0);

    //! Convert anything to a standard string

    //! This auxilliary method can convert any data type to a standard string,
    //! if the << operator has been overloaded for the respective type.
    template<class T> std::string GenStr( const T &value ) {
      std::ostringstream mystream;
      mystream << value << std::ends;
      return mystream.str();
    };

  };
} // end namespace CoupledField
 
#endif
