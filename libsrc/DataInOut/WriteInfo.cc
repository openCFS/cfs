#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <General/environment.hh>
#include <Utils/tools.hh>
#include <cstdarg>
#include "WriteInfo.hh"


namespace CoupledField
{

  WriteInfo::WriteInfo (const Char * name)
  {
#ifdef TRACE
    (*trace) << "Entering WriteInfo::WriteInfo" << std::endl;
#endif

    std::string filename(name);
    filename += ".log";
  
    cfsInfo = new std::ofstream(filename.c_str());

    if (!cfsInfo) 
      Error("Can't open info-file");
  }



  WriteInfo::~WriteInfo ()
  {
#ifdef TRACE
    (*trace) << "Entering WriteInfo::~WriteInfo" << std::endl;
#endif

    if (!cfsInfo)
      delete cfsInfo;
  }
  


  void WriteInfo::PrintHeader()
  {
    std::stringstream header;
  
    header << "=======================================================================" << std::endl
	   << "=======================================================================" << std::endl
	   << "|                                                                     |" << std::endl
 	   << "|                                 CFS++                               |" << std::endl
	   << "|                     (Coupled Field Simulation ++)                   |" << std::endl
	   << "|                                                                     |" << std::endl
	   << "|                                                                     |" << std::endl
	   << "|  Version: 0.9                                                       |" << std::endl
	   << "|  Date:    30-Oct-2003                                               |" << std::endl
	   << "|                                                                     |" << std::endl
	   << "=======================================================================" << std::endl
	   << "=======================================================================" << std::endl << std::endl;

    std::cout << std::endl << header.str();
    *cfsInfo << header.str();
  }
  


  
  void WriteInfo::PrintPiezoMat(MaterialData& material)
  {
    *cfsInfo  << "FULL PIEZO DATAMATRIX OF " << material.GetMaterialName() << ":" << std::endl
	      << std::endl << *(material.GetMatrix()) << std::endl
	      << "density = " << material.GetDensity() << std::endl
	      << "damping coefficient alfa = " << material.GetDampingAlfa() << std::endl
	      << "damping coefficient beta = " << material.GetDampingBeta() << std::endl <<  std::endl;
  }
  


  void WriteInfo::PrintFluidMat(MaterialData& material)
  {
    *cfsInfo << "MATERIAL DATA OF " << material.GetMaterialName() << ":" << std::endl
	     << "compressibility: " << material.GetCompressibility() << std::endl
	     << "density: " << material.GetDensity() << std::endl
	     << "alpha: " << material.GetDampingAlfa() << std::endl
	     << "beta: " << material.GetDampingBeta() << std::endl << std::endl;
  }
  
  void WriteInfo::PrintMagMat(MaterialData& material)
  {
    Double mX, mY, mZ;
    material.GetPermMag(mX, mY, mZ);

    *cfsInfo << "MATERIAL DATA OF " << material.GetMaterialName() << ":" << std::endl
	     << "conductivity:            " << material.GetConductivity() << std::endl
	     << "permeability:            " << material.GetPermeability() << std::endl
	     << "vector of magnetiziation: (" << mX << ", " << mY << ", " << mZ <<")" 
	     << std::endl << std::endl;
  }
  

  void WriteInfo::WriteNonLinIter(Integer iterationCounter, 
				  Double residualErr, Double incrementalErr)
  {
    *cfsInfo << std::endl << " NONLINEAR ITERATION " << iterationCounter 
	     << " ==========================================" << std::endl
	     << " === Residual error          " << residualErr << std::endl
	     << " === Incremental error       " << incrementalErr << std::endl;
  }




  void WriteInfo::PrintCoil(std::string& coilDomain, struct MagEdgePDE::coilDefStruct& coilDef,  AnalysisType& analysistype_)
  {
    *cfsInfo <<  "COIL DESCRIPTION ======================================= " << myEndl
	     <<  "Coil domain " << coilDomain << std::endl
	     <<  "  parameters: current direction = " << coilDef.iDir << std::endl
	     <<  "              current value     = " << coilDef.current << std::endl
	     <<  "              coil area         = " << coilDef.coilArea << std::endl;
    
    if (coilDef.iDir > 3)
      *cfsInfo<< "              coilMidPt         = " << coilDef.coilMidPt << std::endl;
    if (analysistype_==HARMONIC)
      *cfsInfo<< "              current phase     = " << coilDef.currentPhase << std::endl;
    
    *cfsInfo << std::endl << myEndl;
  }




  void WriteInfo::PrintVec(Vector<Double>& vec)
  {
#ifdef TRACE
    (*trace) << "Entering WriteInfo::PrintVec" << std::endl;
#endif

    *cfsInfo << vec << std::endl;
  }
  

  void WriteInfo::PrintF(char * formatChar ...)
  {
#ifdef TRACE
    (*trace) << "entering WriteInfo::PrintF" << std::endl;
#endif
    const Integer maxSize = 100;
    typedef std::string::size_type ST;
    ST actPos=0;
    ST foundPos;
    char charOut[maxSize];
    std::string myStr;
    
    
    std::string formatStr(formatChar);  // conversion to type string: more convenient!
    std::string formatted;              // final output string
    
    
    va_list argList;
    va_start(argList, formatChar);   // init the argument list
    
		
    
    do
      {
	// search for actual position of %-sign
	foundPos = formatStr.find("%",actPos);
	
	// write string before %-sign into formatted string
	formatted += formatStr.substr(actPos, foundPos-actPos);

	// if not already at end of string
	if(foundPos != std::string::npos)
	  {
	    int wsPos = formatStr.find_first_of(" \t",foundPos);
	    std::string subFormatStr = formatStr.substr(foundPos, wsPos-foundPos);
	    std::cout << "foundPos " << foundPos << std::endl;
	    std::cout << "wsPos " << wsPos << std::endl
		      << "subFormatStr " << subFormatStr << myEndl;
	    
	    
	    
	    //	    char formatChar = (formatStr.substr(foundPos+1,foundPos+2)).c_str()[0];
	    char formatChar = subFormatStr[subFormatStr.length()-1];
	    
	    switch (formatChar)
	      {
	      case 'i': 
	      case 'd': 
	      case 'u': 
		int myInt;
		myInt = va_arg(argList, int);
 		sprintf(charOut, subFormatStr.c_str(), myInt);
		break;

	      case 'g': 
	      case 'G': 
	      case 'e': 
	      case 'E': 
	      case 'f': 
		double myDouble;
		myDouble = va_arg(argList, double);
 		sprintf(charOut, subFormatStr.c_str(), myDouble);
		break;
		
	      case 's':
		myStr = va_arg(argList, char *);
 		sprintf(charOut, subFormatStr.c_str(), myStr.c_str());
		break;

	      case 'c': 
		char myChar;
		myChar = va_arg(argList, int); // no char allowed!
 		sprintf(charOut, subFormatStr.c_str(), myChar);
		break;


	      default:
		std::string errStr;
		errStr = "Format character " + formatChar;
		errStr += " not yet defined!";
		Error(errStr.c_str());
		break;
	      }	    

	    formatted += charOut;
	
	    // the percent character and the format character have to be "erased"
	    actPos = foundPos+subFormatStr.length();
	  }
	
	//find() returns npos if nothing is found
      }while(foundPos != std::string::npos);

    
    *cfsInfo << formatted << std::endl << std::endl << std::flush;
    
    va_end(argList);
  }  
    
}





