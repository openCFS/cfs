#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include "Utils/tools.hh"
#include <cstdarg>
#include "WriteInfo.hh"


namespace CoupledField
{

  WriteInfo::WriteInfo (const Char * name)
  {
    ENTER_FCN( "WriteInfo::WriteInfo" );

    std::string filename(name);
    filename += ".info";
  
    cfsInfo = new std::ofstream(filename.c_str());

    if (!cfsInfo) 
      Error("Can't open info-file");
  }



  WriteInfo::~WriteInfo ()
  {
    ENTER_FCN( "WriteInfo::~WriteInfo" );

    if (!cfsInfo)
      delete cfsInfo;
  }
  


  void WriteInfo::PrintHeader()
  {
    ENTER_FCN( "WriteInfo::PrintHeader" );

    std::stringstream header;
    std::string compileDate = __DATE__;
  
    header << "============================================================"
	   << "===========" << std::endl
	   << "============================================================"
	   << "===========" << std::endl
	   << "|                                                                     |" << std::endl
 	   << "|                                 CFS++                               |" << std::endl
	   << "|                     (Coupled Field Simulation ++)                   |" << std::endl
	   << "|                                                                     |" << std::endl
	   << "|                                                                     |" << std::endl
	   << "|  Version: 0.092                                                     |" << std::endl
	   << "|  Date:    " << compileDate
	   << "                                               |\n"
	   << "|                                                                     |" << std::endl
	   << "============================================================"
	   << "===========" << std::endl
	   << "============================================================"
	   << "===========" << std::endl << std::endl;

    std::cout << std::endl << header.str();
    *cfsInfo << header.str();
  }
  

  
  void WriteInfo::PrintPiezoMat(MaterialData& material)
  {
    ENTER_FCN( "WriteInfo::PrintPiezoMat" );

    *cfsInfo  << "FULL PIEZO DATAMATRIX OF " << material.GetMaterialName()
	      << ":" << std::endl
	      << std::endl << *(material.GetMatrix()) << std::endl
	      << "density = " << material.GetDensity() << std::endl
	      << "damping coefficient alfa = " << material.GetDampingAlfa()
	      << std::endl
	      << "damping coefficient beta = " << material.GetDampingBeta()
	      << std::endl <<  std::endl;
  }
  


  void WriteInfo::PrintFluidMat(MaterialData& material)
  {
    ENTER_FCN( "WriteInfo::PrintFluidMat" );

    *cfsInfo << "MATERIAL DATA OF " << material.GetMaterialName() << ":"
	     << std::endl << "compressibility: "
	     << material.GetCompressibility() << std::endl
	     << "density: " << material.GetDensity() << std::endl
	     << "alpha: " << material.GetDampingAlfa() << std::endl
	     << "beta: " << material.GetDampingBeta() << std::endl
	     << std::endl;
  }
  


  void WriteInfo::PrintMagMat(MaterialData& material)
  {
    ENTER_FCN( "WriteInfo::PrintMagMat" );

    Double mX, mY, mZ;
    material.GetPermMag(mX, mY, mZ);

    *cfsInfo << "MATERIAL DATA OF " << material.GetMaterialName()
	     << ":" << std::endl
	     << "conductivity:            " << material.GetConductivity()
	     << std::endl
	     << "permeability:            " << material.GetPermeability()
	     << std::endl
	     << "vector of magnetiziation: (" << mX << ", " << mY << ", "
	     << mZ <<")" 
	     << std::endl << std::endl;
  }
  


  void WriteInfo::WriteNonLinIter(const std::string& pdeName,
				  const Integer iterationCounter,
				  const Double residualErr,
				  const Double incrementalErr,
				  double etaLineSearch)
  {
    ENTER_FCN( "WriteInfo::WriteNonLinIter" );

    std::string pdeNameLong(pdeName);
    
    pdeNameLong += "-PDE: ";
    
    *cfsInfo << std::endl << pdeNameLong << "NONLINEAR ITERATION "
	     << iterationCounter 
	     << " ==========================================" << std::endl
	     << pdeNameLong << "=== Residual error          " << residualErr
	     << std::endl
	     << pdeNameLong << "=== Incremental error       "
	     << incrementalErr << std::endl;

    if (etaLineSearch)
      *cfsInfo << pdeNameLong << "=== eta (line search)       " << etaLineSearch << std::endl;
  }




  void WriteInfo::WriteTimeStep(const std::string& pdeName,
				const Integer timeStep, const Double time)
  {
    ENTER_FCN( "WriteInfo::WriteTimeStep" );

    std::string pdeNameLong(pdeName);

    // write std::out info    
    std::cout << myEndl << pdeName << ": Time step " 
	      << timeStep <<" ======================= " << std::endl;      


    *cla << myEndl << pdeName << ": Time step " 
	 << timeStep <<" ********************************************"
	 << std::endl;


    // write to info-file
    pdeNameLong += "-PDE: ";    
    *cfsInfo << std::endl << std::endl << std::endl 
	     << "**********************************************************"
	     << "**********************" 
	     << std::endl << pdeNameLong << "TIME STEP " << timeStep 
	     << ", time: " << time << std::endl;
  }
  


  void WriteInfo:: WriteResult(std::string pdename, std::string resulttype,
			       std::vector<std::string> subdoms,
			       std::vector<Double> results)
  {
    ENTER_FCN( "WriteInfo::WriteResult" );

    if (subdoms.size() != results.size())
      Error("Problem in WriteResults",__FILE__,__LINE__);
 
    *cfsInfo << std::endl << " PostProcessing Result for PDE " << pdename
	     << ": " << resulttype << " ==========" << std::endl;
    for (Integer i=0; i<subdoms.size(); i++)
      *cfsInfo << " === " << subdoms[i] << " : " << results[i] << " (Ws)"
	       << std::endl << std::endl;
  }



  void WriteInfo::PrintCoil(std::string& coilDomain,
			    struct coilDefStruct& coilDef,
			    AnalysisType& analysistype_)
  {
    ENTER_FCN( "WriteInfo::PrintCoil" );

    *cfsInfo <<  "COIL DESCRIPTION ======================================= "
	     << myEndl <<  "Coil domain: " << coilDomain << std::endl;
    if (coilDef.type == CURRENT)
      *cfsInfo << "Coil type  : current-loaded" << std::endl;
    else if (coilDef.type == VOLTAGE)
      *cfsInfo << "Coil type  : voltage-loaded" << std::endl;
    else if (coilDef.type == MEASUREMENT)
	*cfsInfo << "Coil type  : sensing coil (induced voltage)"
		 << std::endl;
    // <<  "  parameters: current direction = " << coilDef.iDir << std::endl

    *cfsInfo   <<  "              ID                = " << coilDef.ID
	       << std::endl;

    if (coilDef.type == CURRENT)
      {
	*cfsInfo << "              current value     = " << coilDef.current
		 << std::endl;
	if (analysistype_==HARMONIC)
	  *cfsInfo << "              current phase     = " << coilDef.phase
		   << std::endl;
      }

    else if (coilDef.type == VOLTAGE)
      {
	*cfsInfo << "              voltage value     = " << coilDef.voltage
		 << std::endl
		 << "              resistance        = "
		 << coilDef.resistance << std::endl;
	if (analysistype_==HARMONIC)
	  *cfsInfo << "              voltage pahse     = " << coilDef.phase << std::endl;
      }
      
    *cfsInfo << "              coil area         = " << coilDef.coilArea << std::endl;

    if (analysistype_==TRANSIENT)
	  *cfsInfo << "              time function     = " << coilDef.timefnc << std::endl;

    if (coilDef.UIfile != "--not--defined")
      {
	if (coilDef.type == CURRENT)
	  *cfsInfo << "              save voltage in: " << coilDef.UIfile << std::endl;
	else if (coilDef.type == VOLTAGE)
	  *cfsInfo << "              save current in: " << coilDef.UIfile << std::endl;	  
	else if (coilDef.type == MEASUREMENT)
	  *cfsInfo << "              save voltage in: " << coilDef.UIfile << std::endl;
      }

    if (coilDef.Lfile != "--not--defined")
      *cfsInfo << "              save inductance in: " << coilDef.Lfile << std::endl;


//     if (coilDef.iDir > 3)
//       *cfsInfo<< "              coilMidPt         = " << coilDef.coilMidPt << std::endl;
//     if (analysistype_==HARMONIC)
//       *cfsInfo<< "              current phase     = " << coilDef.currentPhase << std::endl;
    
    *cfsInfo << std::endl << myEndl;
  }



  void WriteInfo::PrintVec(Vector<Double>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );
    *cfsInfo << vec << std::endl;
  }
  


  void WriteInfo::PrintVec(std::vector<Integer>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );
    *cfsInfo << vec << std::endl;
  }
  


  void WriteInfo::PrintVec(const char * comment, std::vector<Integer>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );
    *cfsInfo << comment << myEndl << vec << myEndl << myEndl;
  }



  void WriteInfo::PrintVec(const char * comment,
			   std::vector<std::string>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );

    *cfsInfo << comment << myEndl;

    for (int i=0; i< vec.size(); i++)
      *cfsInfo << vec[i] << std::endl;

    *cfsInfo << std::endl;
  }


  
  // prints warning to info-file
  void WriteInfo::Warning(const std::string & Text)
  {
    ENTER_FCN( "WriteInfo::Warning" );
    std::cerr << "\033[31mWARNING:\033[0m " << Text << myEndl << myEndl;

    *cfsInfo << myEndl << myEndl << myEndl
	     << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	     << "!!!!!!!!!!!!!" << myEndl
	     << "                          WARNING " << myEndl
	     << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	     << "!!!!!!!!!!!!!" << myEndl
	     << "WARNING: " << Text << myEndl 
	     << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	     << "!!!!!!!!!!!!!" << myEndl
	     << myEndl << myEndl;
  }


  
    
  // prints error to both std::out and info-file
  void WriteInfo::Error(const std::string & Text,
			const Char * const filename, const Integer numline)
  {
    ENTER_FCN( "WriteInfo::Error" );
    
    std::cerr << std::endl << "\033[31mERROR:\033[0m " << Text;

#ifdef TRACE
    FcnTraceHandler::Dump();
#endif

    if (filename) 
      {
	std::cerr <<" (" << filename <<" ";
	if (numline) 
	  std::cerr << numline;
	std::cerr << ")";
	}
    std::cerr << std::endl << std::endl;
    
    

    *cfsInfo << myEndl << myEndl << myEndl
	     << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
	     << myEndl
	     << "                          ERROR " << myEndl
	     << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
	     << myEndl
	     << " ERROR: " << Text;
    
    if (filename) 
      {
	*cfsInfo <<" (" << filename <<" ";
	if (numline) 
	  *cfsInfo << numline;
	*cfsInfo << ")";
      }
    *cfsInfo << std::endl;  
    
    exit(-1);
  }
  


  void WriteInfo::WriteHomBC(const std::string& pdeName,
			     const std::string& subDom, Integer dof)
  {
    ENTER_FCN( "WriteInfo::WriteHomBC" );

    *cfsInfo << pdeName << "-PDE: Homogenous Dirichlet BC on \"" << subDom
	     << "\"";
    if (dof)
      *cfsInfo << " with DOF number " << dof;
    
    *cfsInfo << myEndl;
  }


  void WriteInfo::WriteInHomBC(const std::string& pdeName,
			       const std::string& subDom, 
			       const Double& val, const std::string & fnc,
			       const Integer& dof)
  {
    ENTER_FCN( "WriteInfo::WriteInHomBC" );

    *cfsInfo << pdeName << "-PDE: Inhomogenous Dirichlet BC on \""
	     << subDom  << "\"";
    if (dof)
      *cfsInfo << " with DOF number " << dof;
    *cfsInfo << ", value = " <<  val << ", FncName: " << fnc; 
    *cfsInfo << myEndl;
  }


  void WriteInfo::WriteLoad(const std::string& pdeName,
			    const std::string& subDom, 
			    Double value, const std::string & fnc,
			    Integer dof)
  {
    ENTER_FCN( "WriteInfo::WriteLoad" );

    *cfsInfo << pdeName << "-PDE: Loads on \"" << subDom << "\"";
    
    if (dof)
      *cfsInfo << " with DOF number " << dof;

    *cfsInfo << ", Value=" << value <<  ", FncName: " << fnc <<  myEndl
	     << myEndl;
  }
  



  void WriteInfo::PrintF(const std::string& pdeName,
			 const char * formatChar ...)
  {
    ENTER_FCN( "WriteInfo::PrintF" );
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
    
    // for classes which are not a pde, this string is ""
    if (pdeName.length())
      *cfsInfo << pdeName << "-PDE:";
    
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

    
    *cfsInfo << formatted << std::endl << std::flush;
    
    va_end(argList);
  }  
    
}
