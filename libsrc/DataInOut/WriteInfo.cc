#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>

#ifdef __sgi
#include <stdarg.h>
#else
#include <cstdarg>
#endif

#include "WriteInfo.hh"

#include "Utils/tools.hh"
#include "Utils/Coil.hh"
#include "MaterialData.hh"
#include "PDE/pdes_header.hh"
#include "Utils/vector.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"


namespace CoupledField
{

  WriteInfo::WriteInfo ()
  {

    ENTER_FCN( "WriteInfo::WriteInfo");
    warningOccured_ = FALSE;
    progressRunning_ = FALSE;
    needAck_ = FALSE;
    
    cfsInfo = NULL;
  }


  WriteInfo::~WriteInfo ()
  {
    ENTER_FCN( "WriteInfo::~WriteInfo" );

    if (!cfsInfo)
      delete cfsInfo;
  }
  
  void WriteInfo::CreateFile(const Char * name)
  {
    ENTER_FCN( "WriteInfo::CreateFile" ) 

    std::string filename(name);
    filename += ".info";
    
    cfsInfo = new std::ofstream(filename.c_str());
    
    if (!cfsInfo) 
      Error("Can't open info-file");  
  }


  void WriteInfo::PrintHeader()
  {
   
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
	   << "|  Version: 0.093                                                     |" << std::endl
	   << "|  Date:    " << compileDate
	   << "                                               |\n"
	   << "|                                                                     |" << std::endl
	   << "============================================================"
	   << "===========" << std::endl
	   << "============================================================"
	   << "===========" << std::endl << std::endl;

    std::cout << std::endl << header.str();

    if (cfsInfo)
      *cfsInfo << header.str();
  }
  

  
  void WriteInfo::PrintPiezoMat(MaterialData& material)
  {
    ENTER_FCN( "WriteInfo::PrintPiezoMat" );
    
    if (cfsInfo)
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

    if (cfsInfo)
      *cfsInfo << "MATERIAL DATA OF " << material.GetMaterialName() << ":"
			   << std::endl << "compressibility: "
			   << material.GetCompressibility() << std::endl
			   << "density: " << material.GetDensity() << std::endl
			   << "alpha: " << material.GetDampingAlfa() << std::endl
			   << "beta: " << material.GetDampingBeta() << std::endl
			   << "BoverA: " << material.GetBoverA() << std::endl
			   << std::endl;
  }
  


  void WriteInfo::PrintMagMat(MaterialData& material)
  {
    ENTER_FCN( "WriteInfo::PrintMagMat" );

    Double mX, mY, mZ;
    material.GetPermMag(mX, mY, mZ);
    Double perm,cond;
    material.GetPermeability(2,2,perm);
    material.GetConductivity(2,2,cond);


    *cfsInfo << "MATERIAL DATA OF " << material.GetMaterialName()
	     << ":" << std::endl
	     << "conductivity:            " << cond
	     << std::endl
	     << "permeability:            " << perm
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
    
    if (cfsInfo) 
      {
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
  }


  void WriteInfo::WriteMultiSequenceStep(const Integer sequenceStep, 
					 const AnalysisType analysis)
  {
    std::string analysisString;

    Enum2String(analysis, analysisString);

    ENTER_FCN( "WriteInfo::WriteMultiSequenceStep" );
 
    // write std::out info 
    std::cout << myEndl 
	      << " ***************************** " << myEndl
	      << " MultiSequenceStep: " << sequenceStep << myEndl
	      << " AnalysisType:      " << analysisString << myEndl
	      << " ***************************** " << myEndl << myEndl;



    *cla <<  myEndl 
	 << " ***************************** " << myEndl
	 << " MultiSequenceStep: " << sequenceStep << myEndl
	 << " AnalysisType:      " << analysisString << myEndl
	 << " ***************************** " << myEndl << myEndl;
    


    
    if (cfsInfo)
      *cfsInfo << myEndl 
	       << myEndl<< " ***************************** " << myEndl
	       << " MultiSequenceStep: " << sequenceStep << myEndl
	       << " AnalysisType:      " << analysisString << myEndl
	       << " ***************************** " << myEndl << myEndl;
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
    
    if (cfsInfo)
      *cfsInfo << std::endl << std::endl << std::endl 
	       << "**********************************************************"
	       << "**********************" 
	       << std::endl << pdeNameLong << "TIME STEP " << timeStep 
	       << ", time: " << time << std::endl;
  }
  

  void WriteInfo::WriteHarmonicStep(const std::string& pdeName,
				    const Integer freqStep, const Double frequency)
  {
    ENTER_FCN( "WriteInfo::WriteHarmonicStep" );

    std::string pdeNameLong(pdeName);
    std::string analysis;

    params->Get( "type", analysis, "analysis" );

    // write std::out info    


    if(analysis != "paramIdent"){ // since the paramIdent driver calls iteratively this step, the output is rather disturbing than helpful

        std::cout << myEndl << pdeName << ": Harmonic step " 
          << freqStep <<" ======================= " << std::endl;      


    *cla << myEndl << pdeName << ": Harmonic step " 
	 << freqStep <<" ********************************************"
	 << std::endl;


    // write to info-file
    pdeNameLong += "-PDE: ";    
    if (cfsInfo)
      *cfsInfo << std::endl << std::endl << std::endl 
	       << "**********************************************************"
	       << "**********************" 
	       << std::endl << pdeNameLong << "HARMONIC STEP " << freqStep 
	       << ", frequency: " << frequency << std::endl;
    }
  }


  void WriteInfo:: WriteResult(std::string pdename, std::string resulttype,
			       StdVector<std::string> subdoms,
			       Vector<Double> results,
			       std::string unit,
			       std::string analysis,
			       Double analysisVal)
  {
    ENTER_FCN( "WriteInfo::WriteResult" );

    if (subdoms.GetSize() != results.GetSize())
      Error("Problem in WriteResults",__FILE__,__LINE__);
 
    if (cfsInfo) {
      *cfsInfo << std::endl << " PostProcessing Result for PDE " << pdename
               << ": " << resulttype << " ==========" << std::endl;

      for ( Integer i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "        === " << analysis << " " << analysisVal 
                 << "; " << subdoms[i] << ": " << results[i] 
                 << " "  << unit << std::endl << std::endl;
      }
  }


  void WriteInfo::PrintCoil( Coil &coil, AnalysisType &analysistype ) {

    ENTER_FCN( "WriteInfo::PrintCoil" );

    if (!cfsInfo)
      return;
    
    *cfsInfo << "COIL DESCRIPTION ======================================= "
	     << myEndl;

    // Basic coil info
    *cfsInfo << "Coil domain: "              << coil.coilName_  << std::endl;
    *cfsInfo << "Coil type: "                << coil.coilTypeS_ << std::endl;
    *cfsInfo << "Cross-Section of winding: " << coil.windingCrossSection_
	     << std::endl;

    // ID tag for 2D coils
    if ( coil.coilType_ == Coil::MEASUREMENT2D ||
	 coil.coilType_ == Coil::VOLTAGE2D ||
	 coil.coilType_ == Coil::CURRENT2D ) {
      *cfsInfo <<  "Coil ID: " << coil.id_ << std::endl;
    }

    // Special things for voltage coils
    if ( coil.coilType_ == Coil::VOLTAGE2D ||
	 coil.coilType_ == Coil::VOLTAGE3D ) {

      *cfsInfo << "Voltage value = " << coil.value_      << std::endl;
      *cfsInfo << "Resistance    = " << coil.resistance_ << std::endl;

      if ( analysistype == HARMONIC ) {
	*cfsInfo << "Voltage phase = " << coil.phase_ << std::endl;
      }
      else if ( analysistype == TRANSIENT ) {
	*cfsInfo << "File describing dynamics = " << coil.dynamicsFile_
		 << std::endl;
      }
    }

    // Special things for current coils
    if ( coil.coilType_ == Coil::CURRENT2D ||
	 coil.coilType_ == Coil::CURRENT3D ) {

      *cfsInfo << "Current value = " << coil.value_      << std::endl;

      if ( analysistype == HARMONIC ) {
	*cfsInfo << "Current phase = " << coil.phase_ << std::endl;
      }
      else if ( analysistype == TRANSIENT ) {
	*cfsInfo << "File describing dynamics = " << coil.dynamicsFile_
		 << std::endl;
      }

    }

    // Special things for measurement coils
    if ( coil.coilType_ == Coil::MEASUREMENT2D ||
	 coil.coilType_ == Coil::MEASUREMENT3D ) {

      if ( coil.saveFileL_ != "none" ) {
	*cfsInfo << " Storing inductivity in file: " << coil.saveFileL_
		 << std::endl;
      }

      if ( coil.saveFileU_ != "none" ) {
	*cfsInfo << " Storing current/voltage in file: " << coil.saveFileU_
		 << std::endl;
      }
    }

    // Special things for 3D current and voltage coils
    if ( coil.coilType_ == Coil::CURRENT3D ) {
      if ( coil.isRotational_ ) {
	*cfsInfo << "Rotational current flow specification." << std::endl;
	*cfsInfo << "Midpoint = ( "
		 << coil.midX_ << " , "
		 << coil.midY_ << " , "
		 << coil.midZ_ << " )" << std::endl;
	*cfsInfo << "Axis = ( "
		 << coil.rotAxisX_ << " , "
		 << coil.rotAxisY_ << " , "
		 << coil.rotAxisZ_ << " )" << std::endl;
      }
      else {
	*cfsInfo << "Direction of current flow: ";
	if ( coil.flowDir_ == Coil::XDIR ) {
	  *cfsInfo << "xDir" << std::endl;
	}
	if ( coil.flowDir_ == Coil::YDIR ) {
	  *cfsInfo << "yDir" << std::endl;
	}
	if ( coil.flowDir_ == Coil::ZDIR ) {
	  *cfsInfo << "zDir" << std::endl;
	}
      }
    }

    *cfsInfo << std::endl << myEndl;
  }


  void WriteInfo::PrintVec(Vector<Double>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );
    if (cfsInfo)
      *cfsInfo << vec << std::endl;
  }
  


  void WriteInfo::PrintVec(StdVector<Integer>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );
    if (cfsInfo)
      *cfsInfo << vec << std::endl;
  }
  


  void WriteInfo::PrintVec(const char * comment, StdVector<Integer>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );
    if (cfsInfo)
      *cfsInfo << comment << myEndl << vec << myEndl << myEndl;
  }



  void WriteInfo::PrintVec(const char * comment,
			   StdVector<std::string>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );

    if (cfsInfo)
      {
	*cfsInfo << comment << myEndl;
	
	for (int i=0; i< vec.GetSize(); i++)
	  *cfsInfo << vec[i] << std::endl;
	
	*cfsInfo << std::endl;
      }
  }

  void WriteInfo::PrintMatrix(std::string &comment, const Matrix<Double> &mat)
  {
    ENTER_FCN( "WriteInfo::PrintMatrix" );

    if (cfsInfo)
      *cfsInfo << comment << myEndl << mat << myEndl << myEndl;
  }

  
  // prints warning to info-file
  void WriteInfo::Warning(const std::string & Text,
			  const Char * const filename, const Integer numline)
  {
    ENTER_FCN( "WriteInfo::Warning" );

    if (progressRunning_) {
      std::cerr <<  "\033[31mWARNING\033[0m " << myEndl << myEndl;
    }

    std::cerr << "\033[31mWARNING:\033[0m " << Text << myEndl << myEndl;
    
    warningOccured_ = TRUE;

    if (filename) {
      std::cerr <<" (" << filename <<" ";
      if (numline) {
	std::cerr << numline;
      }
      std::cerr << ")" << std::endl;
    }else
      std::cerr << std::endl;

    if (cfsInfo) {
      *cfsInfo << myEndl << myEndl << myEndl
	       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	       << "!!!!!!!!!!!!!" << myEndl
	       << "                          WARNING " << myEndl
	       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	       << "!!!!!!!!!!!!!" << myEndl
	       << "WARNING: " << Text;
	
      if (filename) {
	*cfsInfo <<" (" << filename <<" ";
	if (numline) 
	  *cfsInfo << numline;
	*cfsInfo << ")";
      }
    
      *cfsInfo << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	       << "!!!!!!!!!!!!!!" << myEndl
	       << myEndl;
    }
  }
  
  
  // prints error to both std::out and info-file
  void WriteInfo::Error(const std::string & Text,
			const Char * const filename, const Integer numline)
  {
    ENTER_FCN( "WriteInfo::Error" );
    
    if (progressRunning_);
    std::cerr << "\033[31mFAILED\033[0m" << std::endl << std::endl;

    std::cerr << std::endl << "\033[31mERROR:\033[0m " << myEndl;
    std::cerr << Text;

    if (filename) {
      std::cerr <<"\n\nThe error occurred in '" << filename << "'";
      if (numline) {
	std::cerr << " on line " << numline << ".";
      }
    }

#ifdef TRACE
    OutInfo::FcnTraceHandler::Dump();
    std::cerr << "\nSee .trace-file for trace dump of function "
	      << "call tree.";
#endif

    std::cerr << std::endl << std::endl;
    
    if (cfsInfo)
      {
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
      }
    
    // exit program
    exit(-1);
  }
  


  void WriteInfo::WriteHomBC(const std::string& pdeName,
			     const std::string& subDom, Integer dof)
  {
    ENTER_FCN( "WriteInfo::WriteHomBC" );
    
    if (cfsInfo)
      {
	
	*cfsInfo << pdeName << "-PDE: Homogenous Dirichlet BC on \"" << subDom
		 << "\"";
	if (dof)
	  *cfsInfo << " with DOF number " << dof;
	
	*cfsInfo << myEndl;
      }
  }


  void WriteInfo::WriteInHomBC(const std::string& pdeName,
			       const std::string& subDom, 
			       const Double& val, const std::string & fnc,
			       const Integer& dof)
  {
    ENTER_FCN( "WriteInfo::WriteInHomBC" );
    
    if (cfsInfo)
      {
	*cfsInfo << pdeName << "-PDE: Inhomogenous Dirichlet BC on \""
		 << subDom  << "\"";
	if (dof)
	  *cfsInfo << " with DOF number " << dof;
	*cfsInfo << ", value = " <<  val << ", FncName: " << fnc; 
	*cfsInfo << myEndl;
      }
  }


  void WriteInfo::WriteLoad(const std::string& pdeName,
			    const std::string& subDom, 
			    Double value, const std::string & fnc,
			    Integer dof)
  {
    ENTER_FCN( "WriteInfo::WriteLoad" );

    if (cfsInfo)
      {
	*cfsInfo << pdeName << "-PDE: Loads on \"" << subDom << "\"";
	
	if (dof)
	  *cfsInfo << " with DOF number " << dof;
	
	*cfsInfo << ", Value=" << value <<  ", FncName: " << fnc <<  myEndl
		 << myEndl;
      }
  }
  



  void WriteInfo::PrintF( const std::string& pdeName,
			  const char * formatChar ... ) {

    ENTER_FCN( "WriteInfo::PrintF" );

    // list of supported format specifiers
    std::string supported = "idugGeEfsc";

    const Integer maxSize = 100;
    typedef std::string::size_type ST;
    ST actPos=0;
    ST foundPos;
    char charOut[maxSize];
    std::string myStr;
    
    // conversion to type string: more convenient!
    std::string formatStr(formatChar);

    // final output string
    std::string formatted;
    
    // init the argument list
    va_list argList;
    va_start(argList, formatChar);
    
    // for classes which are not a pde, this string is ""
    if ( pdeName.length() ) {
      *cfsInfo << pdeName << "-PDE: ";
    }
    
    do {

      // search for actual position of %-sign
      foundPos = formatStr.find("%",actPos);
	
      // write string before %-sign into formatted string
      formatted += formatStr.substr(actPos, foundPos-actPos);

      // if not already at end of string
      if(foundPos != std::string::npos) {
	int wsPos = formatStr.find_first_of( supported, foundPos );
	std::string subFormatStr = formatStr.substr(foundPos,wsPos-foundPos+1);
	char formatChar = subFormatStr[subFormatStr.length()-1];

	switch (formatChar) {

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
    } while(foundPos != std::string::npos);

    if (cfsInfo)
      //*cfsInfo << formatted << std::endl << std::flush;
      *cfsInfo << formatted << std::flush;
    
    va_end(argList);
  }
    
  void WriteInfo::StartProgress(const std::string &name,
				Boolean needAck)
  {
    ENTER_IFCN( "WriteInfo::StartProgress" );
   
    std::string modifiedName = name + " ...";

    needAck_ = needAck;
    
    std::cerr << "++ " << std::setw(60) << std::left << modifiedName;

    if (needAck)
      {
	warningOccured_ = FALSE;
	progressRunning_ = TRUE;
      }
    else
      std::cerr << std::endl;
  }


  void WriteInfo::FinishProgress(const Boolean success)
  {
    ENTER_IFCN( "WriteInfo::StartProgress" );
    
 
    if (!warningOccured_)
      if (success)
	std::cerr << std::setw(10) << "\033[32mOK\033[0m" << std::endl;
      else
	std::cerr << std::setw(10) << "\033[31mFAILED\033[0m" << std::endl;

    warningOccured_ = FALSE;
    progressRunning_ = FALSE;
 }
  
  
}
