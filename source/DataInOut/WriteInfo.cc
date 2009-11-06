// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdarg>

#include <def_cfs_stats.hh>

#include "WriteInfo.hh"
#include "Utils/tools.hh"
#include "Utils/Coil.hh"
#include "Utils/coordSystem.hh"
#include "Materials/baseMaterial.hh"
#include "coloredConsole.hh"
#ifndef INTEGLIB
#include "Domain/resultInfo.hh"
#include "DataInOut/programOptions.hh"
#endif
#include "MatVec/vector.hh"

#ifdef USE_SCRIPTING 
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

#define PROGRESS_TEXT_WIDTH 62

using std::cout;
using std::endl;

namespace CoupledField {


  WriteInfo::WriteInfo () {

    warningOccured_ = false;
    progressRunning_ = false;
    needAck_ = false;

    cfsInfo = NULL;
  }


  // **************
  //   Destructor
  // **************
  WriteInfo::~WriteInfo () {
    delete cfsInfo;
  }


#ifndef INTEGLIB


  // **************
  //   CreateFile
  // **************
  void WriteInfo::CreateFile() {


    // Check if a file is already open
    if ( cfsInfo != NULL ) {
      EXCEPTION("WriteInfo::CreateFile: cfsInfo already points to a "
               << "file! Cowardly refusing to open another one\n");
    }

    std::string filename = progOpts->GetSimName() + ".info";
    cfsInfo = new std::ofstream(filename.c_str());

    // Check if everything went fine
    if ( cfsInfo == NULL ) {
      EXCEPTION("WriteInfo::CreateFile: Failed to open file '"
               << filename << "' for writing status messages\n");
    }
  }


 // *****************
  //   Print Material Data
  // *****************
  void WriteInfo::PrintMaterial( BaseMaterial* material ) {


    *cfsInfo << *material << endl;
  }


  void WriteInfo::WriteNonLinIter(const std::string& pdeName,
                                  const UInt iterationCounter,
                                  const Double residualErr,
                                  const Double incrementalErr,
                                  double etaLineSearch)
  {

    std::string pdeNameLong(pdeName);
    
    pdeNameLong += "-PDE: ";
    
    if (cfsInfo) 
      {
        *cfsInfo << endl << pdeNameLong << "NONLINEAR ITERATION "
                 << iterationCounter 
                 << " ==========================================\n"
                 << pdeNameLong << "=== Residual error          "
                 << residualErr
                 << endl
                 << pdeNameLong << "=== Incremental error       "
                 << incrementalErr << endl;
    
        if (etaLineSearch)
          *cfsInfo << pdeNameLong << "=== eta (line search)       "
                   << etaLineSearch << endl;
      }
  }

  void WriteInfo::WriteNonLinIter2(const std::string& pdeName,
                                  const UInt iterationCounter,
                                  const Double residualErr,
                                  const Double incrementalErr,
                                  const Double energyErr,
                                  double etaLineSearch)
  {

    std::string pdeNameLong(pdeName);

    pdeNameLong += "-PDE: ";

    if (cfsInfo)
      {
        *cfsInfo << endl << pdeNameLong << "NONLINEAR ITERATION "
                 << iterationCounter
                 << " ==========================================\n"
                 << pdeNameLong << "=== Residual error          "
                 << residualErr
                 << endl
                 << pdeNameLong << "=== Incremental error       "
                 << incrementalErr << endl
                 << pdeNameLong << "=== Energy error       "
                 << energyErr << endl;

        if (etaLineSearch)
          *cfsInfo << pdeNameLong << "=== eta (line search)       "
                   << etaLineSearch << endl;
      }
  }

  void WriteInfo::WriteMultiSequenceStep(const UInt sequenceStep, 
                                         const BasePDE::AnalysisType analysis)
  {
    std::string analysisString = BasePDE::analysisType.ToString(analysis);

 
    // write std::out info 
    cout << endl 
              << " ***************************** " << endl
              << " MultiSequenceStep: " << sequenceStep << endl
              << " AnalysisType:      " << analysisString << endl
              << " ***************************** " << endl << endl;



    *cla <<  endl 
         << " ***************************** " << endl
         << " MultiSequenceStep: " << sequenceStep << endl
         << " AnalysisType:      " << analysisString << endl
         << " ***************************** " << endl << endl;
    


    
    if (cfsInfo)
      *cfsInfo << endl 
               << endl<< " ***************************** " << endl
               << " MultiSequenceStep: " << sequenceStep << endl
               << " AnalysisType:      " << analysisString << endl
               << " ***************************** " << endl << endl;
  }
  
  

  // ***************
  //   WriteResult
  // ***************
  void WriteInfo::WriteResult(std::string pdename, std::string resulttype,
                               StdVector<std::string> & subdoms,
                               Vector<Double> & results,
                               std::string unit,
                               std::string analysis,
                               Double analysisVal)
  {

    if (subdoms.GetSize() != results.GetSize())
      Error("Problem in WriteResults",__FILE__,__LINE__);
 
    if (cfsInfo) {
      *cfsInfo << endl << " PostProcessing Result for PDE " << pdename
               << ": " << resulttype << " ==========" << endl;

      for ( UInt i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "        === " << analysis << " " << analysisVal 
                 << "; " << subdoms[i] << ": " 
                 << results[i] 
                 << " "  << unit << endl << endl;
      }
    }
  }

  void WriteInfo::WriteResult(std::string pdename, std::string resulttype,
                               StdVector<std::string> & subdoms,
                               Vector<Complex> & results,
                               std::string unit,
                               std::string analysis,
                               Double analysisVal)
  {

    if (subdoms.GetSize() != results.GetSize())
      Error("Problem in WriteResults",__FILE__,__LINE__);
 
    if (cfsInfo) {
      *cfsInfo << endl << " PostProcessing Result for PDE " << pdename
               << ": " << resulttype << " ==========" << endl;

      for ( UInt i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "        === " << analysis << " " << analysisVal 
                 << "; " << subdoms[i] << ": "
                 << results[i] 
                 << " "  << unit << endl << endl;
      }
    }
  }


  void WriteInfo::PrintCoil( Coil &coil, BasePDE::AnalysisType &analysistype ) {


    if (!cfsInfo)
      return;
    
    *cfsInfo << "COIL DESCRIPTION ======================================= "
             << endl;

    // Basic coil info
    *cfsInfo << "Coil domain: "              << coil.coilRegionName_ << '\n'
             << "Coil type: "                << coil.coilTypeS_ << '\n'
             << "Cross-Section of winding: " << coil.windingCrossSection_
             << endl;

    // ID tag for 2D coils
    if ( coil.coilType_ == Coil::MEASUREMENT2D ||
         coil.coilType_ == Coil::VOLTAGE2D ||
         coil.coilType_ == Coil::CURRENT2D ) {
      *cfsInfo <<  "Coil ID: " << coil.id_ << endl;
    }

    // Special things for voltage coils
    if ( coil.coilType_ == Coil::VOLTAGE2D ||
         coil.coilType_ == Coil::VOLTAGE3D ) {

      *cfsInfo << "Voltage value = " << coil.value_      << endl;
      *cfsInfo << "Resistance    = " << coil.resistance_ << endl;

      if ( analysistype == BasePDE::HARMONIC) {
        *cfsInfo << "Voltage phase = " << coil.phase_ << endl;
      }
    }

    // Special things for current coils
    if ( coil.coilType_ == Coil::CURRENT2D ||
         coil.coilType_ == Coil::CURRENT3D ) {

      *cfsInfo << "Current value = " << coil.value_      << endl;

      if ( analysistype == BasePDE::HARMONIC) {
        *cfsInfo << "Current phase = " << coil.phase_ << endl;
      }

    }

    // Special things for measurement coils
    if ( coil.coilType_ == Coil::MEASUREMENT2D ||
         coil.coilType_ == Coil::MEASUREMENT3D ) {

      if ( coil.saveFileL_ != "none" ) {
        *cfsInfo << " Storing inductivity in file: " << coil.saveFileL_
                 << endl;
      }

      if ( coil.saveFileU_ != "none" ) {
        *cfsInfo << " Storing current/voltage in file: " << coil.saveFileU_
                 << endl;
      }
    }

    // Special things for 3D current and voltage coils
    if ( coil.coilType_ == Coil::CURRENT3D ) {
      if ( coil.isRotational_ ) {
        *cfsInfo << "Rotational current flow specification." << endl;
        *cfsInfo << "Midpoint = ( "
                 << coil.midX_ << " , "
                 << coil.midY_ << " , "
                 << coil.midZ_ << " )" << endl;
        *cfsInfo << "Axis = ( "
                 << coil.rotAxisX_ << " , "
                 << coil.rotAxisY_ << " , "
                 << coil.rotAxisZ_ << " )" << endl;
      }
      else {
        *cfsInfo << "Direction of current flow: ";
        *cfsInfo << coil.locFlowDir_.ToString() << endl;
        *cfsInfo << " in coordinate system " 
                 << coil.flowCoordSys_->GetName();
      }
    }

    *cfsInfo << endl << endl;
  }

  template <class TYPE>
  void WriteInfo::WriteAcouPower(std::string pdename, 
                                 StdVector<std::string> & subdoms,
                                 Vector<TYPE>& power)
  {
 
    if (cfsInfo) {
      *cfsInfo << endl << " PostProcessing Result for PDE " << pdename
               << ": " << " ==========" << endl;
      *cfsInfo << "   Acoustic Power: \n";
      for ( UInt i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "    Subdomain: " <<  subdoms[i] << " : " <<  power[i] 
                 << " W" << endl;
      }
    }
  }


  void WriteInfo::PrintVec(Vector<Complex>& vec)
  {
    if (cfsInfo)
      *cfsInfo << vec << endl;
  }

  void WriteInfo::PrintVec(Vector<Double>& vec)
  {
    if (cfsInfo)
      *cfsInfo << vec << endl;
  }
  

  void WriteInfo::WriteCombustionNoiseInfo(std::string filename, std::string cplRegion,
                                           UInt sosIdx, UInt src1, UInt src2, UInt src3, 
                                           UInt src4, UInt src5, UInt src6, UInt src7) {


    *cfsInfo << "\nCombustion Noise Info:\n" 
	     << " Name of file: " << filename  
	     << "\n Coupling region: " << cplRegion <<  endl;

    if ( sosIdx > 0 ) {      
      *cfsInfo << " Use variable speed of sound: yes" << endl;
    }
    else {
      *cfsInfo << " Use variable speed of sound: no" << endl;
    }

    *cfsInfo << "\n Use the following source terms" << endl;

    if ( src1 > 0 ) {
      *cfsInfo << "      Reynolds stress tensor " << endl;
    }
    if ( src2 > 0 ) {
      *cfsInfo << "      Fluctuation of momentum flux " << endl;
    }
    if ( src3 > 0 ) {
      *cfsInfo << "      Unsteady reaction rates " << endl;
    }
    if ( src4 > 0 ) {
      *cfsInfo << "      2nd order derivative of density " << endl;
    }
    if ( src5 > 0 ) {
      *cfsInfo << "      Heat release " << endl;
    }    
    if ( src6 > 0 ) {
      *cfsInfo << "      Gas Const. 2nd time derivative " << endl;
    }
    if ( src7 > 0 ) {
      *cfsInfo << "      Shear term " << endl;
    }

    *cfsInfo << endl;
  }


  void WriteInfo::PrintVec(StdVector<Integer>& vec)
  {
    if (cfsInfo)
      *cfsInfo << vec << endl;
  }
  


  void WriteInfo::PrintVec(const char * comment, StdVector<Integer>& vec)
  {
    if (cfsInfo)
      *cfsInfo << comment << endl << vec << endl << endl;
  }



  void WriteInfo::PrintVec(const char * comment,
                           StdVector<std::string>& vec)
  {

    if (cfsInfo)
      {
        *cfsInfo << comment << endl;
        
        for (UInt i=0; i< vec.GetSize(); i++)
          *cfsInfo << vec[i] << endl;
        
        *cfsInfo << endl;
      }
  }

  void WriteInfo::PrintMatrix(std::string &comment, const Matrix<Double> &mat)
  {

    if (cfsInfo)
      *cfsInfo << comment << endl << mat << endl << endl;
  }

#endif //INTEGLIB

  // prints warning to info-file
  void WriteInfo::Warning( const std::string & Text,
                           const char* const filename, const UInt numline ) {

#ifdef INTEGLIB
    std::cerr << "INTEGLIB WARNING: " << Text << endl;
#else
#ifdef USE_SCRIPTING
    if ( messenger ) {
      if ( messenger->IsEvaluating() ) {
        messenger->Warning( Text.c_str(), filename, numline );
      }
    }
#endif

    if ( progressRunning_ == true ) {
      cout << endl;
    }

    std::cerr << "\n "
              << fg_yellow << "WARNING:" << fg_reset << "\n "
              << Text << endl;

    warningOccured_ = true;

    if (filename) {
      std::cerr <<" (" << filename <<" ";
      if (numline) {
        std::cerr << numline;
      }
      std::cerr << ")" << endl;
    }
    else {
      std::cerr << endl;
    }

    if (cfsInfo) {
      *cfsInfo << endl << endl << endl
               << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               << "!!!!!!!!!!!!!" << endl
               << "                          WARNING " << endl
               << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               << "!!!!!!!!!!!!!" << endl
               << "WARNING: " << Text;
        
      if (filename) {
        *cfsInfo <<" (" << filename <<" ";
        if (numline) 
          *cfsInfo << numline;
        *cfsInfo << ")";
      }
    
      *cfsInfo << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               << "!!!!!!!!!!!!!!!" << endl
               << endl;
    }
#endif // INTEGLIB
  }


  // *********
  //   Error
  // *********
  void WriteInfo::Error( const std::string &Text, const char *const filename,
                         const UInt numline ) {


#ifdef INTEGLIB
    std::cerr << "INTEGLIB ERROR: " << Text << endl;
#else
#ifdef USE_SCRIPTING
    if ( messenger ) {
      if ( messenger->IsEvaluating() ) {
        messenger->Error( Text.c_str(), filename, numline );
      }
    }
#endif
    

    // If a progress part is still there, then finish it with
    // a failure
    if ( progressRunning_ ) {
      FinishProgress( false );
    }

    std::cerr << endl << " "
              << fg_red << "ERROR:"
              << fg_reset << "\n" << endl;

    if ( Text != "" ) {
      std::cerr << ' ' << Text;
    }
    else {
      std::cerr << ' ' << "I've got the feeling that something is wrong!\n"
                << " Can't say what, however :(\n\n"
                << " (No error message was provided by the"
                << " programmer)";
    }
    if ( filename ) {
      std::cerr << "\n\n This error message was brought to you by\n "
                << filename << ", line " << numline;
    }

    std::cerr << endl << endl;
    
    if (cfsInfo) {
      *cfsInfo << endl << endl << endl
               << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
               << endl
               << "                          ERROR " << endl
               << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
               << endl
               << " ERROR: " << Text;

      if (filename) {
        *cfsInfo <<" (" << filename <<" ";
        if (numline) 
          *cfsInfo << numline;
        *cfsInfo << ")";
      }
      *cfsInfo << endl;  
    }
#endif // INTEGLIB

    // exit program
    exit(-1);

    //     char* a = NULL;
    //     std::cerr << a[0] << a[1] << a[2] << a[3];
  }


#ifndef INTEGLIB



  void WriteInfo::PrintF( const std::string& pdeName,
                          const char * formatChar ... ) {


    // list of supported format specifiers
    std::string supported = "idugGeEfsc";

    const UInt maxSize = 100;
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
    std::va_list argList;
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
        std::string subFormatStr= formatStr.substr(foundPos,wsPos-foundPos+1);
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
          EXCEPTION( "Format character " << formatChar
                   << " not yet defined!" );
          break;
        }

        formatted += charOut;
        
        // the percent character and the format character have to be "erased"
        actPos = foundPos+subFormatStr.length();
      }
        
      //find() returns npos if nothing is found
    } while(foundPos != std::string::npos);

    if (cfsInfo)
      //*cfsInfo << formatted << endl << std::flush;
    *cfsInfo << formatted << std::flush;
    
    va_end(argList);
  }


  // *****************
  //   StartProgress
  // *****************
  void WriteInfo::StartProgress( const std::string &name, bool needAck ) {

   
    std::string modifiedName = name + " ";

    needAck_ = needAck;

    cout << "++ " << std::setw(PROGRESS_TEXT_WIDTH)
              << std::setfill( '.' ) << std::left
              << modifiedName << std::setfill( ' ' ) << ' '
              << std::flush;

    if ( needAck ) {
      warningOccured_ = false;
      progressRunning_ = true;
    }
    else {
      cout << endl;
    }
  }

#endif


  // ******************
  //   FinishProgress
  // ******************
  void WriteInfo::FinishProgress( const bool success ) {


    bool okay = ( warningOccured_ == false ) && ( success == true );

    if ( okay == true ) {
      cout << std::setw(10) << fg_green << "OK" << fg_reset<< endl;
    }
    else if ( success == false ) {
      cout << std::setw(10) << fg_red << "FAILED" << fg_reset << endl;
    }
    else {
      cout << endl;
    }

    // re-set flags
    warningOccured_  = false;
    progressRunning_ = false;
  }

  // explicit template instantiation for GCC compiler
#ifdef __GNUC__
  template
  void  WriteInfo::WriteAcouPower<Double>(std::string pdename, 
                                          StdVector<std::string> & subdoms,
                                          Vector<Double>& power);
  template 
  void  WriteInfo::WriteAcouPower<Complex>( std::string pdename, 
                                            StdVector<std::string> & subdoms,
                                            Vector<Complex>& power );
#endif

}
