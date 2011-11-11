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
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Materials/BaseMaterial.hh"
#include "ColoredConsole.hh"
#include "Domain/Results/ResultInfo.hh"
#include "DataInOut/ProgramOptions.hh"
#include "MatVec/Vector.hh"

#ifdef USE_SCRIPTING 
#include "DataInOut/Scripting/CFSMessenger.hh"
#endif

using std::cout;
using std::endl;

namespace CoupledField {


  WriteInfo::WriteInfo () {
    cfsInfo = NULL;
  }


  // **************
  //   Destructor
  // **************
  WriteInfo::~WriteInfo () {
    delete cfsInfo;
  }


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

  void WriteInfo::PrintVec(Vector<Double>& vec)
  {
    if (cfsInfo)
      *cfsInfo << vec << endl;
  }
  
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



}
