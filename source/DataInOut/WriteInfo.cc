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


    *cfsInfo << *material << std::endl;
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
        *cfsInfo << std::endl << pdeNameLong << "NONLINEAR ITERATION "
                 << iterationCounter 
                 << " ==========================================\n"
                 << pdeNameLong << "=== Residual error          "
                 << residualErr
                 << std::endl
                 << pdeNameLong << "=== Incremental error       "
                 << incrementalErr << std::endl;
    
        if (etaLineSearch)
          *cfsInfo << pdeNameLong << "=== eta (line search)       "
                   << etaLineSearch << std::endl;
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
        *cfsInfo << std::endl << pdeNameLong << "NONLINEAR ITERATION "
                 << iterationCounter
                 << " ==========================================\n"
                 << pdeNameLong << "=== Residual error          "
                 << residualErr
                 << std::endl
                 << pdeNameLong << "=== Incremental error       "
                 << incrementalErr << std::endl
                 << pdeNameLong << "=== Energy error       "
                 << energyErr << std::endl;

        if (etaLineSearch)
          *cfsInfo << pdeNameLong << "=== eta (line search)       "
                   << etaLineSearch << std::endl;
      }
  }

  void WriteInfo::WriteMultiSequenceStep(const UInt sequenceStep, 
                                         const BasePDE::AnalysisType analysis)
  {
    std::string analysisString = BasePDE::analysisType.ToString(analysis);

 
    // write std::out info 
    std::cout << std::endl 
              << " ***************************** " << std::endl
              << " MultiSequenceStep: " << sequenceStep << std::endl
              << " AnalysisType:      " << analysisString << std::endl
              << " ***************************** " << std::endl << std::endl;



    *cla <<  std::endl 
         << " ***************************** " << std::endl
         << " MultiSequenceStep: " << sequenceStep << std::endl
         << " AnalysisType:      " << analysisString << std::endl
         << " ***************************** " << std::endl << std::endl;
    


    
    if (cfsInfo)
      *cfsInfo << std::endl 
               << std::endl<< " ***************************** " << std::endl
               << " MultiSequenceStep: " << sequenceStep << std::endl
               << " AnalysisType:      " << analysisString << std::endl
               << " ***************************** " << std::endl << std::endl;
  }
  


  void WriteInfo::WriteTimeStep(const std::string& pdeName,
                                const UInt timeStep, const Double time)
  {

    std::string pdeNameLong(pdeName);

    // write std::out info    
    std::cout << std::endl << pdeName << ": Time step " 
              << timeStep <<" ======================= " << std::endl;      


    *cla << std::endl << pdeName << ": Time step " 
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
  


  // *********************
  //   WriteHarmonicStep
  // *********************
  void WriteInfo::WriteHarmonicStep( const std::string& pdeName,
                                     const UInt freqStep,
                                     const Double frequency ) {


    // We do not log something in the case of the paramIdent driver,
    // since the output will be more disturbing than helpful in this case
    std::string analysis;
//     params->Get( "type", analysis, "analysis" );
//     if( analysis == "paramIdent" ) {
//       return;
//     }

    // Compute size of freqStep when written as string
    UInt fw = (Integer)std::log10( (float)freqStep ) + 1;

    // Report 1: Goes to standard output
    std::cout << std::endl << pdeName << ": Harmonic step " 
              << freqStep <<" ======================= " << std::endl;      

    // Report 2: Goes into logfile for algebraic sub-system
    (*cla) << std::endl << ' '
           << std::setw( 79 ) << std::setfill( '*' ) << "\n"
           << " *** " << pdeName
           << ": Harmonic step " << freqStep << ' '
           << std::setw( 80 - pdeName.size() - 21 - fw )
           << std::setfill( '*' ) << "*\n "
           << std::setw( 78 ) << std::setfill( '*' ) << "*"
           << std::endl;

    // Report 3: Goes into CFS++ info-file
    if ( cfsInfo != NULL ) {
      (*cfsInfo) << "\n\n"
                 << std::setw( 80 ) << std::setfill( '*' ) << "*\n"
                 << pdeName << "-PDE: HARMONIC STEP " << freqStep 
                 << ", frequency: " << frequency << std::endl;
    }
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
      *cfsInfo << std::endl << " PostProcessing Result for PDE " << pdename
               << ": " << resulttype << " ==========" << std::endl;

      for ( UInt i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "        === " << analysis << " " << analysisVal 
                 << "; " << subdoms[i] << ": " 
                 << results[i] 
                 << " "  << unit << std::endl << std::endl;
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
      *cfsInfo << std::endl << " PostProcessing Result for PDE " << pdename
               << ": " << resulttype << " ==========" << std::endl;

      for ( UInt i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "        === " << analysis << " " << analysisVal 
                 << "; " << subdoms[i] << ": "
                 << results[i] 
                 << " "  << unit << std::endl << std::endl;
      }
    }
  }


  void WriteInfo::PrintCoil( Coil &coil, BasePDE::AnalysisType &analysistype ) {


    if (!cfsInfo)
      return;
    
    *cfsInfo << "COIL DESCRIPTION ======================================= "
             << std::endl;

    // Basic coil info
    *cfsInfo << "Coil domain: "              << coil.coilRegionName_ << '\n'
             << "Coil type: "                << coil.coilTypeS_ << '\n'
             << "Cross-Section of winding: " << coil.windingCrossSection_
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

      if ( analysistype == BasePDE::HARMONIC) {
        *cfsInfo << "Voltage phase = " << coil.phase_ << std::endl;
      }
    }

    // Special things for current coils
    if ( coil.coilType_ == Coil::CURRENT2D ||
         coil.coilType_ == Coil::CURRENT3D ) {

      *cfsInfo << "Current value = " << coil.value_      << std::endl;

      if ( analysistype == BasePDE::HARMONIC) {
        *cfsInfo << "Current phase = " << coil.phase_ << std::endl;
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
        *cfsInfo << coil.locFlowDir_.ToString() << std::endl;
        *cfsInfo << " in coordinate system " 
                 << coil.flowCoordSys_->GetName();
      }
    }

    *cfsInfo << std::endl << std::endl;
  }

  template <class TYPE>
  void WriteInfo::WriteAcouPower(std::string pdename, 
                                 StdVector<std::string> & subdoms,
                                 Vector<TYPE>& power)
  {
 
    if (cfsInfo) {
      *cfsInfo << std::endl << " PostProcessing Result for PDE " << pdename
               << ": " << " ==========" << std::endl;
      *cfsInfo << "   Acoustic Power: \n";
      for ( UInt i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "    Subdomain: " <<  subdoms[i] << " : " <<  power[i] 
                 << " W" << std::endl;
      }
    }
  }


  void WriteInfo::PrintVec(Vector<Complex>& vec)
  {
    if (cfsInfo)
      *cfsInfo << vec << std::endl;
  }

  void WriteInfo::PrintVec(Vector<Double>& vec)
  {
    if (cfsInfo)
      *cfsInfo << vec << std::endl;
  }
  

  void WriteInfo::WriteCombustionNoiseInfo(std::string filename, std::string cplRegion,
                                           UInt sosIdx, UInt src1, UInt src2, UInt src3, 
                                           UInt src4, UInt src5, UInt src6, UInt src7) {


    *cfsInfo << "\nCombustion Noise Info:\n" 
	     << " Name of file: " << filename  
	     << "\n Coupling region: " << cplRegion <<  std::endl;

    if ( sosIdx > 0 ) {      
      *cfsInfo << " Use variable speed of sound: yes" << std::endl;
    }
    else {
      *cfsInfo << " Use variable speed of sound: no" << std::endl;
    }

    *cfsInfo << "\n Use the following source terms" << std::endl;

    if ( src1 > 0 ) {
      *cfsInfo << "      Reynolds stress tensor " << std::endl;
    }
    if ( src2 > 0 ) {
      *cfsInfo << "      Fluctuation of momentum flux " << std::endl;
    }
    if ( src3 > 0 ) {
      *cfsInfo << "      Unsteady reaction rates " << std::endl;
    }
    if ( src4 > 0 ) {
      *cfsInfo << "      2nd order derivative of density " << std::endl;
    }
    if ( src5 > 0 ) {
      *cfsInfo << "      Heat release " << std::endl;
    }    
    if ( src6 > 0 ) {
      *cfsInfo << "      Gas Const. 2nd time derivative " << std::endl;
    }
    if ( src7 > 0 ) {
      *cfsInfo << "      Shear term " << std::endl;
    }

    *cfsInfo << std::endl;
  }


  void WriteInfo::PrintVec(StdVector<Integer>& vec)
  {
    if (cfsInfo)
      *cfsInfo << vec << std::endl;
  }
  


  void WriteInfo::PrintVec(const char * comment, StdVector<Integer>& vec)
  {
    if (cfsInfo)
      *cfsInfo << comment << std::endl << vec << std::endl << std::endl;
  }



  void WriteInfo::PrintVec(const char * comment,
                           StdVector<std::string>& vec)
  {

    if (cfsInfo)
      {
        *cfsInfo << comment << std::endl;
        
        for (UInt i=0; i< vec.GetSize(); i++)
          *cfsInfo << vec[i] << std::endl;
        
        *cfsInfo << std::endl;
      }
  }

  void WriteInfo::PrintMatrix(std::string &comment, const Matrix<Double> &mat)
  {

    if (cfsInfo)
      *cfsInfo << comment << std::endl << mat << std::endl << std::endl;
  }

#endif //INTEGLIB

  // prints warning to info-file
  void WriteInfo::Warning( const std::string & Text,
                           const char* const filename, const UInt numline ) {

#ifdef INTEGLIB
    std::cerr << "INTEGLIB WARNING: " << Text << std::endl;
#else
#ifdef USE_SCRIPTING
    if ( messenger ) {
      if ( messenger->IsEvaluating() ) {
        messenger->Warning( Text.c_str(), filename, numline );
      }
    }
#endif

    if ( progressRunning_ == true ) {
      std::cout << std::endl;
    }

    std::cerr << "\n "
              << fg_yellow << "WARNING:" << fg_reset << "\n "
              << Text << std::endl;

    warningOccured_ = true;

    if (filename) {
      std::cerr <<" (" << filename <<" ";
      if (numline) {
        std::cerr << numline;
      }
      std::cerr << ")" << std::endl;
    }
    else {
      std::cerr << std::endl;
    }

    if (cfsInfo) {
      *cfsInfo << std::endl << std::endl << std::endl
               << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               << "!!!!!!!!!!!!!" << std::endl
               << "                          WARNING " << std::endl
               << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               << "!!!!!!!!!!!!!" << std::endl
               << "WARNING: " << Text;
        
      if (filename) {
        *cfsInfo <<" (" << filename <<" ";
        if (numline) 
          *cfsInfo << numline;
        *cfsInfo << ")";
      }
    
      *cfsInfo << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               << "!!!!!!!!!!!!!!!" << std::endl
               << std::endl;
    }
#endif // INTEGLIB
  }


  // *********
  //   Error
  // *********
  void WriteInfo::Error( const std::string &Text, const char *const filename,
                         const UInt numline ) {


#ifdef INTEGLIB
    std::cerr << "INTEGLIB ERROR: " << Text << std::endl;
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

    std::cerr << std::endl << " "
              << fg_red << "ERROR:"
              << fg_reset << "\n" << std::endl;

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

    std::cerr << std::endl << std::endl;
    
    if (cfsInfo) {
      *cfsInfo << std::endl << std::endl << std::endl
               << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
               << std::endl
               << "                          ERROR " << std::endl
               << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
               << std::endl
               << " ERROR: " << Text;

      if (filename) {
        *cfsInfo <<" (" << filename <<" ";
        if (numline) 
          *cfsInfo << numline;
        *cfsInfo << ")";
      }
      *cfsInfo << std::endl;  
    }
#endif // INTEGLIB

    // exit program
    exit(-1);

    //     char* a = NULL;
    //     std::cerr << a[0] << a[1] << a[2] << a[3];
  }


#ifndef INTEGLIB


  // **************
  //   WriteHomBC
  // **************
  void WriteInfo::WriteHomDirBC( const std::string& pdeName,
                                 HdBcList& list ) {
                                 
    std::string prefix = pdeName + "-PDE: ";

    if (cfsInfo && list.GetSize() > 0) {
      *cfsInfo << prefix << "Homogeneous Dirichlet BC defined on \n";
      
      for( UInt i = 0; i < list.GetSize(); i++ ) {
        HomDirichletBc const & actBc = *list[i];
        EntityList const & actList = *actBc.entities;
        std::string listType;
        EntityList::Enum2String( actList.GetType(), listType );
        *cfsInfo << prefix << "\t'" << actList.GetName() << "' (" 
                 <<  listType
                 << "), dof = " 
                 << actBc.result->GetDofName(actBc.dof)
                 << std::endl;
      }

      *cfsInfo << std::endl;
    }
  }


  // ****************
  //   WriteInHomBC
  // ****************
  void WriteInfo::WriteInhomDirBC( const std::string& pdeName, 
                                   IdBcList& list ) {
    
    
    std::string prefix = pdeName + "-PDE: ";

    if (cfsInfo && list.GetSize() > 0 ) {
      *cfsInfo << prefix << "Inhomogeneous Dirichlet BC defined on \n";
      
      for( UInt i = 0; i < list.GetSize(); i++ ) {
        InhomDirichletBc const & actBc = *list[i];
        EntityList const & actList = *actBc.entities;
        std::string listType;
        EntityList::Enum2String( actList.GetType(), listType );
        *cfsInfo << prefix << "\t'" << actList.GetName() << "' (" 
                 << listType
                 << "), dof = " 
                 << actBc.result->GetDofName(actBc.dof)
                 << ", value = " << actBc.value
                 << ", phase = " << actBc.phase 
                 << std::endl;
      }

      *cfsInfo << std::endl;
    }
  }

  // *******************
  //   WriteInhomNeuBC
  // *******************
  void WriteInfo::WriteInhomNeuBC( const std::string& pdeName, 
                                   InBcList& list ) {
    
    std::string prefix = pdeName + "-PDE: ";
    
    if (cfsInfo && list.GetSize() > 0 ) {
      *cfsInfo << prefix << "Inhomogeneous Neumann BC defined on \n";
      
      for( UInt i = 0; i < list.GetSize(); i++ ) {
        InhomNeumannBc const & actBc = *list[i];
        EntityList const & actList = *actBc.entities;
        std::string listType;
        EntityList::Enum2String( actList.GetType(), listType );
        *cfsInfo << prefix << "\t'" << actList.GetName() << "' (" 
                 << listType
                 << "), dof = " 
                 << actBc.result->GetDofName(actBc.dof)
                 << ", value = " << actBc.value
                 << ", phase = " << actBc.phase 
                 << std::endl;
      }
    }
    *cfsInfo << std::endl;
  }


  // ********************
  //   WriteConstraints
  // ********************
  void WriteInfo::WriteConstraints(  const std::string& pdeName, 
                                     ConstraintList& list ) {
    
   std::string prefix = pdeName + "-PDE: ";
    
    if (cfsInfo && list.GetSize() > 0 ) {
      *cfsInfo << prefix << "Constraints defined on \n";
      
      for( UInt i = 0; i < list.GetSize(); i++ ) {
        Constraint const & actBc = *list[i];
        EntityList const & masterList = *actBc.masterEntities;
        EntityList const & slaveList = *actBc.masterEntities;
        std::string masterListType, slaveListType;
        EntityList::Enum2String( masterList.GetType(), masterListType );
        EntityList::Enum2String( slaveList.GetType(), slaveListType );
        *cfsInfo << prefix << "\tMaster: '" << masterList.GetName() 
                 << "' (" 
                 <<  masterListType
                 << "), dof = " 
                 << actBc.result->GetDofName(actBc.masterDof)
                 << "\tSlave: '"<< slaveList.GetName() 
                 << "' (" 
                 <<  slaveListType
                 << "), dof = " 
                 << actBc.result->GetDofName(actBc.slaveDof)
          
                 << std::endl;
      }
    }
    *cfsInfo << std::endl;

  }


  void WriteInfo::WriteLoad(const std::string& pdeName, 
                            LoadList& list ) {
    
    std::string prefix = pdeName + "-PDE: ";
    
    if (cfsInfo && list.GetSize() > 0 ) {
      *cfsInfo << prefix << "Load defined on \n";
      
      for( UInt i = 0; i < list.GetSize(); i++ ) {
        LoadBc const & actBc = *list[i];
        EntityList const & actList = *actBc.entities;
        std::string listType;
        EntityList::Enum2String( actList.GetType(), listType );
        *cfsInfo << prefix << "\t'" << actList.GetName() << "' (" 
                 << listType
                 << "), dof = " 
                 << actBc.result->GetDofName(actBc.dof)
                 << ", value = " << actBc.value
                 << ", phase = " << actBc.phase 
                 << std::endl;
      }
    }
    *cfsInfo << std::endl;
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
      //*cfsInfo << formatted << std::endl << std::flush;
    *cfsInfo << formatted << std::flush;
    
    va_end(argList);
  }


  // *****************
  //   StartProgress
  // *****************
  void WriteInfo::StartProgress( const std::string &name, bool needAck ) {

   
    std::string modifiedName = name + " ";

    needAck_ = needAck;

    std::cout << "++ " << std::setw(PROGRESS_TEXT_WIDTH)
              << std::setfill( '.' ) << std::left
              << modifiedName << std::setfill( ' ' ) << ' '
              << std::flush;

    if ( needAck ) {
      warningOccured_ = false;
      progressRunning_ = true;
    }
    else {
      std::cout << std::endl;
    }
  }

#endif


  // ******************
  //   FinishProgress
  // ******************
  void WriteInfo::FinishProgress( const bool success ) {


    bool okay = ( warningOccured_ == false ) && ( success == true );

    if ( okay == true ) {
      std::cout << std::setw(10) << fg_green << "OK" << fg_reset<< std::endl;
    }
    else if ( success == false ) {
      std::cout << std::setw(10) << fg_red << "FAILED" << fg_reset << std::endl;
    }
    else {
      std::cout << std::endl;
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
