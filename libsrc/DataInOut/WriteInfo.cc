#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdarg>

#include "WriteInfo.hh"
#include "Utils/tools.hh"
#include "Utils/Coil.hh"
#include "Materials/baseMaterial.hh"
#include "PDE/pdes_header.hh"
#include "Utils/vector.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

#define CFS_VERSION  "0.2"

#define PROGRESS_TEXT_WIDTH 62

namespace CoupledField {

  WriteInfo::WriteInfo () {

    ENTER_FCN( "WriteInfo::WriteInfo");
    warningOccured_ = FALSE;
    progressRunning_ = FALSE;
    needAck_ = FALSE;

    cfsInfo = NULL;
  }


  // **************
  //   Destructor
  // **************
  WriteInfo::~WriteInfo () {
    ENTER_FCN( "WriteInfo::~WriteInfo" );
    delete cfsInfo;
  }


  // **************
  //   CreateFile
  // **************
  void WriteInfo::CreateFile() {

    ENTER_FCN( "WriteInfo::CreateFile" );

    // Check if a file is already open
    if ( cfsInfo != NULL ) {
      (*error) << "WriteInfo::CreateFile: cfsInfo already points to a "
               << "file! Cowardly refusing to open another one\n";
      CoupledField::Error( __FILE__, __LINE__ );
    }

    std::string filename = commandLine->GetSimName() + ".info";
    cfsInfo = new std::ofstream(filename.c_str());

    // Check if everything went fine
    if ( cfsInfo == NULL ) {
      (*error) << "WriteInfo::CreateFile: Failed to open file '"
               << filename << "' for writing status messages\n";
      CoupledField::Error( __FILE__, __LINE__ );  
    }
  }


  void WriteInfo::PrintHeader()
  {
   
    std::stringstream header;
    std::string compileDate = __DATE__;
    std::string version = CFS_VERSION;
    header << "============================================================"
           << "===========\n"
           << "============================================================"
           << "===========\n"
           << "|                                                           "
           << "          |\n"
           << "|                                 CFS++                     "
           << "          |\n"
           << "|                     (Coupled Field Simulation ++)         "
           << "          |\n"
           << "|                                                           "
           << "          |\n"
           << "|                                                           "
           << "          |\n"
           << "|  Version:  " << version << std::setw(70-11-version.length())
           << "|\n"
           << "|  Compiled: " << compileDate
           << "                                              |\n"
           << "|                                                           "
           << "          |\n"
           << "============================================================"
           << "===========\n"
           << "============================================================"
           << "===========\n";

    std::cout << std::endl << header.str() << std::endl;

    if (cfsInfo) {
      *cfsInfo << header.str();
    }
  }
  

  // *****************
  //   PrintPiezoMat
  // *****************
  void WriteInfo::PrintPiezoMat( BaseMaterial* material, Boolean imag ) {

    ENTER_FCN( "WriteInfo::PrintPiezoMat" );

    Matrix<Double> myMat;

    if ( cfsInfo != NULL ) {
      if ( imag ) {
	material->GetTensor(myMat,PIEZO_TENSOR,IMAG);
	*cfsInfo  << "PIEZO-TENSOR (imaginary part)" << material->GetName()
		  << ":\n\n" << myMat << std::endl
		  << std::endl <<  std::endl;
      }
      else {
	material->GetTensor(myMat,PIEZO_TENSOR,REAL);
	*cfsInfo  << "PIEZO-TENSOR (real part)" << material->GetName()
		  << ":\n\n" << myMat << std::endl
		  << std::endl <<  std::endl;
      }
    }
  }

  // *****************
  //   PrintMechanicMat
  // *****************
  void WriteInfo::PrintMechanicMat( BaseMaterial* material, Boolean imag ) {

    ENTER_FCN( "WriteInfo::PrintMechanicMat" );

    Matrix<Double> myMat;

    if ( cfsInfo != NULL ) {
      if (imag ) {
	material->GetTensor(myMat,MECH_STIFFNESS_TENSOR,IMAG);
	*cfsInfo  << "MECHANIC STIFFNESS TENSOR (imaginary part) OF " << material->GetName()
		  << ":\n\n" << myMat
		  << std::endl <<  std::endl;
      }
      else {
	Double density, alpha, beta;
	material->GetTensor(myMat,MECH_STIFFNESS_TENSOR,REAL);
	material->GetScalar(density,DENSITY,REAL);
	*cfsInfo  << "MECHANIC STIFFNESS TENSOR (real part) OF " << material->GetName()
		  << ":\n\n" << myMat << std::endl
		  << "density = " << density << std::endl;
        
        // Check, if damping was specified
        if( material->IsSet( RAYLEIGH_ALPHA) &&
            material->IsSet( RAYLEIGH_BETA) ) {
          
          
          material->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
          material->GetScalar(beta,RAYLEIGH_BETA,REAL);
          

          *cfsInfo << "damping coefficient alfa = " << alpha
                   << std::endl
                   << "damping coefficient beta = " << beta
		   << std::endl <<  std::endl;
        } else {
          *cfsInfo << "damping: not defined"<< std::endl << std::endl;
        }
      }
    }
  }


  // *****************
  //   PrintAcousticMat
  // *****************
  void WriteInfo::PrintAcousticMat( BaseMaterial* material ) {
    ENTER_FCN( "WriteInfo::PrintAcousticMat" );

    if (cfsInfo) {
      Double compress, density, alpha, beta, BoverA, acousticAlpha, fracExp;
      material->GetScalar(compress,ACOU_BULK_MODULUS,REAL);
      material->GetScalar(density,DENSITY,REAL);
      material->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
      material->GetScalar(beta,RAYLEIGH_BETA,REAL);
      material->GetScalar(BoverA,BOVERA,REAL);
      material->GetScalar(acousticAlpha, ACOU_ALPHA, REAL);    
      material->GetScalar(fracExp, FRACTIONAL_EXPONENT, REAL);

      *cfsInfo << "MATERIAL DATA OF " << material->GetName() << ":"
               << std::endl 
               << "compression modulus [N/m^2] = "
               << compress << std::endl
               << "density [kg/m^3]            = " 
               << density << std::endl
               << "alpha                       = " 
               << alpha << std::endl
               << "beta                        = " 
               << beta << std::endl
               << "B/A                         = " 
               << BoverA << std::endl
               << "Acoustic alpha [Np/m]  = " 
               << acousticAlpha << std::endl
               << "fractional exponent         = " 
               << fracExp << std::endl << std::endl;
    }
  }
  

  // *****************
  //   PrintFlowMat
  // *****************
  void WriteInfo::PrintFlowMat( BaseMaterial* material ) {
    ENTER_FCN( "WriteInfo::PrintFlowMat" );

    if (cfsInfo) {
      Double density, dynViscosity;
      material->GetScalar(density,DENSITY,REAL);
      material->GetScalar(dynViscosity,DYNAMIC_VISCOSITY,REAL);

      *cfsInfo << "MATERIAL DATA OF " << material->GetName() << ":"
               << std::endl 
               << "density [kg/m^3]            = " 
               << density << std::endl
               << "dynamic Viscosity      = " 
               << dynViscosity << std::endl << std::endl;
    }
  }
  

  void WriteInfo::PrintThermicMat( BaseMaterial* material ) {
    ENTER_FCN( "WriteInfo::PrintThermicMat" );

    if (cfsInfo) {
      Double density, capacity, conductivity;
      material->GetScalar(density,DENSITY,REAL);
      material->GetScalar(capacity,HEAT_CAPACITY,REAL);
      material->GetScalar(conductivity,HEAT_CONDUCTIVITY,REAL);

      *cfsInfo << "MATERIAL DATA OF " << material->GetName() << ":"
               << std::endl
               << "density [kg/m^3]             = " 
               << density << std::endl
               << "heat capacity [J/kg K]       = " 
               << capacity << std::endl
               << "thermal conductivity [W/m K] = " 
               << conductivity << std::endl << std::endl;
    }
  }

  void WriteInfo::PrintMagMat(BaseMaterial* material)
  {
    ENTER_FCN( "WriteInfo::PrintMagMat" );

    if (cfsInfo) {
      //      Double mX, mY, mZ;
      //      material.GetPermMag(mX, mY, mZ);
      Double perm,cond;
      material->GetScalar(perm,MAG_PERMEABILITY,REAL);
      material->GetScalar(cond,MAG_CONDUCTIVITY,REAL);

      *cfsInfo << "MATERIAL DATA OF " << material->GetName() << ":" 
	       << std::endl
	       << "conductivity:            " << cond
	       << std::endl
	       << "permeability:            " << perm
	       << std::endl
	//	       << "vector of magnetiziation: (" << mX << ", " << mY << ", "
	//	       << mZ <<")" 
	       << std::endl << std::endl;
    }
  }
  
  // ************************
  //   PrintElectrostaticMat
  // ************************
  void WriteInfo::PrintElectrostaticMat( BaseMaterial* material, Boolean imag ) {

    ENTER_FCN( "WriteInfo::PrintElectrostaticMat" );

    Matrix<Double> myMat;

    material->GetTensor(myMat,ELEC_PERMITTIVITY,REAL);

    if ( cfsInfo != NULL ) {
      if ( imag ) {
	material->GetTensor(myMat,ELEC_PERMITTIVITY,IMAG);
	*cfsInfo  << "PERMITTIVITY TENSOR (imaginray part) OF " << material->GetName()
		  << ":\n\n" << myMat << std::endl
		  << std::endl;
      }
      else {
	material->GetTensor(myMat,ELEC_PERMITTIVITY,REAL);
	*cfsInfo  << "PERMITTIVITY TENSOR (real part) OF " << material->GetName()
		  << ":\n\n" << myMat << std::endl
		  << std::endl;
      }
    }
  }


  void WriteInfo::WriteNonLinIter(const std::string& pdeName,
                                  const UInt iterationCounter,
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


  void WriteInfo::WriteMultiSequenceStep(const UInt sequenceStep, 
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
                                const UInt timeStep, const Double time)
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
  


  // *********************
  //   WriteHarmonicStep
  // *********************
  void WriteInfo::WriteHarmonicStep( const std::string& pdeName,
                                     const UInt freqStep,
                                     const Double frequency ) {

    ENTER_FCN( "WriteInfo::WriteHarmonicStep" );

    // We do not log something in the case of the paramIdent driver,
    // since the output will be more disturbing than helpful in this case
    std::string analysis;
    params->Get( "type", analysis, "analysis" );
    if( analysis == "paramIdent" ) {
      return;
    }

    // Compute size of freqStep when written as string
    UInt fw = (Integer)std::log10( (float)freqStep ) + 1;

    // Report 1: Goes to standard output
    std::cout << myEndl << pdeName << ": Harmonic step " 
              << freqStep <<" ======================= " << std::endl;      

    // Report 2: Goes into logfile for algebraic sub-system
    (*cla) << myEndl << ' '
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
  void WriteInfo:: WriteResult(std::string pdename, std::string resulttype,
                               StdVector<std::string> & subdoms,
                               Vector<Double> & results,
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

      for ( UInt i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "        === " << analysis << " " << analysisVal 
                 << "; " << subdoms[i] << ": " 
                 << results[i] 
                 << " "  << unit << std::endl << std::endl;
      }
    }
  }

  void WriteInfo:: WriteResult(std::string pdename, std::string resulttype,
                               StdVector<std::string> & subdoms,
                               Vector<Complex> & results,
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

      for ( UInt i = 0; i < subdoms.GetSize(); i++ ) {
        *cfsInfo << "        === " << analysis << " " << analysisVal 
                 << "; " << subdoms[i] << ": "
                 << results[i] 
                 << " "  << unit << std::endl << std::endl;
      }
    }
  }


  void WriteInfo::PrintCoil( Coil &coil, AnalysisType &analysistype ) {

    ENTER_FCN( "WriteInfo::PrintCoil" );

    if (!cfsInfo)
      return;
    
    *cfsInfo << "COIL DESCRIPTION ======================================= "
             << myEndl;

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

      if ( analysistype == HARMONIC || analysistype == MULTIHARMONIC) {
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

      if ( analysistype == HARMONIC || analysistype == HARMONIC) {
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


  void WriteInfo::PrintVec(Vector<Complex>& vec)
  {
    ENTER_FCN( "WriteInfo::PrintVec" );
    if (cfsInfo)
      *cfsInfo << vec << std::endl;
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
        
        for (UInt i=0; i< vec.GetSize(); i++)
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
  void WriteInfo::Warning( const std::string & Text,
                           const Char* const filename, const UInt numline ) {

    ENTER_FCN( "WriteInfo::Warning" );

#ifdef TCL_INTERFACE
    if ( messenger->IsEvaluating() ) {
      messenger->Warning( Text.c_str(), filename, numline );
    }
#endif

    if ( progressRunning_ == TRUE ) {
      std::cout << std::endl;
    }

    std::cerr << "\n \033[31mWARNING:\033[0m\n " << Text << myEndl;

    warningOccured_ = TRUE;

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
    
      *cfsInfo << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
               << "!!!!!!!!!!!!!!!" << myEndl
               << myEndl;
    }
  }


  // *********
  //   Error
  // *********
  void WriteInfo::Error( const std::string &Text, const Char *const filename,
                         const UInt numline ) {

    ENTER_FCN( "WriteInfo::Error" );

#ifdef TCL_INTERFACE
    if ( messenger->IsEvaluating() ) {
      messenger->Error( Text.c_str(), filename, numline );
    }
#endif
    

    // If a progress part is still there, then finish it with
    // a failure
    if ( progressRunning_ ) {
      FinishProgress( FALSE );
    }

    std::cerr << std::endl << " \033[31mERROR:\033[0m\n" << myEndl;
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

#ifdef TRACE
    if ( commandLine->GetTraceDepth() > 0 ) {
      OutInfo::FcnTraceHandler::Dump();
      std::cerr << "\n\n See '" << commandLine->GetSimName()
                << ".trace' for trace dump of function call tree.";
    }
#endif

    std::cerr << std::endl << std::endl;
    
    if (cfsInfo) {
      *cfsInfo << myEndl << myEndl << myEndl
               << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
               << myEndl
               << "                          ERROR " << myEndl
               << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
               << myEndl
               << " ERROR: " << Text;

      if (filename) {
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


  // **************
  //   WriteHomBC
  // **************
  void WriteInfo::WriteHomDirBC( const std::string& pdeName,
                                 const std::string& subDom, UInt dof ) {
    ENTER_FCN( "WriteInfo::WriteHomBC" );
    
    if (cfsInfo) {
      *cfsInfo << pdeName << "-PDE: Homogenous Dirichlet BC on \"" << subDom
               << "\"";
      if (dof) {
        *cfsInfo << " with DOF number " << dof;
      }
      *cfsInfo << myEndl;
    }
  }


  // ****************
  //   WriteInHomBC
  // ****************
  void WriteInfo::WriteInhomDirBC( const std::string& pdeName,
                                   const std::string& subDom, 
                                   const std::string& val, const std::string & fnc,
                                   const UInt& dof ) {

    ENTER_FCN( "WriteInfo::WriteInHomBC" );
    
    if (cfsInfo) {
      *cfsInfo << pdeName << "-PDE: Inhomogenous Dirichlet BC on \""
               << subDom  << "\"";
      if (dof) {
        *cfsInfo << " with DOF number " << dof;
      }
      *cfsInfo << ", value = " <<  val << ", FncName: " << fnc; 
      *cfsInfo << myEndl;
    }
  }

  // *******************
  //   WriteInhomNeuBC
  // *******************
  void WriteInfo::WriteInhomNeuBC( const std::string& pdeName,
                                   const std::string& subDom) { 
    
    ENTER_FCN( "WriteInfo::WriteInHomNeuBC" );
    
    if (cfsInfo) {
      *cfsInfo << pdeName << "-PDE: Inhomogenous Neumann BC on \""
               << subDom  << "\"";
      *cfsInfo << myEndl;
    }
  }


  // ********************
  //   WriteConstraints
  // ********************
  void WriteInfo::WriteConstraints( const std::string& pdeName,
                                    const std::string& subDom, UInt dof ) {
    ENTER_FCN( "WriteInfo::WriteConstraints" );
    
    if (cfsInfo) {
      *cfsInfo << pdeName << "-PDE: Constraints on \"" << subDom
               << "\"";
      if (dof) {
        *cfsInfo << " with DOF number " << dof;
      }
      *cfsInfo << myEndl;
    }
  }


  void WriteInfo::WriteLoad(const std::string& pdeName,
                            const std::string& subDom, 
                            std::string value, const std::string & fnc,
                            UInt dof)
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
          (*error) << "Format character " << formatChar
                   << " not yet defined!";
          CoupledField::Error( __FILE__, __LINE__ );
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
  void WriteInfo::StartProgress( const std::string &name, Boolean needAck ) {

    ENTER_IFCN( "WriteInfo::StartProgress" );
   
    std::string modifiedName = name + " ";

    needAck_ = needAck;

    std::cout << "++ " << std::setw(PROGRESS_TEXT_WIDTH)
              << std::setfill( '.' ) << std::left
              << modifiedName << std::setfill( ' ' ) << ' '
              << std::flush;

    if ( needAck ) {
      warningOccured_ = FALSE;
      progressRunning_ = TRUE;
    }
    else {
      std::cout << std::endl;
    }
  }


  // ******************
  //   FinishProgress
  // ******************
  void WriteInfo::FinishProgress( const Boolean success ) {

    ENTER_IFCN( "WriteInfo::StartProgress" );

    bool okay = ( warningOccured_ == FALSE ) && ( success == TRUE );

    if ( okay == true ) {
      std::cout << std::setw(10) << "\033[32mOK\033[0m" << std::endl;
    }
    else if ( success == FALSE ) {
      std::cout << std::setw(10) << "\033[31mFAILED\033[0m" << std::endl;
    }
    else {
      std::cout << std::endl;
    }

    // re-set flags
    warningOccured_  = FALSE;
    progressRunning_ = FALSE;
  }
  
}
