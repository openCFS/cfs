#include <fstream>
#include <iostream>
#include <string>

#include <cmath>

#include "bubbledriver.hh"

#include "Utils/vector.hh"
#include "DataInOut/timefunc.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "Domain/domain.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODEDescr/Gilmoredimlos.hh"
#include "ODEDescr/LinearKellerMiksis.hh"
#include "ODEDescr/RPNNP.hh"


namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  BubbleDriver::BubbleDriver(Domain * adomain, 
                             UInt   stepOffset,
                             Double timeOffset,
                             std::string driverTag,
                             Boolean isPartOfSequence) 
    : SingleDriver(adomain, stepOffset, timeOffset, 
                   driverTag, isPartOfSequence) {
    ENTER_FCN( "BubbleDriver:: BubbleDriver" );

    //Pointer to read timefunctions for pressure and its derivative
    ptTimeFunc_ = ptdomain_->GetTimeFncPointer();

    // vectors for accessing parameters
    StdVector<std::string> keyVec, attrVec, valVec;

    attrVec = "tag";
    valVec = driverTag_;
   

    // Get time stepping information from parameter object
    keyVec = "transient", "numSteps";
    params->Get(keyVec, attrVec, valVec, numstep_);
    
    keyVec = "transient", "firstDt";
    params->Get(keyVec, attrVec, valVec, firstdt_);

    keyVec = "transient", "stepSaveBeg";
    params->Get(keyVec, attrVec, valVec, isavebegin_);

    keyVec = "transient", "stepSaveEnd";
    params->Get(keyVec, attrVec, valVec, isaveend_);

    keyVec = "transient", "stepSaveInc";
    params->Get(keyVec, attrVec, valVec, isaveincr_);  


    // Make consistency check. In fact in the XML case the Schema 
    // should catch this error. But one can never be sure.
    if(isavebegin_ <= 0)
      {
        Error("Value of stepsavebegin must be positive and nonzero!",
	      __FILE__, __LINE__ );
      }
   

    // =============================================
    //  Query ParamHandler for material parameters
    // =============================================

    // First make sure that there is a section 'bubbleDynamic'
    // We rely on params using a validating Schema parser here.
    StdVector<Double> auxVec;
    params->GetList( "initRadius", auxVec, "bubbleDynamic" );
    if ( auxVec.GetSize() != 1 ) {
      Error( "Cannot find initRadius! Assuming that section 'bubbleDynamic' "
             "is missing in xml-file", __FILE__, __LINE__ );
    }


    initRadius_ = auxVec[0];

    params->Get( "initVel", initVel_, "bubbleDynamic" );

    params->Get( "density", density_, "bubbleDynamic" );

    params->Get( "sonicVel", sonicVel_, "bubbleDynamic" );

    params->Get( "pStatic", pStatic_, "bubbleDynamic" );

    params->Get( "pVapour", pVapour_, "bubbleDynamic" );

    params->Get( "surfaceTension", surfaceTension_, "bubbleDynamic" );

    params->Get( "polytrop", polytrop_, "bubbleDynamic" );

    params->Get( "viscosity", viscosity_, "bubbleDynamic" );

    // Only necessary if bubbledriver should compute 
    // the functions of pressure und its derivative
    params->Get( "pressure", pressureAmpl_, "bubbleDriver" );
 
    params->Get( "frequency", frequency_, "bubbleDriver" );




    // **************************************************
    //   Check what type of bubble model should be used
    // **************************************************
    StdVector<std::string> auxVec2;
    params->GetList( "bubbleType", auxVec2, "bubbleDriver" );
    if ( auxVec2.GetSize() != 1 ) {
      Error( "Cannot find BubbleType! Assuming that section 'bubbleDriver' "
             "is missing in xml-file", __FILE__, __LINE__ );
    }
    
    String2Enum( auxVec2[0], bubbleDynType_ );

    // Choice which bubbledynamical method is used and 
    // creation of pointer to choosen class
    switch(bubbleDynType_){

    case KELLERMIKSIS:
//       ptBubble_ = new KellerMiksis(initRadius_, density_, 
//  				   sonicVel_, pStatic_, 
//  				   pVapour_, surfaceTension_, 
//  				   polytrop_, viscosity_);
      // Should be expanded in XML-File, such that one
      // can choose linearised KellerMiksis Model
      // ----------------------------------------------------
      ptBubble_ = new LinearKellerMiksis(initRadius_,density_, 
					 sonicVel_, pStatic_, 
					 pVapour_, surfaceTension_, 
					 polytrop_, viscosity_);
      // -----------------------------------------------------
      break;
    case GILMORE:
//       ptBubble_ = new Gilmore(initRadius_,density_, 
//                               sonicVel_, pStatic_, 
//                		pVapour_, surfaceTension_, 
//                               polytrop_, viscosity_);

      // Should be expanded in XML-File, such that one
      // can choose between dimensionless and
      // non-dimensionless computation
      // ----------------------------------------------------
      ptBubble_ = new Gilmoredimlos(initRadius_,density_, 
				    sonicVel_, pStatic_, 
				    pVapour_, surfaceTension_, 
 				    polytrop_, viscosity_);
      // -----------------------------------------------------

      // Should be expanded in XML-File, such that one
      // can choose RPNNP-Model
      // ----------------------------------------------------
//       ptBubble_ = new RPNNP(initRadius_,density_, 
// 			    sonicVel_, pStatic_, 
// 			    pVapour_, surfaceTension_, 
// 			    polytrop_, viscosity_);
      // -----------------------------------------------------



      break;

    default:
      Error("No bubblemethod specified ",__FILE__,__LINE__);
    }

    // **************************************************
    //  Set initial values
    // **************************************************

    y_.Push_back( initRadius_ );
    y_.Push_back( initVel_ );


    // **************************************************
    //   Open output file for bubbledynamics
    // **************************************************
    // Generate filename
    std::string simName = commandLine->GetSimName();
    Char *auxfile = new Char[ strlen( simName.c_str() ) + 4 ];
    strcpy( auxfile, simName.c_str() );
    strcat( auxfile, ".bl" );


    // Open file
    fp = fopen( auxfile , "w" );
    if ( fp == NULL ) {
      (*error) << "Could not open outputfile '" << auxfile << "' for writing!";
      Error( __FILE__, __LINE__ );
    }

    // Free space for filename
    delete[] auxfile;

    // Write file header
    fprintf( fp, "# Time\t\tPressure\t\tRadius\t\tVelocity\t\tKlammerRHS\t\tRHS\n" );

    // First line with initial values
    fprintf( fp, "0\t0\t%10.6e\t%10.6e\t0\t0\n", y_[0], y_[1] );
    //    fprintf( fp, "0\t%10.6e\t%10.6e\n", y_[0], y_[1] );

    // **************************************************
    //   Generate ODE solver object
    // **************************************************
    // Later here could be a switch option to choose the solver method
    ptODESolver_ = new ODESolver_RKF45;
    // ptODESolver_ = new ODESolver_Rosenbrock;
    // ptODESolver_ = new ODESolver_ExplEuler;
 
  }
 



  // **********************
  //   Default destructor
  // **********************
  BubbleDriver::~BubbleDriver() {
    ENTER_FCN( " BubbleDriver::~ BubbleDriver" );
  }


  // *****************
  //   Solve problem
  // *****************

  // Method to compute bubble behaviour 
  void  BubbleDriver::SolveProblem() {
    ENTER_FCN( " BubbleDriver::SolveProblem" );

    Double  steptime  = 0;
    UInt stepsave  = isavebegin_;
    Double  dt = firstdt_;

    Double KlammerRHS, rhs,ndichte;

    UInt nstep;

    // Two alternatives to be chosen here in method 
    // a) follow predetermined steps and give output
    // b) Initial- and Endtime is given, choose steps 
    // freely, output should be defined in ODESolve

    Double tEnd = 1.0e-3;  // should mayby included in XML-File
    Double h1;             // Startvalue for adaptive stepping
    h1 = 1e-10;

    StdVector<Double> dydtaussen(2);
    dydtaussen[0] = 0.0;
    dydtaussen[1] = 0.0;
    ndichte = 1e8;



    //Dimlos Case ---------------------------

    dt = dt / initRadius_ * (sqrt( pStatic_ / density_));

    y_[0] /= initRadius_;
    y_[1] /= (sqrt( pStatic_/ density_));
    

    steptime   = steptime / initRadius_ * (sqrt(pStatic_/ density_));
    tEnd       = tEnd / initRadius_ * (sqrt(pStatic_/ density_)); 
    h1         = h1 / initRadius_ * (sqrt(pStatic_/ density_)); 

    pressureAmpl_= pressureAmpl_ / pStatic_;
    frequency_   = frequency_ * initRadius_ /(sqrt(pStatic_/ density_)); 

    Double tdim;
    StdVector<Double> ydim(2);
    // ------------------------------------------------------------



    for (nstep = 1; nstep <= numstep_; nstep++){
      tEnd = steptime + dt;

      // a) forward the actual pressure values to BubbleODE object
      //    - therefore read the values from input file 
      //    - or compute them directly
      // b) no forwarding since they will be computed in BubbleODE


      // In case of self-computed pressure---------------------------------
       pressure_ = - pressureAmpl_ * sin( 2 * PI * frequency_ * (steptime+dt));
       dpresdt_  = - pressureAmpl_ * (2 * PI * frequency_ ) * 
	 cos( 2 * PI * frequency_ * (steptime+dt));     
       ptBubble_->SetP(pressure_);
       ptBubble_->SetDpdt(dpresdt_);

      // In case of input file---------------------------------------
      //Double val_tfunc;
      //Double val_tfuncDerv;

      //For testing write names of function directly
      //std::string tfname = "linpres1500.dat";
      //std::string tfnameDerv = "linpresDerv1500.dat";

      //val_tfunc = ptTimeFunc_->TimeFuncAtTime(steptime+dt,tfname);
      //val_tfuncDerv = ptTimeFunc_->TimeFuncAtTime(steptime+dt,tfnameDerv);

      //ptBubble_->SetP(val_tfunc);
      //ptBubble_->SetDpdt(val_tfuncDerv);
       //------------------------------------------------------------

      // Call for ODESolver
      ptODESolver_->Solve( steptime, tEnd, y_, *ptBubble_, h1, 0, dt);
	  
      steptime += dt;


      // writing results of bubbledynamic in output file
      if (nstep == stepsave && (nstep <= isaveend_)) {
	//dimensionless---------------------------------------
	ydim[0] = y_[0] * initRadius_ ;
  	ydim[1] = y_[1] * sqrt( pStatic_/ density_); 
 	tdim = tEnd *  initRadius_ / (sqrt( pStatic_/ density_));



    ptBubble_->CompDeriv(tEnd, y_, dydtaussen);
    dydtaussen[1] = dydtaussen[1] * pStatic_ / ( density_ * initRadius_ );
    KlammerRHS = 6.0*ydim[0]*ydim[1]*ydim[1]+3.0*ydim[0]*ydim[0]*dydtaussen[1];
    // Achtung hier ist zwei Mal die Dichte eingefügt, damit man mit dem gekoppleten Output vergleichen kann.
    rhs = 4.0/3.0*(M_PI)* ndichte * KlammerRHS* 998.0*998.0;


	
	fprintf(fp, "%10.6e\t%10.6e\t%10.6e\t%10.6e\t%16.10e\t%16.10e\n", tdim,pressure_*pStatic_,ydim[0],ydim[1], KlammerRHS, rhs);
	// ---------------------------------------------------

	//fprintf( fp, "%10.6e\t%10.6e\t%10.6e\t%10.6e\n", steptime,
	//		 pressure_ ,y_[0],y_[1]);
	stepsave+=isaveincr_;
         }

    }


    // Close output file
    if ( fclose( fp ) == EOF ) 
      Error( "Could not close file after writing", __FILE__, __LINE__ );
    
  }

} // end of namespace
