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




namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  BubbleDriver::BubbleDriver(Domain * adomain, 
                             UInt stepOffset,
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


    // Make consistency check. In fact in the XML case the Schema should catch
    // this error. But one can never be sure.
    if(isavebegin_ <= 0)
      {
        Error( "Value of stepsavebegin must be positive and nonzero!",
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

    Double initRadius = auxVec[0];
    Double initVel = 0.0;
    params->Get( "initVel", initVel, "bubbleDynamic" );
    Double density = 0.0;
    params->Get( "density", density, "bubbleDynamic" );
    Double sonicVel = 0.0;
    params->Get( "sonicVel", sonicVel, "bubbleDynamic" );
    Double pStatic = 0.0;
    params->Get( "pStatic", pStatic, "bubbleDynamic" );
    Double pVapour = 0.0;
    params->Get( "pVapour", pVapour, "bubbleDynamic" );
    Double surfaceTension = 0.0;
    params->Get( "surfaceTension", surfaceTension, "bubbleDynamic" );
    Double polytrop = 0.0;
    params->Get( "polytrop", polytrop, "bubbleDynamic" );
    Double viscosity = 0.0;
    params->Get( "viscosity", viscosity, "bubbleDynamic" );

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

    //Online necessary if bubbledriver should compute 
    //the functions of pressure und its derivative
    params->Get( "pressure", pressureAmpl_, "bubbleDriver" );
 
    params->Get( "frequency", frequency_, "bubbleDriver" );
        


    // Choice which bubbledynamical method is used and 
    // creation of pointer to choosen class
    switch(bubbleDynType_){
    case KELLERMIKSIS:
      ptBubble_ = new KellerMiksis(initRadius,density, sonicVel, pStatic, 
                                   pVapour, surfaceTension, polytrop, viscosity);
      break;
    case GILMORE:
      ptBubble_ = new Gilmore(initRadius,density, sonicVel, pStatic, 
                              pVapour, surfaceTension, polytrop, viscosity);
      break;
    default:
      Error("No bubblemethod specified ",__FILE__,__LINE__);
    }

    // For test phase only two values are needed.
    // Need to be exchanged if bubbledynamic and acoustic is coupled then
    // y_[0] and y_[1] each will become a vector of length number of elements
    y_.Push_back( initRadius );
    y_.Push_back( initVel );



    //________________________________________________________________

    // Testanfangwerte für Dreikörperproblem
   

    //   y_.Push_back(0.994);
    //   y_.Push_back( 0.0);
    //   y_.Push_back(0.0);
    //   y_.Push_back(-2.001585106);

    //________________________________________________________________


    // ----------------------------------------------------------------
    // Open output file in which time step, pressure, Radius and 
    // Velocity at wall of bubble is written



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
    fprintf( fp, "# Time\t\tPressure\t\tRadius\t\tVelocity\n" );

    // First line with initial values
    fprintf( fp, "0\t0\t%16.10e\t%16.10e\n", y_[0], y_[1] );



    // Generate ODE solver object
    ptODESolver_ = new ODESolver_RKF45;
    
    // Optional solver Explicit Euler 
    // Later here could be a switch option to choose the solver method
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

  // Method to compute bubble behavour, 
  void  BubbleDriver::SolveProblem() {
    ENTER_FCN( " BubbleDriver::SolveProblem" );

    //    UInt level     = 0;
    Double  steptime  = 0;
    UInt stepsave  = isavebegin_;

    Double  dt = firstdt_;

    UInt nstep;

    for (nstep = 1; nstep <= numstep_; nstep++){

      // forward the actual pressure values to BubbleODE object
      // in test phase read the values from input file 
      // or compute them directly

      Double val_tfunc;
      Double val_tfuncDerv;

      //For testing write names of function directly
      std::string tfname = "linpres1500.dat";
      std::string tfnameDerv = "linpresDerv1500.dat";

      val_tfunc = ptTimeFunc_->TimeFuncAtTime(steptime+dt,tfname);
      val_tfuncDerv = ptTimeFunc_->TimeFuncAtTime(steptime+dt,tfnameDerv);
      //    std::cout << "val=" << val_tfunc << std::endl;

      pressure_ = - pressureAmpl_ * sin( 2 * PI * frequency_ * (steptime+dt));
      dpresdt_  = - pressureAmpl_ * (2 * PI * frequency_ ) * 
        cos( 2 * PI * frequency_ * (steptime+dt));

      // for test of dimensionless
      // Double pressureDimlos;
      //Double dpresdtDimlos; 

      //pressureDimlos = pressure_ /1e5;
      //dpresdtDimlos = dpresdt_ * 10e-6 /1e5 /(sqrt(1e5/998));
  
      //StdVector<Double> yDimlos;
      //yDimlos = y_;
      //yDimlos[0] /= 10e-6;
      //yDimlos[1] /= (sqrt(1e5/998));
      //ptBubble_->SetP(pressureDimlos);
      //ptBubble_->SetDpdt(dpresdtDimlos);

      //Double  steptimeDimlos  = 0;
      //Double  dtDimlos = 0;

      //steptimeDimlos = steptime / 10-6 * (sqrt(1e5/998));
      //dtDimlos       = dt / 10-6 * (sqrt(1e5/998)); 

      // In case of self-computed pressure
      //ptBubble_->SetP(pressure_);
      //ptBubble_->SetDpdt(dpresdt_);




      // In case pressure ans its derivative are given
      ptBubble_->SetP(val_tfunc);
      ptBubble_->SetDpdt(val_tfuncDerv);


      ptODESolver_->Solve( steptime, steptime+dt, y_, *ptBubble_, 
                           dt / 3.0, 0, dt);

      //ptODESolver_->Solve( steptime, steptime+dt, yDimlos, *ptBubble_, 
      //                   dt / 3.0, 0, dt);
          
      steptime += dt;

      // for test of dimensionless
      //y_[0] = yDimlos[0] * (sqrt(1e5/998)) ;
      //y_[1] = yDimlos[1] * 1e5 / 10e-6 / 998; 


      // writing results of bubbledynamic in output file


      if (nstep == stepsave && (nstep <= isaveend_)) {
        fprintf( fp, "%e\t%e\t%16.10e\t%16.10e\n", steptime,
                 pressure_ ,y_[0],y_[1]);
        stepsave+=isaveincr_;
      }
    }


    // Close output file
    if ( fclose( fp ) == EOF ) 
      Error( "Could not close file after writing", __FILE__, __LINE__ );
    
  }

} // end of namespace
