#include <fstream>
#include <iostream>
#include <string>


#include <cmath>

#include "bubbledriver.hh"
#include "Utils/vector.hh"

#include "DataInOut/GMV/outGMV.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"






//----------------------------------------------------------------------
//  test phase: later constants will be read from external file
//              xml, mat or own xml file
//  constants are needed for bubbledynamics

//Daten Daehnke (stabile Kavitation mit KellerMiksis)
#define RADIUSINIT 1e-5
#define VELINIT    0
#define DENSITY    998.0  
#define SONICVEL   1.43e03
#define PSTATIC    1e05
#define PVAPOUR    2.33e03
#define SURFACTEN  7.25e-02
#define POLYTROP   4.0/3.0
#define VISCOSITY  1e-03
//----------------------------------------------------------------------

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  BubbleDriver::BubbleDriver(Domain * adomain, 
			     Integer stepOffset,
			     Double timeOffset,
			     std::string driverTag,
			     Boolean isPartOfSequence) 
    : SingleDriver(adomain, stepOffset, timeOffset, 
		   driverTag, isPartOfSequence) {
    ENTER_FCN( "BubbleDriver:: BubbleDriver" );


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
   

    // Choice which bubbledynamical method is used and 
    // creation of pointer to choosen class
    switch(bubbleDyn){
      case KELLERMIKSIS:
	ptBubble_ = new KellerMiksis(RADIUSINIT, DENSITY, 
				     SONICVEL, PSTATIC, PVAPOUR, SURFACTEN,
				     POLYTROP, VISCOSITY);
	break;
      case GILMORE:
 	ptBubble_ = new Gilmore(RADIUSINIT, DENSITY, 
				SONICVEL, PSTATIC, PVAPOUR, SURFACTEN,
				POLYTROP, VISCOSITY);
	break;
    default:
      Error("No bubblemethod specified ",__FILE__,__LINE__);
    }

    // For test phase only two values are needed.
    // Need to be exchanged if bubbledynamic and acoustic is coupled then
    // y_[0] and y_[1] each will become a vector of length number of elements
    y_.Push_back( RADIUSINIT );
    y_.Push_back( VELINIT );


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

    // Output procedure needs to be adapted to later proceeding!!!!
    // Later only values for choosen elements should be printed
    // ----------------------------------------------------------------
    fp = fopen( "bubblevalues", "w" );
    if ( fp == NULL ) {
      //     Char *errmsg = NULL;
      //   NewArray( errmsg, Char, strlen(bubblevalues)+40 );
      //   errmsg++;
      //   sprintf( errmsg, "Cannot open file %s for writing!", bubblevalues );
      Error( "Could not open outputfile", __FILE__, __LINE__ );
    }

    // Write file header
    fprintf( fp, "# Time\t\tPressure\t\tRadius\t\tVelocity\n" );

    // First line with initial values
    // Still some print problems
    //    Double aux = 1.0e5;
    // fprintf( fp, "%e\t%e\t%16.10e\t%16.10e\n", 0, aux, y_[0], y_[1] );
    // std::cout << 0 << " " << aux << " " << y_[0] << std::endl;

    fprintf( fp, "0\t1.0e5\t%16.10e\t%16.10e\n", y_[0]/RADIUSINIT, y_[1] );
 
 


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

  // So far just a method to compute bubble behavour, 
  // later supposed to controll also coupling between 
  // bubbledynamic and acoustic
  void  BubbleDriver::SolveProblem() {
    ENTER_FCN( " BubbleDriver::SolveProblem" );

    //    Integer level     = 0;
    Double  steptime  = 0;
    Integer stepsave  = isavebegin_;

    Double  dt = firstdt_;

    Integer nstep;

    for (nstep = 1; nstep <= numstep_; nstep++){

      // forward the actual pressure values to BubbleODE object
      // in test phase read the values from input file 
      // or compute them directly

      pressure_ = - 1e05 * sin( 2 * PI * 500000 * steptime);
      dpresdt_  = - 1e05 * (2 * PI * 500000 ) * 
	cos( 2 * PI * 5000000 * steptime);

      ptBubble_->SetP(pressure_);
      ptBubble_->SetDpdt(dpresdt_);

      // In case of explicit Euler watch out suggested stepsize
      ptODESolver_->Solve( steptime, steptime+dt, y_, *ptBubble_, dt / 3.0,
			   0, dt);

      // Later the resulting values of y_ need to be parsed to acoustic
	  
      steptime += dt;

      // writing results of bubbledynamic in output file

      //      std::cerr << "bubble" << y_[0] / RADIUSINIT << y_[1] ;
      if (nstep == stepsave && (nstep <= isaveend_)) {
	fprintf( fp, "%e\t%e\t%16.10e\t%16.10e\n", steptime,pressure_ 
		 + PSTATIC,y_[0]/RADIUSINIT,y_[1]);
	stepsave+=isaveincr_;
      }
    }


    // Close output file
    if ( fclose( fp ) == EOF ) {
      //      Char *errmsg = NULL;
      //      NewArray( errmsg, Char, strlen(fname)+40 );
      //      errmsg++;
      //    sprintf( errmsg, "Could not close file %s after writing!", fname );
      Error( "Could not close file after writing", __FILE__, __LINE__ );
    }
  }

} // end of namespace
