#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

#include "Domain/grid.hh"
#include "Domain/GridCFS/interface_gridcfs.hh"
#include "Domain/GridStruct/interface_gridstruct.hh"
#include "timefunc.hh"
#include "WriteInfo.hh"
#include "ParamHandling/BaseParamHandler.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/domain.hh"

namespace CoupledField {

  TimeFunc::TimeFunc() {

    ENTER_FCN( "TimeFunc::TimeFunc" );

    maxnumTF_  = 0;

    timeFncDatFiles_=FALSE;
    std::string nametf;

    // Determine type of analysis and, if adaptivity is to be used
    AnalysisType analysisType;
    std::string  analysis;
    params->Get( "type", analysis, "analysis" );
    String2Enum( analysis, analysisType );

    if (analysisType == TRANSIENT4SLICE) {
      params->GetList( "name", fnc_names_, "transient4Slice", "timeDataFile" );
    } else {
      params->GetList( "name", fnc_names_, "transient", "timeDataFile" );
    }

    if (fnc_names_.GetSize())
      timeFncDatFiles_ = TRUE;

    //read in the time functions
    if (timeFncDatFiles_)  
      {
        Info->StartProgress("Reading in time data function(s)");
        ReadTimeFuncs();
        Info->FinishProgress();
      }

  }

  void TimeFunc :: ReadTimeFuncs()
  {
    ENTER_FCN( "TimeFunc::ReadTimeFuncs" );

    maxnumTF_ =  fnc_names_.GetSize();

    valTF_    =  new std::list<Double>[maxnumTF_];
    timeTF_     =  new std::list<Double>[maxnumTF_]; 
    std::string errMsg;



    // loop over time-fncs
    for ( UInt i = 0; i < maxnumTF_; i++ ) {

      std::ifstream timefile;

      if ((fnc_names_[i])=="spc_dependent_fnc")
        {
          std::cout<<"time-space dependent function '"
                   <<fnc_names_[i].c_str()<<"' chosen."<<std::endl;
        }
      else
        {
          timefile.open( fnc_names_[i].c_str() );
          if ( !timefile ) {
            (*error) << "Failed to open file '" << fnc_names_[i]
                     << "' containing data of time function";
            Error( __FILE__, __LINE__ );
          }

          timefile.clear(); // clear flags

          // we don't trust .eof() =)
          timefile.seekg(0,std::ios::end);
          std::string::size_type pos = 0, pos_end = timefile.tellg(),
            line_end_pos , line_start_pos = 0;

          timefile.seekg(0,std::ios::beg); // start from the beginning
          std::string     buf;
          Double          timeT, valT;

          while( pos <= pos_end ) {         
            buf = "";
            line_start_pos = timefile.tellg();
            std::getline(timefile,buf,'\n');
            line_end_pos = timefile.tellg();
        
            // big choice of signs for comment's
            if (buf.length() != 0 &&
                buf[0] != '#' &&
                buf[0] != '%' && 
                buf[0] != '!') {
              timefile.seekg(line_start_pos); // rewind
              timefile >> timeT >> valT ;                    

              valTF_[i].push_back(valT);
              timeTF_[i].push_back(timeT);
              timefile.ignore(100,'\n');
          
            }

            pos = timefile.tellg();  // and, where we are ?    
          
            if( pos != line_end_pos)
              {
                errMsg  = "The time data file '";
                errMsg += fnc_names_[i];
                errMsg += "' is not correctly formatted.\n";
                errMsg += "Please correct it!";
                Error(errMsg.c_str(), __FILE__, __LINE__);
              }
          
          }
      
          timefile.close();
        }
    } // loop over fncs
  }

  Double TimeFunc::TimeFuncAtTime(const Double time,  const std::string fncname)
  {
    ENTER_FCN( "TimeFunc::TimeFuncAtTime" );
  
    Integer numfnc;
 
    //if name of time function not defined, then a constant time function with value
    //1.0 is assumed 
    if (fncname == "none")
      return 1.0;
      
    //get correct time function
    numfnc = fnc_names_.Find(fncname);

    if (numfnc == -1)
      {
        std::string fncstring = "Time Function '";
        fncstring += fncname;
        fncstring += "' is not defined within 'timeDataFile' ";
        fncstring += "element in section 'Transient'!";
        Error(fncstring.c_str(), __FILE__, __LINE__);
      }

    if (time < *timeTF_[numfnc].begin() )
      Error("Wrong time in TimeFuncAtTime",__FILE__,__LINE__);

    //if time larger as defined in time function, return the last value
    if (time > timeTF_[numfnc].back()) 
      return valTF_[numfnc].back();

    //loop over array with time
    Double                              hlp , hlp1, hlp2, hlp_prev=0;
    Double                              v_prev =0;
    std::list<Double>::const_iterator   p,v,p_prev;

    for ( p=timeTF_[numfnc].begin(), v=valTF_[numfnc].begin(); p!= timeTF_[numfnc].end(); p++,v++)
      {

        hlp_prev = hlp;
        hlp      = *p;

        if (time == hlp) return *v;
        if (time < hlp) break;
     
        v_prev   = *v;
      }
 
    // linear interpolation
    hlp1 = hlp - hlp_prev;
    hlp2 = ((hlp-time)/hlp1)*v_prev+((time-hlp_prev)/hlp1)*(*v);
    return hlp2;
 
  }

  StdVector<Double> TimeFunc::TimeSpcFuncAtTime(const Double time, 
                                                const std::string fncname,
                                                StdVector<UInt> nodes, 
                                                Grid * ptgrid)
  {
    ENTER_FCN( "TimeFunc::TimeSpcFuncAtTime" );
  
    StdVector<Double> fncVal;
    fncVal.Resize(nodes.GetSize());
    Double PI=3.1415926535; 
    //Specific parameter values for space dependent function
      Double y0 = 0;
      Double r = 1.5;
      Double t0 = 0.5;
 
    for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
      Vector<Double> ptCoordNode;
      ptgrid->GetNodeCoordinate(ptCoordNode, nodes[iNode]);
      Double x, y, z;
      x = ptCoordNode[0];
      y = ptCoordNode[1];
      z = ptCoordNode[2];
      // Here an switch instruction could be used to add additional functions.

      //This case is for setting the analytical pressure fluctuations P_ak 
      //from vortex_analytical as problem solution
      Double r_sqr=(x*x+y*y);
      if (r_sqr<=((200./81.)*(200./81.)))
        {
          //  if (VortexFlag==2)
          //                             nodalval[ii]=0;
          //                           else if (VortexFlag==3)
          fncVal[iNode]=cos( (PI/(2*r)) * (y-y0) );
         }
       else
         {
           Double P_ak=0;
           Vector<Double> dTijdi;
           SinglePDE * acouPDE_ = domain->GetSinglePDE("acoustic");
           acouPDE_->VortexAnalytical(P_ak, dTijdi, x, y, time, 1);
           //std::cout<<"After getting P_ak in timefunc.cc, P_ak= "<<P_ak<<std::endl;
           fncVal[iNode]=P_ak;
         } 

      

// For benchmark test with paper
// Finite Element analysis of semi-infinite wave guides with high-order
// boundary treatment (Int. J. Numer. Meth. 2003, 58:1955-1983)
//       if ( (abs(y-y0) <= r) && (0<=time) && (time<=t0) )
//         {
//           fncVal[iNode]=cos( (PI/(2*r)) * (y-y0) );
//         }
//       else
//         {
//           fncVal[iNode]=0;
//         } 
    }
    return  fncVal;
  }
  


  TimeFunc :: ~TimeFunc()
  {
    ENTER_FCN( "TimeFunc::~TimeFunc" );
 
    if (maxnumTF_) 
      {
        UInt i;
        for(i=0; i<maxnumTF_; i++)
          {
            timeTF_[i].clear();
            valTF_[i].clear();
          }

        delete [] valTF_;
        delete [] timeTF_;

      }
  }

  void TimeFunc::Print(std::ostream * outfileDat) const
  {
    ENTER_FCN( "TimeFunc::Print" );
    (*outfileDat) << "------------- Print Time function ----------------" 
                  << std::endl;
    for (UInt i=0; i < maxnumTF_; i++)
      {
        (*outfileDat) << " Number of time function is " << i+1 << std::endl;
        (*outfileDat) << "\t" << "time" <<
          "\t" << "value" << std::endl;
        std::list<Double>::const_iterator   t,v;
        for (t=valTF_[i].begin(),v=valTF_[i].end(); t!=valTF_[i].end(); t++,v++ )
          (*outfileDat) << "\t" << *t << "\t" << *v << std::endl;
      } // 
 
  } // end of TimeFunc::Print


  /****************************************************************************
   Function for setting up the start time vector, it is:
   TASK: - built in an if condition for the dimension,
         - calculate the start time vector for the 3D-Case 
  ****************************************************************************/
  void TimeFunc::SetStartTimeVector(Integer numBcs)
  {
    ENTER_FCN( "TimeFunc::SetStartTimeVector");

    Double alpha;
    Double meshSize;
    Double nL;
    Double distToFocus;
    Double deltaT = 0.0;
    Double x;
    Integer nuOfBcsNodes;
    Integer i;
    Double soundSpeed;
    Double waveLength;
    Integer elementsPerWavelength;
    Double dimx, dimy;
    Integer dim_;
    Integer maxnumelemy_, maxnumelemx_;

    //Get slice data
    std::string dimension;
    std::string pdename = "acoustic";
    StdVector<std::string> keyVec;
    // NEW: Read in the Dimension of the Problem
    dim_ = 2;
    params->Get("type", dimension, "geometry");
    if(dimension == "3d")
      dim_ = 3;
    //std::cerr << "Dimension in the timefunc.cc file: " << dim_ << std::endl;
    keyVec  = pdename, "sliceData", "soundSpeed";
    params->Get( keyVec, soundSpeed);
    //std::cout << "sound" << soundSpeed << std::endl;
    keyVec  = pdename, "sliceData", "dimension", "dimY";
    params->Get( keyVec, dimy);
    keyVec  = pdename, "sliceData", "dimension", "dimX";
    params->Get( keyVec, dimx);
    //std::cout << "dim" << dim << std::endl;
    keyVec  = pdename, "sliceData", "waveLength";
    params->Get( keyVec, waveLength);
    //std::cout << "wave" << waveLength << std::endl;
    keyVec  = pdename, "sliceData", "elementsPerWavelength";
    params->Get( keyVec, elementsPerWavelength);
    //std::cout << "elemperwave" << elementsPerWavelength << std::endl;

    if(numBcs!=0) {
      meshSize = waveLength/elementsPerWavelength;
      // std::cout << "meshsize   " << meshSize << std::endl;

      //Size of array
      nuOfBcsNodes = numBcs;
  
      distToFocus = 0.05; //Kann der einfach so festglegt werden??
      /*
       * Now the switchin for either 3D-Case or 2D-Case
       */  
      switch(dim_){
      case 2: //The 2-Dimensional Case
	// std::cout << "nL " << nL << "   nuOfBcsBNodes " << nuOfBcsNodes <<std::endl;
	nL = (nuOfBcsNodes-1.0)* meshSize;
	alpha = atan(nL/distToFocus) * (180 /PI);
	// std::cout << "alpha: " << alpha << "    ";

	// Init StartTime
	startTime_.Resize(nuOfBcsNodes);
	startTime_[0] = 0.0;
	startTime_[nuOfBcsNodes-1] = 0.0;
	x=nL;

	for(i=2;i<nuOfBcsNodes+1;i++) {
	  Double sinus;
	  sinus = sin(alpha*(PI/180));
	  // std::cout << "sinus " << sinus <<"  ";
	  deltaT += (meshSize*(sin(alpha*(PI/180))))/soundSpeed;
	  //  std::cout << "deltaT " << deltaT <<"  ";
	  startTime_[nuOfBcsNodes -i]=deltaT;
	  x=x-meshSize;
	  if(x>0) {
	    alpha = atan(x/distToFocus) * (180/PI);
	    //std::cout << "alpha neu " << alpha << std::endl;
	  }
	}
	// Double test=0.0;
	// for (i=0;i<nuOfBcsNodes;i++) {
	// startTime_[i]=test;
	// test+=2*2E-8;
	// std::cout << "startTime ["<<i<<"]: " << startTime_[i] << "   " ;
	// std::cout << "startpoint["<<i<<"]: " << startTime_[i]/2E-8 << std::endl;
	// }
	break;

      case 3: //The 3-Dimensional Case
	//std::cerr << "Bin jetzt in set timefunc" << std::endl;
	//first of all calculate the maximum number of elements for
	// x and y direction;
	Integer j;
	maxnumelemy_ = (Integer) (dimy/meshSize);
	maxnumelemx_ = (Integer) (dimx/meshSize);
	//if(nuOfBcsNodes == (maxnumelemy_+1)*(maxnumelemx_+1)){
	//std::cerr << "Welcome to timefunc.cc 3D" << std::endl;
	// The mastertime, is the node, that needs the longest time to reach the focus
	startTime_.Resize((maxnumelemx_+1)*(maxnumelemy_+1));
	Integer abs1 = (Integer) (maxnumelemx_+1) - (maxnumelemx_+1)/2;
	//std::cout << "abs1" << abs1 << std::endl;
	Integer abs2 = (Integer) (maxnumelemy_+1) - (maxnumelemx_+1)/2;
	//Insert due to the testing the Channel Case
	if(maxnumelemy_ ==1 && maxnumelemx_ ==1){
	  abs1 = 1;
	  abs2 = 1;
	}
	Double mastertime = sqrt(((abs1*abs1)+(abs2*abs2))*(meshSize*meshSize)
				 +distToFocus*distToFocus)/soundSpeed;
	//std::cout << "mastertime: " << mastertime << std::endl;
	//Loop over all elements of the Boundary Slice
	Integer k = 0;
	Double temp;
	for(i=0; i<maxnumelemx_+1; i++){
	  for(j=0; j<maxnumelemy_+1; j++){
	    /* The startTime_ vector can't be simply set up linear any more.
	     * The element-numbering is not linear anymore. Clearly spoken, 
	     * it isn't not like it was in 2D that: the bigger the node-number
	     * the smaller the Delta.
	     */		
	    //the inhom. Dirichlet nodes are set to its delta Values
	    if(i<abs1 && j<abs2){
	      temp = sqrt((((i+1)*(i+1))+((j+1)*(j+1)))*(meshSize*meshSize)
			  +distToFocus*distToFocus)/soundSpeed;
	      startTime_[k] = mastertime - temp;
	      //std::cout << "Delta Time at node " << k << " = " <<  startTime_[k] << std::endl;					
	    }
	    else{
	      startTime_[k] = -1;
	    }
	    k++;
	  }
	}
	//}else{
	//std::cerr << "Error in timefunc.cc" << std:: endl;
	//std::cerr << "nuOfBcsNodes was: " << nuOfBcsNodes << std::endl;
	//std::cerr << "Actual BC-Nodes is: " << (maxnumelemy_+1)*(maxnumelemx_+1)<<std::endl;
	//exit(0);
	//}			
	break;
      }
    }
  }


} // end of namespace
