#include <iostream>
#include <fstream>

#include "pierceInt.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField {


  PierceDampInt::PierceDampInt(const Double aFactor, SimpleFlow * flow,
			       bool axi, bool coordUpdate )
    : BaseForm(NULL), 
      dampFactor_(aFactor), flowHandler_(flow)
      
  {
    ENTER_FCN( "PierceDampInt::PierceDampInt" );
    name_ = "PierceDampInt";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    baseType_ = DAMPING;
  }

 
  PierceDampInt::~PierceDampInt()
  {
    ENTER_FCN( "PierceDampInt::~PierceDampInt" );
  }

 

  void PierceDampInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1, 
                                   EntityIterator& ent2  )  {
    ENTER_FCN( "PierceDampInt::CalcElemMatrix" );
   
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet, factor;

    UInt dim = ptCoord_.GetSizeRow();

    Matrix<Double> xiDx;
    Vector<Double> shapeFncAtIp;
    Vector<Double> partVec;
    Vector<Double> coordAtIp;
    Vector<Double> velAtIp(dim);
    Vector<Double> velDerivAtIp(dim);

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes);
    elemMat.Init();
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt);

      // compute velocity and derivative at IP        

      coordAtIp = ptCoord_ * shapeFncAtIp;
      flowHandler_->ComputeActFlow(coordAtIp, velAtIp, velDerivAtIp); 
      //      std::cout << "velAtIp:\n" << velAtIp << std::endl;

      partVec = xiDx * velAtIp;
        
//       std::cout << "velAtIP\n " << velAtIp << std::endl;
//       std::cout << "xiDx:\n " << xiDx << std::endl;
//       std::cout << "partVec:\n" << partVec << std::endl;
//       std::cout << "Ni:\n" << shapeFncAtIp << std::endl;

      if (isaxi_) 
        factor = 2 * PI * intWeights[actIntPt-1] * dampFactor_* jacDet * coordAtIp[0];
      else 
        factor = intWeights[actIntPt-1] * dampFactor_ * jacDet;
        
      for (UInt i=0; i<nrNodes; i++) {
	for (UInt j=0; j<nrNodes; j++) {
	  elemMat[i][j] += shapeFncAtIp[i] * partVec[j] * factor;
	}
      }
      //      std::cout << "ElemMatDamp:\n" << elemMat << std::endl;
    }
  
    //    std::cout << "ElemMatDamp:\n" << elemMat << std::endl;
  }



  //================================== Stiffness part ===========================

  PierceStiffInt::PierceStiffInt(const Double aFactor, SimpleFlow* flow, 
				 bool axi, bool coordUpdate )
    : BaseForm(NULL), 
      stiffFactor_(aFactor), flowHandler_(flow)
      
  {
    ENTER_FCN( "PierceStiffInt::PierceStiffInt" );
    name_ = "PierceDampInt";
    isaxi_ = axi;
    coordUpdate_ = coordUpdate;
    baseType_ = STIFFNESS;
  }

 
  PierceStiffInt::~PierceStiffInt()
  {
    ENTER_FCN( "PierceStiffInt::~PierceStiffInt" );
  }

 

  void PierceStiffInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                   EntityIterator& ent1, 
                                   EntityIterator& ent2  )  {
    ENTER_FCN( "PierceStiffInt::CalcElemMatrix" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const UInt dim = ptCoord_.GetSizeRow();


    Double jacDet, factor;

    Matrix<Double> xiDx;
    Vector<Double> shapeFncAtIp;
    Vector<Double> partVec;
    Vector<Double> coordAtIp;

    Vector<Double> velAtIp(dim);
    Vector<Double> velDerivAtIp(dim);

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes);
    elemMat.Init();
    

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);
      ptelem->GetShFncAtIp(shapeFncAtIp, actIntPt);
        
      // compute velocity and derivative at IP        
      coordAtIp = ptCoord_ * shapeFncAtIp;

      flowHandler_->ComputeActFlow(coordAtIp, velAtIp, velDerivAtIp); 
//       std::cout << "velAtIp:\n" << velAtIp << std::endl;
//       std::cout << "velDerivAtIp:\n" << velDerivAtIp << std::endl;

      //Ni,x*v_x + Ni,y*v_y + Ni,z*v_z
      partVec      = xiDx * velAtIp;
        
      if (isaxi_) 
        factor = 2 * PI * intWeights[actIntPt-1] * stiffFactor_ * 
	  jacDet * coordAtIp[0];
      else 
        factor = intWeights[actIntPt-1] * stiffFactor_ * jacDet;
        
      //divergence of velocity
      Double divVel = 0.0;
      for ( UInt i=0; i<dim; i++) 
	divVel += velDerivAtIp[i];

      for (UInt i=0; i<nrNodes; i++) {
	Double part1 =  partVec[i] + shapeFncAtIp[i] * divVel; 
	for (UInt j=0; j<nrNodes; j++) {
	  elemMat[i][j] += part1 * partVec[j] * factor;
	}
      }
      //      std::cout << "ElemMatDamp:\n" << elemMat << std::endl;
    }
  
    // std::cout << "ElemMatStiff:\n" << elemMat << std::endl;
  }


  //================================= class for handling flow =====================


  SimpleFlow::SimpleFlow() {

    ENTER_FCN( "SimpleFlow::SimpleFlow" );

  }

 
  SimpleFlow::~SimpleFlow()
  {
    ENTER_FCN( "SimpleFlow::~SimpleFlow" );
  }


  // ***********************************************************************
  //   Obtain information for flow data
  // ***********************************************************************
  void SimpleFlow::ReadFlowData( std::string pdeName, std::string regionName,
				 UInt dim) {

    ENTER_FCN( "AcousticPDE::ReadFlowData" );

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // Construct vectors for restricted parameter search
    keyVec = pdeName , "region" , "flowData" , "flowDir";
    attrVec= ""         , "name"   , "";
    valVec = ""         , regionName, "";

    StdVector<std::string> stringVal;
    params->GetList( keyVec, attrVec, valVec, stringVal);
    if ( stringVal[0] == "x" ) 
      flowDir_ = X;
    else if ( stringVal[0] == "y" ) 
      flowDir_ = Y;
    else if ( stringVal[0] == "z" ) {
      if ( dim = 2) {
	Error("Direction of flow in 2D can just be x or y",__FILE__,__LINE__);
      }
      flowDir_ = Z;
    }
    else 
      Error("Direction of flow data not valid",__FILE__,__LINE__);

    StdVector<UInt> ival;
    keyVec = "acoustic" , "region" , "flowData" , "flowOrd";
    params->GetList( keyVec, attrVec, valVec, ival);
    flowOrder_ = ival[0];

    StdVector<Double> val;
    keyVec = "acoustic" , "region" , "flowData" , "flowVel";
    params->GetList( keyVec, attrVec, valVec, val);
    flowVelocity_ = val[0];

    keyVec = "acoustic" , "region" , "flowData" , "flowDia";
    params->GetList( keyVec, attrVec, valVec, val);
    flowRadius_ = val[0] / 2.0;

    keyVec = "acoustic" , "region" , "flowData" , "flowLen";
    params->GetList( keyVec, attrVec, valVec, val);
    flowLength_ = val[0];
 
    keyVec = "acoustic" , "region" , "flowData" , "flowDecay";
    params->GetList( keyVec, attrVec, valVec, val);
    flowDecay_ = val[0];

    flowCenter_.Resize(dim);
    keyVec = "acoustic" , "region" , "flowData" , "flowCenterX";
    params->GetList( keyVec, attrVec, valVec, val);
    flowCenter_[0] = val[0];

    keyVec = "acoustic" , "region" , "flowData" , "flowCenterY";
    params->GetList( keyVec, attrVec, valVec, val);
    flowCenter_[1] = val[0];
    if ( dim == 3) {
      keyVec = "acoustic" , "region" , "flowData" , "flowCenterZ";
      params->GetList( keyVec, attrVec, valVec, val);
      flowCenter_[2] = val[0];
    }

  }


  //! set flow data
  void SimpleFlow::ComputeActFlow( const Vector<Double> coord, 
				   Vector<Double>& velocity, 
				   Vector<Double>& derivVel) {

    ENTER_FCN( "SimpleFlow::ComputeActFlow" );

    //help variables 
    Double len1, len2, xa, pos, der, der2, omg, arg, factor, vel;

    UInt dim = coord.GetSize();
    Double dx = coord[0] - flowCenter_[0];
    Double dy = coord[1] - flowCenter_[1];
    Double dz = 0.0;
    if ( dim == 3) 
      dz = coord[2] - flowCenter_[2];

    //    std::cout << "dx=" << dx << " dy=" << dy << std::endl;

    len1 = flowLength_*0.5;
    len2 = len1 + flowDecay_;

    velocity.Init();
    derivVel.Init();

    if ( flowDir_ == X ) {
       xa  = sqrt(dy*dy + dz*dz);
       pos = abs(dx);
    }
    else if ( flowDir_ == Y ) {
       xa  = sqrt(dx*dx + dz*dz);
       pos = abs(dy);
    }
    else if ( flowDir_ == Z ) {
       xa  = sqrt(dx*dx + dy*dy);
       pos = abs(dz);
    }

    /* factor according to length of flow and position */
    factor = 0.0;
    der2   = 0.0;
    if ( pos <= len2 ) {
      if (pos <= len1) {
	factor = flowVelocity_;
	der2   = 0.0;
      }
      else {
	omg    = 4.0 * std::atan(1.0) / (2.0*flowDecay_);
	arg    = omg * ( pos-len1 );
	factor = flowVelocity_ * std::pow(std::cos(arg),2.0);
	/* der2 = -(*flowmax)*omg*sin(2.0*arg); */
	der2 = 0.0;
      }
    }

    if ( xa < flowRadius_ ) {
      if ( flowOrder_ == 0 ) {
	vel = factor;
	der = 0.0;
      } 
      else {
	arg = xa * xa/ ( flowRadius_ * flowRadius_ );
	vel = std::pow( 1.0 - arg, 1./Double(flowOrder_) );
	der = -2.*factor/ ( Double(flowOrder_) * flowRadius_ 
			    * flowRadius_) * vel / ( 1.0-arg );
	vel *= factor;
      }
      
      if ( flowDir_ == X ) {
	velocity[0] = vel;
	derivVel[0] = der2*vel;
      }
      else if ( flowDir_ == Y ) {
	velocity[1] = vel;
	derivVel[1] = der2*vel;
      }
      else if ( flowDir_ == Z ) {
	velocity[2] = vel;
	if ( dim == 3) 
	  derivVel[2] = der2*vel;
      }
    }
    
  }
  
} // end namespace CoupledField



