#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "wedgeFE.hh"

namespace CoupledField
{

WedgeFE::WedgeFE()
{
  ENTER_FCN( "WedgeFE::WedgeFE" );

  Dim_      = 3;
  NumEdges_ = 9;
  NumFaces_ = 5;
  //    numChilds_ = 8;
  
#ifndef XMLPARAMS
  std::string integtype = "GaussOrder2";
  std::string IntRule;
  if (conf->ifget("IntegRules", IntRule)==TRUE)
    conf->ifget("wedge", integtype, "IntegRules");
#else
  std::string integtype;
  params->Get( "type", integtype, "integRules", "wedge" );
#endif
  
  IntegType = String2EnumIntegrationType(integtype.c_str());
  
  //  isSetAtCenter_=FALSE;
}
  
WedgeFE::~WedgeFE()
{
  ENTER_FCN( "WedgeFE::~WedgeFE" );
}

void WedgeFE::SetIntPoints()
{
  ENTER_IFCN( "WedgeFE::SetIntPoints" );

    switch(IntegType)
      {
      case GaussOrder1:

	Error("Integration type GaussOrder1 is not implemented for wedges",__FILE__,__LINE__);
	break;

      case GaussOrder2:
   // From Finite Element Library. www.mathsoft.cse.clrc.ac.uk/felib3/Docs/html/Level-0/wdg6/wdg6.html.
 // corner nodes
	NumIntPoints_ = 6;
	DegreeInteg_  = 3;

	IntWeights_.Resize(NumIntPoints_);
	// weights are different for the lower and upper integration points
	for(Integer i=0; i<IntWeights_.GetSize(); i++)
	  {
	    IntWeights_[i] = 0.25*1.732050808;
	  }
	
	if (!IntPoints_)
	  IntPoints_ = new Vector<Double>[NumIntPoints_];

	for(Integer i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].Resize(Dim_);

	IntPoints_[0][0] =  1.0;
	IntPoints_[0][1] =  0.0;
	IntPoints_[0][2] =  -1.0;
	
	IntPoints_[1][0] = -0.5;
	IntPoints_[1][1] = -0.5*1.732050808;
	IntPoints_[1][2] = -1.0;

	IntPoints_[2][0] = -0.5;
	IntPoints_[2][1] =  0.5*1.732050808;
	IntPoints_[2][2] = -1.0;
	
	IntPoints_[3][0] =  1.0;
	IntPoints_[3][1] =  0.0;
	IntPoints_[3][2] =  1.0;
	
	IntPoints_[4][0] = -0.5;
	IntPoints_[4][1] = -0.5*1.732050808;
	IntPoints_[4][2] =  1.0;

	IntPoints_[5][0] = -0.5;
	IntPoints_[5][1] =  0.5*1.732050808;
	IntPoints_[5][2] =  1.0;




	break;
      
      default:
	Error("Integration type is not implemented",__FILE__,__LINE__);
      }  
  }


} // end of namespace
