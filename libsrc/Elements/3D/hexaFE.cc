#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "hexaFE.hh"

namespace CoupledField
{

HexaFE::HexaFE()
{
  ENTER_FCN( "HexaFE::HexaFE" );

    Dim_      = 3;
    NumEdges_ = 12;
    NumFaces_ = 6;
    numChilds_ = 8;

#ifndef XMLPARAMS
    std::string integtype="GaussOrder2";
    std::string IntRule;
    if (conf->ifget("IntegRules", IntRule)==TRUE)
      conf->ifget("hexa", integtype, "IntegRules");
#else
    std::string integtype;
    params->Get( "type", integtype, "integRules", "hexa" );
#endif

    IntegType=String2EnumIntegrationType(integtype.c_str());

}

HexaFE::~HexaFE()
{
  ENTER_FCN( "HexaFE::~HexaFE" );
  //  if (IntWeights) delete IntWeights;
}

void HexaFE::SetIntPoints()
{
  ENTER_FCN( "HexaFE::SetIntPoints" );
 
  switch(IntegType) 
    {
    case GaussOrder1:
	  
      NumIntPoints_ = 1;
      DegreeInteg_  = 2;


      IntWeights_.Resize(NumIntPoints_);
	// all weights are 1.0
      for(Integer i=0; i<IntWeights_.GetSize(); i++)
	IntWeights_[i] = 1.0;
      
      if (!IntPoints_)
	IntPoints_ = new Vector<Double>[NumIntPoints_];
      
      for(Integer i=0; i<NumIntPoints_; i++)
	IntPoints_[i].Resize(Dim_);

      IntPoints_[0][0] = 0;
      IntPoints_[0][1] = 0;
      IntPoints_[0][2] = 0;      
      
      break;

    case GaussOrder2:
	  
      NumIntPoints_=8;
      DegreeInteg_=3;


	IntWeights_.Resize(NumIntPoints_);
	// all weights are 1.0
	for(Integer i=0; i<IntWeights_.GetSize(); i++)
	  IntWeights_[i] = 1.0;
      
	if (!IntPoints_)
	  IntPoints_ = new Vector<Double>[NumIntPoints_];

	for(Integer i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].Resize(Dim_);

      IntPoints_[0][0] = -0.57735026919;
      IntPoints_[0][1] = -0.57735026919;
      IntPoints_[0][2] = -0.57735026919;
      
      IntPoints_[1][0] = 0.57735026919;
      IntPoints_[1][1] = -0.57735026919;
      IntPoints_[1][2] = -0.57735026919;

      IntPoints_[2][0] = 0.57735026919;
      IntPoints_[2][1] = 0.57735026919;
      IntPoints_[2][2] = -0.57735026919;

      IntPoints_[3][0] = -0.57735026919;
      IntPoints_[3][1] = 0.57735026919;
      IntPoints_[3][2] = -0.57735026919;

      IntPoints_[4][0] = -0.57735026919;
      IntPoints_[4][1] = -0.57735026919;
      IntPoints_[4][2] = 0.57735026919;
      
      IntPoints_[5][0] = 0.57735026919;
      IntPoints_[5][1] = -0.57735026919;
      IntPoints_[5][2] = 0.57735026919;

      IntPoints_[6][0] = 0.57735026919;
      IntPoints_[6][1] = 0.57735026919;
      IntPoints_[6][2] = 0.57735026919;

      IntPoints_[7][0] = -0.57735026919;
      IntPoints_[7][1] = 0.57735026919;
      IntPoints_[7][2] = 0.57735026919;

      break;

   default:
      Error("Integration type is not implemented",__FILE__,__LINE__);
   }
}

} // end of namespace

