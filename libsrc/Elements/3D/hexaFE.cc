#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include "hexaFE.hh"

namespace CoupledField
{

HexaFE::HexaFE()
{
#ifdef TRACE
    (*trace) << "entering HexaFE::HexaFE" << std::endl;
#endif

    Dim_      = 3;
    NumEdges_ = 12;
    NumFaces_ = 6;
    numChilds_ = 8;

  std::string integtype="GaussOrder2";

    std::string IntRule;
    if (conf->ifget("IntegRules", IntRule)==TRUE)
      conf->ifget("hexa", integtype, "IntegRules");

    IntegType=String2EnumIntegrationType(integtype.c_str());

}

HexaFE::~HexaFE()
{
#ifdef TRACE
    (*trace) << "entering HexaFE::~HexaFE" << std::endl;
#endif
    //  if (IntWeights) delete IntWeights;
}

void HexaFE::SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering HexaFE::SetIntPoints" << std::endl;
#endif
 
  switch(IntegType) 
    {
    case GaussOrder1:
	  
      NumIntPoints_ = 1;
      DegreeInteg_  = 2;


      IntWeights_.resize(NumIntPoints_);
	// all weights are 1.0
      for(Integer i=0; i<IntWeights_.size(); i++)
	IntWeights_[i] = 1.0;
      
      if (!IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];
      
      for(Integer i=0; i<NumIntPoints_; i++)
	IntPoints_[i].resize(Dim_);

      IntPoints_[0][0] = 0;
      IntPoints_[0][1] = 0;
      IntPoints_[0][2] = 0;      
      
      break;

    case GaussOrder2:
	  
      NumIntPoints_=8;
      DegreeInteg_=3;


	IntWeights_.resize(NumIntPoints_);
	// all weights are 1.0
	for(Integer i=0; i<IntWeights_.size(); i++)
	  IntWeights_[i] = 1.0;
      
	if (!IntPoints_)
	  IntPoints_ = new std::vector<Double>[NumIntPoints_];

	for(Integer i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].resize(Dim_);

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

