#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include <Elements/basefe.hh>
#include "rectanglefe.hh"

namespace CoupledField
{

RectangleFE::RectangleFE()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::RectangleFE" << std::endl;
#endif
  
  Dim_ = 2;
  NumEdges_ = 4;
  NumFaces_ = 1;
  numChilds_ = 4;
  
  std::string integtype="GaussOrder2";

  std::string IntRule;
  if (conf->ifget("IntegRules", IntRule)==TRUE)
      conf->ifget("rectangle",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());

}

RectangleFE :: ~RectangleFE()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::~RectangleFE" << std::endl;
#endif
}


void RectangleFE:: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::SetIntPoints" << std::endl;
#endif
 

  switch(IntegType) 
    {
    case GaussOrder1:
      
      NumIntPoints_ = 1;
      DegreeInteg_  = 1;
      
      
      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];
      
      IntWeights_.resize(NumIntPoints_);
      
      for(Integer i=0; i<NumIntPoints_; i++)
	{
	  IntWeights_[i]=1;
	  IntPoints_[i].resize(Dim_);
	}
      
      IntPoints_[0][0] = 0;
      IntPoints_[0][1] = 0;
      break;
      
    case GaussOrder2:

      NumIntPoints_=4;
      DegreeInteg_=2;

      
      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];

      IntWeights_.resize(NumIntPoints_);

      for(Integer i=0; i<NumIntPoints_; i++)
	{	  
	  IntWeights_[i]=1;
	  IntPoints_[i].resize(Dim_);
	}
      
      IntPoints_[0][0] = -0.57735026919;
      IntPoints_[1][0] =  0.57735026919;
      IntPoints_[2][0] =  0.57735026919;
      IntPoints_[3][0] = -0.57735026919;

      IntPoints_[0][1] = -0.57735026919;
      IntPoints_[1][1] = -0.57735026919;
      IntPoints_[2][1] =  0.57735026919;
      IntPoints_[3][1] =  0.57735026919;


  

      break;

    case GaussOrder5:

      NumIntPoints_=9;
      DegreeInteg_=5;

      if( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];
      
      for(Integer i=0; i<NumIntPoints_; i++)
	IntPoints_[i].resize(Dim_);

      IntWeights_.resize(NumIntPoints_);

      IntPoints_[0][0] = -0.774596669241483;
      IntPoints_[1][0] =  0.0;
      IntPoints_[2][0] =  0.774596669241483;
      IntPoints_[3][0] = -0.774596669241483;
      IntPoints_[4][0] = 0.0;
      IntPoints_[5][0] =  0.774596669241483;
      IntPoints_[6][0] = -0.774596669241483;
      IntPoints_[7][0] = 0.0;
      IntPoints_[8][0] =  0.774596669241483; 

      IntPoints_[0][1] = -0.774596669241483;
      IntPoints_[1][1] = -0.774596669241483;
      IntPoints_[2][1] = -0.774596669241483;
      IntPoints_[3][1] = 0.0;
      IntPoints_[4][1] = 0.0;
      IntPoints_[5][1] = 0.0;
      IntPoints_[6][1] =  0.774596669241483;
      IntPoints_[7][1] =  0.774596669241483;
      IntPoints_[8][1] =  0.774596669241483;

      IntWeights_[0]= 0.308642;
      IntWeights_[1]= 0.493827;
      IntWeights_[2]= 0.308642;
      IntWeights_[3]= 0.493827;
      IntWeights_[4]= 0.790123; 
      IntWeights_[5]= 0.493827; 
      IntWeights_[6]= 0.308642; 
      IntWeights_[7]= 0.493827; 
      IntWeights_[8]= 0.308642; 

      break;

    case GaussOrder7:

      Error("Type of integration Gauss with order 7 is incorrect", __FILE__, __LINE__);
      NumIntPoints_=16;
      DegreeInteg_=7;
      if ( !IntPoints_) 
	IntPoints_ = new std::vector<Double>[NumIntPoints_];

      for(Integer i=0; i<NumIntPoints_; i++)
	IntPoints_[i].resize(Dim_);
      
      IntPoints_[0][0] = -0.861136311594053;
      IntPoints_[1][0] =  -0.339981043584856;
      IntPoints_[2][0] =  0.339981043584856;
      IntPoints_[3][0] = 0.861136311594053;
      IntPoints_[4][0] = -0.861136311594053;
      IntPoints_[5][0] = -0.339981043584856;
      IntPoints_[6][0] =   0.339981043584856;
      IntPoints_[7][0] =  0.861136311594053;
      IntPoints_[8][0] =  -0.861136311594053;
      IntPoints_[9][0] =  -0.339981043584856;
      IntPoints_[10][0] =  0.339981043584856;
      IntPoints_[11][0] =  0.861136311594053;
      IntPoints_[12][0] =  -0.861136311594053;
      IntPoints_[13][0] =  -0.339981043584856;
      IntPoints_[14][0] =  0.339981043584856;
      IntPoints_[15][0] =  0.861136311594053;


      IntPoints_[0][1] =  -0.861136311594053;
      IntPoints_[1][1] =   -0.861136311594053;
      IntPoints_[2][1] =  -0.861136311594053;
      IntPoints_[3][1] =  -0.861136311594053;
      IntPoints_[4][1] =  -0.339981043584856;
      IntPoints_[5][1] =  -0.339981043584856;
      IntPoints_[6][1] =  -0.339981043584856;
      IntPoints_[7][1] =  -0.339981043584856;
      IntPoints_[8][1] =  0.339981043584856;
      IntPoints_[9][1] =  0.339981043584856;
      IntPoints_[10][1] =  0.339981043584856;
      IntPoints_[11][1] =  0.339981043584856;
      IntPoints_[12][1] =  0.861136311594053;
      IntPoints_[13][1] =  0.861136311594053;
      IntPoints_[14][1] =  0.861136311594053;
      IntPoints_[15][1] =  0.861136311594053;

      IntWeights_[0]= 0.121003;
      IntWeights_[1]= 0.226852;
      IntWeights_[2]= - 0.226852;
      IntWeights_[3]= - 0.121003;
      IntWeights_[4]= 0.226852;
      IntWeights_[5]= 0.425293;
      IntWeights_[6]= -0.425293;
      IntWeights_[7]= -0.226852;
      IntWeights_[8]= -0.226852;
      IntWeights_[9]= -0.425293;
      IntWeights_[10]= 0.425293;
      IntWeights_[11]= 0.226852;
      IntWeights_[12]= -0.121003;
      IntWeights_[13]= -0.226852;
      IntWeights_[14]= 0.226852;
      IntWeights_[15]= 0.121003;

      break;
 
    default:
      std::cerr << "Integration type " << IntegType 
	   << " is not implemented\n" << std::endl; 
      exit (-1);
    }
}

void RectangleFE::CalcSize(Vector<Double> &size, Array<Double> &coordinates)
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::CalcSize" << std::endl;
#endif

}

Double RectangleFE::CalcDistortion(Matrix<Double> &cornerCoords, Vector<Double> &size, Array<Double> &displacements)
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::CalcSize" << std::endl;
#endif
}

} // end of namespace
