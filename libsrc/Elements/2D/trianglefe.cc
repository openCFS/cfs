#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include <Elements/basefe.hh>
#include "trianglefe.hh"

namespace CoupledField
{

TriangleFE::TriangleFE()
{
#ifdef TRACE
  (*trace) << "entering TriangleFE::TriangleFE" << std::endl;
#endif
  
  Dim_ = 2;
  NumEdges_ = 3;
  NumFaces_ = 1;

  std::string integtype="GaussOrder2";

  std::string IntRule;
  if (conf->ifget("IntegRules", IntRule)==TRUE)
      conf->ifget("triangle",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());

}

TriangleFE :: ~TriangleFE()
{
#ifdef TRACE
  (*trace) << "entering TriangleFE::~TriangleFE" << std::endl;
#endif
}


void TriangleFE:: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering TriangleFE::SetIntPoints" << std::endl;
#endif
 
 switch(IntegType)
    {

    case GaussOrder2:

      NumIntPoints_=3;
      DegreeInteg_=2;

      
      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];

      IntWeights_.resize(NumIntPoints_);

      for(Integer i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].resize(Dim_);
      
      IntPoints_[0][0] = 0.166666666666667;
      IntPoints_[1][0] = 0.666666666666667; 
      IntPoints_[2][0] = 0.166666666666667;
      IntPoints_[0][1] = 0.166666666666667;
      IntPoints_[1][1] = 0.166666666666667;
      IntPoints_[2][1] = 0.666666666666667;

      IntWeights_[0]= 0.166666666666667 ;
      IntWeights_[1]= 0.166666666666667 ;
      IntWeights_[2]= 0.166666666666667 ;

      break;

    case GaussOrder3:

      NumIntPoints_ = 4;
      DegreeInteg_  = 3;

      
      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];

      IntWeights_.resize(NumIntPoints_);

      for(Integer i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].resize(Dim_);
      IntPoints_[0][0] = 1.0/3.0;
      IntPoints_[1][0] = 3.0/5.0;
      IntPoints_[2][0] = 1.0/5.0;
      IntPoints_[3][0] = 1.0/5.0;
      IntPoints_[0][1] = 1.0/3.0;
      IntPoints_[1][1] = 1.0/5.0;
      IntPoints_[2][1] = 3.0/5.0;
      IntPoints_[3][1] = 1.0/5.0;
      
      IntWeights_[0]=-0.5625;
      IntWeights_[1]=0.520833333333333;
      IntWeights_[2]=0.520833333333333;
      IntWeights_[3]=0.520833333333333;

      break;

   case GaussOrder4:

// D. A. Dunavant: High Degree Efficient Symmetrical ...
// Int. J. Numer. Methods in Eng.  Vol 21  S. 1129-1148   1985
// (c) Kaskade

     NumIntPoints_ = 6;
     DegreeInteg_  = 4;

      
      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];

      IntWeights_.resize(NumIntPoints_);

      for(Integer i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].resize(Dim_);

      IntPoints_[0][0] = 4.45948490915965e-01;
      IntPoints_[1][0] = 1.08103018168070e-01 ;
      IntPoints_[2][0] = 4.45948490915965e-01;
      IntPoints_[3][0] = 9.15762135097710e-02; 
      IntPoints_[4][0] = 8.16847572980459e-01; 
      IntPoints_[5][0] = 9.15762135097710e-02;  

      IntPoints_[0][1] = 4.45948490915965e-01;
      IntPoints_[1][1] = 4.45948490915965e-01;
      IntPoints_[2][1] = 1.08103018168070e-01;
      IntPoints_[3][1] = 9.15762135097710e-02; 
      IntPoints_[4][1] = 9.15762135097710e-02; 
      IntPoints_[5][1] = 8.16847572980459e-01;  

      IntWeights_[0]= 1.116907948390055e-01*2;
      IntWeights_[1]= 1.116907948390055e-01*2;
      IntWeights_[2]= 1.116907948390055e-01*2;
      IntWeights_[3]= 5.497587182766100e-02*2;
      IntWeights_[4]= 5.497587182766100e-02*2; 
      IntWeights_[5]= 5.497587182766100e-02*2; 
 
      break;

   case GaussOrder5:
 
// D. A. Dunavant: High Degree Efficient Symmetrical ...
// Int. J. Numer. Methods in Eng.  Vol 21  S. 1129-1148   1985
// (c) Kaskade
 
      NumIntPoints_ = 7;
      DegreeInteg_  = 5;
      
      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];

      IntWeights_.resize(NumIntPoints_);

      for(Integer i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].resize(Dim_);      
 
      IntPoints_[0][0] = 3.333333333333333e-01;
      IntPoints_[1][0] = 4.701420641051151e-01;
      IntPoints_[2][0] = 5.971587178976981e-02;
      IntPoints_[3][0] = 4.701420641051151e-01;
      IntPoints_[4][0] = 1.012865073234563e-01;
      IntPoints_[5][0] = 7.974269853530872e-01;
      IntPoints_[6][0] = 1.012865073234563e-01;      
 

      IntPoints_[0][1] = 3.333333333333333e-01;
      IntPoints_[1][1] = 4.701420641051151e-01;
      IntPoints_[2][1] = 4.701420641051151e-01;
      IntPoints_[3][1] = 5.971587178976981e-02;
      IntPoints_[4][1] = 1.012865073234563e-01;
      IntPoints_[5][1] = 1.012865073234563e-01;
      IntPoints_[6][1] = 7.974269853530872e-01;
 
      IntWeights_[0]= 0.2250300003;
      IntWeights_[1]= 0.132394152788506;
      IntWeights_[2]= 0.132394152788506;
      IntWeights_[3]= 0.132394152788506 ;
      IntWeights_[4]= 0.125939180544827;
      IntWeights_[5]=  0.125939180544827;
      IntWeights_[6]=  0.125939180544827; 

      break;
 
    default:
      std::cerr << "Integration type " << IntegType
           << " is not implemented \n" << std::endl; exit(-1);
    }

}

} // end of namespace
