#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "tetraFE.hh"

namespace CoupledField
{

TetraFE::TetraFE()
{
  ENTER_FCN( "TetraFE::TetraFE" );

    Dim_      = 3;
    NumEdges_ = 6;
    NumFaces_ = 4;
    NumCorners_ = 4;
    numChilds_ = 8;
    MidPoint_.Resize(3);
    MidPoint_[0] = 1./4;
    MidPoint_[1] = 1./4;
    MidPoint_[2] = 1./4;
    
    std::string integtype;
    params->Get( "type", integtype, "integRules", "tetra" );

    IntegType = String2EnumIntegrationType(integtype.c_str());

    //  isSetAtCenter_=false;
  }

TetraFE::~TetraFE()
{
  ENTER_FCN( "TetraFE::~TetraFE" );
}

void TetraFE::SetIntPoints()
{
  ENTER_IFCN( "TetraFE::SetIntPoints" );

    switch(IntegType)
      {
      case GaussOrder1:
	// Thomas Hughes "The finite element method", p. 174
    
	NumIntPoints_ = 1;
	DegreeInteg_  = 1;

	if ( !IntPoints_)
	  IntPoints_ = new Vector<Double>[NumIntPoints_];
      
	IntWeights_.Resize(NumIntPoints_);


	// just one point, but keep the procedure general ... ;O)
	for(UInt i=0; i<NumIntPoints_; i++)
	  {	  
	    IntWeights_[i] = 1;

	    IntPoints_[i].Resize(Dim_);

	    // all integration coords are at 1/4
	    for(UInt j=0; j<Dim_; j++)
	      IntPoints_[i][j] = 1./4;
	  }

	CorrectIntWeights();

	break;

      case GaussOrder2:
	// Thomas Hughes "The finite element method", p. 174   

	NumIntPoints_ = 4;
	DegreeInteg_  = 3;


	IntWeights_.Resize(NumIntPoints_);
	// all weights are 0.25
	for(UInt i=0; i<IntWeights_.GetSize(); i++)
	  IntWeights_[i] = 0.25;


	if (!IntPoints_)
	  IntPoints_ = new Vector<Double>[NumIntPoints_];

	for(UInt i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].Resize(Dim_);
      
	IntPoints_[0][0]=0.5854102;
	IntPoints_[0][1]=0.1381966;
	IntPoints_[0][2]=0.1381966;

	IntPoints_[1][0]=0.1381966;
	IntPoints_[1][1]=0.5854102;
	IntPoints_[1][2]=0.1381966;

	IntPoints_[2][0]=0.1381966;
	IntPoints_[2][1]=0.1381966;
	IntPoints_[2][2]=0.5854102;

	IntPoints_[3][0]=0.1381966;
	IntPoints_[3][1]=0.1381966;
	IntPoints_[3][2]=0.1381966;

	CorrectIntWeights();

	break;
      
      case GaussOrder3:
	// Thomas Hughes "The finite element method", p. 174

	NumIntPoints_ = 5;
	DegreeInteg_  = 4;
     
	if (!IntPoints_)
	  IntPoints_ = new Vector<Double>[NumIntPoints_];
     
	for(UInt i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].Resize(Dim_);
     
	IntWeights_.Resize(NumIntPoints_);

	IntPoints_[0][0]=0.25;
	IntPoints_[0][1]=0.25;
	IntPoints_[0][2]=0.25;

	IntPoints_[1][0]=0.5;
	IntPoints_[1][1]=0.1666667; 
	IntPoints_[1][2]=0.1666667; 

	IntPoints_[2][0]=0.1666667;
	IntPoints_[2][1]=0.5;
	IntPoints_[2][2]=0.1666667;

	IntPoints_[3][2]=0.5;
	IntPoints_[3][0]=0.1666667;
	IntPoints_[3][1]=0.1666667;

	IntPoints_[4][0]=0.1666667;
	IntPoints_[4][1]=0.1666667;
	IntPoints_[4][2]=0.1666667;

	IntWeights_[0]=-0.8; 
	IntWeights_[1]=0.45;
	IntWeights_[2]=0.45;
	IntWeights_[3]=0.45;
	IntWeights_[4]=0.45;

	CorrectIntWeights();

	break;

      case GaussOrder5:
	// A.H.Stroud "Approximate Calculation of Multiple Integrals"
	// Printice Hall 1971

	NumIntPoints_ = 15;
	DegreeInteg_  = 5;

	if (!IntPoints_)
	  IntPoints_ = new Vector<Double>[NumIntPoints_];

	for(UInt i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].Resize(Dim_);

	IntWeights_.Resize(NumIntPoints_);

	IntPoints_[0][0]=0.25;  IntPoints_[0][1]=0.25;  IntPoints_[0][2]=0.25;
	IntPoints_[1][0]=0.09197107805272303;  IntPoints_[1][1]=0.09197107805272303;  IntPoints_[1][2]=0.09197107805272303;
	IntPoints_[2][0]=0.72408676584183096;  IntPoints_[2][1]=0.09197107805272303;  IntPoints_[2][2]=0.09197107805272303;
	IntPoints_[3][0]=0.09197107805272303;  IntPoints_[3][1]=0.72408676584183096;  IntPoints_[3][2]=0.09197107805272303;
	IntPoints_[4][0]=0.09197107805272303;  IntPoints_[4][1]=0.09197107805272303;  IntPoints_[4][2]=0.72408676584183096;
	IntPoints_[5][0]=0.44364916731037080;  IntPoints_[5][1]=0.05635083268962915;  IntPoints_[5][2]=0.05635083268962915;
	IntPoints_[6][0]=0.05635083268962915;  IntPoints_[6][1]=0.44364916731037080;   IntPoints_[6][2]=0.05635083268962915;
	IntPoints_[7][0]=0.05635083268962915;  IntPoints_[7][1]=0.05635083268962915;  IntPoints_[7][2]=0.44364916731037080;
	IntPoints_[8][0]=0.05635083268962915;  IntPoints_[8][1]=0.44364916731037080;  IntPoints_[8][2]=0.44364916731037080;
	IntPoints_[9][0]=0.44364916731037080;  IntPoints_[9][1]=0.05635083268962915;  IntPoints_[9][2]=0.44364916731037080;
	IntPoints_[10][0]=0.44364916731037080; IntPoints_[10][1]=0.44364916731037080; IntPoints_[10][2]=0.05635083268962915;
	IntPoints_[11][0]=0.31979362782962989; IntPoints_[11][1]=0.31979362782962989; IntPoints_[11][2]=0.31979362782962989;
	IntPoints_[12][0]=0.04061911651111023; IntPoints_[12][1]=0.31979362782962989; IntPoints_[12][2]=0.31979362782962989;
	IntPoints_[13][0]=0.31979362782962989; IntPoints_[13][1]=0.04061911651111023; IntPoints_[13][2]=0.31979362782962989;
	IntPoints_[14][0]=0.31979362782962989; IntPoints_[14][1]=0.31979362782962989; IntPoints_[14][2]=0.04061911651111023;

 
	IntWeights_[0]=0.019753086419753086;
	IntWeights_[1]=0.011989513963169772;
	IntWeights_[2]=0.011989513963169772;
	IntWeights_[3]=0.011989513963169772;
	IntWeights_[4]=0.011989513963169772;
	IntWeights_[5]=0.008818342151675485;
	IntWeights_[6]=0.008818342151675485;
	IntWeights_[7]=0.008818342151675485;
	IntWeights_[8]=0.008818342151675485;
	IntWeights_[9]=0.008818342151675485;
	IntWeights_[10]=0.00881834215167548;
	IntWeights_[11]=0.011511367871045397;
	IntWeights_[12]=0.011511367871045397;
	IntWeights_[13]=0.011511367871045397;
	IntWeights_[14]=0.011511367871045397;
 
	CorrectIntWeights();
	
	break;

      default:
	Error("Integration type is not implemented",__FILE__,__LINE__);
      }  
  }



  // Division throug 6 is necessary to use the correct quadrature rule for tets
  // this is necessary, because in traditional IntWeights tables, they sum up to 1
  // and no volume is regarded! (see Hughes, p. 174)
  void TetraFE::CorrectIntWeights()
  { 
    ENTER_IFCN( "TetraFE::CorrectIntWeights" );

    for(UInt i=0; i<IntWeights_.GetSize(); i++)
      IntWeights_[i] /= 6.;
  }


void TetraFE::GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                        const StdVector<UInt> & volConnect,
                                        const Vector<Double> & surfIntPoint,
                                        Vector<Double> & volIntPoint)
{
  ENTER_IFCN( "TetraFE::GetLocalIntPoints4Surface" );
  
  // Try to find out, which vertices are in common with
  // the surface element. Then calculate the product of all four
  // and compare them
  //
  //
  // 4+\
  //  |\ \           zeta 	
  //  | \  \ 	      ^ eta 	
  //  |  \  +3	      |/	
  //  |   \ |	      0--> xi
  //  |    \ \
  //  |     \|     REFERENCE TETRAHEDRAL ELEMENT
  //  +------+
  //  1      2

  StdVector<UInt> commonIndex(3);
  UInt found = 0;
  UInt indexProduct = 0;
  std::string errMsg;
  
  volIntPoint.Resize(3);
  
  // loop over surface connect
  for (UInt iSurf=0; iSurf<3; iSurf++)
    // loop over volume connect
    for (UInt iVol=0; iVol<4; iVol++)
      if (surfConnect[iSurf] == volConnect[iVol])
	{
	  commonIndex[found++] = iVol+1;
	}

  indexProduct =  commonIndex[0] * commonIndex[1] * commonIndex[2];

  //std::cerr << "indexProduct = " << indexProduct << std::endl;
  switch(indexProduct)
    {
    case 8:
      std::cerr << "surface [1,2,4] is common" << std::endl;
      // Surface[1,2,4] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = 0.0;
      volIntPoint[2] = surfIntPoint[1];
      break;

    case 24:
      std::cerr << "surface [2,3,4] is common" << std::endl;
      // Surface[2,3,4] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = surfIntPoint[1];
      volIntPoint[2] = 1.0 - surfIntPoint[0] - surfIntPoint[1];
      break;

    case 12:
      std::cerr << "surface [1,3,4] is common" << std::endl;
      // Surface[1,3,4] is common
      volIntPoint[0] = 0.0;
      volIntPoint[1] = surfIntPoint[0];
      volIntPoint[2] = surfIntPoint[1];
      break;
      
    case 6:
      std::cerr << "surface [1,2,3] is common" << std::endl;
      // Surface[1,2,3] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = surfIntPoint[1];
      volIntPoint[2] = 0.0;
      break;
    default:
      errMsg = "TetraFE::GetLocalIntPoints4Surface: surface and volume element ";
      errMsg = "have not three nodes in common. Check your .mesh-file.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
}


} // end of namespace
