#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "wedgeFE.hh"

namespace CoupledField
{

WedgeFE::WedgeFE()
{
  ENTER_FCN( "WedgeFE::WedgeFE" );

  Dim_        = 3;
  NumEdges_   = 9;
  NumFaces_   = 5;
  NumCorners_ = 6;
  //    numChilds_ = 8;
  MidPoint_ = 0.0, 0.0, 0.0;
  
  
  std::string integtype;
  params->Get( "type", integtype, "integRules", "wedge" );
  
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
	NumIntPoints_ = 6;
	DegreeInteg_  = 2;

	IntWeights_.Resize(NumIntPoints_);
	for(Integer i=0; i<IntWeights_.GetSize(); i++)
	  {
	    IntWeights_[i] = 0.166666666666667;
	  }
	
	if (!IntPoints_)
	  IntPoints_ = new Vector<Double>[NumIntPoints_];

	for(Integer i=0; i<NumIntPoints_; i++)
	  IntPoints_[i].Resize(Dim_);

	IntPoints_[0][0] =  0.166666666666667;
	IntPoints_[0][1] =  0.166666666666667;
	IntPoints_[0][2] = -0.57735026919;

	IntPoints_[1][0] =  0.666666666666667;
	IntPoints_[1][1] =  0.166666666666667;
	IntPoints_[1][2] = -0.57735026919;
	
	IntPoints_[2][0] =  0.166666666666667;
	IntPoints_[2][1] =  0.666666666666667;
	IntPoints_[2][2] = -0.57735026919;	

	IntPoints_[3][0] =  0.166666666666667;
	IntPoints_[3][1] =  0.166666666666667;
	IntPoints_[3][2] =  0.57735026919;

	IntPoints_[4][0] =  0.666666666666667;
	IntPoints_[4][1] =  0.166666666666667;
	IntPoints_[4][2] =  0.57735026919;
	
	IntPoints_[5][0] =  0.166666666666667;
	IntPoints_[5][1] =  0.666666666666667;
	IntPoints_[5][2] =  0.57735026919;

	break;
      
      default:
	Error("Integration type is not implemented",__FILE__,__LINE__);
      }  
  }


void WedgeFE::GetLocalIntPoints4Surface(const StdVector<Integer> & surfConnect,
					const StdVector<Integer> & volConnect,
					const Vector<Double> & surfIntPoint,
					Vector<Double> & volIntPoint)
{
  ENTER_IFCN( "WedgeFE::GetLocalIntPoints4Surface" );
  
  // Try to find out, which vertices are in common with
  // the surface element. Then calculate the product of all four
  // and compare them

  
  // Check if surface element is triangle 
  // or quadrilateral
  if (surfConnect.GetSize() == 3 ||
      surfConnect.GetSize() == 6) 
    {
      // ---- Triangle Surface ---
      StdVector<Integer> commonIndex(3);
      Integer found = 0;
      Integer indexSum = 0;
      std::string errMsg;
      
      volIntPoint.Resize(3);
      
      // loop over surface connect
      for (Integer iSurf=0; iSurf<3; iSurf++)
	// loop over volume connect
	for (Integer iVol=0; iVol<6; iVol++)
	  if (surfConnect[iSurf] == volConnect[iVol])
	    {
	      commonIndex[found++] = iVol+1;
	    }
      
      
      indexSum =  commonIndex[0] + commonIndex[1] + commonIndex[2];
      
      switch(indexSum)
	{
	case 6:
	  // Surface[1,2,3] is common
	  volIntPoint[0] = surfIntPoint[0];
	  volIntPoint[1] = surfIntPoint[1];
	  volIntPoint[2] = -1.0;
	  break;
	  
	case 15:
	  // Surface[4,5,6] is common
	  volIntPoint[0] = surfIntPoint[0];
	  volIntPoint[1] = surfIntPoint[1];
	  volIntPoint[2] = 1.0;
	  break;

	default:
	  errMsg = "WedgeFE::GetLocalIntPoints4Surface: surface and volume element ";
	  errMsg = "have not three nodes in common. Check your .mesh-file.";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	}
      
    } else {
      // ---- Quadrilateral Surface ---
      StdVector<Integer> commonIndex(4);
      Integer found = 0;
      Integer indexSum = 0;
      std::string errMsg;
      
      volIntPoint.Resize(3);
      
      // loop over surface connect
      for (Integer iSurf=0; iSurf<4; iSurf++)
	// loop over volume connect
	for (Integer iVol=0; iVol<6; iVol++)
	  if (surfConnect[iSurf] == volConnect[iVol])
	    {
	      commonIndex[found++] = iVol+1;
	    }
      
      
      indexSum =  commonIndex[0] + commonIndex[1] 
                + commonIndex[2] + commonIndex[3];
      
      // NOTE: Since the line quad-element is defined in the range [-1;+1]
      // we have to calculate (1+surfCoord)/2 in order to get the right
      // position on the wedge element
      
      switch(indexSum)
	{
	case 16:
	  // Surface[2,3,5,6] is common
	  volIntPoint[0] = 0.5 - (surfIntPoint[0] / 2.0);
	  volIntPoint[1] = 0.5 + (surfIntPoint[0] / 2.0);
	  volIntPoint[2] = surfIntPoint[0];
	  break;
	  
	case 14:
	  // Surface[1,3,4,6] is common
	  volIntPoint[0] = 0.0;
	  volIntPoint[1] = 0.5 * (1 + surfIntPoint[0]);
	  volIntPoint[2] = 0.5 * (1 + surfIntPoint[1]);
	  break;

	case 12:
	  // Surface[1,2,4,5] is common
	  volIntPoint[0] = 0.5 * (1 + surfIntPoint[1]);
	  volIntPoint[1] = 0.0;
 	  volIntPoint[2] = 0.5 * (1 + surfIntPoint[0]);
	  break;

	default:
	  errMsg = "WedgeFE::GetLocalIntPoints4Surface: surface and volume element ";
	  errMsg = "have not four nodes in common. Check your .mesh-file.";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	} // switch
    } // if
}


} // end of namespace
