#include <iostream>
#include <fstream>

//#include <General/environment.hh>
#include "hexa1FE.hh"

namespace CoupledField
{

Hexa1FE::Hexa1FE():HexaFE()
{ 
  ENTER_FCN( "Hexa1FE::Hexa1FE" );

  Init();
}

Hexa1FE::~Hexa1FE()
{
  ENTER_FCN( "Hexa1FE::~Hexa1FE" );

}

void Hexa1FE::Init()
{
  ENTER_IFCN( "Hexa1FE::Init" );
  
  Dim_ = 3;
  NumNodes_ = 8;
  NumEdges_ = 12;
  //NumFaces_ = 6;

  // first set integration points and corner coords ...
  SetIntPoints();
  SetCornerCoords();


  // ... then calc shape function values at integration points
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();


}

void Hexa1FE::SetCornerCoords()
{
  ENTER_IFCN( "Hexa1FE::SetCornerCoords" );

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] =  -1;
  LCornerCoords_[1][0] =  -1;
  LCornerCoords_[2][0] =  -1;

  LCornerCoords_[0][1] =  1;
  LCornerCoords_[1][1] =  -1;
  LCornerCoords_[2][1] =  -1;

  LCornerCoords_[0][2] =  1;
  LCornerCoords_[1][2] =  1;
  LCornerCoords_[2][2] =  -1;

  LCornerCoords_[0][3] =  -1;
  LCornerCoords_[1][3] =  1;
  LCornerCoords_[2][3] =  -1;

  LCornerCoords_[0][4] =  -1;
  LCornerCoords_[1][4] =  -1;
  LCornerCoords_[2][4] =  1;

  LCornerCoords_[0][5] =  1;
  LCornerCoords_[1][5] =  -1;
  LCornerCoords_[2][5] =  1;

  LCornerCoords_[0][6] =  1;
  LCornerCoords_[1][6] =  1;
  LCornerCoords_[2][6] =  1;

  LCornerCoords_[0][7] =  -1;
  LCornerCoords_[1][7] =  1;
  LCornerCoords_[2][7] =  1;

}



void Hexa1FE :: CalcShapeFnc(Vector<Double> & Shape, 
			     const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Hexa1FE::CalcShapeFnc" );

  Shape.Resize(NumNodes_);

 
  for( Integer i=0; i<NumNodes_; i++)
    Shape[i] = 0.125 * (1 + LCornerCoords_[0][i] * LCoord[0])
                    * (1 + LCornerCoords_[1][i] * LCoord[1]) * (1 + LCornerCoords_[2][i] * LCoord[2]);

}

void Hexa1FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				     const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Hexa1FE::CalcLocalDerivShapeFnc" );

  LDeriv.Resize(NumNodes_,Dim_);

  for( Integer i=0; i<NumNodes_; i++)
    {
      LDeriv[i][0] = 0.125 * LCornerCoords_[0][i] 
	* (1 + LCornerCoords_[1][i] * LCoord[1] )* (1 + LCornerCoords_[2][i] * LCoord[2]);

      LDeriv[i][1] = 0.125 * (1 + LCornerCoords_[0][i] * LCoord[0] )
	* LCornerCoords_[1][i] * (1 + LCornerCoords_[2][i] * LCoord[2]);

      LDeriv[i][2] = 0.125 * (1 + LCornerCoords_[0][i] * LCoord[0])
	* (1 + LCornerCoords_[1][i] * LCoord[1]) * LCornerCoords_[2][i];
    }
  
}


void Hexa1FE::GetLocalIntPoints4Surface(const StdVector<Integer> & surfConnect,
					const StdVector<Integer> & volConnect,
					const Vector<Double> & surfIntPoint,
					Vector<Double> & volIntPoint)
{
  ENTER_IFCN( "Hexa1FE::GetLocalIntPoints4Surface" );
  
  // Try to find out, which vertices are in common with
  // the surface element. Then calculate the product of all four
  // and compare them
  //
  //                    zeta 
  //                     ^ eta 
  //    8 +-------+ 7    |/
  //     /|      /|      0--> xi 
  //    / |     / |
  // 5 +--+----+6 |   
  //   |  +-- -|- + 3    
  //   | / 4   | /    REFERENCE VOLUME ELEMENT
  //   |/      |/
  // 1 +-------+ 2



  StdVector<Integer> commonIndex(4);
  Integer found = 0;
  Integer indexProduct = 0;
  std::string errMsg;
  
  volIntPoint.Resize(3);
  
  // loop over surface connect
  for (Integer iSurf=0; iSurf<4; iSurf++)
    // loop over volume connect
    for (Integer iVol=0; iVol<8; iVol++)
      if (surfConnect[iSurf] == volConnect[iVol])
	{
	  commonIndex[found++] = iVol+1;
	}

  std::cerr << std::endl << std::endl;
  //std::cerr << "commonIndex = " << std::endl << commonIndex << std::endl << std::endl;
  indexProduct =  commonIndex[0] * commonIndex[1];
  indexProduct *= commonIndex[2] * commonIndex[3];

  //std::cerr << "indexProduct = " << indexProduct << std::endl;
  switch(indexProduct)
    {
    case 24:
      // Surface[1,2,3,4] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = surfIntPoint[1];
      volIntPoint[2] = -1;
      break;

    case 1680:
      // Surface[5,6,7,8] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = surfIntPoint[1];
      volIntPoint[2] = 1;
      break;

    case 252:
      // Surface[2,3,7,6] is common
      volIntPoint[0] = 1;
      volIntPoint[1] = surfIntPoint[0];
      volIntPoint[2] = surfIntPoint[1];
      break;
      
    case 160:
      // Surface[1,5,8,4] is common
      volIntPoint[0] = -1;
      volIntPoint[1] = surfIntPoint[0];
      volIntPoint[2] = surfIntPoint[1];
      break;
      
    case 672:
      // Surface[4,3,7,8] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = 1;
      volIntPoint[2] = surfIntPoint[1];
      break;
      
    case 60:
      // Surface[1,2,6,5] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = -1;
      volIntPoint[2] = surfIntPoint[-1];
      break;
      
    default:
      errMsg = "Quad1FE::GetLocalIntPoints4Surface: surface and volume element ";
      errMsg = "have not two nodes in common. Check your .mesh-file.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
}

} // end of namespace

