#include <iostream>
#include <fstream>

#include "quad1fe.hh"

namespace CoupledField
{


Quad1FE :: Quad1FE() : RectangleFE()
{
  ENTER_FCN( "Quad1FE::Quad1FE" );

  Init();
}
  
Quad1FE :: ~Quad1FE()
{
  ENTER_FCN( "Quad1FE::~Quad1FE" );
}

void Quad1FE :: Init()
{
  ENTER_IFCN( "Quad1FE::Init" );

  NumNodes_ = 4;
  SetIntPoints();
  SetCornerCoords();
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();  
}

void Quad1FE :: SetCornerCoords()
{
  ENTER_IFCN( "Quad1FE::SetCornerCoords" );

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] =  -1;
  LCornerCoords_[1][0] =  -1;
  LCornerCoords_[0][1] = 1;
  LCornerCoords_[1][1] =  -1;
  LCornerCoords_[0][2] = 1;
  LCornerCoords_[1][2] = 1;
  LCornerCoords_[0][3] = -1;
  LCornerCoords_[1][3] = 1;


}

void Quad1FE :: CalcShapeFnc(Vector<Double> & Shape, 
			     const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Quad1FE::CalcShapeFnc" );

  Shape.Resize(NumNodes_);
  
  for( Integer i=0; i<NumNodes_; i++)
    Shape[i] = 0.25 * (1 + LCornerCoords_[0][i] * LCoord[0])
                    * (1 + LCornerCoords_[1][i] * LCoord[1]);

}


void Quad1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Quad1FE::CalcLocalDerivShapeFnc" );

  LDeriv.Resize(NumNodes_,Dim_);

  for( Integer i=0; i<NumNodes_; i++)
    {
      LDeriv[i][0] = 0.25 * LCornerCoords_[0][i] 
                          * (1 + LCornerCoords_[1][i] * LCoord[1] );
      LDeriv[i][1] = 0.25 * (1 + LCornerCoords_[0][i] * LCoord[0] )
                          * LCornerCoords_[1][i];
    }
  
}


Double Quad1FE::CalcMeanStrain(Matrix<Double> &cornerCoords, 
			       Matrix<Double> &displacements)
{
  ENTER_IFCN( "Quad1FE::CalcDistortion" );

  Double factor;
  Double eps1, eps2, eps4, eps5, eps11, eps12, eps21, eps22, eps41, eps42, eps51, eps52;
  Double length1, length2, length11, length12, length21, length22;


  // calculate inital size of element
  length11 = abs(cornerCoords[0][0] - cornerCoords[0][1]);
  length12 = abs(cornerCoords[0][3] - cornerCoords[0][2]);
  length1 = (length11+length12) * 0.5;
  
  length21 = abs(cornerCoords[1][0] - cornerCoords[1][3]);
  length22 = abs(cornerCoords[1][1] - cornerCoords[1][2]);
  length2 = (length21+length22) * 0.5;

   // calculate absolute change of size
  eps11 = displacements[0][0] - displacements[0][1];
  eps12 = displacements[0][3] - displacements[0][2];
  eps1 = (eps11+eps12) * 0.5;

  eps21 = displacements[1][0] - displacements[1][3];
  eps22 = displacements[1][1] - displacements[1][2];
  eps2 = (eps21+eps22) * 0.5;
  
  eps41 = displacements[0][2] - displacements[0][1];
  eps42 = displacements[0][3] - displacements[0][0]; 
  eps4 = (eps41+eps42)*0.5;
    
  eps51 = displacements[1][1] - displacements[1][0];
  eps52 = displacements[1][3] - displacements[1][2];
  eps5= (eps51+eps52)*0.5;  

  factor =  0.2 * ((eps1*eps1/(length1*length1))
		 + (eps2*eps2/(length2*length2)) 
		 + (eps5*eps5/(length1*length2))
		 + (eps4*eps4/(length1*length1)));

  return factor;
}

void Quad1FE::GetLocalIntPoints4Surface(const StdVector<Integer> & surfConnect,
					const StdVector<Integer> & volConnect,
					const Vector<Double> & surfIntPoint,
					Vector<Double> & volIntPoint)
{
  ENTER_IFCN( "Quad1FE::GetLocalIntPoints4Surface" );
  
  // Try to find out, which vertices are in common with
  // the surface element. Then calculate the product of both
  // and compare them
  //
  //      eta
  //       ^
  // 4 +---|---+ 3    
  //   |   |   |      
  //   |   0---|-> xi     REFERENCE VOLUME ELEMENT
  //   |       |
  // 1 +-------+ 2



  StdVector<Integer> commonIndex(2);
  Integer found = 0;
  Integer indexProduct = 0;
  std::string errMsg;
  
  volIntPoint.Resize(2);
  
  // loop over surface connect
  for (Integer iSurf=0; iSurf<2; iSurf++)
    // loop over volume connect
    for (Integer iVol=0; iVol<4; iVol++)
      if (surfConnect[iSurf] == volConnect[iVol])
	{
	  commonIndex[found++] = iVol+1;
	}

  indexProduct= commonIndex[0] * commonIndex[1];
  switch(indexProduct)
    {
    case 2:
      // Edge[1,2] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = 1;
      break;

    case 12:
      // Edge[4,3] is common
      volIntPoint[0] = surfIntPoint[0];
      volIntPoint[1] = 1;
      break;

    case 4:
      // Edge[1,4] is common
      volIntPoint[0] = 1;
      volIntPoint[1] = surfIntPoint[0];
      break;

    case 6:
      // Edge[2,3] is common
      volIntPoint[0] = 1;
      volIntPoint[1] = surfIntPoint[0];
      break;

    default:
      errMsg = "Quad1FE::GetLocalIntPoints4Surface: surface and volume element ";
      errMsg = "have not two nodes in common. Check your .mesh-file.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
}

} // end of namespace


