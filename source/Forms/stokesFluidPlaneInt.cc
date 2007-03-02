#include <iostream>
#include <fstream>

#include "stokesFluidPlaneInt.hh"

namespace CoupledField
{

  StokesFluidPlaneInt::StokesFluidPlaneInt(Double density,
                                           Double dynamicViscosity)
    : StokesFluidInt(density, dynamicViscosity)
  {
    ENTER_FCN( "StokesFluidPlaneInt::StokesFluidPlaneInt" );

    name_ = "StokesFluidPlaneInt";
  }


 
  StokesFluidPlaneInt::~StokesFluidPlaneInt()
  {
    ENTER_FCN( "StokesFluidPlaneInt::~StokesFluidPlaneInt" );
  }



  void StokesFluidPlaneInt::CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 )
  {
    ENTER_FCN( "StokesFluidPlaneInt::CalcElementMatrix" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts= ptelem->GetNumIntPoints();
  
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  

    //Double mu=dynamicViscosity_;

    UInt N;  // DOFs per Node

    // derivation of shape functions after global coordinates 
    Matrix<Double> partElemAMat, partElemATMat, locElemMat;
    Matrix<Double> A1, A1a, A1b;
    Matrix<Double> A2, A2a, A2b, A2n;
    Matrix<Double> vA2, vA2n;
    Matrix<Double>A3_1, A3_2;
    Matrix<Double> A4_1, A4_2;
    Matrix<Double> A5, vvA1_A5; 
    Matrix<Double> xiDxDy;

    Vector<Double>  xi, xiDx, xiDy;

    N = 4; // 4 DOFs per Node

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(numFncs*N); 
    elemMat.Init();
    locElemMat.Resize(numFncs*N); 
    locElemMat.Init();

    xiDx.Resize(numFncs);
    xiDy.Resize(numFncs);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());
        ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, 
                                      ptCoord_, jacDet, ent1.GetElem() );

        for (UInt i=0; i< numFncs; i++) 
          {
            xiDx[i] = xiDxDy[i][0];
            xiDy[i] = xiDxDy[i][1];
          }

        //A1 = xiDxDy * xiDxDyTransp;;

        A1a.DyadicMult(xiDx,xiDx);
        A1b.DyadicMult(xiDy,xiDy);
        A1 = A1a + A1b;
        A1 *= intWeights[actIntPt-1] * jacDet;

        A2a.DyadicMult(xiDx,xiDy);
        A2b.DyadicMult(xiDy,xiDx);
        A2 = A2a - A2b;
        A2 *= intWeights[actIntPt-1] * jacDet;

        A3_1.DyadicMult(xiDy,xi);
        A3_1 *= intWeights[actIntPt-1] * jacDet;

        A3_2.DyadicMult(xi,xiDy);
        A3_2 *= intWeights[actIntPt-1] * jacDet;

        A4_1.DyadicMult(xiDx,xi);
        A4_1 *= -1 * intWeights[actIntPt-1] * jacDet;

        A4_2.DyadicMult(xi,xiDx);
        A4_2 *= -1 * intWeights[actIntPt-1] * jacDet;

        A5.DyadicMult(xi,xi);
        A5 *= intWeights[actIntPt-1] * jacDet;
        A5 += A1*dynamicViscosity_*dynamicViscosity_;
        //****************************************************************
	//vvA1_A5 = A1*dynamicViscosity_*dynamicViscosity_ + A5;
        A2n     = -A2;
        vA2     =  A2 * dynamicViscosity_;
        vA2n    = -A2 * dynamicViscosity_;

        //****************************************************************
	locElemMat.AddSubMatrix(A1       ,  0,          0);
        locElemMat.AddSubMatrix(A1       ,  numFncs,    numFncs);
        locElemMat.AddSubMatrix(A1       ,  2*numFncs,  2*numFncs);
        //
        locElemMat.AddSubMatrix(A2       ,  0,          numFncs);
        locElemMat.AddSubMatrix(A2n      ,  numFncs,    0);
        //
        locElemMat.AddSubMatrix(A3_1     ,  0,          3*numFncs);
        locElemMat.AddSubMatrix(A3_2     ,  3*numFncs,  0);
        //
        locElemMat.AddSubMatrix(A4_1     ,  numFncs,    3*numFncs);
        locElemMat.AddSubMatrix(A4_2     ,  3*numFncs,  numFncs);
        //
        locElemMat.AddSubMatrix(vA2      ,  2*numFncs,  3*numFncs);
        locElemMat.AddSubMatrix(vA2n     ,  3*numFncs,  2*numFncs);
        //
        locElemMat.AddSubMatrix(A5       ,  3*numFncs,  3*numFncs);

//  __                                 __
// |                                    |
// |  xiDx   xiDy  0      0             |
// |  0      0     xiDx   dynVisc*xiDy  |
// |  0      0     xiDy  -dynVisc*xiDx  |
// |  xiDy  -xiDx  0      xi            |
// |__                                __|

//  __                                                                  __
// |                                                                      |
// |  xiDx²*xiDy²   0             0             xiDy*xi                   |
// |  0             xiDx²*xiDy²   0            -xiDx*xi                   |
// |  0             0             xiDx²*xiDy²   0                         |
// |  xiDy*xi      -xiDx*xi       0             mu²*xiDy²+mu²*xiDx²+xi²   |
// |__                                                                  __|



//        partElemAMat.Resize(numFncs*N,numFncs*N); partElemAMat.Init();


//        for (i=0; i<numFncs; i++)
//          {
//            for (j=0; j<numFncs; j++)
//              {
//	        xiDx2_Plus_xiDy2 = (xiDxDy[i][0]*xiDxDy[j][0]) + (xiDxDy[i][1]*xiDxDy[j][1]);
//	        partElemAMat[i][j]                     =  xiDx2_Plus_xiDy2;
//	        partElemAMat[i+1*numFncs][j+1*numFncs] =  xiDx2_Plus_xiDy2;
//	        partElemAMat[i+2*numFncs][j+2*numFncs] =  xiDx2_Plus_xiDy2;
//		      
//	        xi_Mult_xiDy = (xi[i]*xiDxDy[j][1]);
//	        partElemAMat[i+3*numFncs][j]           =  xi_Mult_xiDy;
//
//	        xiDy_Mult_xi = (xiDxDy[i][1]*xi[j]);
//	        partElemAMat[i][j+3*numFncs]           =  xiDy_Mult_xi;
//		      
//	        xi_Mult_xiDx = (xi[i]*xiDxDy[j][0]);
//	        partElemAMat[i+3*numFncs][j+numFncs]   = -xi_Mult_xiDx;
//		      
//	        xiDx_Mult_xi = (xiDxDy[i][0]*xi[j]);
//	        partElemAMat[i+numFncs][j+3*numFncs]   = -xiDx_Mult_xi;
//
//	        mu2_xiDy2_Plus_mu2_xiDx2_Plus_xi2 = (mu*mu*xiDxDy[i][1]*xiDxDy[j][1]) + (mu*mu*xiDxDy[i][0]*xiDxDy[j][0]) + (xi[i]*xi[j]);
//	        partElemAMat[i+3*numFncs][j+3*numFncs] = mu2_xiDy2_Plus_mu2_xiDx2_Plus_xi2;
		      
		           
//            partElemAMat[0][j]           =  xiDxDy[j][0];
//            partElemAMat[0][j+numFncs]   =  xiDxDy[j][1];

//            partElemAMat[1][j+2*numFncs] =  xiDxDy[j][0];
//            partElemAMat[1][j+3*numFncs] =  dynamicViscosity_ * xiDxDy[j][1];

//            partElemAMat[2][j+2*numFncs] =  xiDxDy[j][1];
//            partElemAMat[2][j+3*numFncs] = -dynamicViscosity_ * xiDxDy[j][0];

//            partElemAMat[3][j]           =  xiDxDy[j][1];
//            partElemAMat[3][j+numFncs]   = -xiDxDy[j][0];
//            partElemAMat[3][j+3*numFncs] =  xi[j];


//              }
//	  }

//        partElemAMat *= intWeights[actIntPt-1] * jacDet;
        //partElemAMat.Transpose(partElemATMat);

        // assemble element matrix
        //locElemMat += partElemATMat * partElemAMat;
//        locElemMat += partElemAMat;
      }
      
    ResortElementMatrix(elemMat, locElemMat, numFncs, N);
    //std::cout << "locElemMat:" << std::endl
    //          << locElemMat << std::endl;
    //std::cout << "elemMat:" << std::endl
    //          << elemMat << std::endl;
  }


} // end namespace CoupledField
