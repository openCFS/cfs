// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "pyra2FE.hh"

namespace CoupledField
{

Pyra2FE::Pyra2FE():PyraFE()
{
  ENTER_FCN( "Pyra2FE::Pyra2FE" );
  
  Init();
}



Pyra2FE::~Pyra2FE()
{
  ENTER_FCN( "Pyra2FE::~Pyra2FE" );
}

void Pyra2FE::Init()
{
  ENTER_IFCN( "Pyra2FE::Init" );

  NumNodes_ = 13;
  NumEdges_ = 8;
  
  CommonInit();
}


void Pyra2FE::SetCornerCoords()
{
  ENTER_IFCN( "Pyra2FE::SetCornerCoords" );

  LCornerCoords_.Resize(Dim_,NumNodes_);

  // Original ones
//   LCornerCoords_[0][1] =  1;
//   LCornerCoords_[1][1] =  1;
//   LCornerCoords_[2][1] =  0;

//   LCornerCoords_[0][2] =  -1;
//   LCornerCoords_[1][2] =  1;
//   LCornerCoords_[2][2] =  0;

//   LCornerCoords_[0][3] =  -1;
//   LCornerCoords_[1][3] =  -1;
//   LCornerCoords_[2][3] =  0;

//   LCornerCoords_[0][0] =  1;
//   LCornerCoords_[1][0] =  -1;
//   LCornerCoords_[2][0] =  0;

//   LCornerCoords_[0][4] =  0;
//   LCornerCoords_[1][4] =  0;
//   LCornerCoords_[2][4] =  1;

// //Definition of the quadratic nodes

//   LCornerCoords_[0][6] =  0;
//   LCornerCoords_[1][6] =  1;
//   LCornerCoords_[2][6] =  0;

//   LCornerCoords_[0][7] = -1;
//   LCornerCoords_[1][7] =  0;
//   LCornerCoords_[2][7] =  0;

//   LCornerCoords_[0][8] =  0;
//   LCornerCoords_[1][8] =  -1;
//   LCornerCoords_[2][8] =  0;

//   LCornerCoords_[0][5] =  1;
//   LCornerCoords_[1][5] =  0;
//   LCornerCoords_[2][5] =  0;


//   LCornerCoords_[0][10] =  0.5;
//   LCornerCoords_[1][10] =  0.5;
//   LCornerCoords_[2][10] =  0.5;

//   LCornerCoords_[0][11] = -0.5;
//   LCornerCoords_[1][11] =  0.5;
//   LCornerCoords_[2][11] =  0.5;

//   LCornerCoords_[0][12] = -0.5;
//   LCornerCoords_[1][12] = -0.5;
//   LCornerCoords_[2][12] =  0.5;

//   LCornerCoords_[0][9] =  0.5;
//   LCornerCoords_[1][9] = -0.5;
//   LCornerCoords_[2][9] =  0.5;



  LCornerCoords_[0][1] =  1;
  LCornerCoords_[1][1] =  1;
  LCornerCoords_[2][1] =  0;

  LCornerCoords_[0][0] =  -1;
  LCornerCoords_[1][0] =  1;
  LCornerCoords_[2][0] =  0;

  LCornerCoords_[0][3] =  -1;
  LCornerCoords_[1][3] =  -1;
  LCornerCoords_[2][3] =  0;

  LCornerCoords_[0][2] =  1;
  LCornerCoords_[1][2] =  -1;
  LCornerCoords_[2][2] =  0;

  LCornerCoords_[0][4] =  0;
  LCornerCoords_[1][4] =  0;
  LCornerCoords_[2][4] =  1;

//Definition of the quadratic nodes

  LCornerCoords_[0][5] =  0;
  LCornerCoords_[1][5] =  1;
  LCornerCoords_[2][5] =  0;

  LCornerCoords_[0][8] = -1;
  LCornerCoords_[1][8] =  0;
  LCornerCoords_[2][8] =  0;

  LCornerCoords_[0][7] =  0;
  LCornerCoords_[1][7] =  -1;
  LCornerCoords_[2][7] =  0;

  LCornerCoords_[0][6] =  1;
  LCornerCoords_[1][6] =  0;
  LCornerCoords_[2][6] =  0;


  LCornerCoords_[0][10] =  0.5;
  LCornerCoords_[1][10] =  0.5;
  LCornerCoords_[2][10] =  0.5;

  LCornerCoords_[0][9] = -0.5;
  LCornerCoords_[1][9] =  0.5;
  LCornerCoords_[2][9] =  0.5;

  LCornerCoords_[0][12] = -0.5;
  LCornerCoords_[1][12] = -0.5;
  LCornerCoords_[2][12] =  0.5;

  LCornerCoords_[0][11] =  0.5;
  LCornerCoords_[1][11] = -0.5;
  LCornerCoords_[2][11] =  0.5;



}




// std::ostream& operator<< (std::ostream & outStr, Vector<Double> xOut)
// {
//   for (Integer i=0; i<xOut.size(); i++)
//     outStr <<  " " << xOut[i];
//   return outStr;
// }


void Pyra2FE :: CalcShapeFnc(Vector<Double> & Shape,
			     const Vector<Double> & LCoord,
                             const Elem*, UInt dof,
                             AnsatzFct::FctEntityType fctEntityType )
{
  ENTER_IFCN( "Pyra2FE::CalcShapeFnc" );

  Shape.Resize(NumNodes_);

	//"Pyramidal Elements"
	// F. Zgainski, J.L. Coulomb, Y. Marechal. IEEE Transactions on
        // Magnetics, Vol. 32, No. 3, May 1996


  //                                     zeta 
  //             4                         ^  
  //             +                         |
  //           // \                        0--> eta 
  //       12 +/ \ +11                    / 
  //        9/+ 10+ \                    xi 
  //        / /   \  \
  //     3 +-/--+-\ ---+ 2
  //      / /  7  \   /
  //   8 + /      \  + 6   REFERENCE VOLUME ELEMENT
  //    //        \ /
  // 0 +-----+-----+ 1
  //        5




  // Order is the same as in source paper, just index ist changed to follow 
  // our elements standard numbering
  if (LCoord[2]== 1)
    {
      Shape[1] = 0.25*(LCoord[0]+LCoord[1]-1)*((1+LCoord[0])*
					       (1+LCoord[1])-LCoord[2]);
      Shape[0] = 0.25*(LCoord[0]-LCoord[1]-1)*((1+LCoord[0])*
					       (1-LCoord[1])-LCoord[2]);
      Shape[3] = 0.25*(-LCoord[0]-LCoord[1]-1)*((1-LCoord[0])*
						(1-LCoord[1])-LCoord[2]);
      Shape[2] = 0.25*(-LCoord[0]+LCoord[1]-1)*((1-LCoord[0])*
						(1+LCoord[1])-LCoord[2]);
      Shape[4] = LCoord[2]*(2*LCoord[2]-1);
      Shape[5] = 0.0;
      Shape[8] = 0.0;
      Shape[7] = 0.0;
      Shape[6] = 0.0;
      Shape[10] = 0.0;
      Shape[9] = 0.0;
      Shape[12] = 0.0;
      Shape[11] = 0.0;
    }
  else
    {
      Shape[1] = 0.25*( LCoord[0]+LCoord[1]-1) * ((1+LCoord[0])*(1+LCoord[1])
		-LCoord[2]+(LCoord[0] * LCoord[1] * LCoord[2])/(1-LCoord[2]));
      Shape[0] = 0.25*( LCoord[0]-LCoord[1]-1) * ((1+LCoord[0])*(1-LCoord[1])
		-LCoord[2]-(LCoord[0] * LCoord[1] * LCoord[2])/(1-LCoord[2]));
      Shape[3] = 0.25*(-LCoord[0]-LCoord[1]-1) * ((1-LCoord[0])*(1-LCoord[1])
		-LCoord[2]+(LCoord[0] * LCoord[1] * LCoord[2])/(1-LCoord[2]));
      Shape[2] = 0.25*(-LCoord[0]+LCoord[1]-1) * ((1-LCoord[0])*(1+LCoord[1])
		-LCoord[2]-(LCoord[0] * LCoord[1] * LCoord[2])/(1-LCoord[2]));
      Shape[4] = LCoord[2] * (2 * LCoord[2]-1);
      Shape[5] = (1+LCoord[0]-LCoord[2]) * (1-LCoord[1]-LCoord[2])*
	         (1+LCoord[1]-LCoord[2])/(2*(1-LCoord[2]));
      Shape[8] = (1+LCoord[0]-LCoord[2])*(1-LCoord[0]-LCoord[2])*
                 (1-LCoord[1]-LCoord[2])/(2*(1-LCoord[2]));
      Shape[7] = (1-LCoord[0]-LCoord[2])*(1-LCoord[1]-LCoord[2])*
                 (1+LCoord[1]-LCoord[2])/(2*(1-LCoord[2]));
      Shape[6] = (1+LCoord[0]-LCoord[2])*(1-LCoord[0]-LCoord[2])*
	         (1+LCoord[1]-LCoord[2])/(2*(1-LCoord[2]));
      Shape[10] = LCoord[2]*(1+LCoord[0]-LCoord[2])*
	         (1+LCoord[1]-LCoord[2])/(1-LCoord[2]);
      Shape[9] = LCoord[2]*(1+LCoord[0]-LCoord[2])*
                 (1-LCoord[1]-LCoord[2])/(1-LCoord[2]);
      Shape[12] = LCoord[2]*(1-LCoord[0]-LCoord[2])*
                 (1-LCoord[1]-LCoord[2])/(1-LCoord[2]);
      Shape[11] = LCoord[2]*(1-LCoord[0]-LCoord[2])*
                 (1+LCoord[1]-LCoord[2])/(1-LCoord[2]);
    }

//   //Format for  Maple
//       Shape[1] = 0.25*( LCoord_0+LCoord_1-1) * 
//((1+LCoord_0)*(1+LCoord_1)-LCoord_2+
// 		  (LCoord_0 * LCoord_1 * LCoord_2)/(1-LCoord_2));
//       Shape[0] = 0.25*( LCoord_0-LCoord_1-1) * ((1+LCoord_0)*(1-LCoord_1)-
//           LCoord_2-(LCoord_0 * LCoord_1 * LCoord_2)/(1-LCoord_2));
//       Shape[3] = 0.25*(-LCoord_0-LCoord_1-1) * ((1-LCoord_0)*(1-LCoord_1)-
//           LCoord_2+(LCoord_0 * LCoord_1 * LCoord_2)/(1-LCoord_2));
//       Shape[2] = 0.25*(-LCoord_0+LCoord_1-1) * ((1-LCoord_0)*(1+LCoord_1)-
//           LCoord_2-(LCoord_0 * LCoord_1 * LCoord_2)/(1-LCoord_2));
//       Shape[4] = LCoord_2 * (2 * LCoord_2-1);
//       Shape[5] = (1+LCoord_0-LCoord_2) * (1-LCoord_1-LCoord_2)*
//                  (1+LCoord_1-LCoord_2)/(2*(1-LCoord_2));
//       Shape[8] = (1+LCoord_0-LCoord_2)*(1-LCoord_0-LCoord_2)*
//                  (1-LCoord_1-LCoord_2)/(2*(1-LCoord_2));
//       Shape[7] = (1-LCoord_0-LCoord_2)*(1-LCoord_1-LCoord_2)*
//                  (1+LCoord_1-LCoord_2)/(2*(1-LCoord_2));
//       Shape[6] = (1+LCoord_0-LCoord_2)*(1-LCoord_0-LCoord_2)*
//                  (1+LCoord_1-LCoord_2)/(2*(1-LCoord_2));
//       Shape[10] = LCoord_2*(1+LCoord_0-LCoord_2)*(1+LCoord_1-LCoord_2)/
//           (1-LCoord_2);
//       Shape[9] = LCoord_2*(1+LCoord_0-LCoord_2)*(1-LCoord_1-LCoord_2)/
//           (1-LCoord_2);
//       Shape[12] = LCoord_2*(1-LCoord_0-LCoord_2)*(1-LCoord_1-LCoord_2)/
//           (1-LCoord_2);
//       Shape[11] = LCoord_2*(1-LCoord_0-LCoord_2)*(1+LCoord_1-LCoord_2)/
//           (1-LCoord_2);

//   if (Shape[4] < 0)
//     Error("Local coordinates are not inside pyramidal element!",
//    __FILE__,__LINE__);

}


void Pyra2FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
				       const Vector<Double> & LCoord,
                                       const Elem*, UInt dof,
                                       AnsatzFct::FctEntityType fctEntityType )
{
  ENTER_IFCN( "Pyra2FE::CalcLocalDerivShapeFnc" );

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv.Init();

  // Derivatives for the quadratic case.      
  // Calculated symbolically with Maple.

  LDeriv[4][0] = 0;
  LDeriv[4][1] = 0;
  LDeriv[4][2] = 4*LCoord[2]-1;

  if (LCoord[2]==1)
    {
      LDeriv[1][0] = 0.25*( (1+LCoord[0])*(1+LCoord[1])-LCoord[2]+
			    (LCoord[0]+LCoord[1]-1)*(1+LCoord[1]));
      LDeriv[1][1] = 0.25*( (1+LCoord[0])*(1+LCoord[1])-LCoord[2]+
			    (LCoord[0]+LCoord[1]-1)*(1+LCoord[0]));
      LDeriv[1][2] = 0.25*( -LCoord[0]-LCoord[1]+1);
      LDeriv[0][0] = 0.25*( (1+LCoord[0])*(1-LCoord[1])-LCoord[2]+
			    (LCoord[0]-LCoord[1]-1)*(1-LCoord[1]));
      LDeriv[0][1] = 0.25*( -(1+LCoord[0])*(1-LCoord[1])+LCoord[2]+
			    (LCoord[0]-LCoord[1]-1)*(-1-LCoord[0]));
      LDeriv[0][2] = 0.25*( -LCoord[0]+LCoord[1]+1);
      LDeriv[3][0] = 0.25*( -(1-LCoord[0])*(1-LCoord[1])+LCoord[2]+
			    (-LCoord[0]-LCoord[1]-1)*(-1+LCoord[1]));
      LDeriv[3][1] = 0.25*( -(1-LCoord[0])*(1-LCoord[1])+LCoord[2]+
			    (-LCoord[0]-LCoord[1]-1)*(-1+LCoord[0]));
      LDeriv[3][2] = 0.25*( LCoord[0]+LCoord[1]+1);
      LDeriv[2][0] = 0.25*( -(1-LCoord[0])*(1+LCoord[1])+LCoord[2]+
			    (-LCoord[0]+LCoord[1]-1)*(-1-LCoord[1]));
      LDeriv[2][1] = 0.25*( (1-LCoord[0])*(1+LCoord[1])-LCoord[2]+
			    (-LCoord[0]+LCoord[1]-1)*(1-LCoord[0]));
      LDeriv[2][2] = 0.25*( LCoord[0]-LCoord[1]+1);
      
      LDeriv[6][0] = 0.0;
      LDeriv[6][1] = 0.0;
      LDeriv[6][2] = 0.0;
      
      LDeriv[7][0] = 0.0;
      LDeriv[7][1] = 0.0;
      LDeriv[7][2] = 0.0;
      
      LDeriv[8][0] = 0.0;
      LDeriv[8][1] = 0.0;
      LDeriv[8][2] = 0.0;
      
      LDeriv[5][0] = 0.0;
      LDeriv[5][1] = 0.0;
      LDeriv[5][2] = 0.0;
      
      LDeriv[10][0] = 0.0;
      LDeriv[10][1] = 0.0;
      LDeriv[10][2] = 0.0;
      
      LDeriv[11][0] = 0.0;
      LDeriv[11][1] = 0.0;
      LDeriv[11][2] = 0.0;
      
      LDeriv[12][0] = 0.0;
      LDeriv[12][1] = 0.0;
      LDeriv[12][2] = 0.0;
      
      LDeriv[9][0] = 0.0;
      LDeriv[9][1] = 0.0;
      LDeriv[9][2] = 0.0;
    }
  else
    {
      LDeriv[1][0] = 0.25*((1+LCoord[0])*(1+LCoord[1])-LCoord[2]+
			   LCoord[0]*LCoord[1]*LCoord[2]/(1-LCoord[2])+
			   (LCoord[0]+LCoord[1]-1)*(1+LCoord[1]+
				LCoord[1]*LCoord[2]/(1-LCoord[2])));
      LDeriv[1][1] = 0.25*((1+LCoord[0])*(1+LCoord[1])-LCoord[2]+
			   LCoord[0]*LCoord[1]*LCoord[2]/(1-LCoord[2])+
			   (LCoord[0]+LCoord[1]-1)*(1+LCoord[0]+LCoord[0]*
				LCoord[2]/(1-LCoord[2])));
      LDeriv[1][2] = 0.25*((LCoord[0]+LCoord[1]-1)*(-1+LCoord[0]*LCoord[1]/
			(1-LCoord[2])+LCoord[0]*LCoord[1]*LCoord[2]/
				((1-LCoord[2])*(1-LCoord[2]))));

      LDeriv[0][0] = 0.25*((1+LCoord[0])*(1-LCoord[1])-LCoord[2]-LCoord[0]*
			   LCoord[1]*LCoord[2]/(1-LCoord[2])+
			   (LCoord[0]-LCoord[1]-1)*(1-LCoord[1]-
			LCoord[1]*LCoord[2]/(1-LCoord[2])));
      LDeriv[0][1] = 0.25*(-(1+LCoord[0])*(1-LCoord[1])+LCoord[2]+LCoord[0]*
			   LCoord[1]*LCoord[2]/(1-LCoord[2])+
			   (LCoord[0]-LCoord[1]-1)*(-1-LCoord[0]-
				LCoord[0]*LCoord[2]/(1-LCoord[2])));
      LDeriv[0][2] = 0.25*((LCoord[0]-LCoord[1]-1)*(-1-LCoord[0]*LCoord[1]/
			(1-LCoord[2])-LCoord[0]*LCoord[1]*LCoord[2]/
			((1-LCoord[2])*(1-LCoord[2]))));
      
      LDeriv[3][0] = 0.25*(-(1-LCoord[0])*(1-LCoord[1])+LCoord[2]-LCoord[0]*
			   LCoord[1]*LCoord[2]/(1-LCoord[2])+(-LCoord[0]-
			LCoord[1]-1)*(-1+LCoord[1]+LCoord[1]*LCoord[2]/
			(1-LCoord[2])));
      LDeriv[3][1] = 0.25*(-(1-LCoord[0])*(1-LCoord[1])+LCoord[2]-LCoord[0]*
			   LCoord[1]*LCoord[2]/(1-LCoord[2])+
			   (-LCoord[0]-LCoord[1]-1)*(-1+LCoord[0]+LCoord[0]*
				LCoord[2]/(1-LCoord[2])));
      LDeriv[3][2] = 0.25*((-LCoord[0]-LCoord[1]-1)*(-1+LCoord[0]*LCoord[1]/
				(1-LCoord[2])+LCoord[0]*LCoord[1]*LCoord[2]/
				((1-LCoord[2])*(1-LCoord[2]))));
      
      LDeriv[2][0] = 0.25*(-(1-LCoord[0])*(1+LCoord[1])+LCoord[2]+LCoord[0]*
			   LCoord[1]*LCoord[2]/(1-LCoord[2])+(-LCoord[0]+
				LCoord[1]-1)*(-1-LCoord[1]-LCoord[1]*
				LCoord[2]/(1-LCoord[2])));
      LDeriv[2][1] = 0.25*((1-LCoord[0])*(1+LCoord[1])-LCoord[2]-LCoord[0]*
			   LCoord[1]*LCoord[2]/(1-LCoord[2])+
			   (-LCoord[0]+LCoord[1]-1)*(1-LCoord[0]-LCoord[0]*
				LCoord[2]/(1-LCoord[2])));
      LDeriv[2][2] = 0.25*((-LCoord[0]+LCoord[1]-1)*(-1-LCoord[0]*LCoord[1]/
				(1-LCoord[2])-LCoord[0]*LCoord[1]*LCoord[2]/
				((1-LCoord[2])*(1-LCoord[2]))));

      LDeriv[5][0] = .5*(1-LCoord[1]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             (1-LCoord[2]);
      LDeriv[5][1] = -.5*(1+LCoord[0]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             (1-LCoord[2])+.5*(1+LCoord[0]-LCoord[2])*
	             (1-LCoord[1]-LCoord[2])/(1-LCoord[2]);
      LDeriv[5][2] = -.5*(1-LCoord[1]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             (1-LCoord[2])-.5*(1+LCoord[0]-LCoord[2])*
	             (1+LCoord[1]-LCoord[2])/(1-LCoord[2])-.5*
	             (1+LCoord[0]-LCoord[2])*(1-LCoord[1]-LCoord[2])/
	             (1-LCoord[2])+.5*(1+LCoord[0]-LCoord[2])*
                     (1-LCoord[1]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             ((1-LCoord[2])*(1-LCoord[2]));
      
      LDeriv[8][0] = .5*(1-LCoord[0]-LCoord[2])*(1-LCoord[1]-LCoord[2])/
	             (1-LCoord[2])-.5*(1+LCoord[0]-LCoord[2])*
	             (1-LCoord[1]-LCoord[2])/(1-LCoord[2]);
      LDeriv[8][1] = -.5*(1+LCoord[0]-LCoord[2])*(1-LCoord[0]-LCoord[2])/
	             (1-LCoord[2]);
      LDeriv[8][2] = -.5*(1-LCoord[0]-LCoord[2])*(1-LCoord[1]-LCoord[2])/
	             (1-LCoord[2])-.5*(1+LCoord[0]-LCoord[2])*
	             (1-LCoord[1]-LCoord[2])/(1-LCoord[2])-.5*
	             (1+LCoord[0]-LCoord[2])*(1-LCoord[0]-LCoord[2])/
	             (1-LCoord[2])+.5*(1+LCoord[0]-LCoord[2])*
	             (1-LCoord[0]-LCoord[2])*(1-LCoord[1]-LCoord[2])/
	             ((1-LCoord[2])*(1-LCoord[2]));
      
      LDeriv[7][0] = -.5*(1-LCoord[1]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             (1-LCoord[2]);
      LDeriv[7][1] = -.5*(1-LCoord[0]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             (1-LCoord[2])+.5*(1-LCoord[0]-LCoord[2])*
	             (1-LCoord[1]-LCoord[2])/(1-LCoord[2]);
      LDeriv[7][2] = -.5*(1-LCoord[1]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             (1-LCoord[2])-.5*(1-LCoord[0]-LCoord[2])*
	             (1+LCoord[1]-LCoord[2])/(1-LCoord[2])-.5*(1-LCoord[0]-
			LCoord[2])*(1-LCoord[1]-LCoord[2])/(1-LCoord[2])+
	             .5*(1-LCoord[0]-LCoord[2])*(1-LCoord[1]-LCoord[2])*
	             (1+LCoord[1]-LCoord[2])/((1-LCoord[2])*(1-LCoord[2]));
      
      LDeriv[6][0] = .5*(1-LCoord[0]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             (1-LCoord[2])-.5*(1+LCoord[0]-LCoord[2])*
	             (1+LCoord[1]-LCoord[2])/(1-LCoord[2]);
      LDeriv[6][1] = .5*(1+LCoord[0]-LCoord[2])*(1-LCoord[0]-LCoord[2])/
	             (1-LCoord[2]);
      LDeriv[6][2] = -.5*(1-LCoord[0]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             (1-LCoord[2])-.5*(1+LCoord[0]-LCoord[2])*
	             (1+LCoord[1]-LCoord[2])/(1-LCoord[2])-.5*(1+LCoord[0]-
			LCoord[2])*(1-LCoord[0]-LCoord[2])/
	             (1-LCoord[2])+.5*(1+LCoord[0]-LCoord[2])*
	             (1-LCoord[0]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	             ((1-LCoord[2])*(1-LCoord[2]));
      
      LDeriv[10][0] = LCoord[2]*(1+LCoord[1]-LCoord[2])/(1-LCoord[2]);
      LDeriv[10][1] = LCoord[2]*(1+LCoord[0]-LCoord[2])/(1-LCoord[2]);
      LDeriv[10][2] = (1+LCoord[0]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	              (1-LCoord[2])-LCoord[2]*(1+LCoord[1]-LCoord[2])/
	              (1-LCoord[2])-LCoord[2]*(1+LCoord[0]-LCoord[2])/
	              (1-LCoord[2])+LCoord[2]*(1+LCoord[0]-LCoord[2])*
	              (1+LCoord[1]-LCoord[2])/((1-LCoord[2])*(1-LCoord[2]));
      
      LDeriv[9][0] = LCoord[2]*(1-LCoord[1]-LCoord[2])/(1-LCoord[2]);
      LDeriv[9][1] = -LCoord[2]*(1+LCoord[0]-LCoord[2])/(1-LCoord[2]);
      LDeriv[9][2] = (1+LCoord[0]-LCoord[2])*(1-LCoord[1]-LCoord[2])/
	             (1-LCoord[2])-LCoord[2]*(1-LCoord[1]-LCoord[2])/
	             (1-LCoord[2])-LCoord[2]*(1+LCoord[0]-LCoord[2])/
	             (1-LCoord[2])+LCoord[2]*(1+LCoord[0]-LCoord[2])*
	             (1-LCoord[1]-LCoord[2])/((1-LCoord[2])*(1-LCoord[2]));
      
      LDeriv[12][0] = -LCoord[2]*(1-LCoord[1]-LCoord[2])/(1-LCoord[2]);
      LDeriv[12][1] = -LCoord[2]*(1-LCoord[0]-LCoord[2])/(1-LCoord[2]);
      LDeriv[12][2] = (1-LCoord[0]-LCoord[2])*(1-LCoord[1]-LCoord[2])/
	              (1-LCoord[2])-LCoord[2]*(1-LCoord[1]-LCoord[2])/
	              (1-LCoord[2])-LCoord[2]*(1-LCoord[0]-LCoord[2])/
	              (1-LCoord[2])+LCoord[2]*(1-LCoord[0]-LCoord[2])*
	              (1-LCoord[1]-LCoord[2])/((1-LCoord[2])*(1-LCoord[2]));
      
      LDeriv[11][0] = -LCoord[2]*(1+LCoord[1]-LCoord[2])/(1-LCoord[2]);
      LDeriv[11][1] = LCoord[2]*(1-LCoord[0]-LCoord[2])/(1-LCoord[2]);
      LDeriv[11][2] = (1-LCoord[0]-LCoord[2])*(1+LCoord[1]-LCoord[2])/
	              (1-LCoord[2])-LCoord[2]*(1+LCoord[1]-LCoord[2])/
	              (1-LCoord[2])-LCoord[2]*(1-LCoord[0]-LCoord[2])/
	              (1-LCoord[2])+LCoord[2]*(1-LCoord[0]-LCoord[2])*
	              (1+LCoord[1]-LCoord[2])/((1-LCoord[2])*(1-LCoord[2]));
    }
}




} // end of namespace

