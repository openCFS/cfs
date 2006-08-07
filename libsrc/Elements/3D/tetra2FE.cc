#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "tetra2FE.hh"
//
namespace CoupledField
{

Tetra2FE::Tetra2FE():TetraFE()
{ 
  ENTER_FCN( "Tetra2FE::Tetra2FE" );

  Init();
}


Tetra2FE::~Tetra2FE()
{
  ENTER_FCN( "Tetra2FE::~Tetra2FE" );
  
}

void Tetra2FE::Init()
{
  ENTER_IFCN( "Tetra2FE::Init" );

  NumNodes_ = 10;
  NumEdges_ = 6;

  CommonInit();
}


void Tetra2FE::SetCornerCoords()
{
  ENTER_IFCN( "Tetra2FE::SetCornerCoords" );

  LCornerCoords_.Resize(Dim_,NumNodes_);

  LCornerCoords_[0][0] =  0;
  LCornerCoords_[1][0] =  0;
  LCornerCoords_[2][0] =  0;

  LCornerCoords_[0][1] =  1;
  LCornerCoords_[1][1] =  0;
  LCornerCoords_[2][1] =  0;

  LCornerCoords_[0][2] =  0;
  LCornerCoords_[1][2] =  1;
  LCornerCoords_[2][2] =  0;

  LCornerCoords_[0][3] =  0;
  LCornerCoords_[1][3] =  0;
  LCornerCoords_[2][3] =  1;

//Quadratic nodes

  LCornerCoords_[0][4] =  0.5;
  LCornerCoords_[1][4] =  0;
  LCornerCoords_[2][4] =  0;

  LCornerCoords_[0][5] =  0.5;
  LCornerCoords_[1][5] =  0.5;
  LCornerCoords_[2][5] =  0;

  LCornerCoords_[0][6] =  0;
  LCornerCoords_[1][6] =  0.5;
  LCornerCoords_[2][6] =  0;

  LCornerCoords_[0][7] =  0;
  LCornerCoords_[1][7] =  0;
  LCornerCoords_[2][7] =  0.5;

  LCornerCoords_[0][8] =  0.5;
  LCornerCoords_[1][8] =  0;
  LCornerCoords_[2][8] =  0.5;

  LCornerCoords_[0][9] =  0;
  LCornerCoords_[1][9] =  0.5;
  LCornerCoords_[2][9] =  0.5;

}

/// defines the connection between nodes with "their" edge
void Tetra2FE :: SetEdgeVertices()
{
  ENTER_IFCN( "SetEdgeVertices" );
  const UInt nrNodesPerEdge = 3;

  edgeVertices_.Resize(NumEdges_, nrNodesPerEdge);

  edgeVertices_[0][0] = 0;
  edgeVertices_[0][1] = 4;
  edgeVertices_[0][2] = 1;

  edgeVertices_[1][0] = 0;
  edgeVertices_[1][1] = 6;
  edgeVertices_[1][2] = 2;

  edgeVertices_[2][0] = 0;
  edgeVertices_[2][1] = 7;
  edgeVertices_[2][2] = 3;

  edgeVertices_[3][0] = 1;
  edgeVertices_[3][1] = 5;
  edgeVertices_[3][2] = 2;

  edgeVertices_[4][0] = 1;
  edgeVertices_[4][1] = 8;
  edgeVertices_[4][2] = 3;

  edgeVertices_[5][0] = 2;
  edgeVertices_[5][1] = 9;
  edgeVertices_[5][2] = 3;
}



// std::ostream& operator<< (std::ostream & outStr, Vector<Double> xOut)
// {
//   for (Integer i=0; i<xOut.size(); i++)
//     outStr <<  " " << xOut[i];
//   return outStr;
// }


void Tetra2FE :: CalcShapeFnc(Vector<Double> & Shape,
			     const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Tetra2FE::CalcShapeFnc" );


  Shape.Resize(NumNodes_);

//See Finite Element Procedures :Klaus Juergen Bathe, Prentice Hall,
//page 375 Sec. 5.3.
//Definition of the shape functions from 5 to 10
	
	Shape[4]=4*LCoord[0]*(1 - LCoord[0] - LCoord[1] - LCoord[2]);
	Shape[5]=4*LCoord[0]*LCoord[1];
	Shape[6]=4*LCoord[1]*(1 - LCoord[0] - LCoord[1] - LCoord[2]);
	Shape[7]=4*LCoord[2]*(1 - LCoord[0] - LCoord[1] - LCoord[2]);
	Shape[8]=4*LCoord[0]*LCoord[2];
	Shape[9]=4*LCoord[1]*LCoord[2];
        

//definition of the shape functions from 1 to 4
	
	Shape[0]= 1. - LCoord[0] - LCoord[1] - LCoord[2] - 0.5*Shape[4] -
	          0.5*Shape[6] - 0.5*Shape[7];
	Shape[1]= LCoord[0] - 0.5*Shape[4] - 0.5*Shape[5] - 0.5*Shape[8];
	Shape[2]= LCoord[1] - 0.5*Shape[5] - 0.5*Shape[6] - 0.5*Shape[9];
	Shape[3]= LCoord[2] - 0.5*Shape[8] - 0.5*Shape[7] - 0.5*Shape[9];
	
// 	if (Shape[0] < 0)
// 	  Error("Local coordinates are not inside tetrahedral element!",
//       __FILE__,__LINE__);

}


void Tetra2FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
				       const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Tetra2FE::CalcLocalDerivShapeFnc" );

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv.Init();

//Calculation of the local derivatives from 5 to 10

	LDeriv[4][0]=4. - 8*LCoord[0] - 4.*LCoord[1] - 4.*LCoord[2];
	LDeriv[4][1]=-4.*LCoord[0];
	LDeriv[4][2]=-4.*LCoord[0];
	
	LDeriv[5][0]=4.*LCoord[1];
	LDeriv[5][1]=4.*LCoord[0];
	LDeriv[5][2]=0.0;
	
	LDeriv[6][0]=-4.*LCoord[1];
	LDeriv[6][1]=4. - 8*LCoord[1] - 4.*LCoord[0] - 4.*LCoord[2];
	LDeriv[6][2]=-4.*LCoord[1];
	
	LDeriv[8][0]=4.*LCoord[2];
	LDeriv[8][1]=0.0;
	LDeriv[8][2]=4.*LCoord[0];
	
	LDeriv[9][0]=0.0;
	LDeriv[9][1]=4.*LCoord[2];
	LDeriv[9][2]=4.*LCoord[1];
	
	LDeriv[7][0]=-4.*LCoord[2];
	LDeriv[7][1]=-4.*LCoord[2];
	LDeriv[7][2]=4. - 8.*LCoord[2] - 4.*LCoord[0] - 4.*LCoord[1];

	//Calculation of the local derivatives from 1 to 4 

	LDeriv[0][0]= -1.-0.5*LDeriv[4][0]-0.5*LDeriv[6][0]-0.5*LDeriv[7][0];
	LDeriv[0][1]= -1.-0.5*LDeriv[4][1]-0.5*LDeriv[6][1]-0.5*LDeriv[7][1];
	LDeriv[0][2]= -1.-0.5*LDeriv[4][2]-0.5*LDeriv[6][2]-0.5*LDeriv[7][2];
	
	LDeriv[1][0]= 1.-0.5*LDeriv[4][0]-0.5*LDeriv[5][0]-0.5*LDeriv[8][0]; 
	LDeriv[1][1]= -0.5*LDeriv[4][1]-0.5*LDeriv[5][1]-0.5*LDeriv[8][1];
	LDeriv[1][2]= -0.5*LDeriv[4][2]-0.5*LDeriv[5][2]-0.5*LDeriv[8][2];
	
	LDeriv[2][0]= -0.5*LDeriv[5][0]-0.5*LDeriv[6][0]-0.5*LDeriv[9][0];
	LDeriv[2][1]= 1.-0.5*LDeriv[5][1]-0.5*LDeriv[6][1]-0.5*LDeriv[9][1];
	LDeriv[2][2]= -0.5*LDeriv[5][2]-0.5*LDeriv[6][2]-0.5*LDeriv[9][2];
	
	LDeriv[3][0]= -0.5*LDeriv[9][0]-0.5*LDeriv[7][0]-0.5*LDeriv[8][0];
	LDeriv[3][1]= -0.5*LDeriv[9][1]-0.5*LDeriv[7][1]-0.5*LDeriv[8][1];
	LDeriv[3][2]= 1.-0.5*LDeriv[9][2]-0.5*LDeriv[7][2]-0.5*LDeriv[8][2]; 
}

// see Kaltenbacher: "Numerical Sim. of Mechatronic Sensors and Actuators"
//  p. 25 calculates the edge shape function of a tetrahedral of first order.
void Tetra2FE :: CalcEdgeShapeFnc(Matrix<Double> & edgeShape,
				  const Vector<Double> & LCoord,
				  const Matrix<Double> & cornernodes)
{
  ENTER_IFCN( "Tetra2FE::CalcEdgeShapeFnc" );

Error("Tetra2FE::CalcEdgeShapeFnc not yet implemented!!",__FILE__,__LINE__);

//   edgeShape.Resize(NumEdges_, Dim_);
//   edgeShape.Init();

//   // nodal shape functions of a tet
//   Vector<Double> nodeShape;
//   CalcShapeFnc(nodeShape, LCoord);


//   // local derivates of nodal tet, dimension: nrNodes x Dim_
//   Matrix<Double> xDxi;
//   //  CalcLocalDerivShapeFnc(xDxi, localcoord);
//   GetGlobDerivShFnc(xDxi, LCoord, cornernodes);

//   for (Integer actEdge=0; actEdge<NumEdges_; actEdge++)
//     {
//       Integer node1 = edgeVertices_[actEdge][0];
//       Integer node2 = edgeVertices_[actEdge][1];

//       for (Integer actDim=0; actDim<Dim_; actDim++)
// 	edgeShape[actEdge][actDim] =
// 	  nodeShape[node1] * xDxi[node2][actDim] -
// 	  nodeShape[node2] * xDxi[node1][actDim];
//     }
}



// calculated the Nedelec shape function in an arbitrary point
void Tetra2FE::GetEdgeGlobalDerivShapeFnc(Vector<Matrix<Double>*> &shapeDeriv,
					   const Vector<Double> & lCoord,
					   const Matrix<Double> &cornerCoords)
{
  ENTER_IFCN( "Tetra2FE::GetEdgeGlobalDerivShapeFnc" );

Error("Tetra2FE::GetEdgeGlobalDerivShapeFnc not yet implemented!!",
      __FILE__,__LINE__);

//   shapeDeriv.Resize(NumEdges_);


//   // local derivates of nodal tet, dimension: nrNodes x Dim_
//   Matrix<Double> xDxi;
//   GetGlobDerivShFnc(xDxi, lCoord, cornerCoords);


//   for (Integer actEdge=0; actEdge<NumEdges_; actEdge++)
//     {
//       shapeDeriv[actEdge]->Resize(Dim_,Dim_);

//       Integer node1 = edgeVertices_[actEdge][0];
//       Integer node2 = edgeVertices_[actEdge][1];

//       for (Integer dim1=0; dim1<Dim_; dim1++)
// 	for (Integer dim2=0; dim2<Dim_; dim2++)
// 	  (*shapeDeriv[actEdge]) [dim2][dim1] =
// 	    xDxi[node1][dim1] * xDxi[node2][dim2] -
// 	    xDxi[node1][dim2] * xDxi[node2][dim1];
//     }
}


} // end of namespace

