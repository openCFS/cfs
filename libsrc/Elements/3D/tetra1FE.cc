#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "tetra1FE.hh"

namespace CoupledField
{

Tetra1FE::Tetra1FE():TetraFE()
{ 
#ifdef TRACE
  (*trace) << "entering Tetra1FE::Tetra1FE" << std::endl;
#endif
   Init();
}


Tetra1FE::~Tetra1FE()
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::~Tetra1FE" << std::endl;
#endif
 ; 
}

void Tetra1FE::Init()
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::Init" << std::endl;
#endif
  
  Dim_ = 3;
  NumNodes_ = 4;
  NumEdges_ = 6;

  // first set integration points and corner coords ...
  SetIntPoints();
  SetCornerCoords();

  // ... then calc shape function values at integration points
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();
  SetEdgeVertices();
}


void Tetra1FE::SetCornerCoords()
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::SetCornerCoords" << std::endl;
#endif

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
//   LCornerCoords_[0][0] =  0;
//   LCornerCoords_[1][0] =  0;
//   LCornerCoords_[2][0] =  0;

//   LCornerCoords_[0][1] =  1;
//   LCornerCoords_[1][1] =  0;
//   LCornerCoords_[2][1] =  0;

//   LCornerCoords_[0][2] =  0;
//   LCornerCoords_[1][2] =  1;
//   LCornerCoords_[2][2] =  0;

//   LCornerCoords_[0][3] =  0;
//   LCornerCoords_[1][3] =  0;
//   LCornerCoords_[2][3] =  1;

  LCornerCoords_[0][0] =  1;
  LCornerCoords_[1][0] =  0;
  LCornerCoords_[2][0] =  0;

  LCornerCoords_[0][1] =  0;
  LCornerCoords_[1][1] =  1;
  LCornerCoords_[2][1] =  0;

  LCornerCoords_[0][2] =  0;
  LCornerCoords_[1][2] =  0;
  LCornerCoords_[2][2] =  1;

  LCornerCoords_[0][3] =  0;
  LCornerCoords_[1][3] =  0;
  LCornerCoords_[2][3] =  0;
}

/// defines the connection between nodes with "their" edge 
void Tetra1FE :: SetEdgeVertices()
{
  const Integer nrNodesPerEdge = 2;
  
  edgeVertices_.Resize(NumEdges_, nrNodesPerEdge);

//   edgeVertices_[0][0] = 0;
//   edgeVertices_[0][1] = 1;

//   edgeVertices_[1][0] = 0;
//   edgeVertices_[1][1] = 2;

//   edgeVertices_[2][0] = 0;
//   edgeVertices_[2][1] = 3;

//   edgeVertices_[3][0] = 1;
//   edgeVertices_[3][1] = 2;

//   edgeVertices_[4][0] = 1;
//   edgeVertices_[4][1] = 3;

//   edgeVertices_[5][0] = 2;
//   edgeVertices_[5][1] = 3;

  edgeVertices_[0][0] = 3;
  edgeVertices_[0][1] = 0;

  edgeVertices_[1][0] = 3;
  edgeVertices_[1][1] = 1;

  edgeVertices_[2][0] = 3;
  edgeVertices_[2][1] = 2;

  edgeVertices_[3][0] = 0;
  edgeVertices_[3][1] = 1;

  edgeVertices_[4][0] = 0;
  edgeVertices_[4][1] = 2;

  edgeVertices_[5][0] = 1;
  edgeVertices_[5][1] = 2;
}



// std::ostream& operator<< (std::ostream & outStr, std::vector<Double> xOut)
// {
//   for (Integer i=0; i<xOut.size(); i++)
//     outStr <<  " " << xOut[i];
//   return outStr;
// }


void Tetra1FE :: CalcShapeFnc(std::vector<Double> & Shape, 
			     const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::CalcShapeFnc" << std::endl;
#endif

  Shape.resize(NumNodes_);

  // see Hughes p. 170
  
  Shape[0] = 1.0 - LCoord[0] - LCoord[1] - LCoord[2];

  if (Shape[0] < 0)
    Error("Local coordinates are not inside tetrahedral element!",__FILE__,__LINE__);

  for( Integer i=1; i<NumNodes_; i++)
    Shape[i] = LCoord[i-1];

}


void Tetra1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv.Init();
  
  for( Integer i=0; i<Dim_; i++)
    LDeriv[0][i] = -1.0;

  for( Integer i=1; i < NumNodes_; i++)
    LDeriv[i][i-1] = 1.0;
}
  

// see Kaltenbacher: "Numerical Sim. of Mechatronic Sensors and Actuators" p. 25
// calculates the edge shape function of a tetrahedral of first order.
void Tetra1FE :: CalcEdgeShapeFnc(Matrix<Double> & edgeShape, 
			     const std::vector<Double> & LCoord, const Matrix<Double> & cornernodes)
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::CalcShapeFnc" << std::endl;
#endif

  edgeShape.Resize(NumEdges_, Dim_);
  edgeShape.Init();

  // nodal shape functions of a tet
  std::vector<Double> nodeShape;
  CalcShapeFnc(nodeShape, LCoord);


  // local derivates of nodal tet, dimension: nrNodes x Dim_
  Matrix<Double> xDxi;  
  //  CalcLocalDerivShapeFnc(xDxi, localcoord);
  GetGlobDerivShFnc(xDxi, LCoord, cornernodes);    

  for (Integer actEdge=0; actEdge<NumEdges_; actEdge++)
    {
      Integer node1 = edgeVertices_[actEdge][0];
      Integer node2 = edgeVertices_[actEdge][1];
      
      for (Integer actDim=0; actDim<Dim_; actDim++)
	edgeShape[actEdge][actDim] = 
	  nodeShape[node1] * xDxi[node2][actDim] - 
	  nodeShape[node2] * xDxi[node1][actDim];
    }  
}



// calculated the Nedelec shape function in an arbitrary point
void Tetra1FE :: GetEdgeGlobalDerivShapeFnc(std::vector< Matrix<Double>* > & shapeDeriv, 
					    const std::vector<Double> & lCoord,
					    const Matrix<Double> & cornerCoords)
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::GetEdgeGlobalDerivShapeFnc" << std::endl;
#endif

  shapeDeriv.resize(NumEdges_);


  // local derivates of nodal tet, dimension: nrNodes x Dim_
  Matrix<Double> xDxi;  
  GetGlobDerivShFnc(xDxi, lCoord, cornerCoords);
  
  
  for (Integer actEdge=0; actEdge<NumEdges_; actEdge++)
    {
      shapeDeriv[actEdge]->Resize(Dim_,Dim_);

      Integer node1 = edgeVertices_[actEdge][0];
      Integer node2 = edgeVertices_[actEdge][1];

      for (Integer dim1=0; dim1<Dim_; dim1++)
	for (Integer dim2=0; dim2<Dim_; dim2++)
	  (*shapeDeriv[actEdge]) [dim2][dim1] = 
	    //	  (*shapeDeriv[actEdge]) [dim1][dim2] = 
	    xDxi[node1][dim1] * xDxi[node2][dim2] -
	    xDxi[node1][dim2] * xDxi[node2][dim1];
    }
}


} // end of namespace

