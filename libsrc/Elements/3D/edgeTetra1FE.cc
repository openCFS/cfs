#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "edgeTetra1FE.hh"

namespace CoupledField
{

EdgeTetra1FE::EdgeTetra1FE():Tetra1FE()
{ 
#ifdef TRACE
  (*trace) << "entering EdgeTetra1FE::EdgeTetra1FE" << std::endl;
#endif
   Init();
}


EdgeTetra1FE::~EdgeTetra1FE()
{
#ifdef TRACE
  (*trace) << "entering EdgeTetra1FE::~EdgeTetra1FE" << std::endl;
#endif
 ; 
}

void EdgeTetra1FE::Init()
{
#ifdef TRACE
  (*trace) << "entering EdgeTetra1FE::Init" << std::endl;
#endif
  
  Dim_ = 3;
  NumNodes_ = 4;
  NumEdges_ = 6;
  
  // first set integration points and corner coords ...
  SetIntPoints();
  SetCornerCoords();

  // then calc shape function values at integration points
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();
  SetEdgeVertices();
}


/// defines the connection between nodes with "their" edge 
void EdgeTetra1FE :: SetEdgeVertices()
{
  const Integer nrNodesPerEdge = 2;
  
  edgeVertices_.Resize(NumEdges_, nrNodesPerEdge);

 //  edgeVertices_[0][0] = 0;
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



// see Kaltenbacher: "Numerical Sim. of Mechatronic Sensors and Actuators" p. 25
// calculates the edge shape function of a tetrahedral of first order.
void EdgeTetra1FE :: CalcEdgeShapeFnc(Matrix<Double> & edgeShape, 
			     const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering EdgeTetra1FE::CalcShapeFnc" << std::endl;
#endif

  edgeShape.Resize(NumEdges_, Dim_);


  // nodal shape functions of a tet
  std::vector<Double> nodeShape;
  std::vector<Double> localcoord;
  CalcShapeFnc(nodeShape, LCoord);


  // local derivates of nodal tet, dimension: nrNodes x Dim_
  Matrix<Double> xDxi;  
  GetGlobDerivShFnc(xDxi, localcoord, LCoord);
  
  

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


void EdgeTetra1FE :: CalcEdgeGlobalDerivShapeFnc(std::vector< Matrix<Double> > & shapeDeriv, 
				       const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering EdgeTetra1FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  shapeDeriv.resize(NumEdges_);


  // local derivates of nodal tet, dimension: nrNodes x Dim_
  Matrix<Double> xDxi;  
  CalcLocalDerivShapeFnc(xDxi, LCoord);
  
  
  for (Integer actEdge=0; actEdge<NumEdges_; actEdge++)
    {
      shapeDeriv[actEdge].Resize(Dim_,Dim_);

      Integer node1 = edgeVertices_[actEdge][0];
      Integer node2 = edgeVertices_[actEdge][1];

      for (Integer dim1=0; dim1<Dim_; dim1++)
	for (Integer dim2=0; dim2<Dim_; dim2++)
	  shapeDeriv[actEdge][dim1][dim2] = 
	    xDxi[node1][dim2] * xDxi[node2][dim1] -
	    xDxi[node1][dim2] * xDxi[node2][dim1];
    }
}


} // end of namespace

