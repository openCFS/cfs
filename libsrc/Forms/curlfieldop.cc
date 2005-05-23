#include "Forms/curlfieldop.hh"

#include <string>

#include "Elements/basefe.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"

#include <PDE/StdPDE.hh>

namespace CoupledField
{

CurlEdgeOp::CurlEdgeOp(Grid * ptGrid, 
		       StdPDE * ptPDE,
		       NodeEQN * ptEQN,
		       NodeStoreSol<Double> & aSol,
		       BaseSystem * algsys) 
  : BaseOperator(ptGrid, ptPDE, ptEQN), algsys_(algsys)
{
  ENTER_FCN( "CurlEdgeOp::CurlEdgeOp" );

  sol_ = &aSol;
}


CurlEdgeOp::~CurlEdgeOp()
{
  ENTER_FCN( "CurlEdgeOp::~CurlEdgeOp" );
}


void CurlEdgeOp::CalcElemCurlEdge(Vector<Double> & curlField, 
				  const Elem * ptElement,
				  const Vector<Double> & LCoord)
{
  ENTER_FCN( "CurlEdgeOp::CalcElemCurlEdgeOp" );
  
  Error("CurlEdgeOp::CalcElemCurlEdgeOp: Not working due to EQN-class!",
	__FILE__, __LINE__);
 //  ShortInt dim = ptElement->ptElem->GetDim();
  
//   curlField.Resize(dim);
//   for (Integer i=0; i<dim; i++)
//     curlField[i] = 0;

//   Integer nrEdges = ptElement->ptElem->GetNumEdges();
//   Integer nrNodes = ptElement->ptElem->GetNumNodes();

//   Matrix<Double> CornerCoords; 
//   ptGrid_->GetCoordNodesElemMat(ptElement->connect, CornerCoords, level_);


//   Matrix<Double> curlOnEdges;


//   StdVector<Matrix<Double>* > deriv;
//   deriv.Resize(nrEdges);
//   for (Integer actEdge=0; actEdge < nrEdges; actEdge++)
//     deriv[actEdge] = new Matrix<Double>;
  

//   ptElement->ptElem->GetEdgeGlobalDerivShapeFnc (deriv, LCoord, CornerCoords);

//   curlOnEdges.Resize(dim, nrEdges);
  
//   for (Integer actEdge=0; actEdge < nrEdges; actEdge++)
//     for (Integer actDim=0; actDim < dim; actDim++)
//       curlOnEdges[actDim][actEdge] = 
// 	(*deriv[actEdge])[(actDim+2)%dim][(actDim+1)%dim] -
// 	(*deriv[actEdge])[(actDim+1)%dim][(actDim+2)%dim];


//   StdVector<Double> sol(nrEdges);
//   // global edge index
//   StdVector<Integer> epos(nrEdges);
//   StdVector<Integer> esign(nrEdges);


//   StdVector<Integer> pos(nrNodes);
  
//   for (Integer i=0; i < nrNodes; i++)
//      pos[i] = (*ptMesh2PDENode_)[ptElement->connect[i]-1];
  
//   ptElement->ptElem->GetGlobalEdgeIndices(epos, &pos[0], algsys_);

// #ifdef DEBUG
//   (*debug) << "CurlOP pos \n" << pos << std::endl
// 	   << "epos \n " << epos << std::endl;
  
// #endif

//   for (Integer j=0; j<nrEdges; j++)
//     {
//       esign[j] = epos[j]/abs(epos[j]);
//       epos[j]  = abs(epos[j]);
//       sol[j] = (*sol_)[epos[j]-1] * esign[j];
//     }
  
  
//   // loop over edge curls
//   for( Integer i=0; i<dim; i++ )
//     {
//       curlField[i]=0;
      
//       for( Integer j=0; j < nrEdges; j++ )
// 	curlField[i] += curlOnEdges[i][j] * sol[j];
//     }
  
}




void CurlEdgeOp::CalcElemMagVec(Vector<Double> & magVecPot, 
				const Elem * ptElement,
				const Vector<Double> & lCoord)
{
  ENTER_FCN( "CurlEdgeOp::CalcElemMagVec" );
  
  Error( "CurlEdgeOp::CalcElemMagVec: Not working due to EQN-class",
	 __FILE__, __LINE__);
  // Integer nrEdges = ptElement->ptElem->GetNumEdges();
//   Integer nrNodes = ptElement->ptElem->GetNumNodes();
//   BaseFE * ptElem = ptElement->ptElem;
//   ShortInt dim = ptElem->GetDim();

//   Matrix<Double> cornerCoords; 
//   Matrix<Double> shape;

//   StdVector<Double> sol(nrEdges);
//   // global edge index
//   StdVector<Integer> epos(nrEdges);
//   StdVector<Integer> esign(nrEdges);
//   Vector<Integer> pos(nrNodes);



//   magVecPot.Resize(dim);

//   for (Integer i=0; i<dim; i++)
//     magVecPot[i] = 0;

//   ptGrid_->GetCoordNodesElemMat(ptElement->connect, cornerCoords, level_);
  
//   ptElem->CalcEdgeShapeFnc(shape, lCoord, cornerCoords);

  
//   for (Integer i=0; i < nrNodes; i++)
//      pos[i] = (*ptMesh2PDENode_)[ptElement->connect[i]-1];
  
//   ptElem->GetGlobalEdgeIndices(epos, &pos[0], algsys_);


//   for (Integer j=0; j<nrEdges; j++)
//     {
//       esign[j] = epos[j]/abs(epos[j]);
//       epos[j]  = abs(epos[j]);
//       sol[j] = (*sol_)[epos[j]-1] * esign[j];
//     }
  

  
  
//   // loop over edge curls
//   // magVecPot = sol * shape;
//   for( Integer j=0; j<dim; j++ )  
//     {
//       magVecPot[j]=0;
//       for( Integer i=0; i < nrEdges; i++ )
// 	magVecPot[j] += shape[i][j] * sol[i];
//     }
}


//=========================== CurlNodeOperator-Class=======================================

CurlNodeOp::CurlNodeOp(Grid * ptGrid, 
		       StdPDE * ptPDE,
		       NodeEQN * ptEQN,
		       NodeStoreSol<Double> & aSol) 
  : BaseOperator(ptGrid, ptPDE, ptEQN, FALSE)
{
  ENTER_FCN( "CurlNodeOp::CurlNodeOp" );

  sol_ = &aSol;
}

CurlNodeOp::~CurlNodeOp()
{
  ENTER_FCN( "CurlNodeOp::~CurlNodeOp" );

}

void CurlNodeOp::CalcElemCurlNode(Vector<Double> & B, 
				    const Elem * ptElement,
				    const Vector<Double> & LCoord)
{
  ENTER_FCN( "CurlNodeOp::CalcElemCurlNode" );
  
  ShortInt dim;
  Double solEntry;
  dim = ptElement->ptElem->GetDim();
  if (dim ==2)
    {
      B.Resize(dim);
      B.Init();

      Integer nShFnc = 0;
      nShFnc = ptElement->ptElem->GetNumNodes();
      
      Matrix<Double> CornerCoords; 
      ptPDE_->GetElemCoords(ptElement->connect, CornerCoords);
      
      Matrix<Double> GlobalGradient;
      
      ptElement->ptElem->GetGlobDerivShFnc(GlobalGradient, LCoord, CornerCoords);
      
      if (isaxi_)
	{
	  Vector<Double> ShpFncAtIp;
	  Vector<Double> CoordAtIP;
	  ptElement->ptElem->GetShFnc(ShpFncAtIp,LCoord);
	  CoordAtIP = CornerCoords * ShpFncAtIp;
	  for (Integer i=0; i<nShFnc; i++)
	    GlobalGradient[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
	}
      
      // loop over shape functions
      for( Integer i=0; i<dim; i++ )
	for( Integer j=0; j<nShFnc; j++ )
	  {
	    sol_->Get(ptElement->connect[j]-1,0,solEntry);
	    //const Double solEntry = (*sol_)(ptElement->connect[j],1);
	    //B[i] += GlobalGradient[j][i] * (*sol_)[(*ptMesh2PDENode_) [ptElement->connect[j]-1]-1];
	    B[i] += GlobalGradient[j][i] * solEntry;
	  }
      //account, that we compute the curl!
      Double temp = B[0];
      if (isaxi_)
	{
	  B[0] = -B[1];
	  B[1] = temp;
	}
      else
	{
	  B[0] = B[1];
	  B[1] = -temp;
	}
    }

  else if (dim ==3)
    Error("CalcElemCurlNode for 3D not implemented",__FILE__,__LINE__);
 
}



} // namespace CoupledField
