#include "Forms/elecfieldop.hh"

#include <string>

#include "Elements/basefe.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"

#include <PDE/basePDE.hh>

namespace CoupledField
{

ElecFieldOp::ElecFieldOp(Grid * ptGrid, 
			 BasePDE * ptPDE,
			 NodeEQN * ptEQN,
			 NodeStoreSol<Double> & EPotential,
			 const Integer level,
			 Boolean isaxi)
  : BaseOperator( ptGrid, ptPDE, ptEQN, level, isaxi)
{
  ENTER_FCN( "ElecFieldOp::ElecFieldOp" );
  
  this->EPotential_ = &EPotential;
  
}

ElecFieldOp::~ElecFieldOp()
{
  ENTER_FCN( "ElecFieldOp::~ElecFieldOp" );

}

void ElecFieldOp::CalcElemElecField(Vector<Double> & E, 
				    const Elem * ptElement,
				    const Vector<Double> & LCoord)
{
  ENTER_FCN( "ElecFieldOp::CalcElemElecFieldOp" );
  
  ShortInt dim;
  Double elecEntry;
  dim = ptElement->ptElem->GetDim();
  E.Resize(dim);
  E.Init();

  Integer nShFnc = 0;
  nShFnc = ptElement->ptElem->GetNumNodes();
  
  const StdVector<Integer> & connect = ptElement->connect;
  
  Matrix<Double> CornerCoords; 
  ptPDE_->GetElemCoords(connect, CornerCoords, level_);

  Matrix<Double> GlobalGradient;

  ptElement->ptElem->GetGlobDerivShFnc(GlobalGradient, LCoord, CornerCoords);

  // loop over shape functions
  for( Integer i=0; i<dim; i++ )
    for( Integer j=0; j<nShFnc; j++ )
      {
	//std::cerr << "Longing for connect = " << connect[j] << std::endl;
	EPotential_->Get(connect[j]-1,0,elecEntry);
	//istd:cerr << "elecEntry = " << elecEntry << std::endl;
	//E[i] -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [connect[j]-1]-1];
	E[i] -= GlobalGradient[j][i] * elecEntry;
      }
  
}



void ElecFieldOp::CalcSDElecField(NodeStoreSol<Double> & E,
				  const StdVector<std::string> & SD, 
				  const Vector<Double> & LCoord)
{
  ENTER_FCN( "ElecFieldOp::CalcSDElecField" );
  
  Integer nShFnc = 0;
  ShortInt dim;
  Matrix<Double> CornerCoords;
  Matrix<Double> GlobalGradient;
  
  StdVector<Elem *> SubDomain;
  Double elecEntry;
  Integer maxelem;
  maxelem = ptGrid_->GetMaxnumElem(level_, SD);
  dim = ptGrid_->GetDim();

  E.SetNumSolutions(1);
  E.SetNumNodes(maxelem);
  E.SetNumDofs(dim);
            
  // Iterate over all subdomains
  for( Integer iSD=0; iSD<SD.GetSize(); iSD++)
    {
      ptGrid_->GetElemSD(SubDomain,SD[iSD],level_);

      // Iterate over whole SubDomain
      for( Integer k=0; k<SubDomain.GetSize(); k++) 
	{
	  nShFnc = SubDomain[k]->ptElem->GetNumNodes();
	  
	  ptPDE_->GetElemCoords( SubDomain[k]->connect, CornerCoords, level_ );
	  
	  SubDomain[k]->ptElem->GetGlobDerivShFnc(GlobalGradient, LCoord, CornerCoords);
	  
	  // loop over shape functions
	  for( Integer i=0; i<dim; i++ )
	    for( Integer j=0; j<nShFnc; j++ )
	      {
		// NOT WORKING YET
		//elecEntry = (*EPotential_)(SubDomain[k]->connect[j],1);
		///E(k,i) -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [SubDomain[k]->connect[j]-1]-1];	    
		//E(k,i) -= GlobalGradient[j][i] * elecEntry;
	      }
	  
	}
    }
}



CurlEdgeOp::CurlEdgeOp(Grid * ptGrid, 
		       BasePDE * ptPDE,
		       NodeEQN * ptEQN,
		       NodeStoreSol<Double> & aSol,
		       const Integer level,
		       BaseSystem * algsys) 
  : BaseOperator(ptGrid, ptPDE, ptEQN, level), algsys_(algsys)
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
		       BasePDE * ptPDE,
		       NodeEQN * ptEQN,
		       NodeStoreSol<Double> & aSol,
		       const Integer level)
  : BaseOperator(ptGrid, ptPDE, ptEQN, level, FALSE)
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
      ptPDE_->GetElemCoords(ptElement->connect, CornerCoords, level_);
      
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
