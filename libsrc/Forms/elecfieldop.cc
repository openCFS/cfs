#include "Forms/elecfieldop.hh"

#include <Elements/basefe.hh>
#include <string>
#include <Domain/elem.hh>
#include <Domain/grid.hh>
#include <Matrix/matrix.hh>
#include <General/environment.hh>

namespace CoupledField
{

ElecFieldOp::ElecFieldOp(Grid * ptGrid, 
			 std::vector<Integer> * ptMesh2PDENode,
			 Vector<Double> * EPotential,
			 const Integer level) : BaseOperator( ptGrid, ptMesh2PDENode, level)
{
#ifdef TRACE
  (*trace) << "entering ElecFieldOp::ElecFieldOp" << std::endl;
#endif
  
  this->EPotential_ = EPotential;

}

ElecFieldOp::~ElecFieldOp()
{
#ifdef TRACE
  (*trace) << "entering ElecFieldOp::~ElecFieldOp" << std::endl;
#endif

}

void ElecFieldOp::CalcElemElecField(Vector<Double> & E, 
				    const Elem * ptElement,
				    const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering ElecFieldOp::CalcElemElecFieldOp" << std::endl;
#endif
  
  ShortInt Dim;
  Dim = ptElement->ptElem->GetDim();
  
  E.Resize(Dim);
  
  Integer nShFnc = 0;
  nShFnc = ptElement->ptElem->GetNumNodes();

  Matrix<Double> CornerCoords; 
  ptGrid_->GetCoordNodesElemMat(ptElement->connect, CornerCoords, level_);


  for( Integer i=0; i<Dim; i++)
    E[i] = 0.0;
  
  Matrix<Double> GlobalGradient;

  ptElement->ptElem->GetGlobDerivShFnc(GlobalGradient, LCoord, CornerCoords);

  // loop over shape functions
  for( Integer i=0; i<Dim; i++ )
    for( Integer j=0; j<nShFnc; j++ )
      E[i] -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [ptElement->connect[j]-1]-1];
  
}



void ElecFieldOp::CalcSDElecField(Vector<Double> *  E,
				  const std::vector<std::string> & SD, 
				  const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering ElecFieldOp::CalcSDElecField" << std::endl;
#endif

  
  Integer nShFnc = 0;
  ShortInt Dim;
  Matrix<Double> CornerCoords;
  Matrix<Double> GlobalGradient;
  
  std::vector<Elem *> SubDomain;
  
  Integer maxelem;
  maxelem = ptGrid_->GetMaxnumElem(level_, SD);
  Dim = ptGrid_->GetDim();

  for( Integer i=0; i<Dim; i++)
    E[i].Resize(maxelem);
            
  // Iterate over all subdomains
  for( Integer iSD=0; iSD<SD.size(); iSD++)
    {
      ptGrid_->GetElemSD(SubDomain,SD[iSD],level_);

      // Iterate over whole SubDomain
      for( Integer k=0; k<SubDomain.size(); k++) 
	{
	  nShFnc = SubDomain[k]->ptElem->GetNumNodes();
	  
	  ptGrid_->GetCoordNodesElemMat( SubDomain[k]->connect, CornerCoords, level_ );
	  
	  SubDomain[k]->ptElem->GetGlobDerivShFnc(GlobalGradient, LCoord, CornerCoords);
	  
	  // loop over shape functions
	  for( Integer i=0; i<Dim; i++ )
	    for( Integer j=0; j<nShFnc; j++ )
	      E[i][k] -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [SubDomain[k]->connect[j]-1]-1];	    
	  
	}
    }
}




CurlEdgeOp::CurlEdgeOp(Grid * ptGrid, 
		       std::vector<Integer> * ptMesh2PDENode,
		       Vector<Double> * aSol,
		       const Integer level,
		       BaseSystem * algsys) 
  : BaseOperator(ptGrid, ptMesh2PDENode, level), sol_(aSol), algsys_(algsys)
{
#ifdef TRACE
  (*trace) << "entering CurlEdgeOp::CurlEdgeOp" << std::endl;
#endif
}


CurlEdgeOp::~CurlEdgeOp()
{
#ifdef TRACE
  (*trace) << "entering CurlEdgeOp::~CurlEdgeOp" << std::endl;
#endif
}


void CurlEdgeOp::CalcElemCurlEdge(Vector<Double> & curlField, 
				  const Elem * ptElement,
				  const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering CurlEdgeOp::CalcElemCurlEdgeOp" << std::endl;
#endif
  
  ShortInt dim = ptElement->ptElem->GetDim();
  
  curlField.Resize(dim);
  for (Integer i=0; i<dim; i++)
    curlField[i] = 0;

  Integer nrEdges = ptElement->ptElem->GetNumEdges();
  Integer nrNodes = ptElement->ptElem->GetNumNodes();

  Matrix<Double> CornerCoords; 
  ptGrid_->GetCoordNodesElemMat(ptElement->connect, CornerCoords, level_);


  Matrix<Double> curlOnEdges;


  std::vector<Matrix<Double>* > deriv;
  deriv.resize(nrEdges);
  for (Integer actEdge=0; actEdge < nrEdges; actEdge++)
    deriv[actEdge] = new Matrix<Double>;
  

  ptElement->ptElem->GetEdgeGlobalDerivShapeFnc (deriv, LCoord, CornerCoords);

  curlOnEdges.Resize(dim, nrEdges);
  
  for (Integer actEdge=0; actEdge < nrEdges; actEdge++)
    for (Integer actDim=0; actDim < dim; actDim++)
      curlOnEdges[actDim][actEdge] = 
	(*deriv[actEdge])[(actDim+2)%dim][(actDim+1)%dim] -
	(*deriv[actEdge])[(actDim+1)%dim][(actDim+2)%dim];


  std::vector<Double> sol(nrEdges);
  // global edge index
  std::vector<Integer> epos(nrEdges);


  Vector<Integer> pos(nrNodes);
  
  for (Integer i=0; i < nrNodes; i++)
     pos[i] = (*ptMesh2PDENode_)[ptElement->connect[i]-1];
  
  epos[0] = algsys_->GetNode2Edge(pos[3], pos[0]);
  epos[1] = algsys_->GetNode2Edge(pos[3], pos[1]);
  epos[2] = algsys_->GetNode2Edge(pos[3], pos[2]);
  epos[3] = algsys_->GetNode2Edge(pos[0], pos[1]);
  epos[4] = algsys_->GetNode2Edge(pos[0], pos[2]);
  epos[5] = algsys_->GetNode2Edge(pos[1], pos[2]);

#ifdef DEBUG
  (*debug) << "CurlOP pos \n" << pos << std::endl
	   << "epos \n " << epos << std::endl;
  
#endif

  for (Integer j=0; j<nrEdges; j++)
    {
      epos[j]  = abs(epos[j]);
      sol[j] = (*sol_)[epos[j]-1] ; //* (Double)(epos[j]/abs(epos[j]));
    }
  
  
  // loop over edge curls
  for( Integer i=0; i<dim; i++ )
    {
      curlField[i]=0;
      
      for( Integer j=0; j < nrEdges; j++ )
	curlField[i] += curlOnEdges[i][j] * sol[j];
    }
  
}




} // namespace CoupledField
