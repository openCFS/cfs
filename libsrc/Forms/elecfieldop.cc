#include "Forms/elecfieldop.hh"

#include <PDE/basepde.hh>
#include <Elements/basefe.hh>
#include <string>
#include <Domain/elem.hh>
#include <Domain/grid.hh>
#include <Matrix/matrix.hh>
#include <General/environment.hh>
#include <PDE/basepde.hh>

namespace CoupledField
{

ElecFieldOp::ElecFieldOp(Grid * ptGrid, 
			 BasePDE * ptPDE,
			 std::vector<Integer> * ptMesh2PDENode,
			 Vector<Double> * EPotential,
			 const Integer level) : BaseOperator( ptGrid, ptPDE, ptMesh2PDENode, level)
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
  
  ShortInt dim;
  dim = ptElement->ptElem->GetDim();
  E.Resize(dim);
  E.Init();

  Integer nShFnc = 0;
  nShFnc = ptElement->ptElem->GetNumNodes();

  Matrix<Double> CornerCoords; 
  ptPDE_->GetElemCoords(ptElement->connect, CornerCoords, level_);

  Matrix<Double> GlobalGradient;

  ptElement->ptElem->GetGlobDerivShFnc(GlobalGradient, LCoord, CornerCoords);

  // loop over shape functions
  for( Integer i=0; i<dim; i++ )
    for( Integer j=0; j<nShFnc; j++ )
      E[i] -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [ptElement->connect[j]-1]-1];
  
}



void ElecFieldOp::CalcSDElecField(Array<Double> & E,
				  const std::vector<std::string> & SD, 
				  const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering ElecFieldOp::CalcSDElecField" << std::endl;
#endif

  
  Integer nShFnc = 0;
  ShortInt dim;
  Matrix<Double> CornerCoords;
  Matrix<Double> GlobalGradient;
  
  std::vector<Elem *> SubDomain;
  
  Integer maxelem;
  maxelem = ptGrid_->GetMaxnumElem(level_, SD);
  dim = ptGrid_->GetDim();

  E.reshape(dim, maxelem);
  E.init();
            
  // Iterate over all subdomains
  for( Integer iSD=0; iSD<SD.size(); iSD++)
    {
      ptGrid_->GetElemSD(SubDomain,SD[iSD],level_);

      // Iterate over whole SubDomain
      for( Integer k=0; k<SubDomain.size(); k++) 
	{
	  nShFnc = SubDomain[k]->ptElem->GetNumNodes();
	  
	  ptPDE_->GetElemCoords( SubDomain[k]->connect, CornerCoords, level_ );
	  
	  SubDomain[k]->ptElem->GetGlobDerivShFnc(GlobalGradient, LCoord, CornerCoords);
	  
	  // loop over shape functions
	  for( Integer i=0; i<dim; i++ )
	    for( Integer j=0; j<nShFnc; j++ )
	      E[i][k] -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [SubDomain[k]->connect[j]-1]-1];	    
	  
	}
    }
}



CurlEdgeOp::CurlEdgeOp(Grid * ptGrid, 
		       BasePDE * ptPDE,
		       std::vector<Integer> * ptMesh2PDENode,
		       Vector<Double> * aSol,
		       const Integer level,
		       BaseSystem * algsys) 
  : BaseOperator(ptGrid, ptPDE, ptMesh2PDENode, level), sol_(aSol), algsys_(algsys)
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
  std::vector<Integer> esign(nrEdges);


  Vector<Integer> pos(nrNodes);
  
  for (Integer i=0; i < nrNodes; i++)
     pos[i] = (*ptMesh2PDENode_)[ptElement->connect[i]-1];
  
  ptElement->ptElem->GetGlobalEdgeIndices(epos, &pos[0], algsys_);


#ifdef DEBUG
  (*debug) << "CurlOP pos \n" << pos << std::endl
	   << "epos \n " << epos << std::endl;
  
#endif

  for (Integer j=0; j<nrEdges; j++)
    {
      esign[j] = epos[j]/abs(epos[j]);
      epos[j]  = abs(epos[j]);
      sol[j] = (*sol_)[epos[j]-1] * esign[j];
    }
  
  
  // loop over edge curls
  for( Integer i=0; i<dim; i++ )
    {
      curlField[i]=0;
      
      for( Integer j=0; j < nrEdges; j++ )
	curlField[i] += curlOnEdges[i][j] * sol[j];
    }
  
}




void CurlEdgeOp::CalcElemMagVec(Vector<Double> & magVecPot, 
				const Elem * ptElement,
				const std::vector<Double> & lCoord)
{
#ifdef TRACE
  (*trace) << "entering CurlEdgeOp::CalcElemMagVec" << std::endl;
#endif
  
  Integer nrEdges = ptElement->ptElem->GetNumEdges();
  Integer nrNodes = ptElement->ptElem->GetNumNodes();
  BaseFE * ptElem = ptElement->ptElem;
  ShortInt dim = ptElem->GetDim();

  Matrix<Double> cornerCoords; 
  Matrix<Double> shape;

  std::vector<Double> sol(nrEdges);
  // global edge index
  std::vector<Integer> epos(nrEdges);
  std::vector<Integer> esign(nrEdges);
  Vector<Integer> pos(nrNodes);



  magVecPot.Resize(dim);

  for (Integer i=0; i<dim; i++)
    magVecPot[i] = 0;

  ptGrid_->GetCoordNodesElemMat(ptElement->connect, cornerCoords, level_);
  
  ptElem->CalcEdgeShapeFnc(shape, lCoord, cornerCoords);

  
  for (Integer i=0; i < nrNodes; i++)
     pos[i] = (*ptMesh2PDENode_)[ptElement->connect[i]-1];
  
  ptElem->GetGlobalEdgeIndices(epos, &pos[0], algsys_);


  for (Integer j=0; j<nrEdges; j++)
    {
      esign[j] = epos[j]/abs(epos[j]);
      epos[j]  = abs(epos[j]);
      sol[j] = (*sol_)[epos[j]-1] * esign[j];
    }
  

  
  
  // loop over edge curls
  // magVecPot = sol * shape;
  for( Integer j=0; j<dim; j++ )  
    {
      magVecPot[j]=0;
      for( Integer i=0; i < nrEdges; i++ )
	magVecPot[j] += shape[i][j] * sol[i];
    }
}




} // namespace CoupledField
