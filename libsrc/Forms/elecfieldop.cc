#include "Forms/elecfieldop.hh"

#include <Elements/basefe.hh>
#include <string>
#include <Domain/elem.hh>
#include <Domain/grid.hh>
#include <Matrix/matrix.hh>


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


} // namespace CoupledField
