#include "Forms/gradfieldop.hh"

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

GradientFieldOp::GradientFieldOp(Grid * ptGrid, 
				 BasePDE * ptPDE,
				 NodeEQN * ptEQN,
				 NodeStoreSol<Double> & potential,
				 const SolutionType solType,
				 const Integer level,
				 Boolean isaxi)
  : BaseOperator( ptGrid, ptPDE, ptEQN, level, isaxi)
{
  ENTER_FCN( "GradientFieldOp::GradientFieldOp" );
  
  this->potential_ = &potential;
  solType_ = solType;
}

GradientFieldOp::~GradientFieldOp()
{
  ENTER_FCN( "GradientFieldOp::~GradientFieldOp" );

}

void GradientFieldOp::CalcElemGradField(Vector<Double> & elemField,
					const Elem * ptElement,
					const Vector<Double> & lCoord,
					const Double factor)
{
  ENTER_FCN( "GradientFieldOp::CalcElemGradField" );
  
  ShortInt dim;
  Double potEntry;
  dim = ptElement->ptElem->GetDim();
  elemField.Resize(dim);
  elemField.Init();

  Integer nShFnc = 0;
  nShFnc = ptElement->ptElem->GetNumNodes();
  
  const StdVector<Integer> & connect = ptElement->connect;
  
  Matrix<Double> CornerCoords; 
  ptPDE_->GetElemCoords(connect, CornerCoords, level_);

  Matrix<Double> GlobalGradient;

  ptElement->ptElem->GetGlobDerivShFnc(GlobalGradient, lCoord, CornerCoords);

  // loop over shape functions
  for( Integer i=0; i<dim; i++ )
    for( Integer j=0; j<nShFnc; j++ )
      {
	//std::cerr << "Longing for connect = " << connect[j] << std::endl;
	potential_->Get(solType_, connect[j]-1,0,potEntry);
	//istd:cerr << "elecEntry = " << elecEntry << std::endl;
	//E[i] -= GlobalGradient[j][i] * (*EPotential_)[(*ptMesh2PDENode_) [connect[j]-1]-1];
	elemField[i] -= GlobalGradient[j][i] * potEntry * factor;
      }
  
}



void GradientFieldOp::CalcSDGradField(Vector<Double> & elemField,
				      const StdVector<std::string> & SD, 
				      const Vector<Double> & lCoord,
				      const Vector<Double> & factors)
{
  ENTER_FCN( "GradientFieldOp::CalcSDGradField" );
  
  Error( "GradientFieldOp::CalcSDElecField: Not working yet", __FILE__, __LINE__);
  Integer nShFnc = 0;
  ShortInt dim;
  Matrix<Double> cornerCoords;
  Matrix<Double> globalGradient;
  
  StdVector<Elem *> subDomain;
  Double potEntry;
  Integer maxelem;
  maxelem = ptGrid_->GetMaxnumElem(level_, SD);
  dim = ptGrid_->GetDim();

  elemField.Resize(maxelem * dim);
  //elemField.SetNumSolutions(1);
  //elemField.SetNumNodes(maxelem);
  //elemField.SetNumDofs(dim);
            
  // Iterate over all subdomains
  for( Integer iSD=0; iSD<SD.GetSize(); iSD++)
    {
      ptGrid_->GetElemSD(subDomain,SD[iSD],level_);

      // Iterate over whole SubDomain
      for( Integer k=0; k<subDomain.GetSize(); k++) 
	{
	  nShFnc = subDomain[k]->ptElem->GetNumNodes();
	  
	  ptPDE_->GetElemCoords( subDomain[k]->connect, cornerCoords, level_ );
	  
	  subDomain[k]->ptElem->GetGlobDerivShFnc(globalGradient, lCoord, cornerCoords);
	  
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


} // namespace CoupledField
