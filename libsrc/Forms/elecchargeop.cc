#include "elecchargeop.hh"
#include "Domain/grid.hh"
#include "Elements/basefe.hh"
#include "PDE/basePDE.hh"
#include "Domain/elem.hh"

namespace CoupledField
{


ElecChargeOp::ElecChargeOp(Grid * ptGrid,
	       BasePDE * ptPDE,
	       NodeEQN * ptEQN,
	       const Integer level,
	       Boolean isaxi)
  : BaseOperator(ptGrid, ptPDE, ptEQN, level, isaxi)
{
  ENTER_FCN( "ElecChargeOp::ElecChargeOp" );

}

 
ElecChargeOp::~ElecChargeOp()
{
  ENTER_FCN( "ElecChargeOp::~ElecChargeOp" );

}
  

void ElecChargeOp::CalcElemCharge(Double & charge,
				  const Elem * ptElement,
				  const Vector<Double> & lCoord,
				  const Double & eNormalFluxDensity)
{
  ENTER_IFCN( "ElecChargeOp::CalcElemCharge" );
  
  Double jacDet, chargeAux, elemFluxDensity;
  Vector<Double> shFnc, globCoord;
  Matrix<Double> coordMat;
  BaseFE * ptElemFE = ptElement->ptElem;
  Integer nrIntPts, nrNodes;
  
  
  charge = 0;
  chargeAux = 0;
  jacDet =0;
  nrIntPts= ptElemFE->GetNumIntPoints();
  nrNodes = ptElemFE->GetNumNodes();
  const Vector<Double> & intWeights = ptElemFE->GetIntWeights();   
  
  ptPDE_->GetElemCoords(ptElement->connect, coordMat, level_);
  
  // loop over all integration points
  for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
    {
      ptElemFE->GetShFncAtIp(shFnc, actIntPt);
      jacDet = ptElemFE->CalcJacobianDetAtIp(actIntPt, coordMat);
      
      // loop over all shape functions
      for (Integer actShFnc=1; actShFnc<= shFnc.GetSize(); actShFnc++)
	{
	  chargeAux = shFnc[actShFnc-1] * intWeights[actIntPt-1];
	  chargeAux *=  jacDet * eNormalFluxDensity;
	  if (isaxi_)
	    {
	      globCoord = coordMat * shFnc;
	      chargeAux *= 2 * PI * globCoord[0];
	    }
	  charge += chargeAux;
	}
    }
}
 


  

void ElecChargeOp::CalcElemCharges(Vector<Double> & charges,
				   const StdVector<Elem*> & surfElems,
				   const Vector<Double> & lCoord,
				   const Vector<Double> & eNormalFluxDensity)
{
  ENTER_FCN( "ElecChargeOp::CalcElemCharges" );
  
  Double jacDet, charge, elemFluxDensity;
  Vector<Double> shFnc, globCoord;
  Matrix<Double> coordMat;
  BaseFE * ptElem;
  Integer nrIntPts, nrNodes;
  charges.Resize(surfElems.GetSize());
  
  // loop over all surface elements
  for (Integer iElem=0; iElem<surfElems.GetSize(); iElem++)
    {
      ptElem = surfElems[iElem]->ptElem;
      charge = 0;
      jacDet =0;
      nrIntPts= ptElem->GetNumIntPoints();
      nrNodes = ptElem->GetNumNodes();
      const Vector<Double> & intWeights = ptElem->GetIntWeights();   
      
      ptPDE_->GetElemCoords(surfElems[iElem]->connect, coordMat, level_);
      ptElem->GetShFnc(shFnc, lCoord);

      for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
	{
	  jacDet = ptElem->CalcJacobianDet(lCoord, coordMat);
	  charge = shFnc[actIntPt] *intWeights[actIntPt];
	  charge *= jacDet * eNormalFluxDensity[iElem];
	  if (isaxi_)
	    {
	      globCoord = coordMat * shFnc;
	      charge *= 2 * PI * globCoord[0];
	    }
	  
	  charges[iElem] += charge;
	  

	}
    }
}


} // end of namespace
