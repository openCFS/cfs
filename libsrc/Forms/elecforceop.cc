#include "elecforceop.hh"

#include <PDE/basepde.hh>
#include <Elements/basefe.hh>
#include <string>
#include <Domain/elem.hh>
#include <Domain/grid.hh>
#include <Matrix/matrix.hh>

namespace CoupledField
{
 
ElecForceOp::ElecForceOp(Grid * ptGrid,
			 BasePDE * ptPDE,
			 std::vector<Integer> * ptMesh2PDENode,
			 Vector<Double> * EPotential,
			 Integer level) : BaseOperator(ptGrid, ptPDE, ptMesh2PDENode, level)
{
#ifdef TRACE
  (*trace) << "entering ElecForceOp::ElecForceOp" << std::endl;
#endif

  ElecFieldOp_ = new ElecFieldOp(ptGrid, ptPDE, ptMesh2PDENode, EPotential, level);

}

ElecForceOp::~ElecForceOp()
{
#ifdef TRACE
  (*trace) << "entering ElecForceOp::~ElecForceOp" << std::endl;
#endif

  if (ElecFieldOp_) delete ElecFieldOp_;
}

void ElecForceOp::CalcElemElecForce(Vector<Double> & F,
				    const Elem * ptElement,
				    Double epsilon,
				    const std::vector<ShortInt> & IsBoundaryNode)
{
#ifdef TRACE
  (*trace) << "entering ElecForceOp::CalcElemElecForce" << std::endl;
#endif


  Vector<Double> E;
  std::vector<Double> * Ip;
  Matrix<Double> JInv, dJ_dr, CornerCoords,J;
  ShortInt Dim, NumNodes, NumIntPoints;
  Double DetJ, DetdJ_dr;
  
  Dim = ptElement->ptElem->GetDim();
  NumNodes = ptElement->ptElem->GetNumNodes();
  NumIntPoints = ptElement->ptElem->GetNumIntPoints();
  Ip = ptElement->ptElem->GetIntPoints();
 
  ptPDE_->GetElemCoords(ptElement->connect, CornerCoords, level_);
  
  F.Resize(Dim);
  F.Init();
   
  std::vector<Double> temp;
  temp.resize(Dim);

  //std::cerr << std::endl << std::endl;
  //std::cerr << "----------------------------" << std::endl;

  // Loop over integration points
  for (Integer nIp=1; nIp<NumIntPoints+1; nIp++)
    {
      // Calculate E-Field
      //std::cerr << "numIntPoints = " << NumIntPoints << std::endl;
      //std::cerr << Ip[0].size() << std::endl;
      //std::cerr << "Ip[" << nIp <<"] = " << Ip[nIp][0] << "," << Ip[nIp][1] << std::endl;
      ElecFieldOp_->CalcElemElecField(E, ptElement, Ip[nIp-1]);
      
      //std::cerr << "Element [" << ptElement->ElemNum << "] E = " << E[0] << ", " << E[1] << std::endl;
      
      //std::cerr << "CornerCoords:" << std::endl << CornerCoords << std::endl;
      
      // Calculate J
      ptElement->ptElem->CalcJacobianAtIp(J, nIp, CornerCoords);
      //std::cerr << "J = " << std::endl << J << std::endl;
      // Calculate J-1
      ptElement->ptElem->CalcInvJacobianAtIp(JInv, nIp, CornerCoords);
      //std::cerr << "JInv = " << std::endl << JInv << std::endl;
      
      // Calculate Det(J)
      DetJ = ptElement->ptElem->CalcJacobianDetAtIp(nIp, CornerCoords);
      //std::cerr << "DetJ = " << std::endl <<DetJ << std::endl;
      
      Matrix<Double> SpecCornerCoords;
      SpecCornerCoords.Resize(Dim,NumNodes);
      

      // 
      for (Integer i=0; i<Dim; i++) 
	{
	  // form SpecCornerCoords-Array
	  SpecCornerCoords.Init();
	  
	  for( Integer j=0; j<NumNodes; j++)
	    if (IsBoundaryNode[j] == 1)
	      SpecCornerCoords[i][j] = 1;
	  //std::cerr << "SpecBoundary" << i <<" = " << std::endl << SpecCornerCoords << std::endl;
	  
	  // calculate dJ_dr and Det(dJ_dr)
	  ptElement->ptElem->CalcJacobianAtIp(dJ_dr, nIp, SpecCornerCoords);
	  //std::cerr << "dJ_dr = " << std::endl << dJ_dr << std::endl;
	  DetdJ_dr = dJ_dr.Det();
	  //DetdJ_dr = ptElement->ptElem->CalcJacobianDetAtIp(nIp, SpecCornerCoords);
	  
	  // finally calculate electric force per element

	  F[i] -=  0.5 * (( E * ( JInv * (dJ_dr * E) ) * DetJ 
	  	     -  ( E * E ) * DetdJ_dr * 2.0) * epsilon);

	  // std::cerr << "F[" << i << "] = " << F[i] << std::endl;
	  //std::cerr << "E*E = " << E * E << std::endl;
	}
      //std::cerr << "After Ip " << nIp << " F = " << F << std::endl;
    }

  //std::cerr << std::endl;
  //std::cerr << "F = <<" << F << std::endl;
}


} // end of namespace CoupledField
