#include "elecforceop.hh"
#include <string>
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"
#include "PDE/basePDE.hh"

namespace CoupledField
{
 
ElecForceOp::ElecForceOp(Grid * ptGrid,
			 BasePDE * ptPDE,
			 NodeEQN * ptEQN,
			 NodeStoreSol<Double> & EPotential,
			 Integer level,
			 Boolean isaxi) 
  : BaseOperator(ptGrid, ptPDE, ptEQN, level, isaxi)
{
  ENTER_FCN( "ElecForceOp::ElecForceOp" );

  gradFieldOp_ = new GradientFieldOp(ptGrid, ptPDE, ptEQN, 
				     EPotential, ELEC_POTENTIAL, level, isaxi);

}

ElecForceOp::~ElecForceOp()
{
  ENTER_FCN( "ElecForceOp::~ElecForceOp" );

  if (gradFieldOp_) 
    delete gradFieldOp_;
}


void ElecForceOp::CalcElemElecForce(ElemStoreSol<Double> & F,
				    const Elem * ptElement,
				    Double epsilon,
				    const StdVector<ShortInt> & IsBoundaryNode)
{
  ENTER_FCN( "ElecForceOp::CalcElemElecForce" );

  Vector<Double> E;
  Vector<Double> * Ip;
  Vector<Double> intWeights = ptElement->ptElem->GetIntWeights();
  Matrix<Double> JInv, dJ_dr, CornerCoords,J;
  ShortInt Dim, NumNodes, NumIntPoints;
  Double DetJ, DetdJ_dr;
  Double elecEntry;
  
  Dim = ptElement->ptElem->GetDim();
  NumNodes = ptElement->ptElem->GetNumNodes();
  NumIntPoints = ptElement->ptElem->GetNumIntPoints();
  Ip = ptElement->ptElem->GetIntPoints();
 
  // Get element coordinates for calculation of J
  ptPDE_->GetElemCoords(ptElement->connect, CornerCoords, level_);
  
  F.SetNumSolutions(1);
  F.SetNumElems(IsBoundaryNode.GetSize());
  F.SetSolutionType(ELEC_FORCE_VWP);
  F.SetNumDofs(Dim);
  F.Init(0.0);
  
  // TEST TEST
  Matrix<Double> J_Trans, J_Inv_Trans, J_r_Trans;

  // Loop over integration points
  for (Integer nIp=1; nIp<NumIntPoints+1; nIp++)
    {
      // Calculate E-Field
      gradFieldOp_->CalcElemGradField(E, ptElement, Ip[nIp-1], 1);
            
      // Calculate J 
      ptElement->ptElem->CalcJacobianAtIp(J, nIp, CornerCoords);

      // Calculate J-1
      ptElement->ptElem->CalcInvJacobianAtIp(JInv, nIp, CornerCoords);        
      JInv.Transpose(J_Inv_Trans);
      JInv = J_Inv_Trans;

      // Calculate Det(J)
      DetJ = ptElement->ptElem->CalcJacobianDetAtIp(nIp, CornerCoords);
          
      Matrix<Double> SpecCornerCoords, J_Transposed;
      SpecCornerCoords.Resize(Dim,NumNodes);

      //account for axisymmetric case
      Double factor = 1;
      if (isaxi_)
	{
	  Vector<Double> shpFncAtIp;
	  Vector<Double> coordAtIP;
	  ptElement->ptElem->GetShFncAtIp(shpFncAtIp, nIp);
	  coordAtIP = CornerCoords * shpFncAtIp;
	  factor = 2 * PI * coordAtIP[0];
	}     

      // loop over all boundary nodes
      for (Integer nNode=0; nNode<IsBoundaryNode.GetSize(); nNode++)
	{
	  // loop over all dimension
	  for (Integer i=0; i<Dim; i++) 
	    {
	      // form SpecCornerCoords-Array
	      SpecCornerCoords.Init();
	  
	      if (IsBoundaryNode[nNode] == 1)
		{
		SpecCornerCoords[i][nNode] = 1;
		}
	      else {
		break;
	      }


	      // calculate dJ_dr and Det(dJ_dr)
	      ptElement->ptElem->CalcJacobianAtIp(dJ_dr, nIp, SpecCornerCoords);

	      //std::cerr << "dJ_dr = " << std::endl << dJ_dr << std::endl;
	      DetdJ_dr = CalcDetJDr(J, dJ_dr, i);
	      dJ_dr.Transpose(J_r_Trans);
	      dJ_dr = J_r_Trans;
	      
	      // Force Calculation
	      F(nNode,i) += intWeights[nIp-1] * factor * ( (E * ( JInv * (dJ_dr * E) ) * -DetJ  
							     + ( E * E ) * DetdJ_dr * 0.5) * epsilon);
	    } // loop over dimension
	} // loop over boundary nodes
    } // loop over integration points
  
}

// void ElecForceOp::CalcElemElecForce(Vector<Double> & F,
// 				    const Elem * ptElement,
// 				    Double epsilon,
// 				    const StdVector<ShortInt> & IsBoundaryNode)
// {
// ENTER_FCN( "ElecForceOp::CalcElemElecForce");


//   Vector<Double> E;
//   StdVector<Double> * Ip;
//   StdVector<Double> intWeights = ptElement->ptElem->GetIntWeights();
//   Matrix<Double> JInv, dJ_dr, CornerCoords,J;
//   ShortInt Dim, NumNodes, NumIntPoints;
//   Double DetJ, DetdJ_dr;
  
//   Dim = ptElement->ptElem->GetDim();
//   NumNodes = ptElement->ptElem->GetNumNodes();
//   NumIntPoints = ptElement->ptElem->GetNumIntPoints();
//   Ip = ptElement->ptElem->GetIntPoints();
 
//   ptPDE_->GetElemCoords(ptElement->connect, CornerCoords, level_);
  
//   F.Resize(Dim);
//   F.Init();
   
//   StdVector<Double> temp;
//   temp.Resize(Dim);

//   // TEST TEST
//   Matrix<Double> J_Trans, J_Inv_Trans, J_r_Trans;

  
//   // Loop over integration points
//   for (Integer nIp=1; nIp<NumIntPoints+1; nIp++)
//     {
//       // Calculate E-Field
//       //std::cerr << "numIntPoints = " << NumIntPoints << std::endl;
//       //std::cerr << Ip[0].size() << std::endl;
//       //std::cerr << "Ip[" << nIp <<"] = " << Ip[nIp][0] << "," << Ip[nIp][1] << std::endl;
//       ElecFieldOp_->CalcElemElecField(E, ptElement, Ip[nIp-1]);
      
//       //std::cerr << "Element [" << ptElement->elemNum << "] E = " << E[0] << ", " << E[1] << std::endl; 
//       //std::cerr << "CornerCoords:" << std::endl << CornerCoords << std::endl;
      
//       // Calculate J 
//       ptElement->ptElem->CalcJacobianAtIp(J, nIp, CornerCoords);
//       //std::cerr << "J = " << std::endl << J << std::endl;
      
//       //J.Transpose(J_Trans);
//       //J = J_Trans;

//       // Calculate J-1
//       ptElement->ptElem->CalcInvJacobianAtIp(JInv, nIp, CornerCoords);
//       //std::cerr << "JInv = " << std::endl << JInv << std::endl;
      
//       //JInv.Transpose(J_Inv_Trans);
//       //JInv = J_Inv_Trans;

//       // Calculate Det(J)
//       DetJ = ptElement->ptElem->CalcJacobianDetAtIp(nIp, CornerCoords);
//       //std::cerr << "DetJ = " << std::endl <<DetJ << std::endl;
      
//       Matrix<Double> SpecCornerCoords, J_Transposed;
//       SpecCornerCoords.Resize(Dim,NumNodes);
      
//       // 
//       for (Integer i=0; i<Dim; i++) 
// 	{
// 	  // form SpecCornerCoords-Array
// 	  SpecCornerCoords.Init();
	  
// 	  for( Integer j=0; j<NumNodes; j++)
// 	    if (IsBoundaryNode[j] == 1)
// 	      SpecCornerCoords[i][j] = 1;

// 	  if (nIp == 1)
// 	    std::cerr << "SpecCornerArray: " << std::endl << SpecCornerCoords << std::endl;

// 	  // calculate dJ_dr and Det(dJ_dr)
// 	  ptElement->ptElem->CalcJacobianAtIp(dJ_dr, nIp, SpecCornerCoords);

// 	  //std::cerr << "dJ_dr = " << std::endl << dJ_dr << std::endl;
// 	  DetdJ_dr = CalcDetJDr(J, dJ_dr, i);
// 	  //dJ_dr.Transpose(J_r_Trans);
// 	  //dJ_dr = J_r_Trans;

	 
// 	  F[i] += intWeights[nIp-1] * (( E * ( JInv * (dJ_dr * E) ) * DetJ 
// 	  	                         -  ( E * E ) * DetdJ_dr / 2.0) * epsilon);
	  
// 	  std::cerr << "IP[" << nIp << "], Dim[" << i+1 << "]" << std::endl;
// 	  std::cerr << "J = " <<  std::endl << J << std::endl;
// 	  std::cerr << "dJ_dr = " << std::endl << dJ_dr << std::endl;
// 	  std::cerr << "E * ( JInv * (dJ_dr * E) ) * DetJ = " <<  E * ( JInv * (dJ_dr * E) ) * DetJ << std::endl;
// 	  std::cerr << "( E * E ) * DetdJ_dr * 2.0) = "<<  ( E * E ) * DetdJ_dr * 2.0 << std::endl;

// 	  // std::cerr << "F[" << i << "] = " << F[i] << std::endl;
// 	  //std::cerr << "E*E = " << E * E << std::endl;
// 	}
//       //std::cerr << "After Ip " << nIp << " F = " << F << std::endl;
//     }

//   //std::cerr << std::endl;
//   //std::cerr << "F = <<" << F << std::endl;
// }

Double ElecForceOp::CalcDetJDr(Matrix<Double> &J, Matrix<Double> &dJ_dr, Integer dim)
{
  ENTER_FCN( "ElecForceOp::CalcDetJDr" );
  
  Double det;

  if (J.GetSizeRow() == 2)
    {
      det = dJ_dr[0][0]*J[1][1]+dJ_dr[1][1]*J[0][0]-dJ_dr[0][1]*J[1][0]-dJ_dr[1][0]*J[0][1]; 
    } else {
      
      det = dJ_dr[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
	 +  dJ_dr[0][1] * (J[1][0]*J[2][2] - J[1][2]*J[2][0])
	 +  dJ_dr[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
	 +  dJ_dr[1][0] * (J[0][1]*J[2][2] - J[0][2]*J[2][1])
	 +  dJ_dr[1][1] * (J[0][0]*J[2][2] - J[0][2]*J[2][0])
	 +  dJ_dr[1][2] * (J[0][0]*J[2][1] - J[0][1]*J[2][0])
	 +  dJ_dr[2][0] * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
	 +  dJ_dr[2][1] * (J[0][0]*J[1][2] - J[0][2]*J[1][0])
	 +  dJ_dr[2][2] * (J[0][0]*J[1][1] - J[0][1]*J[1][0]);
     
    }

  return det;
}

} // end of namespace CoupledField
