#include "baseForceOp.hh"

#include <string>

#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"
#include "PDE/StdPDE.hh"
#include "Utils/elemstoresol.hh"
#include "DataInOut/MaterialData.hh"

namespace CoupledField
{
 
  BaseForceOp::BaseForceOp(Grid * ptGrid,
                           StdPDE * ptPDE,
                           NodeEQN * ptEQN,
                           NodeStoreSol<Double> & sol,
                           UInt dim,
                           MaterialData* &matData,
                           StdVector<RegionIdType> & allSubdoms,
                           Boolean isaxi) 
    : BaseOperator(ptGrid, ptPDE, ptEQN, isaxi)
  {
    ENTER_FCN( "BaseForceOp::BaseForceOp" );

    dim_ = dim;
    materialData_ = matData;
    PDEsubdoms_ = allSubdoms;

  }

  BaseForceOp::~BaseForceOp()
  {
    ENTER_FCN( "BaseForceOp::~BaseForceOp" );

  }


  void BaseForceOp::Setup(StdVector<RegionIdType>& neighRegions, 
                          StdVector<UInt>& couplingnodes)
  {
    ENTER_FCN( "BaseForceOp::Setup" );

    couplingNodes_ = couplingnodes;
    neighRegions_  = neighRegions;

    //get the interface elements to the coupling nodes
    ptGrid_->GetElemsNextToNodes(interfaceElems_, couplingnodes, 
                                 neighRegions );

    //get memory
    isBoundaryNode_.Resize(interfaceElems_.GetSize());
    elemNodeToCouplingNode_.Resize(interfaceElems_.GetSize());
  
  
    for (UInt ielem=0; ielem<interfaceElems_.GetSize(); ielem++)
      {
        isBoundaryNode_[ielem].Resize(interfaceElems_[ielem]->connect.GetSize());
        elemNodeToCouplingNode_[ielem].Resize(interfaceElems_[ielem]->connect.GetSize());
      
        // Determine Boundary Nodes
        for (UInt ielemnode=0; ielemnode<isBoundaryNode_[ielem].GetSize(); ielemnode++) {
          for (UInt inodes=0; inodes<couplingnodes.GetSize(); inodes++) {
            if (interfaceElems_[ielem]->connect[ielemnode] == couplingnodes[inodes] ) {
              isBoundaryNode_[ielem][ielemnode] = 1;
              elemNodeToCouplingNode_[ielem][ielemnode] = inodes;
              break;
            } // end if
          }
        }
      } 
  
  }


  void BaseForceOp::CalcNodeForce(Vector<Double> & force, Vector<Double> & totalForce )
  {
    ENTER_FCN( "BaseForceOp::CalcNodeForce" );

    ElemStoreSol<Double> force_temp;
  
    force.Init(0.0);
  
    for (UInt ielem=0; ielem<interfaceElems_.GetSize(); ielem++)
      {
        // Get Material Parameter
        Double matVal;
      
        for (UInt i=0; i<PDEsubdoms_.GetSize(); i++)   {
          if (interfaceElems_[ielem]->regionId == PDEsubdoms_[i]) {
            matVal = GetMatVal(i); 
          }
        }
      
        CalcElemElecForce( force_temp, interfaceElems_[ielem], matVal, isBoundaryNode_[ielem]);
      

        // Add the element force to the according coupling node
        for (UInt ielemnode=0; ielemnode<interfaceElems_[ielem]->connect.GetSize(); ielemnode++)
          for( UInt idim=0; idim<dim_; idim++) {
            force[elemNodeToCouplingNode_[ielem][ielemnode]*dim_+idim] += 
              force_temp(ielemnode,idim);
          }
      }


    totalForce.Resize(dim_);
  
    for (UInt i=0; i<couplingNodes_.GetSize(); i++)
      for (UInt dim=0; dim<dim_; dim++)
        totalForce[dim] += force[i*dim_+dim];

    //std::cerr << "Sum of E-Force:" << std::endl << sum << std::endl;
  
  }



  void BaseForceOp::CalcElemElecForce(ElemStoreSol<Double> & F,
                                      const Elem * ptElement,
                                      Double matVal,
                                      const StdVector<ShortInt> & IsBoundaryNode)
  {
    ENTER_FCN( "BaseForceOp::CalcElemElecForce" );

    Vector<Double> E;
    Vector<Double> * Ip;
    Vector<Double> intWeights = ptElement->ptElem->GetIntWeights();
    Matrix<Double> JInv, dJ_dr, CornerCoords,J;
    UInt Dim, NumNodes, NumIntPoints;
    Double DetJ, DetdJ_dr;
  
    Dim = ptElement->ptElem->GetDim();
    NumNodes = ptElement->ptElem->GetNumNodes();
    NumIntPoints = ptElement->ptElem->GetNumIntPoints();
    Ip = ptElement->ptElem->GetIntPoints();
 
    // Get element coordinates for calculation of J
    ptPDE_->GetElemCoords(ptElement->connect, CornerCoords);
  
    F.SetNumSolutions(1);
    F.SetNumElems(IsBoundaryNode.GetSize());
    F.SetSolutionType(solType_);
    F.SetNumDofs(Dim);
    F.Init(0.0);
  
    // Needed variables
    Matrix<Double> J_Trans, J_Inv_Trans, J_r_Trans;
    Double E_square;
    Vector<Double> dJdrTimesE, JInvTimesdJdr;


    // Loop over integration points
    for (UInt nIp=1; nIp<NumIntPoints+1; nIp++)
      {
        // Calculate E (electric field) / B (magnetic field) -Field
        ComputeField(E, ptElement, Ip[nIp-1]);

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
        for (UInt nNode=0; nNode<IsBoundaryNode.GetSize(); nNode++)
          {
            // loop over all dimension
            for (UInt i=0; i<Dim; i++) 
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
                // F(nNode,i) += ( intWeights[nIp-1] * factor * 
//                                 ( (E * ( JInv * (dJ_dr * E) ) * (-DetJ)  
//                                    + ( E * E ) * DetdJ_dr * 0.5) * matVal) ) * sign_;
                
                E_square = E*E;
                dJdrTimesE = dJ_dr * E;
                JInvTimesdJdr = JInv * dJdrTimesE;
                F(nNode,i) += ( intWeights[nIp-1] * factor * 
                                ( (-DetJ) *(E * JInvTimesdJdr)  
                                + E_square * DetdJ_dr * 0.5) * matVal)  * sign_;


                
              } // loop over dimension
          } // loop over boundary nodes
      } // loop over integration points
  
  }

  // void BaseForceOp::CalcElemElecForce(Vector<Double> & F,
  //                                  const Elem * ptElement,
  //                                  Double epsilon,
  //                                  const StdVector<ShortInt> & IsBoundaryNode)
  // {
  // ENTER_FCN( "BaseForceOp::CalcElemElecForce");


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
 
  //   ptPDE_->GetElemCoords(ptElement->connect, CornerCoords);
  
  //   F.Resize(Dim);
  //   F.Init();
   
  //   StdVector<Double> temp;
  //   temp.Resize(Dim);

  //   // TEST TEST
  //   Matrix<Double> J_Trans, J_Inv_Trans, J_r_Trans;

  
  //   // Loop over integration points
  //   for (UInt nIp=1; nIp<NumIntPoints+1; nIp++)
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
  //       for (UInt i=0; i<Dim; i++) 
  //      {
  //        // form SpecCornerCoords-Array
  //        SpecCornerCoords.Init();
          
  //        for( UInt j=0; j<NumNodes; j++)
  //          if (IsBoundaryNode[j] == 1)
  //            SpecCornerCoords[i][j] = 1;

  //        if (nIp == 1)
  //          std::cerr << "SpecCornerArray: " << std::endl << SpecCornerCoords << std::endl;

  //        // calculate dJ_dr and Det(dJ_dr)
  //        ptElement->ptElem->CalcJacobianAtIp(dJ_dr, nIp, SpecCornerCoords);

  //        //std::cerr << "dJ_dr = " << std::endl << dJ_dr << std::endl;
  //        DetdJ_dr = CalcDetJDr(J, dJ_dr, i);
  //        //dJ_dr.Transpose(J_r_Trans);
  //        //dJ_dr = J_r_Trans;

         
  //        F[i] += intWeights[nIp-1] * (( E * ( JInv * (dJ_dr * E) ) * DetJ 
  //                                       -  ( E * E ) * DetdJ_dr / 2.0) * epsilon);
          
  //        std::cerr << "IP[" << nIp << "], Dim[" << i+1 << "]" << std::endl;
  //        std::cerr << "J = " <<  std::endl << J << std::endl;
  //        std::cerr << "dJ_dr = " << std::endl << dJ_dr << std::endl;
  //        std::cerr << "E * ( JInv * (dJ_dr * E) ) * DetJ = " <<  E * ( JInv * (dJ_dr * E) ) * DetJ << std::endl;
  //        std::cerr << "( E * E ) * DetdJ_dr * 2.0) = "<<  ( E * E ) * DetdJ_dr * 2.0 << std::endl;

  //        // std::cerr << "F[" << i << "] = " << F[i] << std::endl;
  //        //std::cerr << "E*E = " << E * E << std::endl;
  //      }
  //       //std::cerr << "After Ip " << nIp << " F = " << F << std::endl;
  //     }

  //   //std::cerr << std::endl;
  //   //std::cerr << "F = <<" << F << std::endl;
  // }

  Double BaseForceOp::CalcDetJDr(Matrix<Double> &J, Matrix<Double> &dJ_dr, UInt dim)
  {
    ENTER_FCN( "BaseForceOp::CalcDetJDr" );
  
    Double det;

    if (J.GetSizeRow() == 2)
      {
        det = dJ_dr[0][0]*J[1][1]+dJ_dr[1][1]*J[0][0]-dJ_dr[0][1]*J[1][0]-dJ_dr[1][0]*J[0][1]; 
      } else {
      
        det = dJ_dr[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
          -  dJ_dr[0][1] * (J[1][0]*J[2][2] - J[1][2]*J[2][0])
          +  dJ_dr[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
          -  dJ_dr[1][0] * (J[0][1]*J[2][2] - J[0][2]*J[2][1])
          +  dJ_dr[1][1] * (J[0][0]*J[2][2] - J[0][2]*J[2][0])
          -  dJ_dr[1][2] * (J[0][0]*J[2][1] - J[0][1]*J[2][0])
          +  dJ_dr[2][0] * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
          -  dJ_dr[2][1] * (J[0][0]*J[1][2] - J[0][2]*J[1][0])
          +  dJ_dr[2][2] * (J[0][0]*J[1][1] - J[0][1]*J[1][0]);
     
      }

    return det;
  }

} // end of namespace CoupledField
