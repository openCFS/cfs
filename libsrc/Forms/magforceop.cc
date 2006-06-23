#include "magforceop.hh"
#include <string>
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"
#include "PDE/StdPDE.hh"
#include "Materials/baseMaterial.hh"

namespace CoupledField
{
 
  MagLorentzForceOp::MagLorentzForceOp(Grid * ptGrid,
                                       StdPDE * ptPDE,
                                       shared_ptr<EqnMap> eqnMap,
                                       NodeStoreSol<Double> & magPotential,
                                       bool isaxi, bool coordUpdate) 
    : BaseOperator(ptGrid, ptPDE, eqnMap, isaxi, coordUpdate )
  {
    ENTER_FCN( "MagLorentzForceOp::MagLorentzForceOp" );

    curlFieldOp_ = new  CurlNodeOp(ptGrid, ptPDE, eqnMap,magPotential, coordUpdate);
    curlFieldOp_->Set2DType(isaxi);
    if ( ptGrid->GetDim() != 2 ) {
      Error( "Currently MagLorentzForceOp just working for 2D problems",
             __FILE__, __LINE__ );
    }
  }

  MagLorentzForceOp::~MagLorentzForceOp()
  {
    ENTER_FCN( "MagLorentzForceOp::~MagLorentzForceOp" );

    if (curlFieldOp_) 
      delete curlFieldOp_;
  }


  void MagLorentzForceOp::CalcElemMagLorentzForce(Matrix<Double>& F,
                                                  Vector<Double>& Jeddy,
                                                  const Elem * ptElement)
  {
    ENTER_FCN( "MagLorentzForceOp::CalcElemElecForce" );

    Vector<Double> B, B1;
    Vector<Double> * Ip;
    Matrix<Double> CornerCoords;
    Vector<Double> intWeights = ptElement->ptElem->GetIntWeights();

    UInt Dim, NumNodes, NumIntPoints;

    Dim = ptElement->ptElem->GetDim();
    NumNodes = ptElement->ptElem->GetNumNodes();
    NumIntPoints = ptElement->ptElem->GetNumIntPoints();
    Ip = ptElement->ptElem->GetIntPoints();
 
    //resize F and set to zero
    F.Resize(NumNodes,Dim);
    F.Init(0);

    // Get element coordinates
    ptGrid_->GetElemNodesCoord( CornerCoords, ptElement->connect, 
                                coordUpdate_ );

    //just for testing with CAPA
    //   Vector<Double> LCoord;
    //   LCoord.Resize(2);
    //   LCoord[0] = 0;
    //   LCoord[1] = 0;
    //   curlFieldOp_->CalcElemCurlNode(B1, ptElement, LCoord);

    //   Double JeddyAtIp1 = 0;
    //   for (UInt k=0; k<Jeddy.GetSize(); k++)
    //     JeddyAtIp1 += Jeddy[k]*0.25;


    // Loop over integration points
    Double factor = 1.0;
    for (UInt nIp=1; nIp<NumIntPoints+1; nIp++)
      {
        // Calculate B-Field
        curlFieldOp_->CalcElemCurlNode(B, ptElement, Ip[nIp-1]);
      
        Vector<Double> shpFncAtIp;
        ptElement->ptElem->GetShFncAtIp(shpFncAtIp, nIp);
        Double JeddyAtIp =  shpFncAtIp * Jeddy;     
 
        if (isaxi_)
          {
            Vector<Double> coordAtIP;
            ptElement->ptElem->GetShFncAtIp(shpFncAtIp, nIp);
            coordAtIP = CornerCoords * shpFncAtIp;
            factor = 2 * PI * coordAtIP[0];
          }     

        Double jacDet = ptElement->ptElem->CalcJacobianDetAtIp(nIp, CornerCoords);
        JeddyAtIp *= factor * intWeights[nIp-1] * jacDet;

        //compute J x B
        if (B.GetSize() == 2) {
          for (UInt node=0; node<NumNodes; node++) {
            if (isaxi_) {
              F[node][0] += JeddyAtIp*B[1]*shpFncAtIp[node];
              F[node][1] -= JeddyAtIp*B[0]*shpFncAtIp[node];
            }
            else {
              F[node][0] -= JeddyAtIp*B[1]*shpFncAtIp[node];
              F[node][1] += JeddyAtIp*B[0]*shpFncAtIp[node];
            }
          }

          //        just for testing with CAPA
          //       //compute J x B
          //       if (B.GetSize() == 2) {
          //      for (UInt node=0; node<NumNodes; node++) {
          //        if (isaxi_) {
          //          F[node][0] += JeddyAtIp1*B1[1]*shpFncAtIp[node]*factor * intWeights[nIp-1] * jacDet;
          //          F[node][1] -= JeddyAtIp1*B1[0]*shpFncAtIp[node]*factor * intWeights[nIp-1] * jacDet;
          //        }
          //        else {
          //          F[node][0] -= JeddyAtIp1*B1[1]*shpFncAtIp[node]*factor * intWeights[nIp-1] * jacDet;
          //          F[node][1] += JeddyAtIp1*B1[0]*shpFncAtIp[node]*factor * intWeights[nIp-1] * jacDet;
          //        }
          //      }
          //       }

        }

      } // loop over integration points
  }



  //---------------------------- VWP ----------------------------------------------------- 
  MagForceOp::MagForceOp(Grid * ptGrid,
                         StdPDE * ptPDE,
                         shared_ptr<EqnMap> eqnMap,
                         NodeStoreSol<Double> & sol,
                         UInt dim,
                         std::map<RegionIdType,BaseMaterial*>& matData,
                         bool isaxi, bool coordUpdate ) 
    : BaseForceOp(ptGrid, ptPDE,  eqnMap, sol, dim, 
                  matData, isaxi, coordUpdate )
  {
    ENTER_FCN( "MagForceOp::MagForceOp" );

    curlFieldOp_ = new CurlNodeOp(ptGrid, ptPDE, eqnMap, sol);
    curlFieldOp_->Set2DType(isaxi);

    solType_ = MAG_FORCE_VWP;
    sign_    = 1.0;
  }

  MagForceOp::~MagForceOp()
  {
    ENTER_FCN( "MagForceOp::~MagForceOp" );

    if (curlFieldOp_) 
      delete curlFieldOp_;
  }


  void MagForceOp::ComputeField(Vector<Double> & Field, const Elem * ptElement,
                                const Vector<Double> & lCoord)
  {
    ENTER_FCN( "MagForceOp::ComputeField" );

    curlFieldOp_->CalcElemCurlNode(Field, ptElement, lCoord);
  } 

  Double MagForceOp::GetMatVal(RegionIdType actRegion)
  {
    ENTER_FCN( "MagForceOp::GetMatVal" );

    Double reluctivity;
    materials_[actRegion]->GetScalar(reluctivity,MAG_RELUCTIVITY,REAL);

    return reluctivity;
  } 


} // end of namespace CoupledField
