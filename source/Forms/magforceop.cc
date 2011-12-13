// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Elements/basefe.hh"
#include "Forms/curlfieldop.hh"
#include "General/exception.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "magforceop.hh"

namespace CoupledField {
class EqnMap;
class StdPDE;
template <class TYPE> class NodeStoreSol;
}  // namespace CoupledField

namespace CoupledField
{

  template<class TYPE>
  MagLorentzForceOp<TYPE>::MagLorentzForceOp(Grid * ptGrid,
                                       StdPDE * ptPDE,
                                       shared_ptr<EqnMap> eqnMap,
                                       NodeStoreSol<TYPE> & magPotential,
                                       bool isaxi, bool coordUpdate)
    : BaseOperator(ptGrid, ptPDE, eqnMap, isaxi, coordUpdate )
  {

    //WARN( "Only working with Lagrange Elements");
    curlFieldOp_ = new  CurlNodeOp<TYPE>(ptGrid, ptPDE, eqnMap,
                                         magPotential, coordUpdate);
    curlFieldOp_->Set2DType(isaxi);
    if ( ptGrid->GetDim() != 2 ) {
      EXCEPTION("Currently MagLorentzForceOp just working for 2D problems");
    }
  }

  template<class TYPE>
  MagLorentzForceOp<TYPE>::~MagLorentzForceOp()
  {

    if (curlFieldOp_)
      delete curlFieldOp_;
  }

  template<class TYPE>
  void MagLorentzForceOp<TYPE>::CalcElemMagLorentzForce(Matrix<TYPE>& F,
                                                        Vector<TYPE>& Jeddy,
                                                        const EntityIterator& ent)
  {

    Vector<TYPE> B, B1;
    Vector<Double> * Ip;
    Matrix<Double> CornerCoords;
    Vector<Double> intWeights = ent.GetElem()->ptElem->GetIntWeights();

    UInt Dim, NumNodes, NumIntPoints;

    BaseFE *elem = ent.GetElem()->ptElem;
    Dim = elem->GetDim();

    NumNodes = elem->GetNumNodes();
    NumIntPoints = elem->GetNumIntPoints();
    Ip = elem->GetIntPoints();

    //resize F and set to zero
    F.Resize(NumNodes,Dim);
    F.Init();

    // Get element coordinates
    ptGrid_->GetElemNodesCoord( CornerCoords, ent.GetElem()->connect,
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
    TYPE factor = 1.0;
    for (UInt nIp=1; nIp<NumIntPoints+1; nIp++)
      {
        // Calculate B-Field
        curlFieldOp_->CalcElemCurlNode(B, ent, Ip[nIp-1]);

        Vector<Double> shpFncAtIp;
        elem->GetShFncAtIp(shpFncAtIp, nIp, ent.GetElem() );
        TYPE JeddyAtIp =  shpFncAtIp * Jeddy;

        if (isaxi_)
          {
            Vector<Double> coordAtIP;
            elem->GetShFncAtIp(shpFncAtIp, nIp, ent.GetElem());
            coordAtIP = CornerCoords * shpFncAtIp;
            factor = 2 * PI * coordAtIP[0];
          }

        Double jacDet = elem->
          CalcJacobianDetAtIp(nIp, CornerCoords, ent.GetElem() );
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

  // explicit teplate instantiation for gcc compiler
#ifdef __GNUC__
  template class  MagLorentzForceOp<Double>;
  template class  MagLorentzForceOp<Complex>;
#endif

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

    curlFieldOp_ = new CurlNodeOp<Double>(ptGrid, ptPDE, eqnMap, sol);
    curlFieldOp_->Set2DType(isaxi);

    solType_ = MAG_FORCE_VWP;
    sign_    = 1.0;
  }

  MagForceOp::~MagForceOp()
  {

    if (curlFieldOp_)
      delete curlFieldOp_;
  }


  void MagForceOp::ComputeField(Vector<Double> & Field,
                                const EntityIterator& ent,
                                const Vector<Double> & lCoord)
  {

    curlFieldOp_->CalcElemCurlNode(Field, ent,lCoord);
  }

  Double MagForceOp::GetMatVal(RegionIdType actRegion)
  {

    Double reluctivity;
    materials_[actRegion]->GetScalar(reluctivity,MAG_RELUCTIVITY,Global::REAL);

    return reluctivity;
  }


} // end of namespace CoupledField
