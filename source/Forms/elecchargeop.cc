#include "elecchargeop.hh"
#include "Domain/grid.hh"
#include "Elements/basefe.hh"
#include "PDE/StdPDE.hh"
#include "Domain/elem.hh"

namespace CoupledField
{

  template<class TYPE>
  ElecChargeOp<TYPE>::ElecChargeOp(Grid * ptGrid,
                                   StdPDE * ptPDE,
                                   shared_ptr<EqnMap> eqnMap,
                                   shared_ptr<ResultInfo> result,
                                   bool isaxi, bool coordUpdate )
    : BaseOperator(ptGrid, ptPDE, eqnMap, isaxi, coordUpdate ) {
    ENTER_FCN( "ElecChargeOp::ElecChargeOp" );
    
    result_ = result;
    
  }
  
  template <class TYPE>
  ElecChargeOp<TYPE>::~ElecChargeOp() {
    ENTER_FCN( "ElecChargeOp::~ElecChargeOp" );
    
  }
  
  template <class TYPE>
  void ElecChargeOp<TYPE>::CalcElemCharge(TYPE & charge,
                                          const EntityIterator& it,
                                          const Vector<Double> & lCoord,
                                          const TYPE & eNormalFluxDensity) {

    ENTER_IFCN( "ElecChargeOp::CalcElemCharge" );
    
    Double jacDet;
    TYPE chargeAux;
    Vector<Double> shFnc, globCoord;
    Matrix<Double> coordMat;
    BaseFE * ptElemFE = it.GetSurfElem()->ptElem;
    UInt nrIntPts, nrNodes;
    
  
    charge = 0;
    chargeAux = 0;
    jacDet =0;
    ptElemFE->SetAnsatzFct( result_->fctType );
    nrIntPts= ptElemFE->GetNumIntPoints();
    nrNodes = ptElemFE->GetNumNodes();
    const Vector<Double> & intWeights = ptElemFE->GetIntWeights();   
  
    ptGrid_->GetElemNodesCoord( coordMat, it.GetSurfElem()->connect, 
                                coordUpdate_ );
  
    // loop over all integration points
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        ptElemFE->GetShFncAtIp(shFnc, actIntPt, it.GetElem() );
        jacDet = ptElemFE->CalcJacobianDetAtIp(actIntPt, coordMat, 
                                               it.GetSurfElem() );
        chargeAux =  eNormalFluxDensity * jacDet * intWeights[actIntPt-1];
        if (isaxi_)
          {
            globCoord = coordMat * shFnc;
            chargeAux *= 2 * PI * globCoord[0];
          }
        charge += chargeAux;
      }

  }


 


  
  template <class TYPE>
  void ElecChargeOp<TYPE>::
  CalcElemCharges(Vector<TYPE> & charges,
                  const shared_ptr<SurfElemList>  surfElems,
                  const Vector<Double> & lCoord,
                  const Vector<TYPE> & eNormalFluxDensity)
  {
    ENTER_FCN( "ElecChargeOp::CalcElemCharges" );
    
    Warning( "Not used anymore", __FILE__, __LINE__ );
  //   Double jacDet; 
//     TYPE charge,helpNormalFluxDensity;
//     Vector<Double> shFnc, globCoord;
//     Matrix<Double> coordMat;
//     BaseFE * ptElem;
//     UInt nrIntPts, nrNodes;
//     charges.Resize(surfElems.GetSize());
//     charges.Init();
    
  
//     // loop over all surface elements
//     for (UInt iElem=0; iElem<surfElems.GetSize(); iElem++)
//       {
//         ptElem = surfElems[iElem]->ptElem;
//         charge = 0;
//         jacDet =0;
//         nrIntPts= ptElem->GetNumIntPoints();
//         nrNodes = ptElem->GetNumNodes();
//         const Vector<Double> & intWeights = ptElem->GetIntWeights();   
      
//         ptGrid_->GetElemNodesCoord(coordMat, surfElems[iElem]->connect, 
//                                    coordUpdate_ );
//         ptElem->GetShFnc(shFnc, lCoord, surfElems[iElem]);
//         // if (charges.IsComplex()){
//           for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
//             {
//               jacDet = ptElem->CalcJacobianDet(lCoord, coordMat, 
//                                                surfElems[iElem]);
//               charge = shFnc[actIntPt] *intWeights[actIntPt] * jacDet;
//               eNormalFluxDensity.GetEntry(iElem,helpNormalFluxDensity);
//               charge *=  helpNormalFluxDensity;
//               if (isaxi_)
//                 {
//                   globCoord = coordMat * shFnc;
//                   charge *= 2 * PI * globCoord[0];
//                 }
          
//               //          charges[iElem] += charge;
//               charges.AddEntry(iElem,charge);      

//             }
//         }
//         else {
//           for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
//             {
//               jacDet = ptElem->CalcJacobianDet(lCoord, coordMat);
//               charge = shFnc[actIntPt] *intWeights[actIntPt];
//               eNormalFluxDensity.GetEntry(iElem,helpNormalFluxDensityD);
//               charge *= jacDet * helpNormalFluxDensityD;
//               if (isaxi_)
//                 {
//                   globCoord = coordMat * shFnc;
//                   charge *= 2 * PI * globCoord[0];
//                 }
          
//               //          charges[iElem] += charge;
//               charges.AddEntry(iElem,charge);     

//             }


    
  }
  
  // Explicit template instantiation
#ifdef __GNUC__
  template class ElecChargeOp<Double>;
  template class ElecChargeOp<Complex>;
#endif
} // end of namespace
