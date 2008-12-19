// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Forms/acouPowerDensityOp.hh"

#include <string>

#include "Elements/basefe.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"

#include <PDE/StdPDE.hh>

namespace CoupledField {

  template<class TYPE>  
  AcouPowerDensityOp<TYPE>::AcouPowerDensityOp(Grid * ptGrid, 
                                               StdPDE * ptPDE,
                                               shared_ptr<EqnMap> eqnMap,
                                               bool isaxi)
    : BaseOperator(ptGrid, ptPDE, eqnMap, isaxi)
  {
    
    isaxi_ = isaxi;
    //Warning( "Only working with Lagrange Functions", __FILE__, __LINE__ );
  }

  template<class TYPE>
  AcouPowerDensityOp<TYPE>::~AcouPowerDensityOp()
  {
  }


  template<class TYPE>
  void AcouPowerDensityOp<TYPE>::CalcElemPD(Vector<Double> & elemPD,
                                            const EntityIterator& it,
                                            const Double density)
  {

    Elem const * ptElement = it.GetElem();
    const UInt nrIntPts = ptElement->ptElem->GetNumIntPoints();
    const UInt nrNodes  = ptElement->ptElem->GetNumNodes();
    const Vector<Double> &intWeights = ptElement->ptElem->GetIntWeights();  

    // initialize output vector
    elemPD.Resize(nrNodes);
    elemPD.Init(0.0);

    StdVector<UInt>  connect = ptElement->connect;
    Matrix<Double> ptCoord;
    ptGrid_->GetElemNodesCoord(ptCoord, connect, coordUpdate_ );

    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx, xiDxTransp;

    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIp;

    Vector<TYPE> elemSol, elemSolDeriv1;
    shared_ptr<ResultInfo> res = (ptPDE_->GetResultInfos())[0];
    ptPDE_->GetSolVecOfElement(elemSol,it,res);
    ptPDE_->GetDerivSolVecOfElement(elemSolDeriv1,it,res);

//     std::cout << "Number Nodes:\n" << nrNodes << std::endl;
//     std::cout << "Number IntPoints:\n" << nrIntPts << std::endl;
//     std::cout << "elemSol:\n" << elemSol << std::endl;
//     std::cout << "elemSolDeriv1:\n" << elemSolDeriv1 << std::endl;

    Vector<TYPE> solGradAtIp, solDeriv1GradAtIp;
    TYPE solDeriv1AtIp;

    Double N1;
    Vector<Double> N2;
    Double jacDet;
    Double sumJacDet=0.0;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = 0;
      ptElement->ptElem->GetShFncAtIp(ShpFncAtIp,actIntPt,ptElement);
      ptElement->ptElem->GetGlobDerivShFncAtIp(xiDx,actIntPt,
                                               ptCoord,jacDet, ptElement);
        
      if (isaxi_) {
        CoordAtIp = ptCoord * ShpFncAtIp;
        for (UInt i=0; i<nrNodes; i++)
          xiDx[i][0] += ShpFncAtIp[i] / CoordAtIp[0];
          
        jacDet *= 2 * PI * CoordAtIp[0];
      }
        
      xiDx.Transpose(xiDxTransp);
      dimensions_ = xiDx.GetSizeCol();
        
      // compute gradient of solution + 1st derivative at integration point
      solGradAtIp       = xiDxTransp * elemSol;
      solDeriv1GradAtIp = xiDxTransp * elemSolDeriv1;
        
      // get 1st derivartive of solution at integration point
      solDeriv1AtIp     = ShpFncAtIp * elemSolDeriv1;

        
      N1 = ComputeN1( solGradAtIp, solDeriv1GradAtIp );
      N2 = ComputeN2( solGradAtIp, solDeriv1AtIp );

      sumJacDet += jacDet;
      for (UInt i=0; i< nrNodes; i++) {

        elemPD[i] += ShpFncAtIp[i] * N1;

        for (UInt j=0; j<xiDx.GetSizeCol(); j++)
          elemPD[i] -= xiDx[i][j] * N2[j];
      }
      elemPD *= jacDet * intWeights[actIntPt-1];
    }
//     std::cout << "sum of jacDet=" << sumJacDet << std::endl;

    elemPD *=  density;

   // (*debug) << "Coupling vector acoustic-heat conduction is: " 
   //          << std::endl << elemPD << std::endl << std::endl;

  }

  template<class TYPE>  
  Double AcouPowerDensityOp<TYPE>::ComputeN1( Vector<TYPE> solGrtIp,
                                              Vector<TYPE> solD1GrAtIp )
  {

    (*error) << "AcouPowerDensityOp::ComputeN1 "
             << "only defined for TYPE=Complex/Double "; 
    Error( __FILE__, __LINE__ );
    return TYPE();
  }

  template<>
  Double AcouPowerDensityOp<Double>::ComputeN1( Vector<Double> solGrAtIp,
                                                Vector<Double> solD1GrAtIp )
  {

    Double N1 = 0;
    for (UInt k=0; k<dimensions_; k++)
      N1 += solGrAtIp[k] * solD1GrAtIp[k];
    return N1;
  }

  template<>
  Double AcouPowerDensityOp<Complex>::ComputeN1( Vector<Complex> solGrAtIp,
                                                 Vector<Complex> solD1GrAtIp )
  {

    Double N1 = 0;
    for (UInt k=0; k<dimensions_; k++)
      N1 += 0.5 * std::real( std::conj(solGrAtIp[k]) * solD1GrAtIp[k] );
    return N1;
  }


  template<class TYPE>  
  Vector<Double> AcouPowerDensityOp<TYPE>::ComputeN2( Vector<TYPE> solGrAtIp,
                                                      TYPE solD1AtIp )
  {

    (*error) << "AcouPowerDensityOp::ComputeN2 "
             << "only defined for TYPE=Complex/Double "; 
    Error( __FILE__, __LINE__ );
    return TYPE();
  }

  template<>
  Vector<Double> AcouPowerDensityOp<Double>::ComputeN2( Vector<Double> solGrAtIp,
                                                        Double solD1AtIp )
  {

    Vector<Double> N2;
    N2 = solGrAtIp * solD1AtIp;
    return N2;
  }

  template<>
  Vector<Double> AcouPowerDensityOp<Complex>::ComputeN2( Vector<Complex> solGrAtIp,
                                                         Complex solD1AtIp )
  {

    Vector<Double> N2;
    N2.Resize( solGrAtIp.GetSize() );

    for (UInt k=0; k<solGrAtIp.GetSize(); k++) {
      N2[k] = 0.5 * std::real( std::conj(solGrAtIp[k]) * solD1AtIp );
    }   
    return N2;
  }

//   void AcouPowerDensityOp::CalcSubdomPD(Vector<Double> & elemPD,
//                                         const StdVector<RegionIdType> & SD,
//                                         const Double density)
//   {
//     (*error) << "AcouPowerDensityOp::CalcSDElecField: Not implemented yet"; 
//     Error( __FILE__, __LINE__ );
//   }



  // explicit teplate instantiation for gcc compiler
#ifdef __GNUC__
template class AcouPowerDensityOp<Double>;
template class AcouPowerDensityOp<Complex>;
#endif

  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate AcouPowerDensityOp<Double>
#pragma instantiate AcouPowerDensityOp<Complex>
#endif


} // namespace CoupledField
