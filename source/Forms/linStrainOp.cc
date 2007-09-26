// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>

#include "Elements/basefe.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"
#include "Forms/linStrainOp.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

#include <PDE/StdPDE.hh>

namespace CoupledField
{

  template<class TYPE>
  LinStrainOp<TYPE>::LinStrainOp(Grid * ptGrid, 
                                 StdPDE * ptPDE,
                                 shared_ptr<EqnMap> eqnMap,
                                 NodeStoreSol<TYPE> & displacement,
                                 const SolutionType solType,
                                 bool isaxi)
    : BaseOperator(ptGrid, ptPDE, eqnMap, isaxi)
  {
    this->displacement_ = &displacement;
    solType_ = solType;
    ptPDE_ = ptPDE;
 
  }

  template<class TYPE>
  LinStrainOp<TYPE>::~LinStrainOp()
  {

  }

  template <class TYPE>
  void LinStrainOp<TYPE>::CalcElemLinearStrain(CFSVector & strain,
                                               const Elem * ptElement,
                                               Matrix<Double> & lCoord)
  {
  
    UInt dim;
    Vector<TYPE> potEntry(1);
    dim = ptElement->ptElem->GetDim();
 
    Vector<TYPE> & helpElemField = dynamic_cast<Vector<TYPE>&> (strain);   

    helpElemField.Resize(dim);
    helpElemField.Init(0.0);
   
    UInt nShFnc = 0;
    UInt ip=1;
    nShFnc = ptElement->ptElem->GetNumNodes();

  
    Vector<TYPE> actDispl;
    ElemList elemList(domain->GetGrid());
    elemList.SetElement( ptElement );
    EntityIterator it = elemList.GetIterator();
    displacement_->GetElemSolution(actDispl,it);
    
    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ptElement->ptElem, ip, lCoord);

    helpElemField=linBMat * actDispl;
    
  }

  template<class TYPE>
  void LinStrainOp<TYPE>::calcBMat(Matrix<Double> &bMat, BaseFE * ptelem, 
                                   UInt ip, Matrix<Double> & ptCoord){

    UInt numFncs = ptelem->GetNumNodes ();
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = spaceDim;

    UInt actDim, actNode, j, k;
    if (spaceDim==2)
      if (isaxi_)      
        bMat.Resize(4, numFncs * nrDofs);
      else // plane strain
        bMat.Resize(3, numFncs * nrDofs);
    else // 3d
      bMat.Resize(6, numFncs*nrDofs);
    
    bMat.Init();
    // local shape functions derived after global coords
    // (format: numFncs x spaceDim)
    Matrix<Double> xiDx;

  //   if (isSetIntPoint_) 
//       ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord);
//     else
    Warning( "In the following line the NULL has to be replaced (Andreas )",
             __FILE__, __LINE__ );
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, NULL );


    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < numFncs; actNode++)
        bMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

    switch(spaceDim)
      {
      case 2:
        j = 1;
        k = 0;
        
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[spaceDim][actNode * spaceDim + 1] = xiDx[actNode][0];
            bMat[spaceDim][actNode * spaceDim]     = xiDx[actNode][1];
          }

        if (isaxi_)
          {
            UInt idxtheta = 4;
            Vector<Double> ShpFncAtIp;
            Vector<Double> CoordAtIP;

           //  if (isSetIntPoint_) 
//               ptelem->GetShFnc(ShpFncAtIp,intPoint_);
//             else
            Warning( "In the following line the NULL has to be replaced (Andreas )",
                     __FILE__, __LINE__ );
            ptelem->GetShFncAtIp(ShpFncAtIp,ip,NULL);

            CoordAtIP = ptCoord * ShpFncAtIp;

            for (actNode = 0; actNode < numFncs; actNode++)          
              bMat[idxtheta-1][actNode * spaceDim] =
                ShpFncAtIp[actNode] / CoordAtIP[0];
          }

        break;


      case 3:
        UInt actDim=spaceDim;
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][1];
          }

        actDim++;
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][0];
          }

        actDim++;
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][1];
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][0];
          }
        break;
      }

    //    isSetIntPoint_ = false;

  
  } // end calcBMat

  template<class TYPE>
  void LinStrainOp<TYPE>::CalcSDLinearStrain(CFSVector & strain,
                                              const StdVector<RegionIdType> & SD, 
                                              Matrix<Double> & lCoords)
  {
  
     Error( "LinStrainOp::CalcSDLinearStrain: Not working yet", __FILE__, __LINE__);

  }

  // explicit template instantiation for GCC compiler
#ifdef __GNUC__
  // Template instantiation for used vectors
  template class LinStrainOp<Double>;
  template class LinStrainOp<Complex>;
#endif


  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate LinStrainOp<Double>
#pragma instantiate LinStrainOp<Complex>
#endif


} // namespace CoupledField
