#include <string>

#include "Elements/basefe.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"
#include "Forms/linStrainOp.hh"

#include <PDE/StdPDE.hh>

namespace CoupledField
{

  template<class TYPE>
  LinStrainOp<TYPE>::LinStrainOp(Grid * ptGrid, 
                                 StdPDE * ptPDE,
                                 NodeEQN * ptEQN,
                                 NodeStoreSol<TYPE> & displacement,
                                 const SolutionType solType,
                                 Boolean isaxi)
    : BaseOperator(ptGrid, ptPDE, ptEQN, isaxi)
  {
    ENTER_FCN( "LinStrainOp::LinStrainOp" );  
    this->displacement_ = &displacement;
    solType_ = solType;
    ptPDE_ = ptPDE;
 
  }

  template<class TYPE>
  LinStrainOp<TYPE>::~LinStrainOp()
  {
    ENTER_FCN( "LinStrainOp::~LinStrainOp" );

  }

  template <class TYPE>
  void LinStrainOp<TYPE>::CalcElemLinearStrain(CFSVector & strain,
                                               const Elem * ptElement,
                                               Matrix<Double> & lCoord)
  {
    ENTER_FCN( "LinStrainOp::CalcElemLinearStrain" );
  
    UInt dim;
    Vector<TYPE> potEntry(1);
    dim = ptElement->ptElem->GetDim();
 
    Vector<TYPE> & helpElemField = dynamic_cast<Vector<TYPE>&> (strain);   

    helpElemField.Resize(dim);
    helpElemField.Init(0.0);
   
    UInt nShFnc = 0;
    UInt ip=1;
    nShFnc = ptElement->ptElem->GetNumNodes();

    StdVector<UInt> connecth = ptElement->connect;
  
    Vector<TYPE> actDispl;
    displacement_->GetElemSolution(actDispl,connecth);
    
    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ptElement->ptElem, ip, lCoord);

    Matrix<TYPE> linBMatTemp = linBMat;
    helpElemField=linBMatTemp * actDispl;
    
  }

  template<class TYPE>
  void LinStrainOp<TYPE>::calcBMat(Matrix<Double> &bMat, BaseFE * ptelem, 
                                   UInt ip, Matrix<Double> & ptCoord){
    ENTER_FCN("LinStrainOP::calcBMat");

    const UInt nrNodes  = ptelem->GetNumNodes();
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = ptPDE_->getPDE_dofspernode();  

    UInt actDim, actNode, j, k;
    if (spaceDim==2)
      if (isaxi_)      
        bMat.Resize(4, nrNodes * nrDofs);
      else // plane strain
        bMat.Resize(3, nrNodes * nrDofs);
    else // 3d
      bMat.Resize(6, nrNodes*nrDofs);
    
    // local shape functions derived after global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;

  //   if (isSetIntPoint_) 
//       ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord);
//     else
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);


    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < nrNodes; actNode++)
        bMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

    switch(spaceDim)
      {
      case 2:
        j = 1;
        k = 0;
        
        for (actNode = 0; actNode < nrNodes; actNode++)
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
              ptelem->GetShFncAtIp(ShpFncAtIp,ip);

            CoordAtIP = ptCoord * ShpFncAtIp;

            for (actNode = 0; actNode < nrNodes; actNode++)          
              bMat[idxtheta-1][actNode * spaceDim] =
                ShpFncAtIp[actNode] / CoordAtIP[0];
          }

        break;


      case 3:
        UInt actDim=spaceDim;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][1];
          }

        actDim++;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][0];
          }

        actDim++;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][1];
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][0];
          }
        break;
      }

    //    isSetIntPoint_ = FALSE;

  
  } // end calcBMat

  template<class TYPE>
  void LinStrainOp<TYPE>::CalcSDLinearStrain(CFSVector & strain,
                                              const StdVector<RegionIdType> & SD, 
                                              Matrix<Double> & lCoords)
  {
    ENTER_FCN( "LinStrainOp::CalcSDLinearStrain" );
  
     Error( "LinStrainOp::CalcSDLinearStrain: Not working yet", __FILE__, __LINE__);

  }

  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate LinStrainOp<Double>
#pragma instantiate LinStrainOp<Complex>
#endif


} // namespace CoupledField
