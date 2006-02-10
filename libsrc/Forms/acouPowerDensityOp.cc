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
  
  AcouPowerDensityOp::AcouPowerDensityOp(Grid * ptGrid, 
                                         StdPDE * ptPDE,
                                         NodeEQN * ptEQN,
                                         Boolean isaxi)
    : BaseOperator(ptGrid, ptPDE, ptEQN, isaxi)
  {
    ENTER_FCN( "AcouPowerDensityOp::AcouPowerDensityOp" );
    
    isaxi_ = isaxi;

    //this->sol_ = & sol;
    //this->solDeriv1_ = & solDeriv1;
  }


  AcouPowerDensityOp::~AcouPowerDensityOp()
  {
    ENTER_FCN( "AcouPowerDensityOp::~AcouPowerDensityOp" );
  }

  void AcouPowerDensityOp::CalcElemPD(Vector<Double> & elemPD,
                                          const Elem * ptElement,
                                          const Double density)
  {
    ENTER_FCN( "AcouPowerDensityOp::CalcElemEnergy" );

    const UInt nrIntPts = ptElement->ptElem->GetNumIntPoints();
    const UInt nrNodes  = ptElement->ptElem->GetNumNodes();

    // initialize output vector
    elemPD.Resize(nrNodes);
    elemPD.Init(0.0);

    StdVector<UInt> connect = ptElement->connect;
    Matrix<Double> ptCoord;
    ptPDE_->GetElemCoords(connect, ptCoord);


    Vector<Double> elemSol, elemSolDeriv1;
    ptPDE_->GetSolVecOfElement(elemSol,connect);
    ptPDE_->GetDerivSolVecOfElement(elemSolDeriv1,connect);

    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx, xiDxTransp;

    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIp;

    Vector<Double> solGradAtIp, solDeriv1GradAtIp;
    Double solDeriv1AtIp;

    Double N1;
    Vector<Double> N2;
    Double jacDet;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = 0;
      ptElement->ptElem->GetShFncAtIp(ShpFncAtIp,actIntPt);
      ptElement->ptElem->GetGlobDerivShFncAtIp(xiDx,actIntPt,ptCoord,jacDet);
        
      if (isaxi_) {
        CoordAtIp = ptCoord * ShpFncAtIp;
        for (UInt i=0; i<nrNodes; i++)
          xiDx[i][0] += ShpFncAtIp[i] / CoordAtIp[0];
          
        jacDet *= 2 * PI * CoordAtIp[0];
      }
        
      xiDx.Transpose(xiDxTransp);
        
      // compute gradient of solution + 1st derivative at integration point
      solGradAtIp       = xiDxTransp * elemSol;
      solDeriv1GradAtIp = xiDxTransp * elemSolDeriv1;
        
      // get 1st derivartive of solution at integration point
      solDeriv1AtIp     = ShpFncAtIp * elemSolDeriv1;

        
      N1 = 0;
      for (UInt j=0; j<xiDx.GetSizeCol(); j++)
        N1 += solGradAtIp[j] * solDeriv1GradAtIp[j];

      N2 = solGradAtIp * solDeriv1AtIp;

      for (UInt i=0; i< nrNodes; i++) {

        elemPD[i] += ShpFncAtIp[i] * N1;

        for (UInt j=0; j<xiDx.GetSizeCol(); j++)
          elemPD[i] -= xiDx[i][j] * N2[j];
      }
      elemPD *= jacDet * density;  
    }
    
  }

  
  void AcouPowerDensityOp::CalcSubdomPD(Vector<Double> & elemPD,
                                        const StdVector<RegionIdType> & SD,
                                        const Double density)
  {
    ENTER_FCN( "AcouPowerDensityOp::CalcSDEnergy" );

    (*error) << "AcouPowerDensityOp::CalcSDElecField: Not working yet"; 
    Error( __FILE__, __LINE__ );

  }
} // namespace CoupledField
