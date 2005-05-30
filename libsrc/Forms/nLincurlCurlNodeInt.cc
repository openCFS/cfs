#include <iostream>
#include <fstream>

#include "Utils/SmoothSpline.hh"
#include "nLincurlCurlNodeInt.hh"

namespace CoupledField
{

  nLinCurlCurlNode2DInt::nLinCurlCurlNode2DInt(BaseFE * aptelem, Double aVal, Boolean axi)
    : BaseForm(aptelem),startmatVal_ (aVal)
  {
    ENTER_FCN( "nLinCurlCurlNode2DInt::nLinCurlCurlNode2DInt");

    isaxi_      = axi;
    nonLinType_ = NEWTON;
  }


  nLinCurlCurlNode2DInt::nLinCurlCurlNode2DInt(ApproxData *nlinFnc, Double aVal, Boolean axi)
    : BaseForm(),startmatVal_ (aVal)
  {
    ENTER_FCN( "nLinCurlCurlNode2DInt::nLinCurlCurlNode2DInt");

    isaxi_      = axi;
    nonLinType_ = NEWTON;

    //set pointer to nonlinear function
    nlinFnc_ = nlinFnc;
  }

 
  nLinCurlCurlNode2DInt::~nLinCurlCurlNode2DInt()
  {
    ENTER_FCN( "nLinCurlCurlNode2DInt::~nLinCurlCurlNode2DInt");
  }


  void nLinCurlCurlNode2DInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
    ENTER_FCN( "nLinCurlCurlNode2DInt::CalcElementMatrix");
  
    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Matrix<Double> partElemMatAxi;
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIP;
    Vector<Double> drAtIp;

    Double reluctivity, derivReluctivity, Hfield, dHfield;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes); elemMat.Init(0);

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
        jacDet = 0;
        
        ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

        if (isaxi_)
          {
            ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
            CoordAtIP = ptCoord * ShpFncAtIp;
            for (UInt i=0; i<nrNodes; i++)
              xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
            
            jacDet *= 2 * PI * CoordAtIP[0];
          }
  
        xiDx.Transpose(xiDxTransp);
        partElemMat = xiDx * xiDxTransp;

        //compute value for nonlinear reluctivity
        Vector<Double> B(2);
        UInt dim = 2;
        for( UInt i=0; i<dim; i++ )
          for( UInt j=0; j<nrNodes; j++ )
            B[i] += xiDx[j][i] * magPot_[j];

        Double Babs = B.NormL2();

        if (Babs ==0) 
          reluctivity = startmatVal_;
        else {
          Hfield      = nlinFnc_->EvaluateFuncInv(Babs);
          reluctivity = Hfield / Babs;
        }

        partElemMat *= reluctivity;
        
        if (nonLinType_ == NEWTON) {
          if (Babs ==0) 
            derivReluctivity = 0;
          else {          
            //Newton method
            Vector<Double> eB(2); eB = B * 1/Babs;
            dHfield = nlinFnc_->EvaluatePrimeInv(Babs);
            derivReluctivity = (dHfield*Babs - Hfield) / (Babs*Babs);
            for (UInt p=0;  p<nrNodes; p++)
              for (UInt q=0; q<nrNodes; q++) {               
                partElemMat[p][q] +=  derivReluctivity * 
                  (eB[0]*eB[0]*xiDx[p][1]*xiDx[q][1] +
                   eB[1]*eB[1]*xiDx[p][0]*xiDx[q][0] -
                   eB[0]*eB[1]*xiDx[p][1]*xiDx[q][0] -
                   eB[1]*eB[0]*xiDx[p][0]*xiDx[q][1] );
              }
          }
        }
        
    
        partElemMat *= intWeights[actIntPt-1] * jacDet;
        elemMat += partElemMat;
      }
  }


  void nLinCurlCurlNode2DInt::SetNonLinMethod(std::string atype)
  {
    ENTER_FCN( "nLinCurlCurlNode2DInt::SetNonLinMethod");
    
    if (atype == "fixPoint")
      nonLinType_ = FIXEDPOINT;
    
  }


  void nLinCurlCurlNode2DInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
    ENTER_FCN( "nLinCurlCurlNode2DInt::Print");

    (*out)<< "CurlCurlNode2D stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
