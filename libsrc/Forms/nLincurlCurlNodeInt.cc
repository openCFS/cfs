#include <iostream>
#include <fstream>

#include "Utils/SmoothSpline.hh"
#include "nLincurlCurlNodeInt.hh"

namespace CoupledField
{

  nLinCurlCurlNode2DInt::nLinCurlCurlNode2DInt(BaseFE * aptelem, Double aVal, Boolean axi)
    : BaseForm(aptelem),startmatVal_ (aVal)
  {
#ifdef TRACE
    (*trace) << "entering nLinCurlCurlNode2DInt::nLinCurlCurlNode2DInt" << std::endl;
#endif

    isaxi_ = axi;
  }


  nLinCurlCurlNode2DInt::nLinCurlCurlNode2DInt(ApproxData *nlinFnc, Double aVal, Boolean axi)
    : BaseForm(),startmatVal_ (aVal)
  {
#ifdef TRACE
    (*trace) << "entering nLinCurlCurlNode2DInt::nLinCurlCurlNode2DInt" << std::endl;
#endif

    isaxi_ = axi;

    //set pointer to nonlinear function
    nlinFnc_ = nlinFnc;
  }

 
  nLinCurlCurlNode2DInt::~nLinCurlCurlNode2DInt()
  {
#ifdef TRACE
    (*trace) << "entering nLinCurlCurlNode2DInt::~nLinCurlCurlNode2DInt" << std::endl;
#endif
  }


  void nLinCurlCurlNode2DInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
#ifdef TRACE
    (*trace) << "entering nLinCurlCurlNode2DInt::CalcElementMatrix" << std::endl;
#endif
  
    const Integer nrIntPts= ptelem->GetNumIntPoints();
    const Integer nrNodes = ptelem->GetNumNodes();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> xiDx;
    Matrix<Double> xiDxTransp;
    Matrix<Double> partElemMat;
    Matrix<Double> partElemMatAxi;
    std::vector<Double> ShpFncAtIp;
    std::vector<Double> CoordAtIP;
    std::vector<Double> drAtIp;

    Double reluctivity;

    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrNodes); elemMat.Init(0);
    
    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	jacDet = 0;
	
	ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

	if (isaxi_)
	  {
	    ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	    CoordAtIP = ptCoord * ShpFncAtIp;
            for (Integer i=0; i<nrNodes; i++)
		xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
            
            jacDet *= 2 * PI * CoordAtIP[0];
	}
  
	xiDx.Transpose(xiDxTransp);
	partElemMat = xiDx * xiDxTransp;

	//compute value for nonlinear reluctivity
	Vector<Double> B(2);
	Integer dim = 2;
	for( Integer i=0; i<dim; i++ )
	  for( Integer j=0; j<nrNodes; j++ )
	    B[i] += xiDx[j][i] * magPot_[j];

	Double Babs = B.NormL2();

	if (Babs ==0) 
	  reluctivity = startmatVal_;
	else
	  reluctivity = nlinFnc_->EvaluateFuncInv(Babs) / Babs;

	partElemMat *= intWeights[actIntPt-1] * jacDet * reluctivity;
	elemMat += partElemMat;
      }
  

#ifdef TRACE
    (*trace) << "leaving nLinCurlCurlNode2DInt::CalcElemMatrix" << std::endl;
#endif
  }



  void nLinCurlCurlNode2DInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
#ifdef TRACE
    (*trace) << "entering nLinCurlCurlNode2DInt::Print" << std::endl;
#endif
    (*out)<< "CurlCurlNode2D stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
