#include <iostream>
#include <fstream>

#include "curlCurlEdgeInt.hh"

namespace CoupledField
{

  CurlCurlEdgeInt::CurlCurlEdgeInt(BaseFE * aptelem, Double aVal)
    : BaseForm(aptelem), reluctivity_ (aVal)
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlEdgeInt::CurlCurlEdgeInt" << std::endl;
#endif
  }


 
  CurlCurlEdgeInt::~CurlCurlEdgeInt()
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlEdgeInt::~CurlCurlEdgeInt" << std::endl;
#endif
  }



  void CurlCurlEdgeInt::CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat)
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlEdgeInt::CalcElementMatrix" << std::endl;
#endif
  
    const Integer nrIntPts= ptelem->GetNumIntPoints();
    const Integer nrEdges = ptelem->GetNumEdges();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    std::vector< Matrix<Double>* > xiDx;
    xiDx.resize(nrEdges);
    for (Integer i=0; i<nrEdges; i++)
      xiDx[i] = new Matrix<Double>;
    
    Matrix<Double> curl;
    Matrix<Double> curlTransp;
    Matrix<Double> partElemMat;
  



    // set matrix to desired size and set all elements to zero
    elemMat.Resize(nrEdges); 
    elemMat.Init();
    

    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	// calc glob derivs of shape functions and jacobian determinante
	ptelem->GetEdgeGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord);

	CalcEdgeCurl(curl, xiDx);
	
 	curl.Transpose(curlTransp);

 	partElemMat = curlTransp * curl;

	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
	
 	partElemMat *= intWeights[actIntPt-1] * jacDet * reluctivity_;

	elemMat += partElemMat;
      }
    

#ifdef DEBUG 
// 	(*debug) << "CurlCurlEdgeInt: ElemMat " << std::endl
// 		 << elemMat << std::endl
// 		 << "\n reluctivity " << reluctivity_ << std::endl;
#endif


#ifdef TRACE
    (*trace) << "leaving CurlCurlEdgeInt::CalcElemMatrix" << std::endl;
#endif
  }


  // calculates the curl, if the global derivates are already given in shapeDeriv
  void CurlCurlEdgeInt::CalcEdgeCurl(Matrix<Double>& curl, 
				     const std::vector< Matrix<Double>* >& shapeDeriv)
  {
    Integer nrEdges = shapeDeriv.size();
    Integer dim = shapeDeriv[0]->GetSizeRow();
    
    curl.Resize(dim, nrEdges);
    
    for (Integer actEdge=0; actEdge < nrEdges; actEdge++)
      for (Integer actDim=0; actDim < dim; actDim++)
	curl[actDim][actEdge] = 
	  (*shapeDeriv[actEdge])[(actDim+2)%dim][(actDim+1)%dim] -
	  (*shapeDeriv[actEdge])[(actDim+1)%dim][(actDim+2)%dim];
    
  }
  

  void CurlCurlEdgeInt::Print(std::ostream * out, const Matrix<Double> Result) const
  {
#ifdef TRACE
    (*trace) << "entering CurlCurlEdgeInt::Print" << std::endl;
#endif
    (*out)<< "CurlCurlEdge stiffness matrix:" << std::endl << Result;
  }

} // end namespace CoupledField
