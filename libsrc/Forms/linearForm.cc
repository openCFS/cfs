#include <iostream>
#include <fstream>
#include <math.h>

#include "linearForm.hh"


namespace CoupledField
{


  LinearForm::LinearForm(BaseFE * aptelem) : BaseForm(aptelem)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::LinearForm" << std::endl;
#endif
  }


  
  LinearForm ::~LinearForm()
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::~LinearForm" << std::endl;
#endif
  }



  void LinearForm::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & Result)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::CalcElemVector" << std::endl;
#endif
  }
  


  LinearEdgeInt::LinearEdgeInt(BaseFE * aptelem, Double aVal, Integer aDirection,
			       std::vector<Double> * aCoilMidPt) 
    : LinearForm(aptelem), val_(aVal), direction_(aDirection)
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::LinearEdgeInt" << std::endl;
#endif
    
    if(aCoilMidPt)
      coilMidPt_ = new std::vector<Double>(*aCoilMidPt);
  }


  
  LinearEdgeInt ::~LinearEdgeInt()
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::~LinearEdgeInt" << std::endl;
#endif
  }



  void LinearEdgeInt::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::CalcElemVector" << std::endl;
#endif

    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Integer nrEdges  = ptelem->GetNumEdges();
    const Integer dim      = ptelem->GetDim();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Integer xDir = 0;
    const Integer yDir = 1;
    const Integer zDir = 2;
    
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> shapeEdge;
    std::vector<Double> partElemVec;
    partElemVec.resize(nrEdges);


    // set vector to desired size and set all elements to zero    
    elemVec.resize(nrEdges); 
    for (Integer i=0; i<nrEdges; i++)
      elemVec[i] = 0;
    
    std::vector<Double> currentVec(dim, 0);
    // is needed, if coil is curved
    std::vector<Double> globCoord;
    Double len;
    
    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {    
	// define direction of current
	switch(direction_)
	  {	
	  case 1:
	  case 2:
	  case 3:
	    // if current direction goes parallel to a coordinate axis
	    currentVec[direction_-1] = val_;
	    break;

	  case 4: // current in xy-plane ==> z=0
	    ptelem->GetGlobalEdgeIndicesAtIP(globCoord, actIntPt, ptCoord);
	    currentVec[xDir] = -(globCoord[yDir] - (*coilMidPt_)[yDir]);
	    currentVec[yDir] = globCoord[xDir] - (*coilMidPt_)[xDir];
	    currentVec[zDir] = 0;

	    len = L2Norm(currentVec);
	    currentVec *= val_ / len;
	    break;

	  case 5: // current in yz-plane ==> x=0
	    ptelem->GetGlobalEdgeIndicesAtIP(globCoord, actIntPt, ptCoord);
	    currentVec[xDir] = 0;
	    currentVec[yDir] = -(globCoord[zDir] - (*coilMidPt_)[zDir]);
	    currentVec[zDir] = globCoord[yDir] - (*coilMidPt_)[yDir];

	    len = L2Norm(currentVec);
	    currentVec *= val_ / len;
	    break;

	  case 6: // current in xz-plane ==> y=0
	    ptelem->GetGlobalEdgeIndicesAtIP(globCoord, actIntPt, ptCoord);
	    currentVec[xDir] = globCoord[zDir] - (*coilMidPt_)[zDir];
	    currentVec[yDir] = 0;
	    currentVec[zDir] = -(globCoord[xDir] - (*coilMidPt_)[xDir]);

	    len = L2Norm(currentVec);
	    currentVec *= val_ / len;

	    break;
	    
	  default:
	    std::string errMsg = "Selected current direction with number ";
	    errMsg += direction_ + " not supported!";
	    Error(errMsg.c_str(),__FILE__,__LINE__);
	  }
	
	ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord);

	partElemVec = shapeEdge * currentVec;

	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

	for(Integer i=0; i<partElemVec.size(); i++)
	  elemVec[i] += partElemVec[i] * intWeights[actIntPt-1] * jacDet;	
      }
  

#ifdef DEBUG 
	(*debug) << "CalcElemVector:  "  << std::endl
		 << partElemVec << std::endl
		 << "\n jacDet " << jacDet << std::endl;
#endif


#ifdef TRACE
    (*trace) << "leaving linearEdgeInt::CalcElementVector" << std::endl;
#endif

  }




  // ==================================================================
  // nLinMech
  // ==================================================================


  nLinMech_linFormInt::nLinMech_linFormInt(BaseFE * aptelem, MaterialData & mat) 
    : LinearForm(aptelem), matData_(mat)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech_linFormInt::nLinMech_linFormInt" << std::endl;
#endif
  }


  
  nLinMech_linFormInt ::~nLinMech_linFormInt()
  {
#ifdef TRACE
    (*trace) << "entering nLinMech_linFormInt::~nLinMech_linFormInt" << std::endl;
#endif
  }



  void nLinMech_linFormInt::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech_linFormInt::CalcElemVector" << std::endl;
#endif

    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer nrDofs   = getNrDofs();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    std::vector<Double> piolaStressVec;    
    std::vector<Double> partElemVec;

    Matrix<Double> linBMat; 
    Matrix<Double> nonLinBMat; 
    Matrix<Double> dMat;
    Matrix<Double> transpSumB;    // we need transposed of the b-matrices

    // This, as friend defined bilinearform holds the necessary differential operators
    nLinMech3dInt_PiolaStress piolaStressBiform(ptelem, matData_);

    if (!elemDisp_.size_row() || !elemDisp_.size_col()) 
      Error("Undefined displacements! ",__FILE__,__LINE__);

    // set vector to desired size and set all elements to zero    
    partElemVec.resize(nrNodes * nrDofs, 0);
    elemVec.resize(nrNodes*nrDofs, 0);
    


    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {    
	piolaStressBiform.setActElemDispl(elemDisp_);
	piolaStressBiform.calcPiolaStressVec(piolaStressVec, actIntPt, ptCoord);
	piolaStressBiform.calcNonLinBMat(nonLinBMat, actIntPt, ptCoord);
	piolaStressBiform.calcLinBMat(linBMat, actIntPt, ptCoord);

	nonLinBMat += linBMat;

	nonLinBMat.Transpose(transpSumB);
	
	partElemVec = transpSumB * piolaStressVec;

	Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
	partElemVec *= jacDet * intWeights[actIntPt-1];
	
	elemVec +=  partElemVec;

#ifdef DEBUG 
// 	(*debug) << "LINFORM: piolaStressVec: "<< std::endl << partElemVec << std::endl;
// 	(*debug) << "LINFORM: linBMAt: "<< std::endl << linBMat << std::endl;
// 	(*debug) << "LINFORM: nonLinBMAt: "<< std::endl << nonLinBMat << std::endl;
// 	(*debug) << "LINFORM: transpSumB: "<< std::endl << transpSumB << std::endl;
#endif

      }
  

#ifdef DEBUG 
	(*debug) << "CalcElemVector:  "  << std::endl
		 << partElemVec << std::endl;
#endif



#ifdef TRACE
    (*trace) << "leaving nLinMech_linFormInt::CalcElementVector" << std::endl;
#endif

  }


  } // end of namespace
