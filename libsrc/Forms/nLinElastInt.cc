#include <iostream>
#include <fstream>

#include "nLinElastInt.hh"

namespace CoupledField
{


  // returns nonlinear B - matrix (first part) for BDB
  void nLinElastInt::calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
#ifdef TRACE
    (*trace) << "entering nLinElastInt::calcBMat " << std::endl;
#endif
    
    if (!elemDisp_.size_row() || !elemDisp_.size_col()) 
      Error("Undefined displacements! ",__FILE__,__LINE__);

    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer nrDofs   = getNrDofs();    

    Matrix<Double> linBMat;    
    linElastInt::calcBMat( linBMat, ip, ptCoord);

    
    bMat.Resize(getDimD(), nrNodes * nrDofs);

    
    // local shape functions derived after global coords 
    // (format: nrNodes x nrDofs)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);


    // calculate derivates of global displacements at element nodes
    // elemDisp_ holds displacement of the actual element (dimension: nrDofs x nrNodes)
    // xiDx holds derivatives of shape functions after global coords (dimension: nrNodes x nrDofs)
    // displDeriv dimension: nrDofs x nrDofs
    Matrix<Double>  displDeriv;  
    displDeriv = elemDisp_ * xiDx;

    // we need transposed of displacement derivatives
    Matrix<Double> displDerivTransp;
    displDeriv.Transpose(displDerivTransp);  


    Matrix<Double> bMatOneNode;
    bMatOneNode.Resize(linBMat.size_row(), nrDofs);
 
    for(int actNode=0; actNode < nrNodes; actNode++)
      {
	linBMat.GetSubMatrix(bMatOneNode, 0, actNode*nrDofs);
	
	bMatOneNode *= displDerivTransp;
	
	bMat.SetSubMatrix(bMatOneNode, 0, actNode*nrDofs);
      }

#ifdef DEBUG
    *debug << std::endl << "element displacement: " << std::endl 
	   << elemDisp_ << std::endl
	   << "ip-nr = " << ip << std::endl
	   << "Nonlinear BMat: " << std::endl << bMat << std::endl;
#endif
  }
  


  // calculates the D-matrix of a 3d-problem 
  void nLinMech3dInt::calcDMat(Matrix<Double> & dMat)
  {
#ifdef TRACE
  (*trace) << "entering mech3DInt::calcDMat " << std::endl;
#endif
    const Integer nrElems3d = getDimD();
    
    Matrix<Double> * matMatrix =  ptMaterial->GetMatrix();
    
    dMat.Resize(nrElems3d);

    for (Integer i=0; i<nrElems3d; i++)
      for (Integer j=0; j<nrElems3d; j++)
	dMat[i][j] = (*matMatrix)[i][j];	
  }





  // ===================================================================================
  // =================== standard con- and destructors (just for tracing) ==============
  // ===================================================================================

  // base class for nonlinear calculation of elasticity
  nLinElastInt::nLinElastInt(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering nLinElastInt::nLinElastInt" << std::endl;
#endif
  }
 

  nLinElastInt::~nLinElastInt()
  {
#ifdef TRACE
    (*trace) << "entering nLinElastInt::~nLinElastInt" << std::endl;
#endif
  }



  // nonlinear calculation of elasticity in 3d
  nLinMech3dInt::nLinMech3dInt(BaseFE * aptelem, MaterialData & matData) 
    : nLinElastInt(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech3dInt::nLinMech3dInt" << std::endl;
#endif
  }
 

  nLinMech3dInt::~nLinMech3dInt()
  {
#ifdef TRACE
    (*trace) << "entering nLinMech3dInt::~nLinMech3dInt" << std::endl;
#endif
  }

} // end namespace CoupledField
