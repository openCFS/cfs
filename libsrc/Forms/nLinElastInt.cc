#include <iostream>
#include <fstream>

#include "nLinElastInt.hh"

namespace CoupledField
{

  /// calculates Piola-Kirchoff-stresses T (vector notation)
// T = c . V with c the stiffness matrix and V the Green-Lagrange strain vector
// V = B_lin * u  +  1/2 * B_nonLin * u
// see Habil. M. Kaltenbacher Eq. (3.33)
void nLinMech3dInt_PiolaStress::
calcPiolaStressVec(std::vector<Double>& piolaStressVec, Integer ip, Matrix<Double> & ptCoord)
{
#ifdef TRACE
    (*trace) << "entering mech3DInt_PiolaStress::calcPiolaStressTensor" << std::endl;
#endif
    Matrix<Double> dMat;
    calcMaterialDMat(dMat);
    
  // convert displacement of all elem nodes into one vector: 
  // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
  std::vector<Double> displVec;
  //  elemDisp_.ConvertToVec_RowsFirst(displVec);
  elemDisp_.ConvertToVec_AppendCols(displVec);


  // linear differential operator B_lin
  Matrix<Double> linBMat;    
  calcLinBMat( linBMat, ip, ptCoord);


  // nonlinear differential operator B_nonLin
  Matrix<Double> nLinBMat;    
  calcNonLinBMat( nLinBMat, ip, ptCoord);


  // special construction: direct addition of these two terms does not work, 
  // because of temporary vectors ... :O(    
  std::vector<Double> part1(linBMat * displVec );
  std::vector<Double> part2(nLinBMat * 0.5 * displVec);
  std::vector<Double> nonLinStrain;
  
  nonLinStrain = part1 + part2;
  piolaStressVec = dMat * nonLinStrain;
}


void nLinMech3dInt_PiolaStress::
calcMaterialDMat(Matrix<Double> & dMat)
{
  // dirty trick: dimension of d-matrix is set to 6 (as it should be in 3d mechanics)
  // this is done, because the special Piola-Kirchhoff-BDB uses a d-matrix of size 9
  // this means, the calculation of the "standard" (linear mech.) and the nonlinear b-matrix
  // would be of wrong dimension
  setPiolaDimD(6);

  // this is the original material matrix
  // calcDMat without 
  nLinMech3dInt_BNonLin::calcDMat(dMat);

  // set size back
  setPiolaDimD(9);
}



// returns linear B - matrix
void nLinMech3dInt_PiolaStress::
calcLinBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
{
  // dirty trick: dimension of d-matrix is set to 6 (as it should be in 3d mechanics)
  // this is done, because the special Piola-Kirchhoff-BDB uses a d-matrix of size 9
  // this means, the calculation of the "standard" (linear mech.) and the nonlinear b-matrix
  // would be of wrong dimension
  setPiolaDimD(6);
  
  // linear differential operator B_lin
  linElastInt::calcBMat(bMat, ip, ptCoord);
  setPiolaDimD(9);
}




  /// returns nonlinear B - matrix
void nLinMech3dInt_PiolaStress::
calcNonLinBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
{
  // dirty trick: dimension of d-matrix is set to 6 (as it should be in 3d mechanics)
  // this is done, because the special Piola-Kirchhoff-BDB uses a d-matrix of size 9
  // this means, the calculation of the "standard" (linear mech.) and the nonlinear b-matrix
  // would be of wrong dimension
  setPiolaDimD(6);
  
  // linear differential operator B_lin
  nLinElastInt::calcBMat(bMat, ip, ptCoord);
  setPiolaDimD(9);
}

  

  /// conversion of stress vector to stress tensor
void nLinMech3dInt_PiolaStress::
convertStressVecToTensor(Matrix<Double>& stressTensor, std::vector<Double>& piolaStress)
{
#ifdef TRACE
    (*trace) << "entering mech3DInt_PiolaStress::convertStressVecToTensor " << std::endl;
#endif

    Integer indexRow[] = {1, 2, 3, 2, 1, 1}; // first index of tensor notation
    Integer indexCol[] = {1, 2, 3, 3, 3, 2}; // second index of tensor notation

    stressTensor.Resize( getNrDofs() );
    
    // build symmetrical tensor
    for (Integer i=0; i<piolaStress.size(); i++)
      {
	stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[i];	
	stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[i];	
      }
}





// calculates the D-matrix needed for the Piola-Kirchhoff-Stress part
// This matrix is equal to a block diagonal matrix with Piola-Stresses in tensor 
// notation used as diagonal blocks (nrDim times)
// (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void nLinMech3dInt_PiolaStress::calcDMat(Matrix<Double> & dMat, Integer ip, Matrix<Double> & ptCoord)
  {
#ifdef TRACE
    (*trace) << "entering mech3DInt_PiolaStress::calcDMat " << std::endl;
#endif

    const Integer dimD = getDimD();
    const Integer nrDofs = getNrDofs();

    std::vector<Double> piolaStressVec;
    Matrix<Double> stressTensor;
    
    calcPiolaStressVec(piolaStressVec, ip, ptCoord );
    convertStressVecToTensor(stressTensor, piolaStressVec);

    // in "Resize", matrix elements are set to zero
    dMat.Resize(dimD);
    

    for (Integer i=0; i<nrDofs; i++)
      dMat.SetSubMatrix(stressTensor, i*nrDofs, i*nrDofs);
  }





// returns B - matrix for piola stresses
// (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void nLinMech3dInt_PiolaStress::
  calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech3dInt_PiolaStress::calcBMat " << std::endl;
#endif    

    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer nrDofs   = getNrDofs();    


    // in "Resize", matrix elements are set to zero
    bMat.Resize(getDimD(), nrNodes * nrDofs);

    
    // local shape functions derived after global coords 
    // (format: nrNodes x nrDofs)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);


    for(int actNode=0; actNode < nrNodes; actNode++)
      for(int globPos=0; globPos < nrDofs; globPos++)
	for(int actDof=0; actDof < nrDofs; actDof++)
	  bMat[globPos*nrDofs + actDof][actNode * nrDofs + globPos] = xiDx[actNode][actDof];      
  }








  


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
    // xiDx holds derivatives of shape functions after global coords in one IP (dimension: nrNodes x nrDofs)
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
  }
  


  // calculates the D-matrix of a 3d-problem 
  // (see mech3dInt::calcDMat)
  void nLinMech3dInt_BNonLin::calcDMat(Matrix<Double> & dMat)
  {
#ifdef TRACE
    (*trace) << "entering mech3DInt_BNonLin::calcDMat " << std::endl;
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
  nLinMech3dInt_BNonLin::nLinMech3dInt_BNonLin(BaseFE * aptelem, MaterialData & matData) 
    : nLinElastInt(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech3dInt_BNonLin::nLinMech3dInt_BNonLin" << std::endl;
#endif
  }
 

  nLinMech3dInt_BNonLin::~nLinMech3dInt_BNonLin()
  {
#ifdef TRACE
    (*trace) << "entering nLinMech3dInt_BNonLin::~nLinMech3dInt_BNonLin" << std::endl;
#endif
  }




  // nonlinear calculation of elasticity in 3d
  nLinMech3dInt_PiolaStress::nLinMech3dInt_PiolaStress(BaseFE * aptelem, MaterialData & matData) 
    : nLinMech3dInt_BNonLin(aptelem, matData), piolaDimD_(9)

  {
#ifdef TRACE
    (*trace) << "entering nLinMech3dInt_PiolaStress::nLinMech3dInt_PiolaStress" << std::endl;
#endif
    updateDMatInEveryIP_ = 1;
  }
 

  nLinMech3dInt_PiolaStress::~nLinMech3dInt_PiolaStress()
  {
#ifdef TRACE
    (*trace) << "entering nLinMech3dInt_PiolaStress::~nLinMech3dInt_PiolaStress" << std::endl;
#endif
  }

} // end namespace CoupledField
