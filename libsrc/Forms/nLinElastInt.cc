#include <iostream>
#include <fstream>

#include "nLinElastInt.hh"

namespace CoupledField
{



  // =============================================================================
  // base class for nonlinear mechanics  
  // =============================================================================

  // base class for nonlinear calculation of elasticity
  nLinElastInt::nLinElastInt(BaseFE * aptelem, MaterialData & matData) 
    : linElastInt(aptelem, matData)
  {
    ENTER_FCN( "nLinElastInt::nLinElastInt" );
  }


  // base class for nonlinear calculation of elasticity
  nLinElastInt::nLinElastInt(MaterialData & matData) 
    : linElastInt(matData)
  {
    ENTER_FCN( "nLinElastInt::nLinElastInt" );
  }
 

  nLinElastInt::~nLinElastInt()
  {
    ENTER_FCN( "nLinElastInt::~nLinElastInt" );
  }



  // returns nonlinear B - matrix (first part) for BDB
  void nLinElastInt::calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "nLinElastInt::calcBMat" );

    
    if (!elemDisp_.GetSizeRow() || !elemDisp_.GetSizeCol()) 
      Error("Undefined displacements! ",__FILE__,__LINE__);

    const Integer nrNodes = ptelem->GetNumNodes();
    const Integer nrDofs  = getNrDofs();    
    const Integer dimD    = getDimD();

    Matrix<Double> linBMat;    
    linElastInt::calcBMat( linBMat, ip, ptCoord);


    bMat.Resize(dimD, nrNodes * nrDofs);

    
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
    // size_row() = nr of rows!!
    bMatOneNode.Resize(linBMat.GetSizeRow(), nrDofs);

    

    for(int actNode=0; actNode < nrNodes; actNode++)
      {
	linBMat.GetSubMatrix(bMatOneNode, 0, actNode*nrDofs);

	bMatOneNode *= displDerivTransp;	

	bMat.SetSubMatrix(bMatOneNode, 0, actNode*nrDofs);
      }


    if (isaxi_)
      {
	const Integer spaceDim  = ptelem->GetDim();

	Vector<Double> shpFncAtIp;
	ptelem->GetShFncAtIp(shpFncAtIp, ip);

	Vector<Double> coordAtIP;
	coordAtIP = ptCoord * shpFncAtIp;

	Double  l33=0;  // for l33, see Bathe, page 552
	for (Integer actPos=0; actPos < shpFncAtIp.GetSize(); actPos++)
	  l33 += elemDisp_[0][actPos] * shpFncAtIp[actPos] / coordAtIP[0];
	
	
	for (Integer actNode = 0; actNode < nrNodes; actNode++)	     
	  {	    
	    //  (N_a/r) * (u_x/r)
	    bMat[dimD - 1][actNode * spaceDim]     = l33 * shpFncAtIp[actNode] / coordAtIP[0];
	    bMat[dimD - 1][actNode * spaceDim + 1] = 0;
	  }
      }

  }
  



  void nLinElastInt::
  Calc2DMaterialMatrix(Matrix<Double> & dMat, enum orientation2D actOrientation)
  {  
    ENTER_FCN( "nLinElastInt::Calc2DMaterialMatrix" );
    const Integer nrElems2d = 3;
  
    Integer rowPtrXY[] = {1,2,6};  // indices of rows and lines for xy-plane
    Integer rowPtrYZ[] = {2,3,4};  // indices of rows and lines for yz-plane
    Integer rowPtrXZ[] = {1,3,5};  // indices of rows and lines for xz-plane
    Integer * rowPtr;
    Integer i,j;
  
    switch(actOrientation)
      {	
      case xy: 
	{
	  rowPtr = rowPtrXY;
	  break;
	}
      case xz: 
	{
	  rowPtr = rowPtrXZ;
	  break;
	}
      
      case yz: 
	{
	  rowPtr = rowPtrYZ;    
	  break;
	}
      }    
  
    Matrix<Double> const & matMatrix =  *(ptMaterial->GetMatrix());
  
    dMat.Resize(nrElems2d);
  
    for (i=0; i<nrElems2d; i++)
      for (j=0; j<nrElems2d; j++)
	dMat[i][j] = matMatrix[rowPtr[i]-1][rowPtr[j]-1];	
  }




  // =============================================================================
  // 3D nonlinear mechanics (part with nonlinear B-mat)
  // =============================================================================


  nLinMech3dInt_BNonLin::nLinMech3dInt_BNonLin(BaseFE * aptelem, MaterialData & matData) 
    : nLinElastInt(aptelem, matData)
  {
    ENTER_FCN( "nLinMech3dInt_BNonLin::nLinMech3dInt_BNonLin" );

    className = "nLinMech3dInt_BNonLin";
  }


  nLinMech3dInt_BNonLin::nLinMech3dInt_BNonLin(MaterialData & matData) 
    : nLinElastInt(matData)
  {
    ENTER_FCN( "nLinMech3dInt_BNonLin::nLinMech3dInt_BNonLin" );
  }
 

  nLinMech3dInt_BNonLin::~nLinMech3dInt_BNonLin()
  {
    ENTER_FCN( "nLinMech3dInt_BNonLin::~nLinMech3dInt_BNonLin" );
  }


  // calculates the D-matrix of a 3d-problem 
  // (see mech3dInt::calcDMat)
  void nLinMech3dInt_BNonLin::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "nLinMech3dInt_BNonLin::calcDMat" );

    const Integer nrElems3d = getDimD();
    
    Matrix<Double> const & matMatrix = *(ptMaterial->GetMatrix());
    
    dMat.Resize(nrElems3d);

    for (Integer i=0; i<nrElems3d; i++)
      for (Integer j=0; j<nrElems3d; j++)
	dMat[i][j] = matMatrix[i][j];	
  }
  





  // =============================================================================
  // 3D nonlinear mechanics (part of Piola-Kirchhoff stress)
  // =============================================================================

  nLinMechInt_PiolaStress::nLinMechInt_PiolaStress(BaseFE * aptelem, MaterialData & matData) 
    : nLinElastInt(aptelem, matData)

  {
    ENTER_FCN( "nLinMechInt_PiolaStress::nLinMechInt_PiolaStress" );

    updateDMatInEveryIP_ = 1;
    //    setPiolaDimD( getFullPiolaDMatSize() );
    className = "  nLinMechInt_PiolaStress";
  }



  nLinMechInt_PiolaStress::nLinMechInt_PiolaStress(MaterialData & matData) 
    : nLinElastInt(matData)
  {
    ENTER_FCN( "nLinMechInt_PiolaStress::nLinMechInt_PiolaStress" );

    updateDMatInEveryIP_ = 1;
    //    setPiolaDimD( getFullPiolaDMatSize() );
  }
 

  nLinMechInt_PiolaStress::~nLinMechInt_PiolaStress()
  {
    ENTER_FCN( "nLinMechInt_PiolaStress::~nLinMechInt_PiolaStress" );
  }






  /// calculates Piola-Kirchoff-stresses T (vector notation)
  // T = c . V with c the stiffness matrix and V the Green-Lagrange strain vector
  // V = B_lin * u  +  1/2 * B_nonLin * u
  // see Habil. M. Kaltenbacher Eq. (3.33)
  void nLinMechInt_PiolaStress::
  CalcStressVec(Vector<Double>& piolaStressVec, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "mech3DInt_PiolaStress::calcPiolaStressTensor" );

    Matrix<Double> dMat;
    calcMaterialDMat(dMat);
  
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<Double> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);


    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcLinBMat( linBMat, ip, ptCoord);


    // nonlinear differential operator B_nonLin
    Matrix<Double> nLinBMat;    
    calcNonLinBMat( nLinBMat, ip, ptCoord);


    // special construction: direct addition of these two terms does not work, 
    // because of temporary vectors ... :O(    
    Vector<Double> part1(linBMat * displVec );
    Vector<Double> part2(nLinBMat * 0.5 * displVec);
    Vector<Double> nonLinStrain;
  
    nonLinStrain = part1 + part2;
    piolaStressVec = dMat * nonLinStrain;
  }




  void nLinMechInt_PiolaStress::
  calcMaterialDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "mech3DInt_PiolaStress::calcMaterialDMat" );

    const Integer nrElems3d = getMaterialDMatSize();
  
    Matrix<Double> const &  matMatrix = (*ptMaterial->GetMatrix());
  
    dMat.Resize(nrElems3d);
  
    for (Integer i=0; i<nrElems3d; i++)
      for (Integer j=0; j<nrElems3d; j++)
	dMat[i][j] = matMatrix[i][j];	
  }




  // returns linear B - matrix
  void nLinMechInt_PiolaStress::
  calcLinBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "mech3DInt_PiolaStress::calcLinBMat" );

    // dirty trick: dimension of d-matrix is set to 6 (as it should be in 3d mechanics)
    // this is done, because the special Piola-Kirchhoff-BDB uses a d-matrix of size 9
    // this means, the calculation of the "standard" (linear mech.) and the nonlinear b-matrix
    // would be of wrong dimension
    setPiolaDimD(  getMaterialDMatSize() );
  
    // linear differential operator B_lin
    linElastInt::calcBMat(bMat, ip, ptCoord);

    setPiolaDimD( getFullPiolaDMatSize() );
  }





  /// returns nonlinear B - matrix
  void nLinMechInt_PiolaStress::
  calcNonLinBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "nLinMechInt_PiolaStress::calcNonLinBMat" );
    // dirty trick: dimension of d-matrix is set to 6 (as it should be in 3d mechanics)
    // this is done, because the special Piola-Kirchhoff-BDB uses a d-matrix of size 9
    // this means, the calculation of the "standard" (linear mech.) and the nonlinear b-matrix
    // would be of wrong dimension
    setPiolaDimD(  getMaterialDMatSize() );
  
    // linear differential operator B_lin
    nLinElastInt::calcBMat(bMat, ip, ptCoord);
    setPiolaDimD( getFullPiolaDMatSize() );
  }

  

  /// conversion of stress vector to stress tensor
  void nLinMechInt_PiolaStress::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {
    ENTER_FCN( "mech3DInt_PiolaStress::convertStressVecToTensor" );

    Integer indexRow[] = {1, 2, 3, 2, 1, 1}; // first index of tensor notation
    Integer indexCol[] = {1, 2, 3, 3, 3, 2}; // second index of tensor notation

    stressTensor.Resize( getNrDofs() );
    
    // build symmetrical tensor
    for (Integer i=0; i<piolaStress.GetSize(); i++)
      {
	stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[i];	
	stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[i];	
      }

  }





  // calculates the D-matrix needed for the Piola-Kirchhoff-Stress part
  // This matrix is equal to a block diagonal matrix with Piola-Stresses in tensor 
  // notation used as diagonal blocks (nrDim times)
  // (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void nLinMechInt_PiolaStress::calcDMat(Matrix<Double> & dMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "nLinMechInt_PiolaStress::calcDMat");

    const Integer dimD = getDimD();
    const Integer nrDofs = getNrDofs();

    Vector<Double> piolaStressVec;
    Matrix<Double> stressTensor;
    
    CalcStressVec(piolaStressVec, ip, ptCoord );
    convertStressVecToTensor(stressTensor, piolaStressVec);

    // in "Resize", matrix elements are set to zero
    dMat.Resize(dimD);
    

    if (isaxi_)
      // the stressTensor has already the correct shape due to the 
      // method "convertStressVecToTensor"
      dMat = stressTensor;
    else
      for (Integer i=0; i<nrDofs; i++)
	dMat.SetSubMatrix(stressTensor, i*nrDofs, i*nrDofs);    
  }






  // returns B - matrix for piola stresses
  // (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void nLinMechInt_PiolaStress::
  calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "nLinMechInt_PiolaStress::calcBMat" );

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


    if (isaxi_)
      {
	const Integer spaceDim = ptelem->GetDim();
	Vector<Double> shpFncAtIp;
	Vector<Double> coordAtIP;
	
	ptelem->GetShFncAtIp(shpFncAtIp,ip);
	coordAtIP = ptCoord * shpFncAtIp;
	
	for (int actNode = 0; actNode < nrNodes; actNode++)	     
	  bMat[getDimD() -1][actNode * spaceDim] = shpFncAtIp[actNode] / coordAtIP[0];
      }
  }










  // =============================================================================
  // class for nonlinear 3d piola stresses
  // =============================================================================




  nLinMech3dInt_PiolaStress::nLinMech3dInt_PiolaStress(BaseFE * aptelem, MaterialData & matData) 
    :nLinMechInt_PiolaStress(aptelem, matData)

  {
    ENTER_FCN( "nLinMech3dInt_PiolaStress::nLinMech3dInt_PiolaStress" );

    updateDMatInEveryIP_ = 1;
    setPiolaDimD( getFullPiolaDMatSize() );
    className = "  nLinMech3dInt_PiolaStress";
  }



  nLinMech3dInt_PiolaStress::nLinMech3dInt_PiolaStress(MaterialData & matData) 
    :nLinMechInt_PiolaStress(matData)
  {
    ENTER_FCN( "nLinMech3dInt_PiolaStress::nLinMech3dInt_PiolaStress");

    updateDMatInEveryIP_ = 1;
    setPiolaDimD( getFullPiolaDMatSize() );
  }
 

  nLinMech3dInt_PiolaStress::~nLinMech3dInt_PiolaStress()
  {
    ENTER_FCN( "nLinMech3dInt_PiolaStress::~nLinMech3dInt_PiolaStress" );

  }







  // =============================================================================
  // class regarding mechanical prestresses
  // =============================================================================


  // base class for regarding prestress
  PreStressInt::PreStressInt(BaseFE * aptelem, MaterialData & matData, Double aPreStressVal, Directions stressDir) 
    : nLinMechInt_PiolaStress(aptelem, matData), preStressVal_(aPreStressVal), preStressDir_(stressDir)
  {
    ENTER_FCN( "PreStressInt::PreStressInt" );

    updateDMatInEveryIP_ = 1;
    setPiolaDimD( getFullPiolaDMatSize() );
    className = "PreStressInt";
  }
 
  // base class for regarding prestress
  PreStressInt::PreStressInt(MaterialData & matData, Double aPreStressVal, Directions stressDir) 
    : nLinMechInt_PiolaStress(matData), preStressVal_(aPreStressVal), preStressDir_(stressDir)
  {
    ENTER_FCN( "PreStressInt::PreStressInt" );

    updateDMatInEveryIP_ = 1;
    setPiolaDimD( getFullPiolaDMatSize() );
    className = "PreStressInt";
  }


  PreStressInt::~PreStressInt()
  {
    ENTER_FCN( "PreStressInt::~PreStressInt" );

  }


  // calculates the D-matrix needed for regarding pre-stresses
  void PreStressInt::CalcStressVec(Vector<Double>& preStressVec, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "PreStressInt::CalcStressVec" );

    // ???????????????????????????????????????????????
    //    const Integer dimD = nLinMechInt_BNonLin::getDimD();
    const Integer dimD = getDimD();
    const Integer nrDofs = getNrDofs();
    const Integer distributedStress = 1;
    

        
    preStressVec.Resize(dimD);
    preStressVec *= 0;
    
    if (! distributedStress)
      switch(preStressDir_)
	{
	case X:
	  preStressVec[0] = preStressVal_;
	  break;
	
	case Y:
	  preStressVec[1] = preStressVal_;
	  break;
	  
	case Z:
	  preStressVec[2] = preStressVal_;
	  break;
	  
	case radXY:
	  {  
	    Vector<Double> shapeFncIP;
	    ptelem->GetShFncAtIp(shapeFncIP, ip);
	    
	    Vector<Double> coordIP(nrDofs);
	    coordIP = ptCoord * shapeFncIP;
	    
	    
	    preStressVec[0] = coordIP[0];
	    preStressVec[1] = coordIP[1];
	    // in xy-plane, no z-coord has to be regarded
	    preStressVec[2] = 0;
	    
	    
	    // center of radial prestress is (0,0)
	    preStressVec *= (preStressVal_ / preStressVec.NormL2());
	    break;
	  }
	  
	default:
	  Error("The prestress direction mentioned in the config-file is not implemented! ",__FILE__,__LINE__);
	}
    else
      {
	Vector<Double> preStrainVec(dimD);
	
	preStrainVec *= 0;   // initialize vector
	
	    
	Matrix<Double> matData;
	if (isaxi_)
	  //axi
	  CalcAxiMaterialMat(matData, actOrientation);
	else if (getNrDofs() == 3)
	  //3D
	  calcMaterialDMat(matData);
	else
	  //2D-plane strain
	  Calc2DMaterialMatrix(matData, actOrientation);	    

	switch(preStressDir_)
	  {
	  case X:
	    preStrainVec[0] = preStressVal_ / matData[0][0];
	    break;
		
	  case Y:
	    preStrainVec[1] = preStressVal_ / matData[1][1];
	    break;
		
	  case Z:
	    preStrainVec[2] = preStressVal_ / matData[2][2];
	    break;
		
	  case radXY:
	    {  
	      Vector<Double> shapeFncIP;
	      ptelem->GetShFncAtIp(shapeFncIP, ip);
		  
	      Vector<Double> coordIP(nrDofs);
	      coordIP = ptCoord * shapeFncIP;
		  
	      // in xy-plane, no z-coord has to be regarded
	      // center of radial prestress is (0,0)
	      coordIP[2]= 0;
		  
	      for (Integer i=0; i<coordIP.GetSize(); i++)
		preStrainVec[i] = coordIP[i] / coordIP.NormL2() * preStressVal_ / matData[i][i];
	      break;
	    }
		
	  default:
	    Error("The prestress direction mentioned in the config-file is not implemented! ",__FILE__,__LINE__);
	  }

	// calc "distributed" stress 
	preStressVec = matData * preStrainVec;
	    
      }
  }




  // (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void PreStressInt::calcDMat(Matrix<Double> & dMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "PreStressInt::calcDMat");

    const Integer dimD = getFullPiolaDMatSize();
    const Integer nrDofs = getNrDofs();

    Vector<Double> piolaStressVec;
    Matrix<Double> stressTensor;
    
    CalcStressVec(piolaStressVec, ip, ptCoord );
    convertStressVecToTensor(stressTensor, piolaStressVec);

    // in "Resize", matrix elements are set to zero
    dMat.Resize(dimD);
    

    if (isaxi_)
      // the stressTensor has already the correct shape due to the 
      // method "convertStressVecToTensor"
      dMat = stressTensor;
    else
      for (Integer i=0; i<nrDofs; i++)
	dMat.SetSubMatrix(stressTensor, i*nrDofs, i*nrDofs);    
  }


  // returns B - matrix for piola stresses
  // (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void PreStressInt::calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord)
  {
    ENTER_FCN( "PreStressInt::calcBMat" );

    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer nrDofs   = getNrDofs();    

    Integer dimD = getFullPiolaDMatSize();
    // in "Resize", matrix elements are set to zero
    bMat.Resize(dimD, nrNodes * nrDofs);

    
    // local shape functions derived after global coords 
    // (format: nrNodes x nrDofs)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);


    for(int actNode=0; actNode < nrNodes; actNode++)
      for(int globPos=0; globPos < nrDofs; globPos++)
	for(int actDof=0; actDof < nrDofs; actDof++)
	  bMat[globPos*nrDofs + actDof][actNode * nrDofs + globPos] = xiDx[actNode][actDof];      


    if (isaxi_)
      {
	const Integer spaceDim = ptelem->GetDim();
	Vector<Double> shpFncAtIp;
	Vector<Double> coordAtIP;
	
	ptelem->GetShFncAtIp(shpFncAtIp,ip);
	coordAtIP = ptCoord * shpFncAtIp;
	
	for (int actNode = 0; actNode < nrNodes; actNode++)	     
	  bMat[getDimD() -1][actNode * spaceDim] = shpFncAtIp[actNode] / coordAtIP[0];
      }
  }

  /// conversion of stress vector to stress tensor
  void PreStressInt::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {
    ENTER_FCN( "PreStressInt3D::convertStressVecToTensor" );

    Integer indexRow[] = {1, 2, 3, 2, 1, 1}; // first index of tensor notation
    Integer indexCol[] = {1, 2, 3, 3, 3, 2}; // second index of tensor notation

    stressTensor.Resize( getNrDofs() );
    
    // build symmetrical tensor
    for (Integer i=0; i<piolaStress.GetSize(); i++)
      {
	stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[i];	
	stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[i];	
      }

  }

  // class for regarding 3d prestress
  PreStressInt3D::PreStressInt3D(BaseFE * aptelem, MaterialData & matData, Double aPreStressVal, 
				 Directions stressDir) 
    : PreStressInt(aptelem, matData, aPreStressVal, stressDir)
  {
    ENTER_FCN( "PreStressInt3D::PreStressInt3D" );
  }
 
  // class for regarding 3d prestress
  PreStressInt3D::PreStressInt3D(MaterialData & matData, Double aPreStressVal, Directions stressDir) 
    : PreStressInt(matData, aPreStressVal, stressDir)
  {
    ENTER_FCN( "PreStressInt3D::PreStressInt3D" );
  }


  PreStressInt3D::~PreStressInt3D()
  {
    ENTER_FCN( "PreStressInt3D::~PreStressInt3D" );

  }


  // class for regarding 2d prestress in plane strain case
  PreStressIntPlaneStrain::PreStressIntPlaneStrain(BaseFE * aptelem, MaterialData & matData, Double aPreStressVal, 
				 Directions stressDir) 
    : PreStressInt(aptelem, matData, aPreStressVal, stressDir)
  {
    ENTER_FCN( "PreStressIntPlaneStrain::PreStressIntPlaneStrain" );
  }
 
  // class for regarding 2d prestress in plane strain case
  PreStressIntPlaneStrain::PreStressIntPlaneStrain(MaterialData & matData, Double aPreStressVal, 
						   Directions stressDir) 
    : PreStressInt(matData, aPreStressVal, stressDir)
  {
    ENTER_FCN( "PreStressIntPlaneStrain::PreStressIntPlaneStrain" );
  }


  PreStressIntPlaneStrain::~PreStressIntPlaneStrain()
  {
    ENTER_FCN( "PreStressIntPlaneStrain::~PreStressIntPlaneStrain" );

  }


  // class for regarding 2d axi prestress 
  PreStressIntAxi::PreStressIntAxi(BaseFE * aptelem, MaterialData & matData, Double aPreStressVal, 
				   Directions stressDir) 
    : PreStressInt(aptelem, matData, aPreStressVal, stressDir)
  {
    ENTER_FCN( "PreStressIntAxi::PreStressIntAxi" );
   isaxi_ = TRUE;
  }
 
  // class for regarding 2d axi prestress 
  PreStressIntAxi::PreStressIntAxi(MaterialData & matData, Double aPreStressVal, 
				   Directions stressDir) 
    : PreStressInt(matData, aPreStressVal, stressDir)
  {
    ENTER_FCN( "PreStressIntAxi::PreStressIntAxi" );
   isaxi_ = TRUE;
  }


  PreStressIntAxi::~PreStressIntAxi()
  {
    ENTER_FCN( "PreStressIntAxi::~PreStressIntAxi" );

  }

  /// conversion of stress vector to stress tensor
  void PreStressIntAxi::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {
    ENTER_FCN( "PreStressIntAxi::convertStressVecToTensor" );

    // indizes see Bathe: "Finite Element Procedures", Sec. 6.3, 
    // page 553: "2. Piola-Kirchhoff stress matrix & vector"
    const Integer indexSize = 9;
    Integer indexRow[]   = {1, 2, 1, 2, 3, 4, 3, 4, 5}; // first index of tensor notation
    Integer indexCol[]   = {1, 2, 2, 1, 3, 4, 4, 3, 5}; // second index of tensor notation
    Integer indexPiola[] = {1, 2, 3, 3, 1, 2, 3, 3, 4}; // index in piola stress vec

    const Integer axiTensSize = 5;
    
    stressTensor.Resize( axiTensSize );
    stressTensor.Init();
    
    // build symmetrical tensor
    for (Integer i=0; i<indexSize; i++)
      {
	stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[ indexPiola[i] -1 ];	
	stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[ indexPiola[i] -1 ];	
      }

  }

  // =============================================================================
  // 2D nonlinear plane strain mechanics
  // =============================================================================

  
  // calculated the D-matrix for the plain strain state
  void  nLinMechPlaneStrainInt_BNonLin::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "mechPlainStrainNLinIntInt::calcDMat" );

    Calc2DMaterialMatrix(dMat, actOrientation);
  
  }

  




  void nLinMechPlaneStrainInt_PiolaStress::calcMaterialDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "mechPlaneStraneInt_PiolaStress::calcMaterialDMat" );

    Calc2DMaterialMatrix(dMat, actOrientation);
  }




  /// conversion of stress vector to stress tensor
  void nLinMechPlaneStrainInt_PiolaStress::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {
    ENTER_FCN( "mechMechPlanseStrainInt_PiolaStress::convertStressVecToTensor" );

    Integer indexRow[] = {1, 2, 1}; // first index of tensor notation
    Integer indexCol[] = {1, 2, 2}; // second index of tensor notation

    stressTensor.Resize( getNrDofs() );
    
    // build symmetrical tensor
    for (Integer i=0; i<piolaStress.GetSize(); i++)
      {
	stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[i];	
	stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[i];	
      }

  }







  // =============================================================================
  // 2D nonlinear axisymmetric mechanics
  // =============================================================================

  
  // calculated the D-matrix for the axi state
  void  nLinMechAxiInt_BNonLin::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "mechPlainStrainNLinIntInt::calcDMat" );

    CalcAxiMaterialMat(dMat, actOrientation);
  }



  nLinMechAxiInt_PiolaStress::nLinMechAxiInt_PiolaStress(BaseFE * aptelem, MaterialData & matData) 
    : nLinMechInt_PiolaStress(aptelem, matData)
  {
    ENTER_FCN( "nLinMechAxiInt_PiolaStress::nLinMechAxiInt_PiolaStress" );

    isaxi_ = TRUE;
    setPiolaDimD( getFullPiolaDMatSize() );
    className = "nLinMechAxiInt_PiolaStress";    
  }



  // nonlinear calculation of elasticity for axi 
  nLinMechAxiInt_PiolaStress::nLinMechAxiInt_PiolaStress(MaterialData & matData) 
    : nLinMechInt_PiolaStress(matData)

  {
    ENTER_FCN( "nLinMechAxiInt_PiolaStress::nLinMechAxiInt_PiolaStress" );

    isaxi_ = TRUE;
    setPiolaDimD( getFullPiolaDMatSize() );
  }
 

  nLinMechAxiInt_PiolaStress::~nLinMechAxiInt_PiolaStress()
  {
    ENTER_FCN( "nLinMechAxiInt_PiolaStress::~nLinMechAxiInt_PiolaStress" );

  }




  void nLinMechAxiInt_PiolaStress::calcMaterialDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "nLinMechAxiInt_PiolaStress::calcMaterialDMat" );

    CalcAxiMaterialMat(dMat, actOrientation);
  }


  /// conversion of stress vector to stress tensor
  void nLinMechAxiInt_PiolaStress::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {
    ENTER_FCN( "nLinMechAxiInt_PiolaStress::convertStressVecToTensor" );

    // indizes see Bathe: "Finite Element Procedures", Sec. 6.3, 
    // page 553: "2. Piola-Kirchhoff stress matrix & vector"
    const Integer indexSize = 9;
    Integer indexRow[]   = {1, 2, 1, 2, 3, 4, 3, 4, 5}; // first index of tensor notation
    Integer indexCol[]   = {1, 2, 2, 1, 3, 4, 4, 3, 5}; // second index of tensor notation
    Integer indexPiola[] = {1, 2, 3, 3, 1, 2, 3, 3, 4}; // index in piola stress vec

    const Integer axiTensSize = 5;
    
    stressTensor.Resize( axiTensSize );
    stressTensor.Init();
    
    // build symmetrical tensor
    for (Integer i=0; i<indexSize; i++)
      {
	stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[ indexPiola[i] -1 ];	
	stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[ indexPiola[i] -1 ];	
      }

  }











  // ===================================================================================
  // nonlinear calculation of elasticity in plane strain state 
  // ===================================================================================

  nLinMechPlaneStrainInt_BNonLin::nLinMechPlaneStrainInt_BNonLin(BaseFE * aptelem, MaterialData & matData) 
    : nLinElastInt(aptelem, matData)
  {
    ENTER_FCN( "nLinMechPlaneStrainInt_BNonLin::nLinMechPlaneStrainInt_BNonLin" );

    className = "nLinMechPlaneStrainInt_BNonLin";
  }


  nLinMechPlaneStrainInt_BNonLin::nLinMechPlaneStrainInt_BNonLin(MaterialData & matData) 
    : nLinElastInt(matData)
  {
    ENTER_FCN( "nLinMechPlaneStrainInt_BNonLin::nLinMechPlaneStrainInt_BNonLin" );

  }
 

  nLinMechPlaneStrainInt_BNonLin::~nLinMechPlaneStrainInt_BNonLin()
  {
    ENTER_FCN( "nLinMechPlaneStrainInt_BNonLin::~nLinMechPlaneStrainInt_BNonLin" );

  }






  nLinMechPlaneStrainInt_PiolaStress::nLinMechPlaneStrainInt_PiolaStress(BaseFE * aptelem, MaterialData & matData) 
    : nLinMechInt_PiolaStress(aptelem, matData)
  {
    ENTER_FCN( "nLinMechPlaneStrainInt_PiolaStress::nLinMechPlaneStrainInt_PiolaStress" );

    setPiolaDimD( getFullPiolaDMatSize() );

    className = "nLinMechPlaneStrainInt_PiolaStress";    
  }



  // nonlinear calculation of elasticity in 3d
  nLinMechPlaneStrainInt_PiolaStress::nLinMechPlaneStrainInt_PiolaStress(MaterialData & matData) 
    : nLinMechInt_PiolaStress(matData)

  {
    ENTER_FCN( "nLinMechPlaneStrainInt_PiolaStress::nLinMechPlaneStrainInt_PiolaStress" );

    setPiolaDimD( getFullPiolaDMatSize() );
  }
 

  nLinMechPlaneStrainInt_PiolaStress::~nLinMechPlaneStrainInt_PiolaStress()
  {
    ENTER_FCN( "nLinMechPlaneStrainInt_PiolaStress::~nLinMechPlaneStrainInt_PiolaStress" );
  }





  // ===================================================================================
  // axisymmetric nonlinear calculation of elasticity
  // ===================================================================================

  nLinMechAxiInt_BNonLin::nLinMechAxiInt_BNonLin(BaseFE * aptelem, MaterialData & matData) 
    : nLinElastInt(aptelem, matData)
  {
    ENTER_FCN( "nLinMechAxiInt_BNonLin::nLinMechAxiInt_BNonLin" );

    isaxi_ = TRUE;
    className = "nLinMechAxiInt_BNonLin";    
  }


  nLinMechAxiInt_BNonLin::nLinMechAxiInt_BNonLin(MaterialData & matData) 
    : nLinElastInt(matData)
  {
    ENTER_FCN( "nLinMechAxiInt_BNonLin::nLinMechAxiInt_BNonLin" );

    isaxi_ = TRUE;
  }
 

  nLinMechAxiInt_BNonLin::~nLinMechAxiInt_BNonLin()
  {
    ENTER_FCN( "nLinMechAxiInt_BNonLin::~nLinMechAxiInt_BNonLin" );

  }








} // end namespace CoupledField
