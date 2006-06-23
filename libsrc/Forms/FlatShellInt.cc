#include <iostream>
#include <fstream>

#include "FlatShellInt.hh"
#define PI 3.141592654

namespace CoupledField {

  void FlatShellInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                        EntityIterator& ent1, 
                                        EntityIterator& ent2 ) {

    ENTER_FCN( "FlatShellInt::CalcElementMatrix" );
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );    

    const UInt nrNodes  = ptelem->GetNumNodes();   
    const UInt nrDofs   = getNrDofs();  
    double jacDet=0.0;
        

    Matrix<Double> bMat; //B matrix
    Matrix<Double> dMat; //D matrix
    Matrix<Double> dB; //D*B matrix
    Matrix<Double> bTrans; //Transpose of B
    Matrix<Double> partElemMat; //Partial Element Matrix
    Matrix<Double> TransMat; //Transformation matrix
    Matrix<Double> ShellCoord; //Shell Coordinates

    //Element matrix in the case of Flat Shells with linear elements is 24x24
    elemMat.Resize(nrNodes * nrDofs );
    elemMat.Init();
       
    //B matrix in the case of flat shells 8x24 
    bMat.Resize( getDimD(), nrNodes * nrDofs );
    bMat.Init();
    
    //D matrix in the case of flat shells 8x8
    dMat.Resize(getDimD()); 
    dMat.Init();
      
    //D*B matrix in the case of flat shells 8x24
    dB.Resize( getDimD(), nrNodes * nrDofs );
    dB.Init();
     
    //The transpose of B in the case of flat shells is 24x8
    bTrans.Resize( nrNodes * nrDofs , getDimD() );
    bTrans.Init();
      
    //The partial Element Matrix with the size of 24x24 
    partElemMat.Resize( nrNodes * nrDofs );
    partElemMat.Init();
         
    // 1. rotate 3D-coordintates to 2D (x-y-reference plane) on the element level,4 nodes x 6 coordinates = 24
    //and stor the full Element transformation matrix in TransMat
	
    FlatShellInt::CoordTrans(ptCoord_, TransMat, ShellCoord );
    
    // 2. calculate partial B- and D-Matrices and assemble into big one
	
    // If the material parameters are constant within the element
    // we can compute the D matrix once and for all
    if ( updateDMatInEveryIP_ == false ) {
      FlatShellInt::calcDMat( dMat );
    }

    for (UInt i =0 ; i<parts_.GetSize(); i++ ) {
	
      //Set the reduced integration          
      if (reducedPart_[i] == true ) {
        ptelem->SetReducedIntegration();
      }
      const UInt nrIntPts = ptelem->GetNumIntPoints();
      const Vector<Double> & intWeights = ptelem->GetIntWeights();
      //Main Cycle	
      for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

        // Check if D matrix must be re-determined for the current integration point
        if ( updateDMatInEveryIP_ == true ) {
          FlatShellInt::calcDMat(dMat);
        }
        //The calcBMat function calculates the B matrix either membrane shear or bending part
        if (parts_[i] == MEMBRANE || parts_[i] == BENDING || parts_[i] == SHEAR) {
	
          FlatShellInt::calcBMat(bMat, actIntPt, ShellCoord, parts_[i]);
          //Calculates the product D*B
          dB=dMat*bMat; 
          //Find the transpose of the bMat and Store it in the bTrans   
          bMat.Transpose(bTrans);
          //3. Calculate B^t*B*D and store it in the partElemMat
          partElemMat = bTrans * dB;
        
        } else if (parts_[i] == COUPLED1) {
          //Calculates B bending part
          FlatShellInt::calcBMat(bMat, actIntPt, ShellCoord, parts_[1]);
          //Calculates the product D*B
          dB=dMat*bMat;
          //Calculates B membrane part
          FlatShellInt::calcBMat(bMat, actIntPt, ShellCoord, parts_[0]); 
          //Find the transpose of the bMat and Store it in the bTrans   
          bMat.Transpose(bTrans);

          partElemMat = bTrans * dB;
		
        } else if (parts_[i] == COUPLED2) {
          //Calculates B membrane part
          FlatShellInt::calcBMat(bMat, actIntPt, ShellCoord, parts_[0]);
          //Calculates the product D*B
          dB=dMat*bMat;
          //Calculates B bending part
          FlatShellInt::calcBMat(bMat, actIntPt, ShellCoord, parts_[1]); 
          //Find the transpose of the bMat and Store it in the bTrans   
          bMat.Transpose(bTrans);
          //3. Calculate B^t*B*D and store it in the partElemMat
          partElemMat = bTrans * dB;
        }
	
        jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ShellCoord );

        // Perform a safety check
        if ( jacDet < 0.0 ) {
          (*error) << "FlatShellInt::CalcElementMatrix: Encountered "
                   << "negative Jacobian determinant!";
          Error( __FILE__, __LINE__ );
        }    
    
        // Sum up part element matrices
        elemMat += partElemMat * jacDet * intWeights[actIntPt-1];

      }
     
      if (reducedPart_[i] == true ) 
        ptelem->SetStandardIntegration();
    }
	
    // 4. Do the back-Rotation and return
    FlatShellInt::LocaltoGlob(elemMat,TransMat );
        
  }

  // For Coordinate Transformations

  void FlatShellInt::CoordTrans( const Matrix<Double> &ptCoord, Matrix<Double> &TransMat, Matrix<Double> &ShellCoord )
  {
    ENTER_FCN( "FlatShellInt::CoordTrans");
   
    const UInt nrNodes  = ptelem->GetNumNodes();   
    const UInt nrDofs   = getNrDofs();
    const UInt lambdaDim = 3; //dimension of the submatrix lambda, part of the transformation matrix
    const UInt TMatDim = getNrDofs();; //dimension of the nodal transformation matrix
    const UInt ElementTransMatDim = nrNodes*nrDofs;//dimension of the element transformation matrix
    Matrix<Double> lambdaMat; //size 3x3 submatrix
    Matrix<Double> TMat;//size 6x6
    Matrix<Double> NewCoord, temp; //matrices of size 3x4
    Vector<Double> Vx, Vr, Vz; //vectors of length 3
    const UInt row = ptCoord.GetSizeRow();// size 3
    const UInt col = ptCoord.GetSizeCol();// size 4
    UInt i=0,j=0,k=0;
    
    TransMat.Resize(ElementTransMatDim);//Resize the Element Transformation Matrix 24x24 and initialize with zeroes
    TransMat.Init();
    TMat.Resize(TMatDim);//Resize the whole node transformation matrix 6x6
    TMat.Init();
    lambdaMat.Resize(lambdaDim);//Resize the lambda submatrix with dimension 3x3
    lambdaMat.Init();
    NewCoord.Resize( row, col ); //Resize NewCoord Matrix 3x4 and initialize with zeroes
    NewCoord.Init();
    temp.Resize( row, col ); //Resize the temporary Matrix 3x4 and initialize with zeroes
    temp.Init();
    ShellCoord.Resize( row - 1, col ); //Resize the shell coordinates 2D 2x4 and initialize with zeroes
    ShellCoord.Init();
    Vx.Resize(row); //Resize the vectors of directional cosines
    Vr.Resize(row);
    Vz.Resize(row);
    Double length=0.0;

    //Vector Vx which represents the edge ij,which coincides with the local axis
 
    Vx[0]= ptCoord[0][1] - ptCoord[0][0] ;
    Vx[1]= ptCoord[1][1] - ptCoord[1][0] ;
    Vx[2]= ptCoord[2][1] - ptCoord[2][0] ;

    // 1/Length of the vector Vx
    length = 1.0 / sqrt( Vx[0] * Vx[0] + Vx[1] * Vx[1] + Vx[2] * Vx[2] );
 
    //directional cosines lambda_x for local x-direction
    //   cos(x,Z), cos(x,Y), cos(x,Z)

    lambdaMat[0][0] = Vx[0] * length;
    lambdaMat[0][1] = Vx[1] * length;
    lambdaMat[0][2] = Vx[2] * length;

    //Reference Vector Vr which FlatShellIntrepresents the edge kj,which is used for cross product with Vx

    Vr[0] = ptCoord[0][2] - ptCoord[0][0] ;
    Vr[1] = ptCoord[1][2] - ptCoord[1][0] ;
    Vr[2] = ptCoord[2][2] - ptCoord[2][0] ;

    //Cross product of vector Vx and Vr,which defines the direction of the vector Vz

    Vz[0] = Vx[1] * Vr[2] - Vx[2] * Vr[1] ;
    Vz[1] = Vx[2] * Vr[0] - Vx[0] * Vr[2] ;
    Vz[2] = Vx[0] * Vr[1] - Vx[1] * Vr[0] ;

    // 1/Length of the vector Vz

    length = 1.0 / sqrt( Vz[0] * Vz[0] + Vz[1] * Vz[1] + Vz[2] * Vz[2] );

    //directional cosines lambda_z for local z-direction
    // cos(z,X), cos(z,Y), cos(z,Z)

    lambdaMat[2][0] = Vz[0] * length;
    lambdaMat[2][1] = Vz[1] * length;
    lambdaMat[2][2] = Vz[2] * length;

    //Directional cosins lambda_y, obtained from the cross product of lamda_x and lambda_z
    // cos(y,X), cos(y,Y), cos(y,Z)

    lambdaMat[1][0] = lambdaMat[0][2] * lambdaMat[2][1] - lambdaMat[0][1] * lambdaMat[2][2];
    lambdaMat[1][1] = lambdaMat[0][0] * lambdaMat[2][2] - lambdaMat[0][2] * lambdaMat[2][0];
    lambdaMat[1][2] = lambdaMat[0][1] * lambdaMat[2][0] - lambdaMat[0][0] * lambdaMat[2][1];
    
    // transform geometry from real(global) coordinate to standard(local) coordinate
    // the temp matrix is (X-x0,Y-y0,Z-z0)T
    // seepage 35 (4.37) from MSc. Thesis of Yue Xiong "Finite Element Simulation of Plate and Shell Elements"

    for( i = 0; i < row; i++ )
      for(  j = 0; j < col; j++ )
        temp[i][j] = ptCoord[i][j] - ptCoord[i][0]; //transform

    //The local 3D coordiantes New Coord

    NewCoord = lambdaMat * temp;

    //The new coordinates of the Shell are 1 dimenstion less.Only x and y dimension is transfered.The z dimension is zero.
    for( i = 0; i < row - 1; i++)
      for( j = 0; j < col; j++)
        ShellCoord[i][j] = NewCoord[i][j];//ptCoord[i][j];

    //The Transformation nodal T matrix is 6x6
    for(UInt k=0;k < TMatDim; k+=3)
      for ( i=0; i < lambdaDim; i++ )
        for ( j=0; j< lambdaDim; j++ )
          TMat[i+k][j+k]=lambdaMat[i][j];
     
    //The transformation Element matrix which is 24x24
    for( k=0;k < ElementTransMatDim ; k+=6)
      for ( i=0; i < TMatDim; i++ )
        for ( j=0; j< TMatDim; j++ )
          TransMat[i+k][j+k]=TMat[i][j];
  
  }

  void FlatShellInt::LocaltoGlob( Matrix<Double> &ElemMat, const Matrix<Double> &TransMat )
  {

    ENTER_FCN( "FlatShellInt::LocaltoGlob" );

    int i, j;
    const Integer row = ElemMat.GetSizeRow(); //equals 24;


    Matrix<Double> StiffTrans;
    Matrix<Double> TFTrans;
    Matrix<Double> KTF;
    
    StiffTrans.Resize(row);
    StiffTrans.Init();
    TFTrans.Resize(row);
    TFTrans.Init();
    KTF.Resize(row);
    KTF.Init();

    //Construction of the big Element Transformation matrix 24x24
    for (i=0 ; i < row; i++)
      for (j=0 ; j < row ; j++)
        StiffTrans[i][j] = TransMat[i][j];

    //adding the Drilling degree of freedom,which can be defined in the beginning of the file
    const Double K = penaltyDof_;
    UInt dofspernode_=6;
    
    for( i = dofspernode_ - 1; i < row; i+= dofspernode_ )
      for( j = dofspernode_ - 1; j < row; j+= dofspernode_ )
        ElemMat[i][j] = K /*0.001 * Max*/;

    //Multiplying the stiffness Matrix with the Transformation matrix
    
    KTF = ElemMat * StiffTrans;

    //Calculates the transpose of the Transformation Matrix

    StiffTrans.Transpose(TFTrans);

    //Calculates the global stiffness matrix 

    ElemMat = TFTrans * KTF;

    // to check if the transform matrix is orthogonal
    //    KTF = TFTrans * StiffTrans;

  }

  void FlatShellInt::calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ShellCoord, SubPartType part) {

    ENTER_FCN( "linElastInt::calcBMat" );

    const UInt nrNodes  = ptelem->GetNumNodes(); //4 nodes for linear element
    const UInt spaceDim = getDimD(); //8 
    const UInt nrDofs   = getNrDofs(); //6 

    UInt actDim, actNode;

    // local shape functions derived after global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx; //derivatives of the shape functions
    Vector<Double> shapeFnc;  //shape functions

    bMat.Resize(spaceDim, nrDofs * nrNodes);// 8x24
    ptelem->GetShFncAtIp(shapeFnc,ip);
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ShellCoord);

    //Filling of the big B matrix

    if (part == MEMBRANE){
      //Membrane part
      for(actDim = 0; actDim < spaceDim-6; actDim++)
        for(actNode = 0; actNode < nrNodes; actNode++)
          bMat[actDim][actNode * nrDofs + actDim] = xiDx[actNode][actDim];

      actDim = spaceDim - 6;
  	
      for(actNode = 0; actNode < nrNodes; actNode++){
        bMat[actDim][actNode * nrDofs]     = xiDx[actNode][1];
        bMat[actDim][actNode * nrDofs + 1] = xiDx[actNode][0];
      }
    
    }

    if (part == BENDING){
      //Bending Part
      for( actDim = 3; actDim < spaceDim - 3; actDim++ )
        for(actNode = 0; actNode < nrNodes; actNode++ )
          bMat[actDim][actNode * nrDofs + actDim] = xiDx[actNode][actDim-3];

      actDim = spaceDim - 3;

      for( actNode = 0; actNode < nrNodes; actNode++ ){
        bMat[actDim][actNode * nrDofs + 3] = xiDx[actNode][1];
        bMat[actDim][actNode * nrDofs + 4] = xiDx[actNode][0];
      }
      
    }

    if (part == SHEAR) {
      //Shear Part
      for( actDim = 6; actDim < spaceDim; actDim++ )
        for( actNode = 0; actNode < nrNodes; actNode++ )
          bMat[actDim][actNode * nrDofs + 2] = xiDx[actNode][actDim-6];

      for( actDim = 6; actDim < spaceDim; actDim++ )
        for( actNode = 0; actNode < nrNodes; actNode++ )
          bMat[actDim][actNode * nrDofs + actDim - 3] = -shapeFnc[actNode];
    } 
 
    //std::cout << "The B TransMatrix is\n" << bMat << std::endl;

  }

  // calculates the summed D-matrix
  void FlatShellInt::calcDMat( Matrix<Double> &dMat){

    ENTER_FCN( "linElastInt::calcDMat" );

    const UInt SizeOfD = 8;
    const UInt nrLayers  = z_.GetSize() - 1;  
    double kappa = 5.0/6.0; 

    Matrix<Double> matMatrix; 
    
    ptMaterial->GetTensor(matMatrix,MECH_STIFFNESS_TENSOR,matDataType_,FULL);
    
    Matrix<Double> RtMat, RsMat, RtMatInv;
    Matrix<Double> A, B, D, K, Ck, Qk, Buffer;
    
    dMat.Resize(SizeOfD);
    dMat.Init();
    
    RtMat.Resize(3);
    RtMat.Init();
    RsMat.Resize(3);
    RsMat.Init();
    RtMatInv.Resize(3);
    RtMatInv.Init();
    
    Ck.Resize(3);
    Ck.Init();
    Buffer.Resize(3);
    Buffer.Init();
    Qk.Resize(3);
    Qk.Init();
        
    A.Resize(3);
    A.Init();
    B.Resize(3);
    B.Init();
    D.Resize(3);
    D.Init();
    K.Resize(2);
    K.Init();
    
    for (UInt k=1 ; k <= nrLayers ; k++)
      {
        //Construction of the Transformation Matrices Rt and Rs 
        //see equation (2.56) and (2.57) Pefort 'Finite Element Modelling of Piezoelectric Active Structures'
        //std::cout << "Layer number\n" << k << std::endl;
        //std::cout << "angle\n" << orAngle_[k] << std::endl;
	
	RtMat[0][0] = cos(orAngle_[k])*cos(orAngle_[k]); 
	RtMat[0][1] = sin(orAngle_[k])*sin(orAngle_[k]);
	RtMat[0][2] = 2*sin(orAngle_[k])*cos(orAngle_[k]);
	RtMat[1][0] = sin(orAngle_[k])*sin(orAngle_[k]);
	RtMat[1][1] = cos(orAngle_[k])*cos(orAngle_[k]);
	RtMat[1][2] = -2*sin(orAngle_[k])*cos(orAngle_[k]);
	RtMat[2][0] = -sin(orAngle_[k])*cos(orAngle_[k]);
	RtMat[2][1] = sin(orAngle_[k])*cos(orAngle_[k]);
	RtMat[2][2] = cos(orAngle_[k])*cos(orAngle_[k]) - sin(orAngle_[k])*sin(orAngle_[k]);
	
        //std::cout << "Rt Matrix is\n" << RtMat << std::endl;
	
	RsMat[0][0] = cos(orAngle_[k])*cos(orAngle_[k]); 
	RsMat[0][1] = sin(orAngle_[k])*sin(orAngle_[k]);
	RsMat[0][2] = sin(orAngle_[k])*cos(orAngle_[k]);
	RsMat[1][0] = sin(orAngle_[k])*sin(orAngle_[k]);
	RsMat[1][1] = cos(orAngle_[k])*cos(orAngle_[k]);
	RsMat[1][2] = -sin(orAngle_[k])*cos(orAngle_[k]);
	RsMat[2][0] = -2*sin(orAngle_[k])*cos(orAngle_[k]);
	RsMat[2][1] = 2*sin(orAngle_[k])*cos(orAngle_[k]);
	RsMat[2][2] = cos(orAngle_[k])*cos(orAngle_[k]) - sin(orAngle_[k])*sin(orAngle_[k]);
	
        //std::cout << "Rs Matrix is\n" << RsMat << std::endl;
        //Transpose of RsMat
	RtMat.Invert(RtMatInv);    
	
        //Matrix Ck
	Ck[0][0] = 2*matMatrix[0][1]*matMatrix[5][5]/matMatrix[0][0] + 2*matMatrix[5][5]; 
	Ck[0][1] = 2*matMatrix[0][1]*matMatrix[5][5]/matMatrix[0][0];
	Ck[0][2] = 0.0;
	Ck[1][0] = Ck[0][1];
	Ck[1][1] = Ck[0][0];
	Ck[1][2] = 0.0;
	Ck[2][0] = 0.0;
	Ck[2][1] = 0.0;
	Ck[2][2] = matMatrix[5][5];
	
        //std::cout << "The Ck Matrix is\n" << Ck << std::endl;
        //Equation 2.61 page 31
	Buffer = Ck*RsMat;
	Qk = RtMatInv*Buffer; 
	
        //std::cout << "The Qk Matrix is\n" << Qk << std::endl;
	
        //Equations (2.65) (2.66) (2.67) and ()
	A += (z_[k] - z_[k-1])*Qk;
	B += 0.5*(z_[k]*z_[k] - z_[k-1]*z_[k-1])*Qk;
	D += 1.0/3*(z_[k]*z_[k]*z_[k] - z_[k-1]*z_[k-1]*z_[k-1])*Qk;
			
	K[0][0] += kappa*(z_[k] - z_[k-1])*matMatrix[5][5];
	K[1][1] += kappa*(z_[k] - z_[k-1])*matMatrix[5][5];
	
        //std::cout << "The A Matrix is\n" << A << std::endl;
        //std::cout << "The B Matrix is\n" << B << std::endl;
        //std::cout << "The D Matrix is\n" << D << std::endl;
        //std::cout << "The K Matrix is\n" << K << std::endl;
      }
	

    //The material matrix for the laminate case is created
    //The whole D matrix is consisted of the separate matrices for membrane, bending and shear part 
    //Matrix part A
    dMat[0][0] = A[0][0]; 
    dMat[0][1] = A[0][1];
    dMat[0][2] = A[0][2];
    dMat[1][0] = A[1][0];
    dMat[1][1] = A[1][1];
    dMat[1][2] = A[1][2];
    dMat[2][0] = A[2][0];
    dMat[2][1] = A[2][1];
    dMat[2][2] = A[2][2];

    //Matrix part B1
    dMat[0][3] = B[0][0]; 
    dMat[0][4] = B[0][1];
    dMat[0][5] = B[0][2];
    dMat[1][3] = B[1][0];
    dMat[1][4] = B[1][1];
    dMat[1][5] = B[1][2];
    dMat[2][3] = B[2][0];
    dMat[2][4] = B[2][1];
    dMat[2][5] = B[2][2];

    //Matrix Part B2
    dMat[3][0] = B[0][0];
    dMat[3][1] = B[0][1];
    dMat[3][2] = B[0][2];
    dMat[4][0] = B[1][0];
    dMat[4][1] = B[1][1];
    dMat[4][2] = B[1][2];
    dMat[5][0] = B[2][0];
    dMat[5][1] = B[2][1];
    dMat[5][2] = B[2][2];
        
    //Matrix part D	
    dMat[3][3] = D[0][0]; 
    dMat[3][4] = D[0][1];
    dMat[3][5] = D[0][2];
    dMat[4][3] = D[1][0];
    dMat[4][4] = D[1][1];
    dMat[4][5] = D[1][2];
    dMat[5][3] = D[2][0];
    dMat[5][4] = D[2][1];
    dMat[5][5] = D[2][2];


    //Matrix part K	
    dMat[6][6] = K[0][0];
    dMat[6][7] = K[0][1];
    dMat[7][6] = K[1][0];
    dMat[7][7] = K[1][1];

    //std::cout << "The D Matrix is\n" << dMat << std::endl;  
  }
  
     
  // ***************
  //   Constructor
  // ***************
  FlatShellInt::FlatShellInt( BaseMaterial * matData ) 
    : BaseForm(matData), updateDMatInEveryIP_(0) {
    ENTER_FCN( "FlatShellInt::FlatShellInt" );

    name_ = "FlatShellInt";
    baseType_ = STIFFNESS;
    
    parts_.Push_back(MEMBRANE);
    parts_.Push_back(BENDING);
    parts_.Push_back(SHEAR);
    parts_.Push_back(COUPLED1);
    parts_.Push_back(COUPLED2);

    reducedPart_.Push_back(false);
    reducedPart_.Push_back(false);
    reducedPart_.Push_back(false);
    reducedPart_.Push_back(false);
    reducedPart_.Push_back(false);
    

    //Thickness of the laminas' layers
    z_.Push_back(-0.125);
    z_.Push_back(+0.125);

   orAngle_.Push_back(0.0);
   orAngle_.Push_back(0.0); 

  }

  // **************
  //   Destructor
  // **************
  FlatShellInt::~FlatShellInt() {
    ENTER_FCN( "FlatShellInt::~FlatShellInt" );
  }
  
} // end namespace CoupledField
