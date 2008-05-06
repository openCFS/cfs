// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "FlatShellInt.hh"
#define PI 3.141592654

namespace CoupledField {


  void FlatShellInt::CoordTrans( const Matrix<Double> &ptCoord, 
				 Matrix<Double> &TransMat, 
				 Matrix<Double> &ShellCoord )
  {
   
    //std::cout << "FlatShellInt::CoordTrans\n" << std::endl;
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrNodes = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrDofs   = getNrDofs();
    const UInt lambdaDim = 3; //dimension of the submatrix lambda, part of the transformation matrix
    const UInt TMatDim = getNrDofs(); //dimension of the nodal transformation matrix
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
    for ( i=0; i < lambdaDim; i++ )
      for ( j=0; j< lambdaDim; j++ )
        TMat[i+k][j+k]=lambdaMat[i][j];

    k = 3;
    for ( i=0; i < TMatDim-lambdaDim; i++ )
      for ( j=0; j< TMatDim-lambdaDim; j++ )
        TMat[i+k][j+k]=lambdaMat[i][j];
     
    //The transformation Element matrix which is 24x24
    for( k=0;k < ElementTransMatDim ; k+=nrDofs)
      for ( i=0; i < TMatDim; i++ )
        for ( j=0; j< TMatDim; j++ )
          TransMat[i+k][j+k]=TMat[i][j];
  
  }

  void FlatShellInt::LocaltoGlob( Matrix<Double> &ElemMat, const Matrix<Double> &TransMat )
  {

    //std::cout << "FlatShellInt::LocaltoGlob\n" << std::endl;

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
    UInt dofspernode_ = getNrDofs();
    
    if ( hasDrillDof_ == true ) {
      for( i = dofspernode_ - 1; i < row;
           i+= dofspernode_ )
        for( j = dofspernode_ - 1; j < row; j+= dofspernode_ )
          ElemMat[i][j] = K /*0.001 * Max*/;
    }

    //Multiplying the stiffness Matrix with the Transformation matrix
    
    KTF = ElemMat * StiffTrans;

    //Calculates the transpose of the Transformation Matrix

    StiffTrans.Transpose(TFTrans);

    //Calculates the global stiffness matrix 

    ElemMat = TFTrans * KTF;

    // to check if the transform matrix is orthogonal
    //    KTF = TFTrans * StiffTrans;
    //std::cerr << "ElemMat = \n" << ElemMat << std::endl;
  }

void FlatShellInt::LocaltoGlobPiezo( Matrix<Double> &ElemMat, const Matrix<Double> &TransMat )
  {

    //std::cout << "FlatShellInt::LocaltoGlob\n" << std::endl;

    int i, j;
    const Integer row = ElemMat.GetSizeRow(); //equals to nrDofs*nrNodes for linear quadrilateral is 24
    const Integer col = ElemMat.GetSizeCol(); //equals to the number of piezoelectric Layers


    Matrix<Double> TransMatInv, ElemMatLocal;
        
    TransMatInv.Resize(row);
    TransMatInv.Init();
    
    ElemMatLocal.Resize(row,col);
    ElemMatLocal.Init();
    
    //adding the Drilling degree of freedom,which can be defined in the beginning of the file
    //const Double K = penaltyDof_;
    //UInt dofspernode_=6;
    
    //for( i = dofspernode_ - 1; i < row; i+= dofspernode_ )
    //  for( j = dofspernode_ - 1; j < row; j+= dofspernode_ )
    //    ElemMat[i][j] = K /*0.001 * Max*/;
    
    //ElemMatLocal is initially equal to ElemMat
    ElemMatLocal = ElemMat;
    
    //Calculates the inverse of the Transformation Matrix
    TransMat.Transpose(TransMatInv);

    //Calculates the global stiffness matrix 

    ElemMat = TransMatInv * ElemMatLocal;

    // to check if the transform matrix is orthogonal
    //    KTF = TFTrans * StiffTrans;

  }   

   // *************************************
  //   Constructor for Composite Material
  // *************************************
  FlatShellInt::FlatShellInt( Composite * composite, bool hasDrillDof ) 
    : BaseForm(NULL) {

    name_ = "FlatShellInt";
    hasDrillDof_ = hasDrillDof;

    // Set flag
    isComposite_ = true;
    isOrthotropic_ = true;
    
    // iterate over all composite materials anc check
    // if all of them are orthotropic
    for( UInt i = 0; i< composite->materials.GetSize(); i++ ){
      if( composite->materials[i]->GetSymmetryType() !=
          BaseMaterial::ORTHOTROPIC ) {
        isOrthotropic_ = false;
        break;
      }
    }
    
    composite_=composite;
  }
    
  // **********************************
  //   Constructor for Normal Material
  // **********************************
  
  FlatShellInt::FlatShellInt( BaseMaterial * matData, bool hasDrillDof ) 
    : BaseForm(matData) {
    
    hasDrillDof_ = hasDrillDof;

    // Set flag
    isComposite_ = false;
  } 
  
  // **************
  //   Destructor
  // **************
  FlatShellInt::~FlatShellInt() {
  }
  
} // end namespace CoupledField
