#include <iostream>
#include <fstream>

#include "FlatShellPiezoInt.hh"
#include "FlatShellStiffInt.hh"

namespace CoupledField {

  void FlatShellPiezoInt::CalcElementMatrix(  Matrix<Double>& elemMat,
                                              EntityIterator& ent1, 
                                              EntityIterator& ent2 ) {

    ENTER_FCN( "FlatShellPiezoInt::CalcElementMatrix" );
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    ptelem->SetAnsatzFct( ansatzFct1_ );    
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    

    // pass ansatz fct to stiffness integrator
    ptr_StiffInt_->SetAnsatzFct( ansatzFct1_, ansatzFct2_ );
    ptr_StiffInt_->ExtractElemInfo( ent1 );
    //    ptr_StiffInt_->it2_ = ent2;

    const UInt nrDofs   = getNrDofs();  
    //Number of Integration points
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    //Number of integration weights
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    //std::cout << "The Number of Dof is\n" << nrDofs << std::endl;
    const UInt nrLayers  = composite_->thickness.GetSize();
    //std::cout << "The Number of Layers is\n" << nrLayers << std::endl;
    double jacDet=0.0;
    
    Matrix<Double> bMat; //B matrix
    Matrix<Double> eMat; //E matrix
    Matrix<Double> bTrans; //Transpose of B
    Matrix<Double> partElemMat; //Partial Element Matrix
    Matrix<Double> TransMat; //Transformation matrix
    Matrix<Double> ShellCoord; //Shell Coordinates
    
    //Element matrix in the case of Flat Shells with linear elements is 24x24
    elemMat.Resize(numFncs * nrDofs, nrLayers);
    elemMat.Init();
       
    //B matrix in the case of flat shells 8x24 
    bMat.Resize( getDimD(), numFncs * nrDofs );
    bMat.Init();
     
    //The transpose of B in the case of flat shells is 24x8
    bTrans.Resize( numFncs * nrDofs , getDimD() );
    bTrans.Init();
      
    //The partial Element Matrix with the size of 24x24 
    partElemMat.Resize(numFncs * nrDofs, nrLayers);
    partElemMat.Init();
    
         
    // 1. rotate 3D-coordintates to 2D (x-y-reference plane) on the element level,4 nodes x 6 coordinates = 24
    //and stor the full Element transformation matrix in TransMat
	
    FlatShellInt::CoordTrans(ptCoord_, TransMat, ShellCoord );
    
    // 2. calculate partial B- and E-Matrices and assemble into big one
	
      
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
          
      //std::cout << "Integration point\n" << actIntPt << std::endl;
      //FlatShellPiezoInt::calcBMat(bMat, actIntPt, ShellCoord);
      ptr_StiffInt_->calcBMat(bMat, actIntPt, ptelem, ShellCoord, FlatShellStiffInt::PIEZO);
      //Find the transpose of the bMat and Store it in the bTrans   
      bMat.Transpose(bTrans);
      //Calculate E matrix
      FlatShellPiezoInt::calcEMat(eMat);
      //3. Calculate B^t*E and store it in the partElemMat
      partElemMat = bTrans * eMat;
	  
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ShellCoord,
                                            ent1.GetElem() );
      // Perform a safety check
      if ( jacDet < 0.0 ) {
        (*error) << "FlatShellStiffInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!";
        Error( __FILE__, __LINE__ );
      }
	  
      // Sum up part element matrices
      elemMat += partElemMat * jacDet * intWeights[actIntPt-1];

    }
     
      
    //std::cout << "The Element Matrix is :"<< elemMat << std::endl << std::endl;

    std::cerr << "Size of piezo-coupling matrix before back-transformation: " 
              << elemMat.GetSizeRow() << " x " << elemMat.GetSizeCol() << std::endl;
    
    //Do the back-Rotation and return
    FlatShellInt::LocaltoGlobPiezo(elemMat,TransMat );
  }
  
  /*void FlatShellPiezoInt::calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ShellCoord) {

  ENTER_FCN( "FlatShellPiezoInt::calcBMat" );

  const UInt numFncs  = ptelem->GetNumNodes(); //4 nodes for linear element
  const UInt spaceDim = getDimD(); //8 
  const UInt nrDofs   = getNrDofs(); //6 

  UInt actDim, actNode;

  // local shape functions derived after global coords
  // (format: numFncs x spaceDim)
  Matrix<Double> xiDx; //derivatives of the shape functions
  Vector<Double> shapeFnc;  //shape functions

  bMat.Resize(spaceDim, nrDofs * numFncs);// 8x24
  ptelem->GetShFncAtIp(shapeFnc,ip);
  ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ShellCoord);
    
  //Filling of the big B matrix

  //Membrane part
  for(actDim = 0; actDim < spaceDim-6; actDim++)
  for(actNode = 0; actNode < numFncs; actNode++)
  bMat[actDim][actNode * nrDofs + actDim] = xiDx[actNode][actDim];

  actDim = spaceDim - 6;
  	
  for(actNode = 0; actNode < numFncs; actNode++){
  bMat[actDim][actNode * nrDofs]     = xiDx[actNode][1];
  bMat[actDim][actNode * nrDofs + 1] = xiDx[actNode][0];
  }
    
  //Bending Part
  for( actDim = 3; actDim < spaceDim - 3; actDim++ )
  for(actNode = 0; actNode < numFncs; actNode++ )
  bMat[actDim][actNode * nrDofs + actDim] = xiDx[actNode][actDim-3];

  actDim = spaceDim - 3;

  for( actNode = 0; actNode < numFncs; actNode++ ){
  bMat[actDim][actNode * nrDofs + 3] = xiDx[actNode][1];
  bMat[actDim][actNode * nrDofs + 4] = xiDx[actNode][0];
  }
      
  //Shear Part
  for( actDim = 6; actDim < spaceDim; actDim++ )
  for( actNode = 0; actNode < numFncs; actNode++ )
  bMat[actDim][actNode * nrDofs + 2] = xiDx[actNode][actDim-6];

  for( actDim = 6; actDim < spaceDim; actDim++ )
  for( actNode = 0; actNode < numFncs; actNode++ )
  bMat[actDim][actNode * nrDofs + actDim - 3] = -shapeFnc[actNode];
    
 
  //std::cout << "The B TransMatrix is\n" << bMat << std::endl;

  }*/
  
  void FlatShellPiezoInt::calcEMat( Matrix<Double> &eMat ){

    ENTER_FCN( "FlatShellPiezoInt::calcEMat" );

    const UInt SizeOfD = 8;
    const UInt nrLayers  = composite_->thickness.GetSize();
    Double thickness_all = 0.0;
    //Piezoelectric coefficients e31 and e32
    Double z_mk =0.0; 
    
    Matrix<Double> RtMat, RsMat, RtMatInv,eTensor; 
    Vector<Double> eVec, Epsillon_k;
    
    eMat.Resize(SizeOfD, nrLayers);
    eMat.Init();
    RtMat.Resize(3);
    RtMat.Init();
    RsMat.Resize(3);
    RsMat.Init();
    RtMatInv.Resize(3);
    RtMatInv.Init();
    Epsillon_k.Resize(3);
    Epsillon_k.Init();
    eVec.Resize(3);
    eVec.Init();
			  
    //Initialisation of the composite structure variables
    z_.Resize(nrLayers + 1);
    z_.Init();
    orAngle_.Resize(nrLayers + 1);
    orAngle_.Init();

    for (UInt k=0; k < nrLayers; k++) {
      thickness_all += composite_->thickness[k];
    }
    z_[0]= - (thickness_all/2.0);
    orAngle_[0]= 0.0;
    
    for(UInt i=1; i <= nrLayers; i++) {
      z_[i] = z_[i-1] + composite_->thickness[i-1];
      orAngle_[i] = composite_->orientation[i-1]*PI/180;
    }
     
    //std::cout << "z coordinates are \n" << z_ << std::endl;
    //std::cout << "orientations are are \n" << orAngle_ << std::endl;
    
                         
    for (UInt k=1 ; k <= nrLayers ; k++) {
    	
      composite_->materials[k-1]->GetTensor(eTensor,PIEZO_TENSOR,REAL,FULL);
   	
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
        
      //z_mk
      z_mk = (z_[k-1] + z_[k])/2.0;
      //std::cout << "z_mk is\n" << z_mk << std::endl;
	
      //Inverse of RsMat
      RtMat.Invert(RtMatInv);
      //eVec initialization
      eVec[0]=eTensor[2][0];
      eVec[1]=eTensor[2][1];
      eVec[2]=0.0;
      //Epsillon_k
      Epsillon_k = RtMatInv*eVec;
      //eMat Formation
      eMat[0][k-1] = Epsillon_k[0];
      eMat[1][k-1] = Epsillon_k[1];
      eMat[2][k-1] = Epsillon_k[2];
	
      eMat[3][k-1] = Epsillon_k[0]*z_mk;
      eMat[4][k-1] = Epsillon_k[1]*z_mk;
      eMat[5][k-1] = Epsillon_k[2]*z_mk;
        
      eMat[6][k-1] = 0.0;
      eMat[7][k-1] = 0.0;
    }
    //std::cout << "The E matrix is\n" << eMat << std::endl;
    
  }
  // *************************************
  //   Constructor for Composite Material
  // *************************************
  FlatShellPiezoInt::FlatShellPiezoInt( Composite * composite ) 
    : FlatShellInt(composite)  {
    ENTER_FCN( "FlatShellPiezoInt::FlatShellPiezoInt" );
    
    name_ = "FlatShellPiezoInt";

    baseType_ = STIFFNESS;
    ptr_StiffInt_ = new FlatShellStiffInt(composite);
    
  }

  // **************
  //   Destructor
  // **************
  FlatShellPiezoInt::~FlatShellPiezoInt() {
    ENTER_FCN( "FlatShellPiezoInt::~FlatShellPiezoInt" );
  }
  
} // end namespace CoupledField
