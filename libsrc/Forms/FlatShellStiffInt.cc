#include <iostream>
#include <fstream>

#include "FlatShellStiffInt.hh"
#define PI 3.141592654

namespace CoupledField {

  void FlatShellStiffInt::CalcElementMatrix(Matrix<Double>& elemMat,
					    EntityIterator& ent1, 
					    EntityIterator& ent2 ) {
    
    
    ENTER_FCN( "FlatShellStiffInt::CalcElementMatrix" );
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    //std::cout << "FlatShellStiffInt::CalcElementMatrix\n" << std::endl;
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
    
    //std::cout << "The Number of nodes is\n" << numFncs << std::endl;   
    const UInt nrDofs   = getNrDofs();  
    //std::cout << "The Number of Dof is\n" << nrDofs << std::endl;
    double jacDet=0.0;
        

    Matrix<Double> bMat; //B matrix
    Matrix<Double> dMat; //D matrix
    Matrix<Double> dB; //D*B matrix
    Matrix<Double> bTrans; //Transpose of B
    Matrix<Double> partElemMat; //Partial Element Matrix
    Matrix<Double> TransMat; //Transformation matrix
    Matrix<Double> ShellCoord; //Shell Coordinates

    //Element matrix in the case of Flat Shells with linear elements is 24x24
    elemMat.Resize(numFncs * nrDofs);
    elemMat.Init();
       
    //B matrix in the case of flat shells 8x24 
    bMat.Resize( getDimD(), numFncs * nrDofs );
    bMat.Init();
    
    //D matrix in the case of flat shells 8x8
    dMat.Resize(getDimD()); 
    dMat.Init();
      
    //D*B matrix in the case of flat shells 8x24
    dB.Resize( getDimD(), numFncs * nrDofs );
    dB.Init();
     
    //The transpose of B in the case of flat shells is 24x8
    bTrans.Resize( numFncs * nrDofs , getDimD() );
    bTrans.Init();
      
    //The partial Element Matrix with the size of 24x24 
    partElemMat.Resize(numFncs * nrDofs);
    partElemMat.Init();
         
    // 1. rotate 3D-coordintates to 2D (x-y-reference plane) on the element level,4 nodes x 6 coordinates = 24
    //and stor the full Element transformation matrix in TransMat
	
    FlatShellInt::CoordTrans(ptCoord_, TransMat, ShellCoord );
    
    // 2. calculate partial B- and D-Matrices and assemble into big one
	
    // If the material parameters are constant within the element
    // we can compute the D matrix once and for all
    if ( isComposite_ == true ) {
      //std::cout << "calcDMatComposite\n" << std::endl;
      FlatShellStiffInt::calcDMatComposite( dMat );
      
    }
    else {
      //std::cout << "calcDMat\n" << std::endl;
      FlatShellStiffInt::calcDMat(dMat);
       
    }

    for (UInt i =0 ; i<parts_.GetSize(); i++ ) {
	
      //Set the reduced integration          
      if (reducedPart_[i] == true ) {
	//std::cout << "Set Reduced Integration\n" << std::endl;  
        ptelem->SetReducedIntegration();
      }
      const UInt nrIntPts = ptelem->GetNumIntPoints();
      const Vector<Double> & intWeights = ptelem->GetIntWeights();
      //Main Cycle	
      for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {
	//std::cout << "Integration point\n" << actIntPt << std::endl;
        //The calcBMat function calculates the B matrix either membrane shear or bending part
        if (parts_[i] == MEMBRANE || parts_[i] == BENDING || parts_[i] == SHEAR) {
	  //std::cout <<"MEMBRANE BENDING OR SHEAR\n" << std::endl;
          //Calculates the B matrix 
          FlatShellStiffInt::calcBMat(bMat, actIntPt, ptelem, ShellCoord, parts_[i]);
          //Calculates the product D*B
          dB=dMat*bMat; 
          //Find the transpose of the bMat and Store it in the bTrans   
          bMat.Transpose(bTrans);
          //3. Calculate B^t*B*D and store it in the partElemMat
          partElemMat = bTrans * dB;
        
        } else if (parts_[i] == COUPLED1) {
	  //std::cout <<"COUPLED 1\n" << std::endl;
          //Calculates B bending part
          FlatShellStiffInt::calcBMat(bMat, actIntPt, ptelem, ShellCoord, parts_[1]);
          //Calculates the product D*B
          dB=dMat*bMat;
          //Calculates B membrane part
          FlatShellStiffInt::calcBMat(bMat, actIntPt, ptelem, ShellCoord, parts_[0]); 
          //Find the transpose of the bMat and Store it in the bTrans   
          bMat.Transpose(bTrans);

          partElemMat = bTrans * dB;
		
        } else if (parts_[i] == COUPLED2) {
          //std::cout <<"COUPLED 2\n" << std::endl;
          //Calculates B membrane part
          FlatShellStiffInt::calcBMat(bMat, actIntPt, ptelem, ShellCoord, parts_[0]);
          //Calculates the product D*B
          //Thickness of the laminas' layers
          dB=dMat*bMat;
          //Calculates B bending part
          FlatShellStiffInt::calcBMat(bMat, actIntPt, ptelem, ShellCoord, parts_[1]); 
          //Find the transpose of the bMat and Store it in the bTrans   
          bMat.Transpose(bTrans);
          //3. Calculate B^t*B*D and store it in the partElemMat
          partElemMat = bTrans * dB;
        }
	
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
     
      if (reducedPart_[i] == true ) 
        ptelem->SetStandardIntegration();
    }
    //std::cout << "The Element Matrix is :"<< elemMat << std::endl;
	
    // 4. Do the back-Rotation and return
    FlatShellInt::LocaltoGlob(elemMat,TransMat );
        
  }
  void FlatShellStiffInt::calcBMat(Matrix<Double> & bMat, Integer ip, BaseFE* elem,
                                   Matrix<Double> & ShellCoord, SubPartType part) {

    ENTER_FCN( "FlatShellStiffInt::calcBMat" );

    const UInt numFncs = elem->GetNumFncs( ansatzFct1_ );
    elem->SetAnsatzFct( ansatzFct1_ );
    const UInt spaceDim = getDimD(); //8 
    const UInt nrDofs   = getNrDofs(); //6 

    UInt actDim, actNode;

    // local shape functions derived after global coords
    // (format: numFncs x spaceDim)
    Matrix<Double> xiDx; //derivatives of the shape functions
    Vector<Double> shapeFnc;  //shape functions

    bMat.Resize(spaceDim, nrDofs * numFncs);// 8x24
    bMat.Init();
    elem->GetShFncAtIp(shapeFnc,ip, it1_.GetElem());
    elem->GetGlobDerivShFncAtIp(xiDx, ip, ShellCoord, it1_.GetElem());
    
    //Filling of the big B matrix

    if (part == MEMBRANE || part == PIEZO){
      //Membrane part
      for(actDim = 0; actDim < spaceDim-6; actDim++)
        for(actNode = 0; actNode < numFncs; actNode++)
          bMat[actDim][actNode * nrDofs + actDim] = xiDx[actNode][actDim];

      actDim = spaceDim - 6;
  	
      for(actNode = 0; actNode < numFncs; actNode++){
        bMat[actDim][actNode * nrDofs]     = xiDx[actNode][1];
        bMat[actDim][actNode * nrDofs + 1] = xiDx[actNode][0];
      }
    
    }

    if (part == BENDING || part == PIEZO){
      //Bending Part
      for( actDim = 3; actDim < spaceDim - 3; actDim++ )
        for(actNode = 0; actNode < numFncs; actNode++ )
          bMat[actDim][actNode * nrDofs + actDim] = xiDx[actNode][actDim-3];

      actDim = spaceDim - 3;

      for( actNode = 0; actNode < numFncs; actNode++ ){
        bMat[actDim][actNode * nrDofs + 3] = xiDx[actNode][1];
        bMat[actDim][actNode * nrDofs + 4] = xiDx[actNode][0];
      }
      
    }

    if (part == SHEAR || part == PIEZO) {
      //Shear Part
      for( actDim = 6; actDim < spaceDim; actDim++ )
        for( actNode = 0; actNode < numFncs; actNode++ )
          bMat[actDim][actNode * nrDofs + 2] = xiDx[actNode][actDim-6];

      for( actDim = 6; actDim < spaceDim; actDim++ )
        for( actNode = 0; actNode < numFncs; actNode++ )
          bMat[actDim][actNode * nrDofs + actDim - 3] = -shapeFnc[actNode];
    } 
 
    //std::cout << "The B TransMatrix is\n" << bMat << std::endl;

  }
 
  // calculates the summed D-matrix
  void FlatShellStiffInt::calcDMatComposite( Matrix<Double> &dMat ){

    ENTER_FCN( "FlatShellIStiffnt::calcDMatComposite" );

    const UInt SizeOfD = 8;
    const UInt nrLayers  = composite_->thickness.GetSize();
    Double kappa = 5.0/6.0;
    Double thickness_all = 0.0; 
    
    Matrix<Double> A, B, D, K, Q, Q_, Qs, Qs_, Ck, Qk, Cs, Buffer, Buffer_shear;
    Matrix<Double> matMatrix;
    Matrix<Double> RtMat, RsMat, RkMat, RtMatInv, RkMatInv;
    Double Ex, Ey, Ez, NUxy, NUyz, NUxz, Gyz, Gzx, Gxy,  NUyx, NUzy, NUzx;
    dMat.Resize(SizeOfD);
    dMat.Init();
    
    
    //Initialisation of the composite structure variables
    z_.Resize(nrLayers + 1);
    z_.Init();
    orAngle_.Resize(nrLayers + 1);
    orAngle_.Init();

    for (UInt k=0; k < nrLayers; k++) {
      thickness_all += composite_->thickness[k];
    }

    z_[0]= - (thickness_all/2.0);
    //z_[0] = composite_->zStart;
    orAngle_[0]= 0.0;
    
    for(UInt i=1; i <= nrLayers; i++) {
      z_[i] = z_[i-1] + composite_->thickness[i-1];
      orAngle_[i] = composite_->orientation[i-1]*PI/180;
    }
     
    //std::cout << "z coordinates are \n" << z_ << std::endl;
    //std::cout << "orientations are are \n" << orAngle_ << std::endl;
        
    //Check for orthotropic material
    if (isOrthotropic_ == true) {
      Q_.Resize(3);
      Q_.Init();
      Q.Resize(3);
      Q.Init();
      Qs_.Resize(2);
      Qs_.Init();
      Qs.Resize(2);
      Qs.Init();
                    
      A.Resize(3);
      A.Init();
      B.Resize(3);
      B.Init();
      D.Resize(3);
      D.Init();
      K.Resize(2);
      K.Init();
                    
      for (UInt k=1 ; k <= nrLayers ; k++){
	//Reading orthotropic material parametres
	composite_->materials[k-1]->GetScalar(Ex,MECH_EMODULUS_X,REAL);
	composite_->materials[k-1]->GetScalar(Ey,MECH_EMODULUS_Y,REAL);
	composite_->materials[k-1]->GetScalar(Ez,MECH_EMODULUS_Z,REAL);
	composite_->materials[k-1]->GetScalar(NUxy,MECH_POISSON_XY,REAL);
	composite_->materials[k-1]->GetScalar(NUyz,MECH_POISSON_YZ,REAL);
	composite_->materials[k-1]->GetScalar(NUxz,MECH_POISSON_XZ,REAL);
	composite_->materials[k-1]->GetScalar(Gyz,MECH_GMODULUS_YZ,REAL);
	composite_->materials[k-1]->GetScalar(Gzx,MECH_GMODULUS_ZX,REAL);
	composite_->materials[k-1]->GetScalar(Gxy,MECH_GMODULUS_XY,REAL);
                        
	//see Finite Element Analysis of Composite Laminates O.O.Ochoa and J.N.Reddy page 10 and 12
        NUyx=(Ey/Ex)*NUxy;
        NUzy=(Ez/Ey)*NUyz;
        NUzx=(Ez/Ex)*NUxz;
	Q_[0][0] = Ex/(1 - NUxy*NUyx); 
	Q_[0][1] = NUxy*Ey/(1 - NUxy*NUyx);
	Q_[0][2] = 0.0;
	Q_[1][0] = NUxy*Ey/(1 - NUxy*NUyx);
	Q_[1][1] = Ey/(1 - NUxy*NUyx);
	Q_[1][2] = 0.0;
	Q_[2][0] = 0.0;
	Q_[2][1] = 0.0;
	Q_[2][2] = Gxy;
	//std::cout << "Q_ Matrix is\n" << Q_ << std::endl;
        
        //std::cout << "dMat = \n" << dMat << std::endl;
	//Matrix Qs_ taken from Reddy page 10
	Qs_[0][0] = Gyz; //Gyz; 
	Qs_[0][1] = 0.0;
	Qs_[1][0] = 0.0;
	Qs_[1][1] = Gzx; //Gzx;
	//std::cout << "Qs_ Matrix is\n" << Qs_ << std::endl;
	
        
	//Matrix Q taken from Reddy page 22
	Q[0][0] = Q_[0][0]*cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]) + 2.0*(Q_[0][1] + 2.0*Q_[2][2])*sin(orAngle_[k])*sin(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]) + Q_[1][1]*sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k]);
 
	Q[0][1] = (Q_[0][0] + Q_[1][1] - 4.0*Q_[2][2])*sin(orAngle_[k])*sin(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]) + Q_[0][1]*(sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k]) + cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]));

	Q[0][2] = (Q_[0][0] - Q_[0][1] - 2.0*Q_[2][2])*sin(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]) + (Q_[0][1] - Q_[1][1] + 2.0*Q_[2][2])*sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k])*cos(orAngle_[k]); 

	Q[1][0] = Q[0][1];

	Q[1][1] = Q_[0][0]*sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k]) + 2.0*(Q_[0][1] + 2.0*Q_[2][2])*sin(orAngle_[k])*sin(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]) + Q_[1][1]*cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]);

	Q[1][2] = (Q_[0][0] - Q_[0][1] - 2.0*Q_[2][2])*sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k])*cos(orAngle_[k]) + (Q_[0][1] - Q_[1][1] + 2.0*Q_[2][2])*sin(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]);

	Q[2][0] = Q[0][2];

	Q[2][1] = Q[1][2];

	Q[2][2] = (Q_[0][0] + Q_[1][1] - 2.0*Q_[0][1] - 2.0*Q_[2][2]) *sin(orAngle_[k])*sin(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]) + Q_[2][2]*(sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k])*sin(orAngle_[k]) + cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k])*cos(orAngle_[k]) );
	
	//std::cout << "Q matrix\n" << Q << std::endl;

	//Matrix Qs
        Qs[0][0] = Qs_[1][1]*cos(orAngle_[k])*cos(orAngle_[k]) +Qs_[0][0]*sin(orAngle_[k])*sin(orAngle_[k]);
	Qs[0][1] = (Qs_[1][1] - Qs_[0][0])*cos(orAngle_[k])*sin(orAngle_[k]);
	Qs[1][0] = (Qs_[1][1] - Qs_[0][0])*cos(orAngle_[k])*sin(orAngle_[k]);
	Qs[1][1] = Qs_[0][0]*cos(orAngle_[k])*cos(orAngle_[k]) +Qs_[1][1]*sin(orAngle_[k])*sin(orAngle_[k]);
	//std::cout << "Qs Matrix is\n" << Qs << std::endl;
	                    
	//Equations (2.65) (2.66) (2.67) and ()
	A += (z_[k] - z_[k-1])*Q;
	B += 0.5*(z_[k]*z_[k] - z_[k-1]*z_[k-1])*Q;
	D += 1.0/3.0*(z_[k]*z_[k]*z_[k] - z_[k-1]*z_[k-1]*z_[k-1])*Q;
	K += kappa*(z_[k] - z_[k-1])*Qs;
	// std::cout << "The A Matrix is\n" << A << std::endl;
	// std::cout << "The B Matrix is\n" << B << std::endl;
	// std::cout << "The D Matrix is\n" << D << std::endl;
	// std::cout << "The K Matrix is\n" << K << std::endl;
      }
    }
    else {
                             
      dMat.Resize(SizeOfD);
      dMat.Init();
      RtMat.Resize(3);
      RtMat.Init();
      RsMat.Resize(3);
      RsMat.Init();
      RkMat.Resize(2);
      RkMat.Init();
      RkMatInv.Resize(2);
      RkMatInv.Init();
      RtMatInv.Resize(3);
      RtMatInv.Init();
    
      Ck.Resize(3);
      Ck.Init();
      Cs.Resize(2);
      Cs.Init();
      Buffer.Resize(3);
      Buffer.Init();
      Buffer_shear.Resize(2);
      Buffer_shear.Init();
      Qk.Resize(3);
      Qk.Init();
      Qs.Resize(2);
      Qs.Init();
                          
      A.Resize(3);
      A.Init();
      B.Resize(3);
      B.Init();
      D.Resize(3);
      D.Init();
      K.Resize(2);
      K.Init();
                          
      for (UInt k=1 ; k <= nrLayers ; k++) {
	//Construction of the Transformation Matrices Rt and Rs 
	//see equation (2.56) and (2.57) Pefort 'Finite Element Modelling of Piezoelectric Active Structures'
	//std::cout << "Layer number\n" << k << std::endl;
	
	//Extracting material matrix
	composite_->materials[k-1]->GetTensor(matMatrix,MECH_STIFFNESS_TENSOR,REAL,FULL);
	//std::cout << "The Material Matrix is\n" << matMatrix << std::endl;
                              
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
	//See Finite Element Analysis of Composite Laminates, O.O.Ochoa and J.N.Reddy at page 12
	RkMat[0][0] = cos(orAngle_[k]);
	RkMat[0][1] = -sin(orAngle_[k]);
	RkMat[1][0] = sin(orAngle_[k]);
	RkMat[1][1] = cos(orAngle_[k]);
	//std::cout << "Rk Matrix is\n" << RkMat << std::endl;
	//Transpose of RsMat
	RtMat.Invert(RtMatInv);
	//Transpose of RkMat
	RkMat.Invert(RkMatInv);         
	//Matrix Ck
	Ck[0][0] = 2.0*matMatrix[0][1]*matMatrix[5][5]/matMatrix[0][0] + 2.0*matMatrix[5][5]; 
	Ck[0][1] = 2.0*matMatrix[0][1]*matMatrix[5][5]/matMatrix[0][0];
	Ck[0][2] = 0.0;
	Ck[1][0] = Ck[0][1];
	Ck[1][1] = Ck[0][0];
	Ck[1][2] = 0.0;
	Ck[2][0] = 0.0;
	Ck[2][1] = 0.0;
	Ck[2][2] = matMatrix[5][5];
	//Matrix Cs
	Cs[0][0] = matMatrix[5][5];
	Cs[0][1] = 0.0;
	Cs[1][0] = 0.0;
	Cs[1][1] = matMatrix[5][5];
	//std::cout << "The Ck Matrix is\n" << Ck << std::endl;
	//std::cout << "The Cs Matrix is\n" << Cs << std::endl;
	//Equation 2.61 page 31
	Buffer = Ck*RsMat;
	Qk = RtMatInv*Buffer;
	Buffer_shear = Cs*RkMat;
	Qs =RkMatInv*Buffer_shear;
	//std::cout << "The Qk Matrix is\n" << Qk << std::endl;
	//std::cout << "The Qs Matrix is\n" << Qs << std::endl;
	//Equations (2.65) (2.66) (2.67) and ()
	A += (z_[k] - z_[k-1])*Qk;
	B += 0.5*(z_[k]*z_[k] - z_[k-1]*z_[k-1])*Qk;
	D += 1.0/3.0*(z_[k]*z_[k]*z_[k] - z_[k-1]*z_[k-1]*z_[k-1])*Qk;
	K += kappa*(z_[k] - z_[k-1])*Qs;
	//std::cout << "The A Matrix is\n" << A << std::endl;
	//std::cout << "The B Matrix is\n" << B << std::endl;
	//std::cout << "The D Matrix is\n" << D << std::endl;
	//std::cout << "The K Matrix is\n" << K << std::endl;
      }
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

    //std::cerr << "The D Composite  Matrix is\n" << dMat << std::endl;  
  }
  
  void FlatShellStiffInt::calcDMat( Matrix<Double> &dMat){

    ENTER_FCN( "FlatShellStiffInt::calcDMat" );

    const UInt nrElems2d = 8;
    double t = thickness_;
    double k = 5.0/6.0; 

    Matrix<Double> matMatrix;
    ptMaterial->GetTensor(matMatrix,MECH_STIFFNESS_TENSOR,REAL,FULL);
    //Matrix<Double> const  & matMatrix = *( ptMaterial->GetMatrix());
    
    //std::cout << "The Material Matrix is\n" << matMatrix << std::endl;
    
    dMat.Resize(nrElems2d);
    dMat.Init();

    //The whole D matrix is consisted of the separate matrices for membrane, bending and shear part 

    dMat[0][0] = t*2.0*(matMatrix[0][1]*matMatrix[5][5]/matMatrix[0][0]+matMatrix[5][5]); 
    dMat[0][1] = t*2.0*(matMatrix[0][1]*matMatrix[5][5]/matMatrix[0][0]);
    dMat[0][2] = 0.0;
    dMat[1][0] = dMat[0][1];
    dMat[1][1] = dMat[0][0];
    dMat[1][2] = 0.0;
    dMat[2][0] = 0.0;
    dMat[2][1] = 0.0;
    dMat[2][2] = t*matMatrix[5][5];
	
    dMat[3][3] = (t*t*t)/6.0*((matMatrix[0][1]*matMatrix[5][5])/matMatrix[0][0]+matMatrix[5][5]); 
    dMat[3][4] = (t*t*t)/6.0*(matMatrix[0][1]*matMatrix[5][5]/matMatrix[0][0]);
    dMat[3][5] = 0.0;
    dMat[4][3] = dMat[3][4];
    dMat[4][4] = dMat[3][3];
    dMat[4][5] = 0.0;
    dMat[5][3] = 0.0;
    dMat[5][4] = 0.0;
    dMat[5][5] = (t*t*t)/12.0*matMatrix[5][5];
	
    dMat[6][6] = k*t*matMatrix[5][5];
    dMat[6][7] = 0.0;
    dMat[7][6] = 0.0;
    dMat[7][7] = k*t*matMatrix[5][5];

    //std::cout << "The D TransMatrix is\n" << dMat << std::endl;  
  }
  
  // *************************************
  //   Constructor for Composite Material
  // *************************************
  FlatShellStiffInt::FlatShellStiffInt( Composite * composite ) 
    : FlatShellInt(composite) {
    ENTER_FCN( "FlatShellStiffInt::FlatShellStiffInt" );

    name_ = "FlatShellStiffInt";

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
              
  }
  // **********************************
  //   Constructor for Normal Material
  // **********************************
  
  FlatShellStiffInt::FlatShellStiffInt( BaseMaterial * matData ) 
    : FlatShellInt(matData) {
    ENTER_FCN( "FlatShellStiffInt::FlatShellStiffInt" );

    name_ = "FlatShellStiffInt";
    
    baseType_ = STIFFNESS;
   
    parts_.Push_back(MEMBRANE);
    parts_.Push_back(BENDING);
    parts_.Push_back(SHEAR);
    
    reducedPart_.Push_back(false);
    reducedPart_.Push_back(false);
    reducedPart_.Push_back(false);
   
  }

  // **************
  //   Destructor
  // **************
  FlatShellStiffInt::~FlatShellStiffInt() {
    ENTER_FCN( "FlatShellStiffInt::~FlatShellStiffInt" );
  }
  
} // end namespace CoupledField
