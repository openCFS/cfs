#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "linPiezoInt.hh"

namespace CoupledField {

  // ========================================================================
  // ========================= linPiezoInt - Part ===========================
  // ========================================================================
 

  // determine the matrix B of the BDB operator
  void linPiezoInt::calcBMat( Matrix<Double> &bMat, Integer ip,
                              Matrix<Double> &ptCoord ) {

    ENTER_FCN( "linPiezoInt::calcBMat" );

    // obtain info on problem sizes
    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer spaceDim = ptelem->GetDim();
    const Integer nrDofs   = getNrDofs();
    const Integer offset   = spaceDim + 1;

    // set correct size of matrix B
    bMat.Resize( getDimD(), nrNodes * nrDofs );
    
    // derivatives of local shape functions with respect to global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    
    if (isSetIntPoint_) 
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord);
    else
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);

    // auxiliary variables
    Integer actDim, actNode;

    // treat mechanical part (same as in linElastInt)
    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < nrNodes; actNode++)
        bMat[actDim][actNode * offset + actDim] = xiDx[actNode][actDim];    
    switch(spaceDim)
      {
      case 2:
        for ( actNode = 0; actNode < nrNodes; actNode++ )
          {
            bMat[spaceDim][actNode * offset + 1] = xiDx[actNode][0];
            bMat[spaceDim][actNode * offset + 0] = xiDx[actNode][1];
          }
        break;


      case 3:
        actDim = spaceDim;
        for ( actNode = 0; actNode < nrNodes; actNode++ )
          {
            bMat[actDim][actNode * offset + 1] = xiDx[actNode][2];
            bMat[actDim][actNode * offset + 2] = xiDx[actNode][1];
          }

        actDim++;
        for ( actNode = 0; actNode < nrNodes; actNode++ )
          {
            bMat[actDim][actNode * offset    ] = xiDx[actNode][2];
            bMat[actDim][actNode * offset + 2] = xiDx[actNode][0];
          }

        actDim++;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * offset    ] = xiDx[actNode][1];
            bMat[actDim][actNode * offset + 1] = xiDx[actNode][0];
          }
        break;
      }

    // treat electrical part
    for( actDim = 0; actDim < spaceDim; actDim++ )
      for( actNode = 0; actNode < nrNodes; actNode++ )
        bMat[2*spaceDim+actDim][(actNode+1)*offset-1] = 
          xiDx[actNode][actDim];

#ifdef DEBUG
    //    (*debug) << std::endl << " Matrix bMat is " << bMat.GetSizeRow()
    //       << " x " << bMat.GetSizeCol() << std::endl;
    //    (*debug) << bMat << std::endl;
#endif

  }


  // determines the material matrix D containing the tensors of mechanical
  // modulus, electrical permittivity and piezoelectric coupling
  // for the 2D plane strain case
  void linPiezoInt::CalcPlaneStrainMaterialMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linPiezoInt::CalcPlaneStrainMaterialMat" );

    Integer rowPtrXY[]={2,3,5,8,9};
    Integer rowPtrYZ[]={2,3,4,8,9};
    Integer rowPtrXZ[]={1,3,5,7,9};
    Integer * rowPtr;

    switch(actOrientation)
      {
      case xy:
        {
          rowPtr=rowPtrXY;
          break;
        }
      case yz:
        {
          rowPtr=rowPtrYZ;
          break;
        }
      case xz:
        {
          rowPtr=rowPtrXZ;
          break;
        }
      default:  //if no orientation was specified
        {
          rowPtr=rowPtrYZ;
          break;
        }
      }

    // set the material matrix
    Integer sizeofD=getDimD();
    dMat.Resize(sizeofD);
    dMat.Init(0);

    Matrix<Double>*matMatrix; 

    if (piezoMatType_ == realMaterialParameter)
      matMatrix = ptMaterial->GetMatrix();
    else if (piezoMatType_ == imagMaterialParameter)
      matMatrix = ptMaterial->GetMatrixC();

    // No damping
    if ( isDamping_ == false ) {

      // Copy entries from material matrix object into D matrix
      for ( Integer i = 0; i < sizeofD; i++ ) {
        for ( Integer j = 0; j < sizeofD; j++ ) {
          dMat[i][j] = (*matMatrix)[rowPtr[i]-1][rowPtr[j]-1];
        }
      }

      // multiply values of permittivity with -1 to obtain the correct
      // bilinearform
      for ( Integer i = sizeofD - 2; i < sizeofD; i++ ) {
        for ( Integer j = sizeofD - 2; j < sizeofD; j++ ) {
          dMat[i][j] *= -1.0;
        }
      }
    }

    // The damping case (Rayleigh currently). Here just the mechanical part
    // damps
    else {

      // Copy entries from mechanical part of material matrix object
      // into D matrix and scale with damping parameter
      Matrix<Double> *matMatrix = NULL;

      if ( piezoMatType_ == realMaterialParameter) {
	matMatrix = ptMaterial->GetMatrix();
      }
      else if ( piezoMatType_ == imagMaterialParameter ) {
	matMatrix = ptMaterial->GetMatrixC();
      }
      else {
	Error( "piezoMatType neither real nor imaginary. I'm lost!",
	       __FILE__, __LINE__ );
      }

      Double beta = ptMaterial->GetDampingBeta();

      for( Integer i = 0; i < sizeofD - 2; i++ ) {
        for ( Integer j = 0; j < sizeofD - 2; j++ ) {
          dMat[i][j] = (*matMatrix)[i][j] * beta * factorDamp_;
        }
      }
    }

    // If desired, try to equilibrate the material parameters
    std::string equilibrate;
    params->Get( "equilibrate", equilibrate, "piezo" );
    if ( equilibrate == "yes" ) {

#ifdef DEBUG_LINPIEZOINT
      (*debug) << "CalcPlaneStrainMaterialMat uses equilibrated material "
	       << "parameters" << std::endl;
#endif

      // Scale mechanical parameters by 10^{-10}
      for( Integer i = 0; i < sizeofD - 2; i++ ) {
        for ( Integer j = 0; j < sizeofD - 2; j++ ) {
          dMat[i][j] *= 1e-10;
        }
      }

      // Scale electrical parameters by 10^{+10}
      for( Integer i = sizeofD - 2; i < sizeofD; i++ ) {
        for ( Integer j = sizeofD - 2; j < sizeofD; j++ ) {
          dMat[i][j] *= 1e+10;
        }
      }

    }
  }


  // determine the material matrix D containing the tensors of mechanical
  // modulus, electrical permittivity and piezoelectric coupling
  // for the axisymmetric case
  void linPiezoInt::CalcAxiMaterialMat(Matrix<Double> & dMat) {

    ENTER_FCN( "linPiezoInt::CalcAxiMaterialMat" );

    Integer rowPtrXY[]={1,2,6,3,7,8};
    Integer rowPtrYZ[]={2,3,4,1,8,9};
    Integer rowPtrXZ[]={1,3,5,2,7,9};
    Integer * rowPtr;   //alte Version

    switch(actOrientation) {

    case xy:
      {
        rowPtr = rowPtrXY;
        break;
      }
    case yz:
      {
        rowPtr = rowPtrYZ;
        break;
      }
    case xz:
      {
        rowPtr = rowPtrXZ;
        break;
      }
    default:    //if no orientation was specified
      {
        rowPtr = rowPtrYZ;
        break;
      }
    }

    // Set the material matrix
    Integer sizeofD=getDimD();
    dMat.Resize(sizeofD);
    dMat.Init(0);


    Matrix<Double> * matMatrix;
    
    if (piezoMatType_ == realMaterialParameter)
      matMatrix = ptMaterial->GetMatrix();
    else if (piezoMatType_ == imagMaterialParameter)
      matMatrix = ptMaterial->GetMatrixC();


    // No damping
    if ( isDamping_ == false ) {

      // Copy entries from material matrix object into D matrix
      for ( Integer i = 0; i < sizeofD; i++ ) {
        for( Integer j = 0; j < sizeofD; j++ ) {
          dMat[i][j] = (*matMatrix)[rowPtr[i]-1][rowPtr[j]-1];
        }
      }

      // Multiply values of permittivity with -1 to obtain the
      // correct bilinearform
      for ( Integer i = sizeofD - 2; i < sizeofD; i++ ) {
        for ( Integer j = sizeofD - 2; j < sizeofD; j++ ) {
          dMat[i][j] *= -1.0;
        }
      }
    }
    
    // The damping case (Rayleigh currently). Here just the mechanical part
    // damps
    else {
      
      // Copy entries from mechanical part of material matrix object
      // into D matrix and multiply with damping parameter


      Matrix<Double> * matMatrix;
    
      if (piezoMatType_ == realMaterialParameter)
	matMatrix = ptMaterial->GetMatrix();
      else if (piezoMatType_ == imagMaterialParameter)
	matMatrix = ptMaterial->GetMatrixC();

      Double beta = ptMaterial->GetDampingBeta();
      for( Integer i = 0; i < sizeofD - 2; i++ ) {
        for ( Integer j = 0; j < sizeofD - 2; j++ ) {
          dMat[i][j] = (*matMatrix)[i][j] * beta;
        }
      }
    }

    // If desired, try to equilibrate the material parameters
    std::string equilibrate;
    params->Get( "equilibrate", equilibrate, "piezo" );
    if ( equilibrate == "yes" ) {

#ifdef DEBUG_LINPIEZOINT
      (*debug) << "CalcAxiMaterialMat uses equilibrated material parameters"
	       << std::endl;
#endif

      // Scale mechanical parameters by 10^{-10}
      for( Integer i = 0; i < sizeofD - 2; i++ ) {
        for ( Integer j = 0; j < sizeofD - 2; j++ ) {
          dMat[i][j] *= 1e-10;
        }
      }

      // Scale electrical parameters by 10^{+10}
      for( Integer i = sizeofD - 2; i < sizeofD; i++ ) {
        for ( Integer j = sizeofD - 2; j < sizeofD; j++ ) {
          dMat[i][j] *= 1e+10;
        }
      }

    }
  }
  

  // determine the material matrix D containing the tensors of mechanical
  // modulus, electrical permittivity and piezoelectric coupling
  // in the 3D case
  void linPiezoInt::Calc3DMaterialMat( Matrix<Double> &dMat ) {

    ENTER_FCN( "linPiezoInt::Calc3DMaterialMat" );

    // Resize and initialise matrix object
    const Integer sizeofD = getDimD();
    dMat.Resize( sizeofD );
    dMat.Init( 0 );

    // Standard case: No damping
    if ( isDamping_ == false ) {

      // Copy entries from material matrix object into D matrix
      Matrix<Double> * matMatrix;
    
      if (piezoMatType_ == realMaterialParameter)
	matMatrix = ptMaterial->GetMatrix();
      else if (piezoMatType_ == imagMaterialParameter)
	matMatrix = ptMaterial->GetMatrixC();
      
      for( Integer i = 0; i < sizeofD; i++ ) {
        for ( Integer j = 0; j < sizeofD; j++ ) {
          dMat[i][j] = (*matMatrix)[i][j];
        }
      }

      // Multiply values of permittivity with -1 to obtain the correct
      // bilinear form
      for ( Integer i = sizeofD - 3; i < sizeofD; i++ ) {
        for ( Integer j = sizeofD - 3; j <sizeofD; j++ ) {
          dMat[i][j] *= -1.0;
        }
      }
    }

    // The damping case (Rayleigh currently). Here just the mechanical part
    // damps
    else {

      // Copy entries from mechanical part of material matrix object
      // into D matrix and multiply with damping parameter

      Matrix<Double> * matMatrix;
    
      if (piezoMatType_ == realMaterialParameter)
	matMatrix = ptMaterial->GetMatrix();
      else if (piezoMatType_ == imagMaterialParameter)
	matMatrix = ptMaterial->GetMatrixC();

      Double beta = ptMaterial->GetDampingBeta();
      for( Integer i = 0; i < sizeofD - 3; i++ ) {
        for ( Integer j = 0; j < sizeofD - 3; j++ ) {
          dMat[i][j] = (*matMatrix)[i][j] * beta;
        }
      }
    }

    // If desired, try to equilibrate the material parameters
    std::string equilibrate;
    params->Get( "equilibrate", equilibrate, "piezo" );
    if ( equilibrate == "yes" ) {

#ifdef DEBUG_LINPIEZOINT
      (*debug) << "Calc3DMaterialMat uses equilibrated material parameters"
	       << std::endl;
#endif

      // Scale mechanical parameters by 10^{-10}
      for( Integer i = 0; i < sizeofD - 3; i++ ) {
        for ( Integer j = 0; j < sizeofD - 3; j++ ) {
          dMat[i][j] *= 1e-10;
        }
      }

      // Scale electrical parameters by 10^{+10}
      for( Integer i = sizeofD - 3; i < sizeofD; i++ ) {
        for ( Integer j = sizeofD - 3; j < sizeofD; j++ ) {
          dMat[i][j] *= 1e+10;
        }
      }

    }

#ifdef DEBUG_LINPIEZOINT
    (*debug) << std::endl << " Matrix dMat:\n" << dMat << std::endl;
#endif

  }


  // calculates of stresses T (vector notation)
  // T = c . S - e^T E with c tensor of mechanical moduli and e the
  // piezoelectric tensor
  // S = Bmech * u 
  // E = Belec V
  // see Habil. M. Kaltenbacher 
  void linPiezoInt::CalcStressVec( Vector<Double>& stressElecVec, Integer ip, 
                                   Matrix<Double> & ptCoord ) {

    ENTER_FCN( "linPiezoInt::CalcStressVec" );
      
    Matrix<Double> dMat;
    calcDMat(dMat);
 
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, VNode1, uNode2X, uNode2Y,VNode2,  ...)
    Vector<Double> solVec;
    elemSol_->ConvertToVec_AppendCols(solVec);

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    Vector<Double> linStrainElec(linBMat * solVec );

    // | c Bmech u - e^T Belec V |
    // | e Bmech u + eps Belec V |
    //    Vector<Double> stressElecVec = dMat * linStrainElec;
    stressElecVec = dMat * linStrainElec;

  }

  void linPiezoInt::CalcStressVec( Vector<Complex>& stressElecVec, Integer ip, 
                                   Matrix<Double> & ptCoord ) {

    ENTER_FCN( "linPiezoInt::CalcStressVec" );
     
    Matrix<Double> dMat;
    calcDMat(dMat);
 
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, VNode1, uNode2X, uNode2Y,VNode2,  ...)
    Vector<Complex> solVec;
    elemSol_->ConvertToVec_AppendCols(solVec);
  

    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    Vector<Complex> linStrainElec; // (linBMat * solVec );
    linStrainElec.Resize(linBMat.GetSizeRow());

    linBMat.MatVecMult_DC(solVec,linStrainElec);

    // | c Bmech u - e^T Belec V |
    // | e Bmech u + eps Belec V |
   
    stressElecVec.Resize(dMat.GetSizeRow());
    dMat.MatVecMult_DC(linStrainElec, stressElecVec);

  }

  


  // ========================================================================
  // ======================== linPiezoAxiInt - Part =========================
  // ========================================================================

  // determine the matrix B of the BDB operator
  void piezoAxiInt::calcBMat( Matrix<Double> &bMat, Integer ip,
                              Matrix<Double> &ptCoord ) {

    ENTER_FCN( "piezoAxiInt::calcBMat" );

    // obtain info on problem sizes
    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer spaceDim = ptelem->GetDim();
    const Integer nrDofs   = getNrDofs();
    const Integer offset   = spaceDim + 1;

    // set correct size of matrix B
    bMat.Resize( getDimD(), nrNodes * nrDofs );
    Integer totalnodes=nrNodes * nrDofs;

    // derivatives of local shape functions with respect to global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;

    if (isSetIntPoint_) 
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord);
    else
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);

    // auxiliary variables
    Integer actDim, actNode, idxtheta = getDimD();
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIp;

    if (isSetIntPoint_) 
      ptelem->GetShFnc(ShpFncAtIp,intPoint_);
    else
      ptelem->GetShFncAtIp(ShpFncAtIp,ip);

    CoordAtIp = ptCoord * ShpFncAtIp;

    // treat mechanical part (same as in linElastInt)
    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < nrNodes; actNode++)
        bMat[actDim][actNode * (spaceDim+1) + actDim] = xiDx[actNode][actDim];

    for (actNode = 0; actNode<nrNodes;actNode++){
      bMat[2][actNode*offset] = xiDx[actNode][1];
      bMat[2][actNode*offset+1] = xiDx[actNode][0];
      bMat[3][actNode*offset] = ShpFncAtIp[actNode]/CoordAtIp[0];
    }

    // treat electrical part
    for( actDim = 0; actDim < spaceDim; actDim++ )
      for( actNode = 0; actNode < nrNodes; actNode++ )
        bMat[2*spaceDim+actDim][(actNode+1)*offset-1] = xiDx[actNode][actDim];

#ifdef DEBUG_LINPIEZOINT
    (*debug) << std::endl << " Matrix bMat is " << bMat.GetSizeRow()
	     << " x " << bMat.GetSizeCol() << std::endl;
    (*debug) << bMat << std::endl;
#endif

  }


  // calculates the D-matrix of a axisymmetric-problem 
  void piezoAxiInt::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linPiezoAxiInt::calcDMat" );
    CalcAxiMaterialMat(dMat);
  }


  // calculates the D-matrix of a axisymmetric-problem 
  void piezoAxiInt::calcDMaterialMatWithComplexDamping( Matrix<Complex> &dMat,
							Double &beta,
							Double &omega ) {

    ENTER_FCN( "linPiezoAxiInt::calcDMatWithComplexDamping" );

    Integer rowPtrXY[]={1,2,6,3,7,8};
    Integer rowPtrYZ[]={2,3,4,1,8,9};
    Integer rowPtrXZ[]={1,3,5,2,7,9};
    Integer * rowPtr;   //alte Version

    switch(actOrientation) {

    case xy:
      {
        rowPtr = rowPtrXY;
        break;
      }
    case yz:
      {
        rowPtr = rowPtrYZ;
        break;
      }
    case xz:
      {
        rowPtr = rowPtrXZ;
        break;
      }
    default:    //if no orientation was specified
      {
        rowPtr = rowPtrYZ;
        break;
      }
    }

    // Set the material matrix
    Integer sizeofD=getDimD();
    dMat.Resize(sizeofD);
    dMat.Init(0);


    // The damping case. The matrix is multiplied by (1+\omega_l*\beta_l*j)
    // Only mech part is damped 
    Complex im=Complex(0,1.0);
    // Copy entries from mechanical part of material matrix object
    // into D matrix and multiply with damping parameters
    Matrix<Double> * matMatrix;

    if (piezoMatType_ == realMaterialParameter)
      matMatrix = ptMaterial->GetMatrix();
    else if (piezoMatType_ == imagMaterialParameter)
      matMatrix = ptMaterial->GetMatrixC();

    for ( Integer i = 0; i < sizeofD; i++ ) {
      for( Integer j = 0; j < sizeofD; j++ ) {
	dMat[i][j] = (*matMatrix)[rowPtr[i]-1][rowPtr[j]-1] *
	  ( 1.0 + beta * omega * im );
      }
    }

  } // end calcDMatWithComplexDamping



  // ========================================================================
  // ======================== linPiezoPlaneStrainInt - Part =================
  // ========================================================================


  // calculates the D-matrix of a axisymmetric-problem 
  void piezoPlainStrainInt::calcDMat(Matrix<Double> & dMat) {
    ENTER_FCN( "piezoPlainStrainInt::calcDMat" );
    CalcPlaneStrainMaterialMat(dMat);
  }

  void piezoPlainStrainInt::
  calcDMaterialMatWithComplexDamping( Matrix<Complex> &dMat, Double &beta,
				      Double &omega ) {

    ENTER_FCN( "piezoPlainStrainInt::calcDMaterialMatWithComplexDamping" );

    Integer rowPtrXY[]={2,3,5,8,9};
    Integer rowPtrYZ[]={2,3,4,8,9};
    Integer rowPtrXZ[]={1,3,5,7,9};
    Integer * rowPtr;

    switch(actOrientation) {

      case xy:
	rowPtr = rowPtrXY;
	break;

      case yz:
	rowPtr = rowPtrYZ;
	break;

      case xz:
	rowPtr = rowPtrXZ;
	break;

	// No orientation was specified
      default:
	rowPtr = rowPtrYZ;
	break;
    }

    // set the material matrix
    Integer sizeofD=getDimD();
    dMat.Resize(sizeofD);
    dMat.Init(0);
       
    // Copy entries from mechanical part of material matrix object
    // into D matrix and scale with damping parameter (1+j*omega*beta)
    Matrix<Double> * matMatrix;

    if (piezoMatType_ == realMaterialParameter)
      matMatrix = ptMaterial->GetMatrix();
    else if (piezoMatType_ == imagMaterialParameter)
      matMatrix = ptMaterial->GetMatrixC();

    Complex imag=Complex(0,1);
    for( Integer i = 0; i < sizeofD; i++ ) {
      for ( Integer j = 0; j < sizeofD; j++ ) {
	dMat[i][j] = (*matMatrix)[i][j] * beta * omega * imag;
      }
    }
  }


  // determine the matrix B of the BDB operator
  void piezoPlainStrainInt::calcBMat( Matrix<Double> &bMat, Integer ip,
                                      Matrix<Double> &ptCoord ) {

    ENTER_FCN( "piezoPlainStrainInt::calcBMat" );

    // obtain info on problem sizes
    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer spaceDim = ptelem->GetDim();
    const Integer nrDofs   = getNrDofs();
    const Integer offset   = spaceDim + 1;

    // set correct size of matrix B
    bMat.Resize( getDimD(), nrNodes * nrDofs );

    Integer totalnodes=nrNodes * nrDofs;

    Matrix<Double> xiDx;

    if (isSetIntPoint_) 
      ptelem->GetGlobDerivShFnc(xiDx, intPoint_, ptCoord);
    else
      ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);

    // auxiliary variables
    Integer actDim, actNode;

    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < nrNodes; actNode++){
        bMat[actDim][actNode * offset + actDim] = xiDx[actNode][actDim];
      }


    for ( actNode = 0; actNode < nrNodes; actNode++ )
      {
        bMat[spaceDim][actNode * offset + 1] = xiDx[actNode][0];
        bMat[spaceDim][actNode * offset + 0] = xiDx[actNode][1];
      }

    // treat electrical part
    for( actDim = 0; actDim < spaceDim; actDim++ ) {
      for( actNode = 0; actNode < nrNodes; actNode++ ) {
        bMat[2*spaceDim+actDim-1][(actNode+1)*offset-1] =
          xiDx[actNode][actDim];
      }
    }

#ifdef DEBUG
    //    (*debug) << std::endl << " Matrix bMat is " << bMat.GetSizeRow()
    //       << " x " << bMat.GetSizeCol() << std::endl;
    //    (*debug) << bMat << std::endl;
#endif

  }



  // ========================================================================
  // ======================== linPiezo3DInt - Part ==========================
  // ========================================================================


  // constructor
  linPiezo3DInt::linPiezo3DInt(BaseFE * aptelem, MaterialData & matData) 
    : linPiezoInt(aptelem, matData) {
    ENTER_FCN( "linPiezo3DInt::linPiezo3DInt" );
  }
 
  // destructor
  linPiezo3DInt::~linPiezo3DInt() {
    ENTER_FCN( "linPiezo3DInt::~linPiezo3DInt" );
  }

  // calculates the D-matrix of a 3d-problem 
  void linPiezo3DInt::calcDMat( Matrix<Double> & dMat ) {
    ENTER_FCN( "linPiezo3DInt::calcDMat" );
    Calc3DMaterialMat(dMat);
  }

  // reimplemented  method of Calc3DMaterialMat. 
  // It is needed to damp the mechanical part in the following way:
  // K=(1+j*\beta_l* \omega_l)K_dd
  // it is needed in the piezoParamIdent Driver during the calculation of the
  // Jacbian matrix of the parameter to solution map.
  void linPiezo3DInt::calcDMaterialMatWithComplexDamping(Matrix<Complex> &dMat,
							 Double &beta,
							 Double &omega ) {

    ENTER_FCN( "linPiezoInt::Calc3DMaterialMatWithComplexDamping" );

    // Resize and initialise matrix object
    const Integer sizeofD = getDimD();
    dMat.Resize( sizeofD );
    dMat.Init( 0 );

    // Copy entries from mechanical part of material matrix object
    // into D matrix and multiply with damping parameter beta_l, omega_l
    // and the imaginary j
    Matrix<Double> * matMatrix;

    if (piezoMatType_ == realMaterialParameter) {
      matMatrix = ptMaterial->GetMatrix();
    }
    else if (piezoMatType_ == imagMaterialParameter) {
      matMatrix = ptMaterial->GetMatrixC();
    }

    Complex imag = Complex( 0, 1.0 );
    for( Integer i = 0; i < sizeofD - 3; i++ ) {
      for ( Integer j = 0; j < sizeofD - 3; j++ ) {
	dMat[i][j] = (*matMatrix)[i][j] * ( 1.0 + beta * imag * omega );
      }
    }

  } // end Calc3dMaterialMatWithComplexDamping

} // end namespace CoupledField
