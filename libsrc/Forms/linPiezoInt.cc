#include <iostream>
#include <fstream>

#include "linPiezoInt.hh"

namespace CoupledField
{

  // ========================================================================
  // ========================= linPiezoInt - Part ===========================
  // ========================================================================
 

  // determine the matrix B of the BDB operator
  void linPiezoInt::calcBMat( Matrix<Double> &bMat, Integer ip,
			      Matrix<Double> &ptCoord )
  {
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
    //	     << " x " << bMat.GetSizeCol() << std::endl;
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
      default:	//if no orientation was specified
	{
	  rowPtr=rowPtrYZ;
	  break;
	}
      }


    Matrix<Double>*matMatrix = ptMaterial->GetMatrix();
    if ( isDamping_ == false ) {
      // Copy entries from material matrix object into D matrix
 
      Integer sizeofD=getDimD();
      dMat.Resize(sizeofD);
      dMat.Init(0);

      for (int i=0; i< sizeofD; i++)
	for (int j=0;j< sizeofD; j++)
	  dMat[i][j]=(*matMatrix)[rowPtr[i]-1][rowPtr[j]-1];


      //multiply values of permittivity with -1 to obtain the correct bilinearform
      //hardcoded for 2D!!!
      for (Integer i = sizeofD-2; i <sizeofD; i++)
	for (Integer j = sizeofD-2; j <sizeofD; j++)
	  dMat[i][j] *= -1.0;
    }
    else {
      // The damping case (Rayleigh currently). Here just the mechanical part
      // damps

      Integer sizeofD = getDimD()-2;
      dMat.Resize(sizeofD);
      dMat.Init(0);

      // Copy entries from mechanical part of material matrix object
      // into D matrix and multiply with damping parameter
      Matrix<Double> * matMatrix = ptMaterial->GetMatrix();
      Double beta = ptMaterial->GetDampingBeta();
      for( Integer i = 0; i < sizeofD; i++ ) {
	for ( Integer j = 0; j < sizeofD; j++ ) {
	  dMat[i][j] = (*matMatrix)[i][j] * beta;
	}
      }
    }
    

  }


  // determine the material matrix D containing the tensors of mechanical
  // modulus, electrical permittivity and piezoelectric coupling
  //for the axisymmetric case
  void linPiezoInt::CalcAxiMaterialMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linPiezoInt::CalcAxiMaterialMat" );

    Integer rowPtrXY[]={1,2,6,3,7,8};
    Integer rowPtrYZ[]={2,3,4,1,8,9};
    Integer rowPtrXZ[]={1,3,5,2,7,9};
    Integer * rowPtr; 	//alte Version

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
      default:	//if no orientation was specified
	{
	  rowPtr=rowPtrYZ;
	  break;
	}
      }

    Matrix<Double> * matMatrix = ptMaterial->GetMatrix();
    

    if ( isDamping_ == false ) {
      // Copy entries from material matrix object into D matrix
      
      Integer sizeofD=getDimD();
      dMat.Resize(sizeofD);
      dMat.Init(0);

      for (int i=0; i< sizeofD; i++)
	for(int j=0; j<sizeofD; j++){
	  dMat[i][j]=(*matMatrix)[rowPtr[i]-1][rowPtr[j]-1];
	  
	  //multiply values of permittivity with -1 to obtain the correct bilinearform
	  //hardcoded for 2D!!!
	  for (Integer i = sizeofD-2; i <sizeofD; i++)
	    for (Integer j = sizeofD-2; j <sizeofD; j++)
	      dMat[i][j] *= -1.0;

	}
    }
    
    else {
      // The damping case (Rayleigh currently). Here just the mechanical part
      // damps
      
      Integer sizeofD = getDimD()-2;
      dMat.Resize(sizeofD);
      dMat.Init(0);
      
      // Copy entries from mechanical part of material matrix object
      // into D matrix and multiply with damping parameter
      Matrix<Double> * matMatrix = ptMaterial->GetMatrix();
      Double beta = ptMaterial->GetDampingBeta();
      for( Integer i = 0; i < sizeofD; i++ ) 
	for ( Integer j = 0; j < sizeofD; j++ ) 
	  dMat[i][j] = (*matMatrix)[i][j] * beta;
    }
  }
  

    // determine the material matrix D containing the tensors of mechanical
    // modulus, electrical permittivity and piezoelectric coupling
    void linPiezoInt::Calc3DMaterialMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linPiezoInt::Calc3DMaterialMat" );

    // Resize and initialise matrix object
    const Integer sizeofD = getDimD();
    dMat.Resize( sizeofD );
    dMat.Init( 0 );

     // Standard case: No damping
    if ( isDamping_ == false ) {
      // Copy entries from material matrix object into D matrix
      Matrix<Double> * matMatrix = ptMaterial->GetMatrix();
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
      Matrix<Double> * matMatrix = ptMaterial->GetMatrix();
      Double beta = ptMaterial->GetDampingBeta();
      for( Integer i = 0; i < sizeofD - 3; i++ ) {
	for ( Integer j = 0; j < sizeofD - 3; j++ ) {
	  dMat[i][j] = (*matMatrix)[i][j] * beta;
	}
      }
    }
  }


  // ========================================================================
  // ======================== linPiezoAxiInt - Part ==========================
  // ========================================================================

  // determine the matrix B of the BDB operator
  void piezoAxiInt::calcBMat( Matrix<Double> &bMat, Integer ip,
			      Matrix<Double> &ptCoord )
  {
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

    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord);

    // auxiliary variables
    Integer actDim, actNode, idxtheta = getDimD();
    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIp;

    ptelem->GetShFncAtIp(ShpFncAtIp, ip);
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

#ifdef DEBUG
    //    (*debug) << std::endl << " Matrix bMat is " << bMat.GetSizeRow()
    //	     << " x " << bMat.GetSizeCol() << std::endl;
    //    (*debug) << bMat << std::endl;
#endif

  }


  // calculates the D-matrix of a axisymmetric-problem 
  void piezoAxiInt::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linPiezoAxiInt::calcDMat" );
    CalcAxiMaterialMat(dMat);
  }


  // ========================================================================
  // ======================== linPiezoPlaneStrainInt - Part =================
  // ========================================================================


  // calculates the D-matrix of a axisymmetric-problem 
  void piezoPlainStrainInt::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "piezoPlainStrainInt::calcDMat" );
    CalcPlaneStrainMaterialMat(dMat);
  }


  // determine the matrix B of the BDB operator
  void piezoPlainStrainInt::calcBMat( Matrix<Double> &bMat, Integer ip, Matrix<Double> &ptCoord )
  {
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
    for( actDim = 0; actDim < spaceDim; actDim++ )
      for( actNode = 0; actNode < nrNodes; actNode++ ){
	bMat[2*spaceDim+actDim-1][(actNode+1)*offset-1] = xiDx[actNode][actDim];
      }

#ifdef DEBUG
    //    (*debug) << std::endl << " Matrix bMat is " << bMat.GetSizeRow()
    //	     << " x " << bMat.GetSizeCol() << std::endl;
    //    (*debug) << bMat << std::endl;
#endif

  }



  // ========================================================================
  // ======================== linPiezo3DInt - Part ==========================
  // ========================================================================


  // constructor
  linPiezo3DInt::linPiezo3DInt(BaseFE * aptelem, MaterialData & matData) 
    : linPiezoInt(aptelem, matData)
  {
    ENTER_FCN( "linPiezo3DInt::linPiezo3DInt" );
  }
 
  // destructor
  linPiezo3DInt::~linPiezo3DInt()
  {
    ENTER_FCN( "linPiezo3DInt::~linPiezo3DInt" );
  }

  // calculates the D-matrix of a 3d-problem 
  void linPiezo3DInt::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linPiezo3DInt::calcDMat" );
    Calc3DMaterialMat(dMat);
  }


} // end namespace CoupledField
