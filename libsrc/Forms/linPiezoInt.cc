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

  // determine the material matrix D containing the tensors of mechanical
  // modulus, electrical permittivity and piezoelectric coupling
  void linPiezoInt::Calc3DMaterialMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "linPiezoInt::calcDMat" );

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
