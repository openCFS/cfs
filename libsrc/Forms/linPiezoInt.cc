#include <iostream>
#include <fstream>

#include "linPiezoInt.hh"

namespace CoupledField
{

  // ==========================================================================
  // =========================== linPiezoInt - Part ===========================
  // ==========================================================================
 

  // determine the matrix B of the BDB operator
  void linPiezoInt::calcBMat( Matrix<Double> &bMat, Integer ip,
			      Matrix<Double> &ptCoord )
  {
#ifdef TRACE
    (*trace) << "entering linPiezoInt::calcBMat " << std::endl;
#endif

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
	bMat[2*spaceDim+actDim][(actNode+1)*offset-1] = xiDx[actNode][actDim];

#ifdef DEBUG
    (*debug) << std::endl << " Matrix bMat is " << bMat.size_row() << " x " <<
      bMat.size_col() << std::endl;
    (*debug) << bMat << std::endl;
#endif

  }


  // ==========================================================================
  // ========================== linPiezo3DInt - Part ==========================
  // ==========================================================================


  // constructor
  linPiezo3DInt::linPiezo3DInt(BaseFE * aptelem, MaterialData & matData) 
    : linPiezoInt(aptelem, matData)
  {
#ifdef TRACE
    (*trace) << "entering linPiezo3DInt::linPiezo3DInt" << std::endl;
#endif
  }

 
  // destructor
  linPiezo3DInt::~linPiezo3DInt()
  {
#ifdef TRACE
    (*trace) << "entering linPiezo3DInt::~linPiezo3DInt" << std::endl;
#endif
  }


  // determine the material matrix D containing the tensors of mechanical
  // modulus, electrical permittivity and piezoelectric coupling
  void linPiezo3DInt::calcDMat(Matrix<Double> & dMat)
  {

#ifdef TRACE
    (*trace) << "entering linPiezoInt::calcDMat " << std::endl;
#endif
    

    // get size of D and resize matrix object
    const Integer sizeofD = getDimD();
    dMat.Resize(sizeofD);
    
    // Copy entries from material matrix object into D matrix
    Matrix<Double> * matMatrix = ptMaterial->GetMatrix();
    for( Integer i = 0; i < sizeofD; i++ )
      for ( Integer j = 0; j < sizeofD; j++ )
	dMat[i][j] = (*matMatrix)[i][j];	

  }


} // end namespace CoupledField
