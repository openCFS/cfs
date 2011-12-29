// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <string>

#include "Domain/elem.hh"
#include "Domain/resultInfo.hh"
#include "Elements/basefe.hh"
#include "Forms/gradfieldop.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/nodestoresol.hh"
#include "nLinPiezoCoupling.hh"


namespace CoupledField {


  // ============
  //   Constructor of nonlinear piezo coupling
  // ============
class EqnMap;
class Grid;
class StdPDE;

  nLinPiezoCoupling::nLinPiezoCoupling( BaseMaterial* matData,
                                        BaseMaterial* matDataMech,
                                        BaseMaterial* matDataElec,                                       
                                        SubTensorType type) 
    : linPiezoCoupling(matData, type){
    
    name_ = "nLinPiezoCoupling";
    isSolDependent_ = true;

    isHysteresis_ = true;

    matDataMech_ = matDataMech;
    matDataElec_ = matDataElec;
    
    if ( type == FULL ) {
      numDofsA_  = 3;
      numDofsB_  = 1;
      matDimRow_ = 6;
      matDimCol_ = 3;
    }
    else if ( type == PLANE_STRAIN || type == PLANE_STRESS ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 3;
      matDimCol_ = 2;
    }
    else if ( type == AXI ) {
      numDofsA_  = 2;
      numDofsB_  = 1;
      matDimRow_ = 4;
      matDimCol_ = 2;
      isaxi_     = true;
    }
  }

    //! Destructor
  nLinPiezoCoupling::~nLinPiezoCoupling() {
  }



  void nLinPiezoCoupling::Set4NonLinMaterial(Grid* ptGrid, 
                                             StdPDE* ptPDE,
                                             shared_ptr<EqnMap> eqnMap,
                                             shared_ptr<ResultInfo> result) 
  {

    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE, eqnMap, *sol2_,
                                             result->fctType, isaxi_, coordUpdate_);

    if (  isHysteresis_ ) {
      // get direction of polarization
      std::string str;
      matDataElec_->GetScalar(str, P_DIRECTION);
      Directions dir;
      String2Enum(str,dir);
      dirP_ = dir;

      // get maximum of polarization
      matDataElec_->GetScalar(Psat_, Y_SATURATION, Global::REAL);

    }
  }

  // ====================
  //   calcElementMatrix
  // ====================
  void nLinPiezoCoupling::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2 ) {


    ent1_ = ent1;
      
    // get displacements of element
    sol1_->GetElemSolutionAsMatrix( elemDispl_, ent1);

    // get elecPot of element
    sol2_->GetElemSolution( elemPot_, ent2);

    Vector<Double> LCoord, Efield;

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 ); 
    ExtractElemInfo( ent2 ); 


    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const UInt nrNodes  = ptelem->GetNumNodes();   
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    double jacDet;

    Matrix<Double> aMat;
    Matrix<Double> bMat;
    Matrix<Double> dMat;
    Matrix<Double> dbMat;
    Double aux;


    elemMat.Resize( nrNodes * getNumDofsA(), nrNodes * getNumDofsB() );
    elemMat.Init();


    // **************************************************
    //  Material matrix independent of integration point
    // **************************************************

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

      // Setup nonlinear material matrix once and for all
      calcDMat( dMat,actIntPt,ptCoord_ );


      // Setup the A matrix for current integration point
       calcAMat( aMat, actIntPt, ptCoord_ );

       // Setup the B matrix for current integration point
       CalcBMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
       jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_,ent1.GetElem() );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        EXCEPTION("ADBInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!");
      }

      // Special things must be done in the axi-symmetric case
      // We need to additionally scale with 2 pi r.
      //
      // NOTE: We assume here that computation is in they-z plane
      //       with z being the axis of symmetry and that y is
      //       represented by the first co-ordinate, thus
      //       2 pi r = "2 pi x"
      if ( isaxi_ ) {
        Vector<Double> ShpFncAtIp;
        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );
        Double aux = 0.0;
        
        for ( UInt i = 0; i < nrNodes; i++ ) {
          aux += ptCoord_[0][i] * ShpFncAtIp[i];
        }
        
        jacDet *= 2.0 * PI * aux;
      }

      // Compute the matrix product D * B and store as intermediate matrix
      dbMat.Resize( dMat.GetNumRows(), bMat.GetNumCols() );
      dMat.Mult( bMat, dbMat );

      // We now compute A * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix
      for ( UInt i = 0; i < aMat.GetNumRows(); i++ ) {
        for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {

          // Compute entry (i,j) of A * D * B
          aux = 0.0;
          for ( UInt k = 0; k < aMat.GetNumCols(); k++ ) {
            aux += aMat[i][k] * dbMat[k][j];
          }

          // Scale result and add to corresponding entry
          // of element matrix
          elemMat[i][j] += aux * jacDet * intWeights[actIntPt-1];
        }
      }
    }    

  }


  void nLinPiezoCoupling::calcDMat(Matrix<Double> & dMat, 
                                   UInt ip, 
                                   Matrix<Double> & ptCoord)
  {

    Matrix<Double> matMatrix;
    ptMaterial->GetTensor(matMatrix,PIEZO_TENSOR,matDataType_,subTensorType_);
    matMatrix.Transpose(dMat);

    Matrix<Double> stiffMat;
    matDataMech_->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR,Global::REAL);

    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord,  it1_.GetElem());

    Vector<Double> LCoord, Efield;

    Vector<Double> TempBu;

    // calc E-Field:
    //compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);

    EfieldOp_->CalcElemGradField( Efield, ent1_, LCoord, 1);


    if ( isHysteresis_ ) {
      // coupling tensor depends on electric polarization
      UInt nrEl = ent1_.GetElem()->elemNum;
      Double actP = matDataElec_->ComputeScalarHystVal( nrEl, Efield[dirP_] );
      Double scaleFactor = actP / Psat_;
      dMat *= scaleFactor;
      //     std::cout << " scaleFactor= " <<  scaleFactor << std::endl;
    }

    //    std::cout << "\n dMat:\n " << dMat << std::endl;
  }
}
