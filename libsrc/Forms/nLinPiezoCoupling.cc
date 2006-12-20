#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
//#include "nLinPiezoCoupling.hh"
//#include "linPiezoCoupling.hh"
#include "Utils/nodestoresol.hh"
#include "Forms/forms_header.hh"


namespace CoupledField {


  // ============
  //   Constructor of nonlinear piezo coupling
  // ============
  nLinPiezoCoupling::nLinPiezoCoupling(ApproxData *nlinFnc, 
                                       BaseMaterial* matData,
                                       BaseMaterial* matDataMech,
                                       BaseMaterial* matDataElec,                                       
                                       SubTensorType type) 
    : linPiezoCoupling(matData, type){
    ENTER_FCN( "nLinPiezoCoupling::nLinPiezoCoupling" );
    
    name_ = "nLinPiezoCoupling";
    isSolDependent_ = true;
    nLinFnc_ = nlinFnc;

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
    ENTER_FCN( "nLinPiezoCoupling::~nLinPiezoCoupling" );
  }



  void nLinPiezoCoupling::Set4NonLinMaterial(Grid* ptGrid, 
                                             StdPDE* ptPDE,
                                             shared_ptr<EqnMap> eqnMap,
                                             shared_ptr<ResultDof> result) 
  {
    ENTER_FCN( "nLinPiezoCoupling::Set4NonLinMaterial" );

    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE, 
                                             eqnMap, *sol2_, 
                                             ELEC_POTENTIAL, 
                                             result, isaxi_, 
                                             coordUpdate_);
  }

  // ====================
  //   calcElementMatrix
  // ====================
  void nLinPiezoCoupling::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1, 
                                             EntityIterator& ent2 ) {
    ENTER_FCN( "nLinPiezoCoupling::CalcElementMatrix" );


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
       calcBMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
       jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_,ent1.GetElem() );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        (*error) << "ADBInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!";
        Error( __FILE__, __LINE__ );
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
      dbMat.Resize( dMat.GetSizeRow(), bMat.GetSizeCol() );
      dMat.Mult( bMat, dbMat );

      // We now compute A * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix
      for ( UInt i = 0; i < aMat.GetSizeRow(); i++ ) {
        for ( UInt j = 0; j < dbMat.GetSizeCol(); j++ ) {

          // Compute entry (i,j) of A * D * B
          aux = 0.0;
          for ( UInt k = 0; k < aMat.GetSizeCol(); k++ ) {
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
    ENTER_FCN( "nLinPiezoCoupling::calcDMat" );

    Matrix<Double> matMatrix;
    ptMaterial->GetTensor(matMatrix,PIEZO_TENSOR,matDataType_,subTensorType_);
    matMatrix.Transpose(dMat);

    Matrix<Double> stiffMat;
    matDataMech_->GetTensor(stiffMat,MECH_STIFFNESS_TENSOR,REAL);

    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord,  it1_.GetElem());

    Vector<Double> LCoord, Efield;

    Vector<Double> TempBu;

    // calc E-Field:
    //compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);

    EfieldOp_->CalcElemGradField( Efield, ent1_, LCoord, 1);


    MechStressStrain<Double> * mechStrainOp = 
      new MechStressStrain<Double>(matDataMech_, subTensorType_);

    mechStrainOp->SetActElemSol(elemDispl_);
    mechStrainOp->SetIntPoint(LCoord);

    mechStrainOp->CalcStrainVec(TempBu,1,ent1_);

    
     Vector<Double> stressVec(stiffMat.GetSizeRow());
     stressVec.Init();
     stressVec = stiffMat * TempBu;	

    std::string nonLinDepend;
    ptMaterial->GetScalar(nonLinDepend, NONLIN_DEPENDENCY);
    
    Integer nonLinCoeff;
    ptMaterial->GetScalar(nonLinCoeff, NONLIN_COEFFICIENT);
    
    std::string nonLinApproxType;
    ptMaterial->GetScalar(nonLinApproxType,NONLIN_APPROXIMATION_TYPE);

    
    if(nonLinDepend=="mechStress"){
      
      Double stressSum=0.0;
      for(UInt i=0;i<stressVec.GetSize();i++)
        stressSum+=std::abs(stressVec[i]);	
      
        if (nonLinApproxType=="polynomial"){
          if(nonLinCoeff==33)
            dMat[2][2]+=0.01*(stressSum + 0.02*stressSum*stressSum)*dMat[2][2];
          else
            Error("The nonlinear parameter dependency is not known", __FILE__, __LINE__);
        }
        else if (nonLinApproxType=="smoothSplines")
         {
           Double approxCoeff;

           approxCoeff = nLinFnc_->EvaluateFunc(stressSum);
           if(nonLinCoeff==33){
             dMat[2][2] = 10.0+approxCoeff;
           }
         }
       else
         Error("The nonlinear approximation type is not known", __FILE__, __LINE__);
     }

    else if(nonLinDepend=="elecField"){
      Double eFieldSum=0.0;
      for(UInt i=0; i<Efield.GetSize();i++)
        eFieldSum+=std::abs(Efield[i]);
      
      if (nonLinApproxType=="polynomial"){
        if(nonLinCoeff==33)
          dMat[2][2]+=0.01*(eFieldSum + 0.02*eFieldSum*eFieldSum)*dMat[2][2];
        else
          Error("The nonlinear parameter dependency is not known", __FILE__, __LINE__);
      }
      else if (nonLinApproxType=="smoothSplines")
        {
          Double approxCoeff;
          
          approxCoeff = nLinFnc_->EvaluateFunc(eFieldSum);

          if(nonLinCoeff==33){
            dMat[2][2] = approxCoeff;
          }
        }
      else
        Error("The nonlinear approximation type is not known", __FILE__, __LINE__);
    }
     else{
       std::cout<<"The data dependency you have chosen is " << nonLinDepend <<std::endl;
       std::cout<<"(If this is true, check if there is a blank in font of choice) " <<std::endl;
       Error("Nonlinear dependency not implemented here", __FILE__, __LINE__);
     }
  }
}
