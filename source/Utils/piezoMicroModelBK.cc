// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "piezoMicroModelBK.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField{


  PiezoMicroModelBK::PiezoMicroModelBK( UInt numElemSD, BaseMaterial* piezoMat, 
                                        BaseMaterial* mechMat, 
                                        BaseMaterial* elecMat,
                                        SubTensorType tensorType, 
                                        Double dt)  {

    // number of elements belonging to this material
    numEl_ = numElemSD;

    // save pointers to material objects
    piezoMat_ = piezoMat;
    mechMat_  = mechMat;
    elecMat_  = elecMat;
    tensorType_ =  tensorType;

    deltaT_ = dt;

    if ( tensorType_ == FULL ) {
      dim_  = 3;
      dimS_ = 6;
    }
    else {
      dim_  = 2;
      dimS_ = 3;
    }

    InitSwitchingSystem();
  }

  // destructor
  PiezoMicroModelBK::~PiezoMicroModelBK(){
  }


  void PiezoMicroModelBK::InitSwitchingSystem() {

    // set number of domain types
    numDomain_ = 14;

    // define the rotation anles
    rotationAngles_.Resize(dim_,numDomain_);
    rotationAngles_[0][0] = 0;
    rotationAngles_[1][0] = -PI/2.0;
    rotationAngles_[2][0] = 0;
    rotationAngles_[0][1] = PI/2.0;
    rotationAngles_[1][1] = 0;
    rotationAngles_[2][1] = 0;
    rotationAngles_[0][2] = 0;
    rotationAngles_[1][2] = 0;
    rotationAngles_[2][2] = 0;
    rotationAngles_[0][3] = 0;
    rotationAngles_[1][3] = PI/2.0;
    rotationAngles_[2][3] = 0;
    rotationAngles_[0][4] = -PI/2.0;
    rotationAngles_[1][4] = 0;
    rotationAngles_[2][4] = 0;
    rotationAngles_[0][5] = PI;
    rotationAngles_[1][5] = 0;
    rotationAngles_[2][5] = 0;
    rotationAngles_[0][6] = 0.9557;
    rotationAngles_[1][6] = 0;
    rotationAngles_[2][6] = PI/4.0;
    rotationAngles_[0][7] = 0.9557;
    rotationAngles_[1][7] = 0;
    rotationAngles_[2][7] = PI/4.0 + PI/2.0;
    rotationAngles_[0][8] = 0.9557;
    rotationAngles_[1][8] = 0;
    rotationAngles_[2][8] = PI/4.0 + PI;
    rotationAngles_[0][9] = 0.9557;
    rotationAngles_[1][9] = 0;
    rotationAngles_[2][9] = PI/4.0 + 3.0*PI/2.0;
    rotationAngles_[0][10] = PI - 0.95570;
    rotationAngles_[1][10] = 0;
    rotationAngles_[2][10] = PI/4.0;
    rotationAngles_[0][11] = PI - 0.9557;
    rotationAngles_[1][11] = 0;
    rotationAngles_[2][11] = PI/4.0 + PI/2.0;
    rotationAngles_[0][12] = PI - 0.9557;
    rotationAngles_[1][12] = 0;
    rotationAngles_[2][12] = PI/4.0 + PI;
    rotationAngles_[0][13] = PI - 0.9557;
    rotationAngles_[1][13] = 0;
    rotationAngles_[2][13] = PI/4.0 + 3*PI/2.0;


    // get the model parameters
    piezoMat_->GetScalar( sponP0_,  SPON_POLARIZATION, Global::REAL);
    piezoMat_->GetScalar( sponS0_,  SPON_STRAIN, Global::REAL);
    piezoMat_->GetScalar( E0_,  EFIELD0, Global::REAL);
    piezoMat_->GetScalar( sigma0_,  STRESS0, Global::REAL);
    piezoMat_->GetScalar( d0Couple_,  DCOUPLE0, Global::REAL);
    piezoMat_->GetScalar( rateConst_,  RATE_CONSTANT, Global::REAL);
    piezoMat_->GetScalar( viscoPlasticIdx_,  VISCO_PLASTIC_INDEX, Global::REAL);
    piezoMat_->GetScalar( saturationIdx_,  SATURATION_INDEX, Global::REAL);
    piezoMat_->GetScalar( volFracInit_, VOLUME_FRAC_INIT, Global::REAL);
    piezoMat_->GetScalar( scaleDriveForceElec_, SCALE_FORCE_ELEC, Global::REAL);
    piezoMat_->GetScalar( scaleDriveForceMech_, SCALE_FORCE_MECH, Global::REAL);
    piezoMat_->GetScalar( scaleDriveForceCouple_, SCALE_FORCE_COUPLE, Global::REAL);
    piezoMat_->GetScalar( meanTemp_, MEAN_TEMPERATURE, Global::REAL);

    // get the material tensor
    Matrix<Double> cTensor;
    mechMat_->GetTensor( cTensor, MECH_STIFFNESS_TENSOR, Global::REAL, tensorType_ );
    Matrix<Double> eTensor;
    piezoMat_->GetTensor( eTensor, PIEZO_TENSOR, Global::REAL, tensorType_ );
    Matrix<Double> epsTensor_cStrain, eTensorTrans;
    elecMat_->GetTensor( epsTensor_cStrain, ELEC_PERMITTIVITY, Global::REAL, tensorType_ );

    cTensor.Invert(sTensorOrig_);
    dTensorOrig_   = eTensor*sTensorOrig_;
    eTensor.Transpose(eTensorTrans);
    epsTensorOrig_ = epsTensor_cStrain + dTensorOrig_*eTensorTrans;

//      std::cout << "sTensorOrig: \n " << sTensorOrig_ << std::endl;
//    std::cout << "dTensorOrig: \n " << dTensorOrig_ << std::endl;
//      std::cout << "epsTensorOrig: \n " << epsTensorOrig_ << std::endl;


    // set the effective tensors (just for correct dimensions)
    effMechStiffTensor_      = sTensorOrig_;
    effMechComplianceTensor_ = sTensorOrig_;
    effPiezoTensor_          = dTensorOrig_;
    effPermittivityTensor_   = epsTensorOrig_;

    //initialize the properties of the material tensors, spontenous polarization
    //and strain in the different directions
    Matrix<Double> baseS(dim_,dim_);
    baseS.Init();
    baseS[0][0] = -0.5*sponS0_;
    baseS[1][1] = -0.5*sponS0_;
    baseS[2][2] = sponS0_;

    Vector<Double> baseP(dim_);
    baseP.Init();
    baseP[2] = sponP0_;

    //define the array dimensions
    Ps_.resize(extents[numDomain_]);
    Ss_.resize(extents[numDomain_]);
    dTensor_.resize(extents[numDomain_]);
    sTensor_.resize(extents[numDomain_]);
    epsTensor_.resize(extents[numDomain_]);

    for ( UInt i=0; i<numDomain_; i++ ) {
      Vector<Double> angle(dim_);
      Matrix<Double> rotMat(dim_,dim_);

      //get angle for direction
      for ( UInt j=0; j<dim_; j++ )
        angle[j] = rotationAngles_[j][i];

      //std::cout << "Angle: \n " << angle << std::endl << std::endl;
      //get rotation matrix
      ComputeRotationMatrix( angle, rotMat );
      //std::cout << "rotMatrix:\n " << rotMat << std::endl;

      //spontenous polarization
      Ps_[i].Resize(dim_);
      Ps_[i] = rotMat * baseP;
      //std::cout << "Ps_: " << i << "\n" << Ps_[i] << "\n" << std::endl;

      //spontenous strain
      Ss_[i].Resize(dimS_);
      Matrix<Double> tmpS;
      RotateMatrix( baseS, tmpS, rotMat );
      //std::cout << "Ss_: " << i << "\n" << tmpS << std::endl;
      ConvertToVoigtNotation( tmpS, Ss_[i] );
      //std::cout << "Ss_: " << i << "\n" << Ss_[i] << "\n" << std::endl;

      //coupling tensor
      RotateMatrix( dTensorOrig_, dTensor_[i], rotMat );
      //std::cout << "dTensor_: " << i << "\n" << dTensor_[i] << "\n" << std::endl;

      //complicance tensor
      RotateMatrix( sTensorOrig_, sTensor_[i], rotMat );
      //std::cout << "sTensor_: " << i << "\n" << sTensor_[i] << std::endl;

      //dielectric tensor
      RotateMatrix( epsTensorOrig_, epsTensor_[i], rotMat );
      // std::cout << "epsTensor_: " << i << "\n" << epsTensor_[i] << std::endl;
    }

    //compute the delta of the switching system
    deltaPs_.resize(extents[numDomain_][numDomain_]);
    deltaSponS_.resize(extents[numDomain_][numDomain_]);; 
    deltaTensorCoupl_.resize(extents[numDomain_][numDomain_]);; 

    for ( UInt i=0; i<numDomain_; i++ ) {
      for ( UInt j=0; j<numDomain_; j++ ) {
        // std::cout << "\n i:" << i << "  j:" << j << std::endl;
        deltaPs_[i][j] = Ps_[j] - Ps_[i];
        deltaSponS_[i][j] = Ss_[j] - Ss_[i];
        //std::cout << "deltaS:\n" << deltaSponS_[i][j] << std::endl << std::endl;
        deltaTensorCoupl_[i][j] = dTensor_[j] - dTensor_[i];
        //std::cout << "deltaCoup:\n" << deltaTensorCoupl_[i][j] << std::endl;
      }
    }

    //compute the switching values
    switchingVal_.Resize(numDomain_,numDomain_);

    for ( UInt i=0; i<numDomain_; i++ ) {
      for ( UInt j=0; j<numDomain_; j++ ) {
        Double vecInner, val;
        Ps_[i].Inner( Ps_[j],  vecInner );
        val = vecInner / ( Ps_[i].NormL2() * Ps_[j].NormL2() );
        if ( val > 1.0 ) 
          val = 1.0;

        if ( val < -1.0 ) 
          val = -1.0;
             
        switchingVal_[i][j] = 2.0 * std::acos( val ) / PI;
      }
    }

    //std::cout << "switchingVal:\n" << switchingVal_ << std::endl;

    // allocate for driving forces and transformation rates
    driveForcesElec_.resize(extents[numEl_][numDomain_][numDomain_]); 
    driveForcesMech_.resize(extents[numEl_][numDomain_][numDomain_]);
    driveForcesCouple_.resize(extents[numEl_][numDomain_][numDomain_]);
    transformationRates_.resize(extents[numEl_][numDomain_][numDomain_]);

    // allocate for volume fractions 
    Double startVolFrac = 1.0/Double(numDomain_);
    volFracAct_.Resize(numEl_,numDomain_);
    volFracPrev_.Resize(numEl_,numDomain_);  
    volFracAct_.InitValue( startVolFrac );
    volFracPrev_.InitValue( startVolFrac );

    // allocate for electric polarization
    effElecPolAct_.resize(extents[numEl_]);
    effElecPolPrev_.resize(extents[numEl_]);
    for ( UInt i=0; i<numEl_; i++ ) {
      effElecPolAct_[i].Resize(dim_);
      effElecPolPrev_[i].Resize(dim_);
      effElecPolPrev_[i].Init();
    }

    //allocate for irreversibel strain
    effStrainIrrAct_.resize(extents[numEl_]);
    effStrainIrrPrev_.resize(extents[numEl_]);
    for ( UInt i=0; i<numEl_; i++ ) {
      effStrainIrrAct_[i].Resize(dimS_);
      effStrainIrrPrev_[i].Resize(dimS_);
      effStrainIrrPrev_[i].Init();
    }

  }

  void PiezoMicroModelBK::GetEffectiveTensors( Matrix<Double>& matMechC,
                                               Matrix<Double>& matMechS,
                                               Matrix<Double>& matElec,
                                               Matrix<Double>& matPiezo,
                                               Vector<Double>& stress, 
                                               Vector<Double>& elecField,
                                               UInt elemIdx,
                                               bool recompute,
                                               bool previous ) {

    ComputeEffectiveTensors( stress, elecField, elemIdx, 
                             recompute, previous );

    matMechC = effMechStiffTensor_;
    matMechS = effMechComplianceTensor_;
    matElec  = effPermittivityTensor_;  //eps-tensor a constant mechanical stress!!
    matPiezo = effPiezoTensor_; // d-Tensor

    // sTensorOrig_.Invert(matMechC);
    //    matMechS = sTensorOrig_;
//     matPiezo = dTensorOrig_;
//     if ( previous ) {
//       Double scale = effElecPolPrev_[elemIdx][2] / sponP0_;
//       matPiezo *= scale;
//         }
//     else {
//       Double scale = effElecPolAct_[elemIdx][2] / sponP0_;
//       matPiezo *= scale;
//     }
    //matElec   = epsTensorOrig_;
  }


  void PiezoMicroModelBK::ComputeEffectiveTensors( Vector<Double>& stress, 
                                                   Vector<Double>& elecField, 
                                                   UInt elemIdx,
                                                   bool recompute,
                                                   bool previous ) {

    //  std::cout << "Efield:\n" << elecField << "\n" << std::endl;

    // just to be sure, that everything goes right
    // programmer may not be aware of this and does quite wrong things!!
    if ( previous )
      recompute = false;

    if ( recompute ) {
      // 1) compute the volume fractions
      ComputeVolumeFractionsExplicit( stress, elecField, elemIdx );
    }

    // 2) compute effective tensors
    effMechComplianceTensor_.Init();
    effPiezoTensor_.Init();
    effPermittivityTensor_.Init();

    for ( UInt i=0; i<numDomain_; i++ ) {
      Matrix<Double>& matd   = dTensor_[i];
      Matrix<Double>& mats   = sTensor_[i];
      Matrix<Double>& mateps = epsTensor_[i];

      if ( previous ) {
        effMechComplianceTensor_ += mats * volFracPrev_[elemIdx][i];
        effPiezoTensor_          += matd * volFracPrev_[elemIdx][i];
        effPermittivityTensor_   += mateps * volFracPrev_[elemIdx][i];
      }
      else {
        effMechComplianceTensor_ += mats * volFracAct_[elemIdx][i];
        effPiezoTensor_          += matd * volFracAct_[elemIdx][i];
        effPermittivityTensor_   += mateps * volFracAct_[elemIdx][i];
      }
    }

    effMechComplianceTensor_.Invert( effMechStiffTensor_ ); 
  }


  void PiezoMicroModelBK::GetEffectiveIrreversibleValues( Vector<Double>& effPirr,
                                                          Vector<Double>& effSirr,
                                                          UInt elemIdx,
                                                          bool recompute,
                                                          bool previous) {

    if (  previous ) {
      effPirr = effElecPolPrev_[elemIdx];
      effSirr = effStrainIrrPrev_[elemIdx];
    }
    else {
      if ( recompute ) {
        //recompute 
        ComputeEffectiveIrreversibleValues( elemIdx );
      }

      effPirr = effElecPolAct_[elemIdx];
      effSirr = effStrainIrrAct_[elemIdx];
    }
    //   effPirr.Init();
    //      effSirr.Init(); 
  }

  void PiezoMicroModelBK::ComputeEffectiveIrreversibleValues( UInt elemIdx ) {

    effElecPolAct_[elemIdx].Init();
    effStrainIrrAct_[elemIdx].Init();

    for ( UInt i=0; i<numDomain_; i++ ) {
      effElecPolAct_[elemIdx]   += Ps_[i] * volFracAct_[elemIdx][i]; 
      //      effStrainIrrAct_[elemIdx] += Ss_[i] * volFracAct_[elemIdx][i]; 
    }

    Double cx, cy, cz, len, lenxy;
    cx  = effElecPolAct_[elemIdx][0];
    cy  = effElecPolAct_[elemIdx][1];
    cz  = effElecPolAct_[elemIdx][2];
    len = effElecPolAct_[elemIdx].NormL2();

    Double threshold =  sponP0_*1e-6;

    if ( std::abs(cx) < threshold ) 
      cx = 0;
    if ( std::abs(cy) < threshold ) 
      cy = 0;
    if ( std::abs(cz) < threshold ) 
      cz = 0;

    Vector<Double> angle(dim_);
    Matrix<Double> rotMat(dim_,dim_);
    Matrix<Double> tmpS;

    lenxy = std::sqrt(cx*cx + cy*cy);
    angle[0] = 0.0;
    angle[1] = std::atan2(cz,lenxy); //
    angle[2] = std::atan2(cy, cx); //

    //get rotation matrix
    ComputeRotationMatrix( angle, rotMat );

    Matrix<Double> baseS(dim_,dim_);
    baseS.Init();
    baseS[0][0] = sponS0_;
    baseS[1][1] = -0.5*sponS0_;
    baseS[2][2] = -0.5*sponS0_;
    RotateMatrix( baseS, tmpS, rotMat );

 //    std::cout << "Pirr: \n " << effElecPolAct_[elemIdx]  << std::endl;
//     std::cout << "Angles: \n " << angle  << std::endl;
//     std::cout << "Rotmat: \n " << rotMat << std::endl;
    ConvertToVoigtNotation( tmpS, effStrainIrrAct_[elemIdx] );
    effStrainIrrAct_[elemIdx] *= len;
    effStrainIrrAct_[elemIdx] /= sponP0_;

//     std::cout << "Pirr:\n " << effElecPolAct_[elemIdx]  << "\n" << std::endl;
//     std::cout << "Sirr:\n " << effStrainIrrAct_[elemIdx]  << "\n" << std::endl;

  }



  void PiezoMicroModelBK::ComputeVolumeFractionsExplicit( Vector<Double>& stress, 
                                                          Vector<Double>& elecField, 
                                                          UInt elemIdx ) {

    //    std::cout << "Efield:\n" << elecField << std::endl;
    // 1) Compute the driving forces
    ComputeDrivingForces( stress, elecField, elemIdx);

    // 2) Compute transformation rates
    ComputeTransformationRates( elemIdx );

    // 3) compute volume fractions
    for ( UInt i=0; i<numDomain_; i++ ) {
      Double deltaVolFrac = 0.0;
      for ( UInt j=0; j<numDomain_; j++ ) {
        deltaVolFrac +=   transformationRates_[elemIdx][j][i] 
                        - transformationRates_[elemIdx][i][j] ;
      }
      volFracAct_[elemIdx][i] =  volFracPrev_[elemIdx][i] +  deltaT_ * deltaVolFrac;
      //std::cout << "VolFrac: \n" << volFracAct_[elemIdx][i] << std::endl;
    }
    //  std::cout << std::endl;

  }
  

  void PiezoMicroModelBK::ComputeDrivingForces( Vector<Double>& stress, 
                                                Vector<Double>& elecField,
                                                UInt elemIdx ) {

    Double partElec, partMech, partCouple;

    for ( UInt i=0; i<numDomain_; i++ ) {
      for ( UInt j=0; j<numDomain_; j++ ) {
        // electric part
        Vector<Double>& dPs = deltaPs_[i][j];
        dPs.Inner( elecField, partElec );
        driveForcesElec_[elemIdx][i][j] = partElec * scaleDriveForceElec_;

        //mechanical part
        Vector<Double>& dSs = deltaSponS_[i][j];
        dSs.Inner( stress, partMech );
        driveForcesMech_[elemIdx][i][j] = partMech * scaleDriveForceMech_;

        //coupling part
        Vector<Double> vecHelp;
        vecHelp = deltaTensorCoupl_[i][j] * stress;
        vecHelp.Inner( elecField, partCouple );
        driveForcesCouple_[elemIdx][i][j] = partCouple * scaleDriveForceCouple_;
        // std::cout << "eForce: " << driveForcesElec_[elemIdx][i][j] << std::endl;
      }
    }
  }


  void PiezoMicroModelBK::ComputeTransformationRates( UInt elemIdx ) {
 
    Double drivingForceElecStar, drivingForceMechStar, drivingForceCoupleStar;  
    Double transRateElec, transRateMech, transRateCouple;
    Double exp1, val;

    for ( UInt i=0; i<numDomain_; i++ ) {
      for ( UInt j=0; j<numDomain_; j++ ) {
        drivingForceElecStar   = switchingVal_[i][j] * E0_ * sponP0_;
        drivingForceMechStar   = switchingVal_[i][j] * sigma0_ * sponS0_;
        drivingForceCoupleStar = switchingVal_[i][j] * E0_ * sigma0_ * d0Couple_;

        //limit the driving forces
        if ( driveForcesElec_[elemIdx][i][j] > 1.1*drivingForceElecStar ) 
          driveForcesElec_[elemIdx][i][j] = 1.1*drivingForceElecStar;

        //compute transformation rates
        if ( abs( switchingVal_[i][j] ) < 1e-10 ) {
           transformationRates_[elemIdx][i][j] = 0; 
        }
        else {
          //electric part
          exp1 = 1.0 - ( driveForcesElec_[elemIdx][i][j] / drivingForceElecStar );
          val  = std::exp(exp1);
          transRateElec = rateConst_ * std::pow( val, viscoPlasticIdx_ )
                                     * std::pow( volFracPrev_[elemIdx][i], saturationIdx_); 
          //std::cout << "rateElec:" << transRateElec << std::endl;

          //mechanical part
          exp1 = 1.0 - ( driveForcesMech_[elemIdx][i][j] / drivingForceMechStar );
          val  = std::exp(exp1);
          transRateMech = rateConst_ * std::pow( val, viscoPlasticIdx_ )
                                     * std::pow( volFracPrev_[elemIdx][i], saturationIdx_); 
          //std::cout << "rateMech:" << transRateMech << std::endl;

          //coupling part part
          exp1 = 1.0 - ( driveForcesCouple_[elemIdx][i][j] / drivingForceCoupleStar );
          val  = std::exp(exp1);
          transRateCouple = rateConst_ * std::pow( val, viscoPlasticIdx_ )
                                     * std::pow( volFracPrev_[elemIdx][i], saturationIdx_); 
          //std::cout << "rateCouple:" << transRateCouple << std::endl;
          //overall rate
          transformationRates_[elemIdx][i][j] = transRateElec + transRateMech + transRateCouple;
//           std::cout << "i:" << i << "  j:" << j << std::endl;
//           std::cout << "rateElec:" << transRateElec << std::endl;
        }
      }
    }

    //limit transformation rates, so that volume fractions are not smaller than 0
    Double part;
    for ( UInt i=0; i<numDomain_; i++ ) {
      Double tmp = 0.0;
      for ( UInt j=0; j<numDomain_; j++ ) {
        tmp += transformationRates_[elemIdx][i][j];
      }

      if ( (volFracPrev_[elemIdx][i] - tmp*deltaT_) < 0 ) {
        for ( UInt j=0; j<numDomain_; j++ ) {
          part = volFracPrev_[elemIdx][i] / ( tmp * deltaT_ ); 
          transformationRates_[elemIdx][i][j] *= part;
        }
      } 
    }

  }

  void PiezoMicroModelBK::ComputeEffectiveCouplingTensor(Matrix<Double>& dMatEff, 
                                                         Vector<Double>& elecFieldAct,
                                                         Vector<Double>& elecFieldPrev,
                                                         UInt elemIdx) {
 
   dMatEff = dTensor_[0];  //dTensor in x-Direction;

   Vector<Double> diffE, diffP, Pact;
   Double scaleVal, scaleDiffVal, addVal;

   diffP = effElecPolAct_[elemIdx] - effElecPolPrev_[elemIdx];;
   scaleDiffVal = diffP.NormL2() / sponP0_;
   scaleVal = effElecPolAct_[elemIdx].NormL2() / sponP0_; 
   dMatEff *= scaleVal;

   diffE    = elecFieldAct - elecFieldPrev;

   if ( diffE.NormL2() > 1e3 ) {
     dMatEff[0][0] +=  scaleDiffVal / diffE.NormL2();
     dMatEff[0][1] -=  0.5*scaleDiffVal / diffE.NormL2();
     dMatEff[0][2] -=  0.5*scaleDiffVal / diffE.NormL2();
   }

   //now, we have to rotate the tensor
   Double cx, cy, cz, len, lenxy;
    cx  = effElecPolAct_[elemIdx][0];
    cy  = effElecPolAct_[elemIdx][1];
    cz  = effElecPolAct_[elemIdx][2];
    len = effElecPolAct_[elemIdx].NormL2();

    Double threshold =  sponP0_*1e-6;

    if ( std::abs(cx) < threshold ) 
      cx = 0;
    if ( std::abs(cy) < threshold ) 
      cy = 0;
    if ( std::abs(cz) < threshold ) 
      cz = 0;

    Vector<Double> angle(dim_);
    Matrix<Double> rotMat(dim_,dim_);
    Matrix<Double> tmpS;

    lenxy = std::sqrt(cx*cx + cy*cy);
    angle[0] = 0.0;
    angle[1] = std::atan2(cz,lenxy); //
    angle[2] = std::atan2(cy, cx); //

    //get rotation matrix
    ComputeRotationMatrix( angle, rotMat );

    RotateMatrix(dMatEff, tmpS, rotMat );
    dMatEff = tmpS;

 }


  void PiezoMicroModelBK::ComputeRotationMatrix( Vector<Double>& angles, 
                                                 Matrix<Double>& rotMat ) {

    Matrix<Double> rotX(3,3), rotY(3,3), rotZ(3,3);
    rotX.Init();
    rotX[0][0] = 1.0;
    rotX[1][1] = std::cos(angles[0]);
    rotX[1][2] = -std::sin(angles[0]);
    rotX[2][1] = std::sin(angles[0]);
    rotX[2][2] = std::cos(angles[0]);

    rotY.Init();
    rotY[0][0] = std::cos(angles[1]);
    rotY[0][2] = -std::sin(angles[1]);
    rotY[1][1] = 1.0;
    rotY[2][0] = std::sin(angles[1]);
    rotY[2][2] = std::cos(angles[1]);

    rotZ.Init();
    rotZ[0][0] = std::cos(angles[2]);
    rotZ[0][1] = -std::sin(angles[2]);
    rotZ[1][0] = std::sin(angles[2]);
    rotZ[1][1] = std::cos(angles[2]);
    rotZ[2][2] = 1.0;

    Matrix<Double> help;
    help = rotX*rotY;
    rotMat = rotZ*help;
  }

  void PiezoMicroModelBK::RotateMatrix( Matrix<Double>& inMat, 
                                        Matrix<Double>& outMat,
                                        Matrix<Double>& trans) {
    
    UInt row = inMat.GetNumRows();
    UInt col = inMat.GetNumCols();

    outMat.Resize(row,col);

    if ( row == dim_ && col == dim_ ) {
      //3x3 Matrix to be rotated
      for (UInt i=0; i<dim_; i++) {
        for (UInt j=0; j<dim_; j++) {
          Double tmp = 0.0;
          for (UInt m=0; m<dim_; m++ ) {
            for ( UInt n=0; n<dim_; n++ ) {
              tmp += trans[i][m] * trans[j][n] * inMat[m][n];
            }
          }
          outMat[i][j] = tmp;
        }
      }
    }
    else {
      Matrix<UInt> ten_to_mat(3,3);
      ten_to_mat[0][0]=0;
      ten_to_mat[1][1]=1;
      ten_to_mat[2][2]=2;
      ten_to_mat[1][2]=3;
      ten_to_mat[0][2]=4;
      ten_to_mat[0][1]=5;
      ten_to_mat[2][1]=3;
      ten_to_mat[2][0]=4;
      ten_to_mat[1][0]=5;
      
      Vector<UInt> mat_to_ten_first(6); 
      mat_to_ten_first[0]=0;
      mat_to_ten_first[1]=1;
      mat_to_ten_first[2]=2;
      mat_to_ten_first[3]=1;
      mat_to_ten_first[4]=0;
      mat_to_ten_first[5]=0;
      
      Vector<UInt> mat_to_ten_sec(6);
      mat_to_ten_sec[0]=0;
      mat_to_ten_sec[1]=1;
      mat_to_ten_sec[2]=2;
      mat_to_ten_sec[3]=2;
      mat_to_ten_sec[4]=2;
      mat_to_ten_sec[5]=1;

      if ( row == 3 ) {
        //3x6 Matrix to be rotated
        
        arrayD3 ten, trans_ten;
        ten.resize(extents[3][3][3]);
        trans_ten.resize(extents[3][3][3]);
        
        for (UInt i=0; i<dim_; i++) {
          for (UInt j=0; j<dim_; j++) {
            for (UInt k=0; k<dim_; k++) {
              ten[i][j][k] = inMat[i][ten_to_mat[j][k]];
            }
          }
        }

        for (UInt ii=0; ii<dim_; ii++) {
          for (UInt jj=0; jj<dim_; jj++) {
            for (UInt kk=0; kk<dim_; kk++) {
              Double tmp = 0.0;
              for (UInt mm=0; mm<dim_; mm++) {
                for (UInt nn=0; nn<dim_; nn++) {
                  for (UInt oo=0; oo<dim_; oo++) {
                    tmp += trans[ii][mm]*trans[jj][nn]*trans[kk][oo]*ten[mm][nn][oo];
                  }
                }
              }
              trans_ten[ii][jj][kk]=tmp;
            }
          }
        }

        for (UInt i=0; i<dim_; i++) {
          for (UInt j=0; j<6; j++) {
            outMat[i][j] = trans_ten[i][mat_to_ten_first[j]][mat_to_ten_sec[j]]; 
          }
        }
      }
      else {
        //6x6 Matrix to be rotated
        arrayD4 ten, trans_ten;
        ten.resize(extents[3][3][3][3]);
        trans_ten.resize(extents[3][3][3][3]);
        
        for (UInt i=0; i<dim_; i++) {
          for (UInt j=0; j<dim_; j++) {
            for (UInt k=0; k<dim_; k++) {
              for (UInt l=0; l<dim_; l++) {
                ten[i][j][k][l] = inMat[ ten_to_mat[i][j] ][ ten_to_mat[k][l] ]; 
              }
            }
          }
        }
        
        for (UInt i=0; i<dim_; i++) {
          for (UInt j=0; j<dim_; j++) {
            for (UInt k=0; k<dim_; k++) {
              for (UInt l=0; l<dim_; l++) {
                Double tmp = 0.0;
                for (UInt m=0; m<dim_; m++) {
                  for (UInt n=0; n<dim_; n++) {
                    for (UInt o=0; o<dim_; o++) {
                      for (UInt p=0; p<dim_; p++) {
                        tmp += trans[i][m]*trans[j][n]*trans[k][o]*trans[l][p]*ten[m][n][o][p];
                      }
                    }
                  }
                }
                trans_ten[i][j][k][l]=tmp;
              }
            }
          }
        }
        
        
        for (UInt i=0; i<6; i++) {
          for (UInt j=0; j<6; j++) {
            outMat[i][j]=trans_ten[ mat_to_ten_first[i] ][mat_to_ten_sec[i] ][mat_to_ten_first[j] ][mat_to_ten_sec[j] ]; 
          }
        }
        
      }
    }
  }
  
    //! 
  void PiezoMicroModelBK::ConvertToVoigtNotation( Matrix<Double>& inMat, 
                                                  Vector<Double>& outVec ) {

    outVec[0] = inMat[0][0];
    outVec[1] = inMat[1][1];
    outVec[2] = inMat[2][2];
    outVec[3] = 2.0*inMat[1][2];
    outVec[4] = 2.0*inMat[0][2];
    outVec[5] = 2.0*inMat[0][1];

  }

} // end of namespace
