#include "PiezoMicroModelHF.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField{


  PiezoMicroModelHF::PiezoMicroModelHF( UInt numElemSD, BaseMaterial* piezoMat, 
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

    if ( tensorType_ == FULL )
      dim_ = 3;
    else 
      dim_ = 2;

    InitSwitchingSystem();
  }

  // destructor
  PiezoMicroModelHF::~PiezoMicroModelHF(){
  }


  void PiezoMicroModelHF::InitSwitchingSystem() {

    // set number of domain types
    numDomain_ = 2;

    // set number of switching systems
    numSwitching_ = numDomain_ * (numDomain_ -1);

    // define the switching systems
    switchingSystems_.Resize(numSwitching_, 2);
    switchingSystems_[0][0] = 1;
    switchingSystems_[0][1] = 2;
    switchingSystems_[1][0] = 2;
    switchingSystems_[1][1] = 1;
 
    // define the connectivity matrix for the switching systems
    connectSystemMatrix_.Resize(numDomain_,numSwitching_);
    connectSystemMatrix_[0][0] = -1.0;
    connectSystemMatrix_[0][1] =  1.0;
    connectSystemMatrix_[1][0] =  1.0;
    connectSystemMatrix_[1][1] = -1.0;

    // define the rotation angles
    rotationAngles_.Resize(3,numDomain_);
    rotationAngles_.Init();
    if ( dim_ == 2) 
      rotationAngles_[2][1] = 180;
    else
      rotationAngles_[1][1] = 180;

    // allocate for driving forces, transformation rates and volume fractions
    driveForces_.Resize(numSwitching_, numEl_);
    transformationRates_.Resize(numSwitching_, numEl_);
    volFracAct_.Resize(numSwitching_, numEl_);
    volFracPrev_.Resize(numSwitching_, numEl_);

    Double startVolFrac = 1.0/Double(numDomain_);
    volFracAct_.InitValue( startVolFrac );
    volFracPrev_.InitValue( startVolFrac );

    // allocate for electric polarization
    effElecPolAct_.Resize(dim_,numEl_);
    effElecPolPrev_.Resize(dim_,numEl_);
    effElecPolPrev_.Init();

    //allocate for irreversibel strain
    if (dim_ == 2 ) {
      effStrainIrrAct_.Resize(3,numEl_);
      effStrainIrrPrev_.Resize(3,numEl_);
      effStrainIrrPrev_.Init();
    }
    else {
      effStrainIrrAct_.Resize(6,numEl_);
      effStrainIrrPrev_.Resize(6,numEl_);
      effStrainIrrPrev_.Init();
    }

    // get the model parameters
    piezoMat_->GetScalar( sponP_, SPON_POLARIZATION, Global::REAL);
    piezoMat_->GetScalar( sponS_, SPON_STRAIN, Global::REAL);
    piezoMat_->GetScalar( driveForce90_,  DRIVING_FORCE_90, Global::REAL);
    piezoMat_->GetScalar( driveForce180_,  DRIVING_FORCE_180, Global::REAL);
    piezoMat_->GetScalar( rateConst_,  RATE_CONSTANT, Global::REAL);
    piezoMat_->GetScalar( viscoPlasticIdx_,  VISCO_PLASTIC_INDEX, Global::REAL);
    viscoPlasticIdx_ -= 1.0;

    piezoMat_->GetScalar( saturationIdx_,  SATURATION_INDEX, Global::REAL);
    saturationIdx_ = 1.0 / saturationIdx_;

    piezoMat_->GetScalar( volFracInit_,   VOLUME_FRAC_INIT, Global::REAL);
    piezoMat_->GetScalar( meanTemp_, MEAN_TEMPERATURE, Global::REAL);

    //! compute value of polarization for each switching system
    switchPolarizationVal_.Resize(numSwitching_);
    switchPolarizationVal_[0] = 2.0 * sponP_;
    switchPolarizationVal_[1] = 2.0 * sponP_;

    //! compute value of strain for each switching system
    switchStrainVal_.Resize(numSwitching_);
    switchStrainVal_[0] = 0.0;
    switchStrainVal_[1] = 0.0;

    // allocate for Schmid tensor
    SchmidTensor_.Resize(dim_, dim_);
    if (dim_ == 2 ) 
      SchmidTensorAsVector_.Resize(3);
    else
      SchmidTensorAsVector_.Resize(6);

    // compute switching directions
    switchingDirection_.Resize(dim_, numSwitching_);
    switchingDirection_.Init();
    if (dim_ == 2) {
      switchingDirection_[1][0] = -1;
      switchingDirection_[1][1] =  1;
    }
    else {
      switchingDirection_[2][0] = -1;
      switchingDirection_[2][1] =  1;
    }

    // compute switching normals
    switchingNormal_.Resize(dim_, numSwitching_);
    switchingNormal_.Init();
    if (dim_ == 2) {
      switchingNormal_[0][0] =  1;
      switchingNormal_[0][1] = -1;
    }
    else {
      switchingNormal_[0][0] =  1;
      switchingNormal_[0][1] = -1;
    }

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

    // init the delta coupling tensor (we just need the correct dimensions)
    deltaCouplingTensor_ = dTensorOrig_;

    // rotate the material tensors
    RotateMaterialTensors();

    // set the effective tensors (just for correct dimensions)
    effMechStiffTensor_      = sTensorOrig_;
    effMechComplianceTensor_ = sTensorOrig_;
    effPiezoTensor_          = dTensorOrig_;
    effPermittivityTensor_   = epsTensorOrig_;

  }

  void PiezoMicroModelHF::GetEffectiveTensors( Matrix<Double>& matMech,
                                               Matrix<Double>& matElec,
                                               Matrix<Double>& matPiezo,
                                               Vector<Double>& stress, 
                                               Vector<Double>& elecField,
                                               UInt elemIdx,
                                               bool recompute,
                                               bool previous ) {

    ComputeEffectiveTensors( stress, elecField, elemIdx, 
                             recompute, previous );

    matMech  = effMechStiffTensor_;
    matElec  = effPermittivityTensor_;
    matPiezo = effPiezoTensor_;
  }

  void PiezoMicroModelHF::GetEffectiveIrreversibleValues( Vector<Double>& effPirr,
                                                          Vector<Double>& effSirr,
                                                          UInt elemIdx,
                                                          bool recompute,
                                                          bool previous) {

    UInt dimS = effStrainIrrAct_.GetNumRows();
    Vector<Double> partPirr( effPirr.GetSize() );
    Vector<Double> partSirr( effSirr.GetSize() );

    if (  previous ) {
      for ( UInt i=0; i<dim_; i++)
        effPirr[i] = effElecPolPrev_[i][elemIdx];
      for ( UInt i=0; i<dimS; i++)
        effSirr[i] = effStrainIrrPrev_[i][elemIdx];
    }
    else {
      if ( recompute ) {
        //recompute 
        ComputeEffectiveIrreversibleValues( elemIdx );
      }

      for ( UInt i=0; i<dim_; i++)
        effPirr[i] = effElecPolAct_[i][elemIdx];
      for ( UInt i=0; i<dimS; i++)
        effSirr[i] = effStrainIrrAct_[i][elemIdx];
    }

  }

  void PiezoMicroModelHF::ComputeEffectiveIrreversibleValues( UInt elemIdx ) {

    UInt dimS = effStrainIrrAct_.GetNumRows();
    Vector<Double> partPirr(dim_);
    Vector<Double> partSirr(dimS);
    Double valP, valS;

    partPirr.Init( 0.0 );
    partSirr.Init( 0.0 );

    for ( UInt sys=0; sys< numSwitching_; sys++ ) {
      valP =  switchPolarizationVal_[sys] * transformationRates_[sys][elemIdx];
      for (UInt i=0; i<dim_; i++) {
        partPirr[i] +=  valP * switchingDirection_[i][sys];
      }

      valS = switchStrainVal_[sys] * transformationRates_[sys][elemIdx]; 
      for (UInt i=0; i<dimS; i++) {
         partSirr[i] =  valS * SchmidTensorAsVector_[i];
      }
    }

    for (UInt i=0; i<dim_; i++) {
      effElecPolAct_[i][elemIdx] = effElecPolPrev_[i][elemIdx] 
                                   + deltaT_* partPirr[i];
    }

    for (UInt i=0; i<dimS; i++) {
      effStrainIrrAct_[i][elemIdx] = effStrainIrrPrev_[i][elemIdx] 
                                  + deltaT_* partSirr[i];
    }
//       partSirr.Init( 0.0 );
//       ComputeSchmidTensor(i);
//       for (UInt j=0; j<SchmidTensorAsVector_.GetSize(); j++) {
//         partSirr[j] =  SchmidTensorAsVector_[j] * switchStrainVal_[i];
//       }
//         }

//       if ( previous ) {
//         effPirr += partPirr * volFracPrev_[i][elemIdx];
//         effSirr += partSirr * volFracPrev_[i][elemIdx];
//       }
//       else {
//         effPirr += partPirr * volFracAct_[i][elemIdx];
//         effSirr += partSirr * volFracAct_[i][elemIdx];
//       }
//     }

  }


  void PiezoMicroModelHF::ComputeEffectiveTensors( Vector<Double>& stress, 
                                                   Vector<Double>& elecField, 
                                                   UInt elemIdx,
                                                   bool recompute,
                                                   bool previous ) {

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
      Matrix<Double>& matd   = dTensors4SwitchingSystem_[i];
      Matrix<Double>& mats   = sTensors4SwitchingSystem_[i];
      Matrix<Double>& mateps = epsTensors4SwitchingSystem_[i];

      if ( previous ) {
        effMechComplianceTensor_ += mats * volFracPrev_[i][elemIdx];
        effPiezoTensor_          += matd * volFracPrev_[i][elemIdx];
        effPermittivityTensor_   += mateps * volFracPrev_[i][elemIdx];
      }
      else {
        effMechComplianceTensor_ += mats * volFracAct_[i][elemIdx];
        effPiezoTensor_          += matd * volFracAct_[i][elemIdx];
        effPermittivityTensor_   += mateps * volFracAct_[i][elemIdx];
      }
    }
    effMechComplianceTensor_.Invert( effMechStiffTensor_ ); 
  }


  void PiezoMicroModelHF::ComputeVolumeFractionsExplicit( Vector<Double>& stress, 
                                                          Vector<Double>& elecField, 
                                                          UInt elemIdx ) {

    // 1) Compute the driving forces
    ComputeDrivingForces( stress, elecField, elemIdx);

    // 2) Compute transformation rates
    ComputeTransformationRates( elemIdx );

    // 3) compute volume fractions
    for ( UInt sys=0; sys<numDomain_; sys++ ) {
      volFracAct_[sys][elemIdx] = volFracPrev_[sys][elemIdx];
//       std::cout << "volFracStart:, sys: " << sys << ",  " 
//                 << volFracAct_[sys][elemIdx] << std::endl; 
      
      for ( UInt k=0; k< numSwitching_; k++ ) {
        volFracAct_[sys][elemIdx] += deltaT_ * connectSystemMatrix_[sys][k] 
                                     * transformationRates_[k][elemIdx];
      }
      if ( volFracAct_[sys][elemIdx] < 0.0 ) 
        volFracAct_[sys][elemIdx] = 0.0;

//       std::cout << "volFrac:, sys: " << sys << ",  " 
//                 << volFracAct_[sys][elemIdx] << std::endl; 
      
    }
  }

  void PiezoMicroModelHF::ComputeDrivingForces( Vector<Double>& stress, 
                                                Vector<Double>& elecField,
                                                UInt elemIdx ) {

    Double partElec, partMech, partCouple;

    //    std::cout << "\n ComputeDrivingForces: Ein = " << elecField << std::endl;

    for ( UInt sys=0; sys< numSwitching_; sys++ ) {
      // electric part
      partElec = 0.0;
      for ( UInt i=0; i<dim_; i++ ) {
        partElec += elecField[i] * switchingDirection_[i][sys];
      }
      partElec *=  switchPolarizationVal_[sys];

      // mechanical part
      partMech = 0.0;
      for ( UInt i=0; i<SchmidTensorAsVector_.GetSize(); i++ ) {
        partMech += SchmidTensorAsVector_[i]*stress[i];
      }
      partMech *= switchStrainVal_[sys];

      // coupling part
      partCouple = 0.0;
      ComputeDeltaCouplingTensor( switchingSystems_[sys][0], switchingSystems_[sys][1] );
      Vector<Double> help;
      help = deltaCouplingTensor_ * elecField;
      stress.Inner( help, partCouple);

      driveForces_[sys][elemIdx] = partElec; // + partMech + partCouple;

//       std::cout << "drivingForce, sys:" << sys << ", "
//                 << driveForces_[sys][elemIdx] << std::endl;
    }
  }


  void PiezoMicroModelHF::ComputeTransformationRates( UInt elemIdx ) {
 
    Double frac1, frac2;
    for ( UInt sys=0; sys< numSwitching_; sys++ ) {
      frac1 = driveForces_[sys][elemIdx] / driveForce180_;
      frac2 = volFracPrev_[sys][elemIdx] / volFracInit_;

      transformationRates_[sys][elemIdx] = rateConst_  
                                           * std::abs( std::pow( frac1, viscoPlasticIdx_ ) )
                                           * frac1 * std::pow( frac2, saturationIdx_ );

      if ( transformationRates_[sys][elemIdx] < 0.0 ) 
        transformationRates_[sys][elemIdx] = 0.0;

//       std::cout << "df, sys: " <<  sys << ", "
//                 << transformationRates_[sys][elemIdx] << std::endl;
    }
  }


  void PiezoMicroModelHF::RotateMaterialTensors() {
  
    // allocate
    dTensors4SwitchingSystem_   = new Matrix<Double>[numDomain_];
    sTensors4SwitchingSystem_   = new Matrix<Double>[numDomain_];
    epsTensors4SwitchingSystem_ = new Matrix<Double>[numDomain_];


    // first domain is standard
    Matrix<Double>& dMat = dTensors4SwitchingSystem_[0];
    dMat = dTensorOrig_;
    Matrix<Double>& sMat = sTensors4SwitchingSystem_[0];
    sMat = sTensorOrig_;
    Matrix<Double>& epsMat = epsTensors4SwitchingSystem_[0];
    epsMat = epsTensorOrig_;


    // now go to the others
    Vector<Double> rotAngle(3);  
    for ( UInt i=1; i<numDomain_; i++ ) {
      // get rotationAngles
      for ( UInt k=0; k<3; k++)
        rotAngle[k] = rotationAngles_[k][i];

//       std::cout << "Rotationangles:\n" << rotAngle << std::endl;
      Matrix<Double> eTensorOld, cTensorOld, epsTensorOld;
      piezoMat_->GetTensor( eTensorOld, PIEZO_TENSOR, Global::REAL, tensorType_ );
      mechMat_->GetTensor( cTensorOld, MECH_STIFFNESS_TENSOR, Global::REAL, tensorType_ );
      elecMat_->GetTensor( epsTensorOld, ELEC_PERMITTIVITY, Global::REAL, tensorType_ );
 //      std::cout << "eTensorOld\n" << eTensorOld << std::endl;
//       std::cout << "cTensorOld\n" << cTensorOld << std::endl;
//       std::cout << "epsTensorOld\n" << epsTensorOld << std::endl;

      // perform rotation of material tensors
      piezoMat_->RotateTensorByRotationAngles( rotAngle, PIEZO_TENSOR, false );
      mechMat_->RotateTensorByRotationAngles( rotAngle, MECH_STIFFNESS_TENSOR, false );
      elecMat_->RotateTensorByRotationAngles( rotAngle, ELEC_PERMITTIVITY, false );

      // get the rotated material tensor
      Matrix<Double> cTensor;
      mechMat_->GetTensor( cTensor, MECH_STIFFNESS_TENSOR, Global::REAL, tensorType_ );
      Matrix<Double> eTensor;
      piezoMat_->GetTensor( eTensor, PIEZO_TENSOR, Global::REAL, tensorType_ );
       Matrix<Double> epsTensor_cStrain, eTensorTrans;
      elecMat_->GetTensor( epsTensor_cStrain, ELEC_PERMITTIVITY, Global::REAL, tensorType_ );

//       std::cout << "eTensor\n" << eTensor << std::endl;
//       std::cout << "cTensor\n" << cTensor << std::endl;
//       std::cout << "epsTensor\n" << epsTensor_cStrain << std::endl;

      Matrix<Double>& matd   = dTensors4SwitchingSystem_[i];
      Matrix<Double>& mats   = sTensors4SwitchingSystem_[i];
      Matrix<Double>& mateps = epsTensors4SwitchingSystem_[i];

      cTensor.Invert( mats );
      matd = eTensor*mats;
      eTensor.Transpose(eTensorTrans);
      mateps = epsTensor_cStrain + matd*eTensorTrans;
    }

    // rotate the tensors back;
    rotAngle.Init( 0.0 );
    piezoMat_->RotateTensorByRotationAngles( rotAngle, PIEZO_TENSOR, false );
    mechMat_->RotateTensorByRotationAngles( rotAngle, MECH_STIFFNESS_TENSOR, false );
    elecMat_->RotateTensorByRotationAngles( rotAngle, ELEC_PERMITTIVITY, false );
  }


  void PiezoMicroModelHF::ComputeSchmidTensor( UInt system ) {

    UInt idx = system - 1;

    for ( UInt i=0; i<dim_; i++ ) {
      for ( UInt j=0; j<dim_; j++ ) {
        SchmidTensor_[i][j] = 0.5 * ( switchingDirection_[i][idx] * switchingNormal_[j][idx]
                                     + switchingDirection_[j][idx] * switchingNormal_[i][idx] );
      }
    }
    if ( dim_ == 2 ) {
      SchmidTensorAsVector_[0] = SchmidTensor_[0][0];
      SchmidTensorAsVector_[1] = SchmidTensor_[1][1];
      SchmidTensorAsVector_[2] = 2.0*SchmidTensor_[0][1];
    }
    else if ( dim_ == 3 ) {
      SchmidTensorAsVector_[0] = SchmidTensor_[0][0];
      SchmidTensorAsVector_[1] = SchmidTensor_[1][1];
      SchmidTensorAsVector_[2] = SchmidTensor_[2][2];
      SchmidTensorAsVector_[3] = 2.0*SchmidTensor_[2][1];
      SchmidTensorAsVector_[4] = 2.0*SchmidTensor_[0][2];
      SchmidTensorAsVector_[5] = 2.0*SchmidTensor_[0][1];
    }
  }

  void PiezoMicroModelHF::ComputeDeltaCouplingTensor( UInt system1, UInt system2 ) {

    Matrix<Double>& mat1 = dTensors4SwitchingSystem_[system1-1];
    Matrix<Double>& mat2 = dTensors4SwitchingSystem_[system2-1];
    Matrix<Double> help;
    help = mat2 - mat1;
    help.Transpose(  deltaCouplingTensor_ );
  }

} // end of namespace
