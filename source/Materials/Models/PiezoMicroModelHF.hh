#ifndef PIEZO_MICRO_MODELL_HF_HH
#define PIEZO_MICRO_MODELL_HF_HH

#include <string>

#include "Domain/Mesh/Grid.hh"
#include "Utils/ToolsFull.hh"
#include "MatVec/Vector.hh"
#include "General/Environment.hh"
#include "Materials/BaseMaterial.hh"

namespace CoupledField {

  // forward class declarations
  

  //! class for Huber-Fleck micro-piezoelectric model
  class PiezoMicroModelHF {

  public:

    //! Default constructor
    PiezoMicroModelHF( UInt numElemSD, BaseMaterial* piezoMat, 
                       BaseMaterial* mechMat, 
                       BaseMaterial* elecMat, 
                       SubTensorType tensorType, 
                       Double dt);

    //! Destructor
    virtual ~PiezoMicroModelHF();
    
    //! Init the switching system
    void InitSwitchingSystem();

    //! get the effective tensor
    void GetEffectiveTensors( Matrix<Double>& matMech,
                              Matrix<Double>& matElec,
                              Matrix<Double>& matPiezo,
                              Vector<Double>& stress, 
                              Vector<Double>& elecField,
                              UInt elemIdx, 
                              bool recompute,
                              bool previous );

    //!
    void GetEffectiveIrreversibleValues( Vector<Double>& Pirr,
                                         Vector<Double>& Sirr,
                                         UInt elemIdx,
                                         bool recompute, 
                                         bool previous);

    //!
    void ComputeEffectiveIrreversibleValues( UInt elemIdx );
 
    //!
    void SetPreviousVolFrac() {
      volFracPrev_ =  volFracAct_;
    }

    //!
    void SetPreviousIrreversibleValues() {
      effElecPolPrev_ = effElecPolAct_;
      effStrainIrrPrev_ = effStrainIrrAct_;
    }

    //!
    void ComputeEffectiveTensors( Vector<Double>& stress, Vector<Double>& elecField, 
                                  UInt elemIdx, bool recompute, bool previous );

    //! computes the volume fractions in an explicit way
    void ComputeVolumeFractionsExplicit( Vector<Double>& stress, 
                                         Vector<Double>& elecField, 
                                         UInt elemIdx );


    //!
    void ComputeDrivingForces( Vector<Double>& stress, Vector<Double>& elecField,
                               UInt elemIdx );

    //! rotate the material tensors towards the domain systems
    void RotateMaterialTensors();

    //! compute the Schmid tensor
    void ComputeSchmidTensor( UInt system );

    //! computes difference of two coupling tensors
    void ComputeDeltaCouplingTensor( UInt system1, UInt system2 );


    //! 
    void ComputeTransformationRates( UInt elemIdx );

  protected:

    SubTensorType tensorType_;  //! type of tensor (FULL, axi, ..)
    UInt dim_;  //! dimension of the problem (2D or 3D)

    UInt numEl_;  //! number of finite elements 
    UInt numDomain_; //!  number of domain types
    UInt numSwitching_; //! number of switching systems

    Double deltaT_;  //! time step

    Double sponP_;           //! spontaneous polarization
    Double sponS_;           //! spontaneous strain
    Double driveForce90_;    //! driving force for 90 degree switching
    Double driveForce180_;   //! driving force for 180 degree switching
    Double rateConst_;       //! rate constant
    Double viscoPlasticIdx_; //! exponent for viscoplastic model
    Double saturationIdx_;    //! stauration exponenent for viscoplastic model
    Double volFracInit_;     //! initial value for volume fraction
    Double meanTemp_;        //! mean temperature

    Vector<Double> switchPolarizationVal_; //! value of polarization for each switching system
    Vector<Double> switchStrainVal_;  //! value of strain for each switching system

    Matrix<Double> switchingDirection_;
    Matrix<Double> switchingNormal_;
    Matrix<UInt>   switchingSystems_;

    Matrix<Double> dTensorOrig_;
    Matrix<Double> sTensorOrig_;
    Matrix<Double> epsTensorOrig_;

    Matrix<Double>* dTensors4SwitchingSystem_;
    Matrix<Double>* sTensors4SwitchingSystem_;
    Matrix<Double>* epsTensors4SwitchingSystem_;

    Matrix<Double> rotationAngles_;

    Matrix<Double> SchmidTensor_;
    Vector<Double> SchmidTensorAsVector_;
    Matrix<Double> deltaCouplingTensor_;
    Matrix<Double> connectSystemMatrix_;

    Matrix<Double> driveForces_;
    Matrix<Double> transformationRates_;
    Matrix<Double> volFracAct_;
    Matrix<Double> volFracPrev_;

    Matrix<Double> effMechStiffTensor_;
    Matrix<Double> effMechComplianceTensor_;
    Matrix<Double> effPiezoTensor_;
    Matrix<Double> effPermittivityTensor_;

    Matrix<Double> effElecPolAct_;
    Matrix<Double> effElecPolPrev_;
    Matrix<Double> effStrainIrrAct_;
    Matrix<Double> effStrainIrrPrev_;

    //pointers to the material objects
    BaseMaterial* piezoMat_;
    BaseMaterial* mechMat_;
    BaseMaterial* elecMat_;
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class PiezoMicroModelBK
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif
