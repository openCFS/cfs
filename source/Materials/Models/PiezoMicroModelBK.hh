#ifndef PIEZO_MICRO_MODELL_BK_HH
#define PIEZO_MICRO_MODELL_BK_HH

#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/multi_array.hpp>
#pragma GCC diagnostic pop

#include "Domain/Mesh/Grid.hh"
#include "Utils/ToolsFull.hh"
#include "MatVec/Vector.hh"
#include "General/Environment.hh"
#include "Materials/BaseMaterial.hh"
#include "ODESolve/BaseODESolver.hh" 
#include "ODESolve/BaseODEProblem.hh"

typedef boost::multi_array<Matrix<double>, 1> array1_Matrix;
typedef boost::multi_array<Matrix<double>, 2> array2_Matrix;
typedef boost::multi_array<Vector<double>, 1> array1_Vector;
typedef boost::multi_array<Vector<double>, 2> array2_Vector;
typedef boost::multi_array<Double, 3> arrayD3;
typedef boost::multi_array<Double, 4> arrayD4;

namespace CoupledField {

  // forward class declarations
  

  //! class for Belov-Kreher micro-piezoelectric model
  class PiezoMicroModelBK {

  public:

    //! Default constructor
    PiezoMicroModelBK( UInt numElemSD, BaseMaterial* piezoMat, 
                       BaseMaterial* mechMat, 
                       BaseMaterial* elecMat, 
                       SubTensorType tensorType, 
                       Double dt);

    //! Destructor
    virtual ~PiezoMicroModelBK();
    
    //! Init the switching system
    void InitSwitchingSystem();

    //! get the effective tensor
    void GetEffectiveTensors( Matrix<Double>& matMechC,
                              Matrix<Double>& matMechS,
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
    void ComputeEffectiveCouplingTensor(Matrix<Double>& dmatEff, 
                                        Vector<Double>& elecFieldAct,
                                        Vector<Double>& elecFieldPrev,
                                        UInt elemIdx);

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

    //! computes the volume fractions in an explicit way
    void ComputeVolumeFractionsImplicit( Vector<Double>& stress, 
                                         Vector<Double>& elecField, 
                                         UInt elemIdx );

    //!
    void ComputeDrivingForces( Vector<Double>& stress, Vector<Double>& elecField,
                               UInt elemIdx );

    //! 
    void ComputeTransformationRates( UInt elemIdx );

    //! 
    void ComputeRotationMatrix( Vector<Double>& angles, 
                                Matrix<Double>& rotMat );

    //!
    void RotateMatrix( Matrix<Double>& inMat, 
                       Matrix<Double>& outMat,
                       Matrix<Double>& rotMat);

    //! 
    void ConvertToVoigtNotation( Matrix<Double>& inMat, 
                                 Vector<Double>& outVec );

  protected:

    SubTensorType tensorType_;  //! type of tensor (FULL, axi, ..)
    UInt dim_;  //! dimension of the problem (2D or 3D)
    UInt dimS_;  //! dimension of strain vector in Voigt notation

    UInt numEl_;  //! number of finite elements 
    UInt numDomain_; //!  number of domain types
    UInt numSwitching_; //! number of switching systems

    Double deltaT_;  //! time step

    Double sponP0_;           //! spontaneous polarization
    Double E0_;               //!  charachteristic electric field intensity
    Double sponS0_;           //! spontaneous strain
    Double sigma0_;           //!  charachteristic stress
    Double d0Couple_;         //! charachteristic coupling coefficient
    Double driveForce90_;    //! driving force for 90 degree switching
    Double driveForce180_;   //! driving force for 180 degree switching
    Double rateConst_;       //! rate constant
    Double viscoPlasticIdx_; //! exponent for viscoplastic model
    Double saturationIdx_;    //! stauration exponenent for viscoplastic model
    Double volFracInit_;     //! initial value for volume fraction
    Double meanTemp_;        //! mean temperature

    Double scaleDriveForceElec_;
    Double scaleDriveForceMech_;
    Double scaleDriveForceCouple_;

    Matrix<Double> switchingVal_;

    array2_Vector deltaPs_;
    array2_Vector deltaSponS_; 
    array2_Matrix deltaTensorCoupl_; 

    array1_Vector Ps_;
    array1_Vector Ss_;
    array1_Matrix dTensor_;
    array1_Matrix cTensor_;
    array1_Matrix sTensor_;
    array1_Matrix epsTensor_;

    Matrix<Double> driveForces_;

    Matrix<Double> volFracAct_;
    Matrix<Double> volFracPrev_;

    array1_Vector effElecPolAct_;
    array1_Vector effElecPolPrev_;
    array1_Vector effStrainIrrAct_;
    array1_Vector effStrainIrrPrev_;

    arrayD3 driveForcesElec_;
    arrayD3 driveForcesMech_;
    arrayD3 driveForcesCouple_;
    arrayD3 transformationRates_;

    Matrix<Double> dTensorOrig_;
    Matrix<Double> cTensorOrig_;
    Matrix<Double> sTensorOrig_;
    Matrix<Double> epsTensorOrig_;
    Matrix<Double> effMechStiffTensor_;
    Matrix<Double> effMechComplianceTensor_;
    Matrix<Double> effPiezoTensor_;
    Matrix<Double> effPermittivityTensor_;

    Matrix<Double> rotationAngles_;

    //pointers to the material objects
    BaseMaterial* piezoMat_;
    BaseMaterial* mechMat_;
    BaseMaterial* elecMat_;

    //! Pointer to ODE for piezo switching
    BaseODEProblem *ptODEPiezo_;

    //! Pointer to the solver class for ode's
    BaseODESolver *ptODESolver_;

    StdVector<Double> yInitOut_;

    bool explicit_;
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
