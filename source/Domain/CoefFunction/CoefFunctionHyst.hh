#ifndef FILE_COEFFUNCTION_HYST_HH
#define FILE_COEFFUNCTION_HYST_HH

#include <boost/tr1/type_traits.hpp>
#include "General/Environment.hh"
#include "CoefFunction.hh"
#include "Materials/BaseMaterial.hh"
#include "FeBasis/FeFunctions.hh"
#include "Utils/Timer.hh"
//#include "Materials/Models/Preisach.hh"

namespace CoupledField {

// forward class declaration
class ApproxData;
class BaseBOperator;
class Grid;
class FeFunctions;
class CoefFunctionHyst;

class CoefFunctionHystMatTensor : public CoefFunction{
public:
		
  // note: coupling tensor has to be in non-transposed form (i.e. 3x6 or 2x3, 2x4)
	CoefFunctionHystMatTensor(PtrCoefFct fieldTensor,PtrCoefFct elasticityTensor,
	                          PtrCoefFct couplingTensor, PtrCoefFct CoefFncHyst,
	                          std::string tensorName, bool transposed)
  :CoefFunction(){

    dimType_ = TENSOR;
	  isAnalytic_ = false;
	  isComplex_  = false;
    dependType_ = SOLUTION;
    
	  fieldTensor_ = fieldTensor;
	  elastTensor_ = elasticityTensor;
	  couplTensor_ = couplingTensor;

    UInt couplRows, couplCols;
    if(couplTensor_ != NULL){
      couplTensor_->GetTensorSize(couplRows,couplCols);
      if(couplCols < couplRows){
        EXCEPTION("Coupling tensor has to be in non-transposed form");
      }
    }
 
	  hystCoefFunction_ = CoefFncHyst;
	
    tensorName_ = tensorName;
    transposed_ = transposed;
	}
  
	virtual ~CoefFunctionHystMatTensor(){};

	void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
		EXCEPTION("CoefFunctionHystMatTensor only returns tensors")
	}

	void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
		EXCEPTION("CoefFunctionHystMatTensor only returns tensors")
	}

  virtual UInt GetVecSize() const {
		EXCEPTION("GetVecSize not used here")
  }
  
	void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
    
    UInt numCols, numRows;
    // GetTensorSize will return the size of the actual tensor that shall
    // be returned, i.e. if it is the transposed of the couplTensor_
    // numCols = numRows(coupTensor_)
    GetTensorSize(numRows,numCols);
    Matrix<Double> tmp;
    
	  outputTensor.Resize(numRows,numCols);
	  outputTensor.Init();
    
    if(transposed_){
      tmp = Matrix<Double>(numCols,numRows);
    } else {
      tmp = Matrix<Double>(numRows,numCols);
    }
    tmp.Init();

    int deltaForm_ = hystCoefFunction_->GetDeltaForm();
    bool deltaFormActive_ = hystCoefFunction_->deltaMatActive();
    bool rotate = true;
    int timelevel = 0; // get current value
    
    if(tensorName_ == "Permittivity"){
    
      if(deltaFormActive_ && (deltaForm_ != 0) ) {
//        std::cout << "Compute DeltaMatrix" << std::endl;
//        std::cout << "deltaForm_ = " << deltaForm_ << std::endl;
        if(deltaForm_ == -1){
//          std::cout << "Compute deltaMat wrt lastTS value" << std::endl;
        } else if (deltaForm_ != 0){
//          std::cout << "Compute deltaMat wrt to the " << deltaForm_ << " last iteration value" << std::endl;
        }
      } else if (deltaForm_ != 0){
//        std::cout << "deltaFormSet_ == true, but deltaMatNotActive!" << std::endl;
      } else {
//        std::cout << "deltaFormSet_ == false" << std::endl;
      }
  //    
      // TODO: implemement actual computation of deltaMatrix
      fieldTensor_->GetTensor(tmp,lpm);
      
    } else if (tensorName_ == "CouplingMechToElec"){
//      std::cout << "ComputeMechToElec" << std::endl;
     
      Matrix<Double> rotatedCouplTensor;
//      std::cout << "couplTensor:" << couplTensor_->ToString() << std::endl;
      if(rotate){
        hystCoefFunction_->ScaleAndRotateCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel);
      }
//      std::cout << "rotated couplTensor:" << rotatedCouplTensor.ToString() << std::endl; 
      tmp = rotatedCouplTensor;
      
    } else if (tensorName_ == "CouplingElecToMech"){ 
//      std::cout << "ComputeElecToMech" << std::endl;
      
      Matrix<Double> rotatedCouplTensor;
//      std::cout << "couplTensor:" << couplTensor_->ToString() << std::endl;
      if(rotate){
        hystCoefFunction_->ScaleAndRotateCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel);
      }
//      std::cout << "rotated couplTensor:" << rotatedCouplTensor.ToString() << std::endl; 
      tmp = rotatedCouplTensor;
      
      //couplTensor_->GetTensor(tmp,lpm);
      //std::cout << "couplTensor:" << couplTensor_->ToString() << std::endl;
    }
//    std::cout << "outTensor = " << outputTensor.ToString() << std::endl;
//    
//    Vector<Double> curPol;
//    std::cout << "Call GetVector of pol" << std::endl;
//    hystCoefFunction_->GetVector(curPol,lpm);
//    std::cout << "Retrieved value: " << std::endl;
//    std::cout << curPol.ToString() << std::endl;
//    Vector<Double> curPol2;
//    std::cout << "Get current Output of hystOperator" << std::endl;
//    curPol2 = hystCoefFunction_->GetOutputOfHysteresisOperator(lpm,0,false);
//    std::cout << "Retrieved value: " << std::endl;
//    std::cout << curPol2.ToString() << std::endl;
//    std::cout << "Diff between getVector and GetOutput" << std::endl;
//    Vector<Double> diff = curPol;
//    diff -= curPol2;
//    std::cout << diff.ToString() << std::endl;
//    
//    std::cout << "Get lastIt Output of hystOperator" << std::endl;
//    curPol2 = hystCoefFunction_->GetOutputOfHysteresisOperator(lpm,1,false);
//    std::cout << "Retrieved value: " << std::endl;
//    std::cout << curPol2.ToString() << std::endl;
//    
//    std::cout << "Get lastTS Output of hystOperator" << std::endl;
//    curPol2 = hystCoefFunction_->GetOutputOfHysteresisOperator(lpm,2,false);
//    std::cout << "Retrieved value: " << std::endl;
//    std::cout << curPol2.ToString() << std::endl;
//    
    
    if(transposed_){
//      std::cout << "Perform transpose" << std::endl;
      tmp.Transpose(outputTensor);
//      std::cout << "Transposed tensor: " << tmp.ToString() << std::endl;
    } else {
      outputTensor = tmp;
    }
    
	}

  //! Return row and columns size of tensor if coefficient function is a tensor
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    if( (tensorName_ == "Permittivity") || (tensorName_ == "Reluctivity")){
      fieldTensor_->GetTensorSize(numRows,numCols);
    } else if ( (tensorName_ == "CouplingMechToElec") || (tensorName_ == "CouplingElecToMech") 
                || (tensorName_ == "CouplingMechToMag") || (tensorName_ == "ComputeMagToMech") ){
      couplTensor_->GetTensorSize(numRows,numCols);
    } else {
      EXCEPTION("Tensor type unknown");
    }
    if(transposed_){
      UInt tmp = numRows;
      numRows = numCols;
      numCols = tmp;
    }
  }

private:

  PtrCoefFct fieldTensor_;
  PtrCoefFct elastTensor_;
  PtrCoefFct couplTensor_;
  
  PtrCoefFct hystCoefFunction_;
  
  std::string tensorName_;
  bool transposed_;
};

class CoefFunctionHystRHSLoad : public CoefFunction{
  
public:
  
  CoefFunctionHystRHSLoad(PtrCoefFct fieldTensor,PtrCoefFct elasticityTensor,
	                          PtrCoefFct couplingTensor, PtrCoefFct CoefFncHyst,
	                          std::string vectorName, bool transposed, bool onBoundary)
  :CoefFunction(){

    dimType_ = VECTOR;
	  isAnalytic_ = false;
	  isComplex_  = false;
    dependType_ = SOLUTION;
    
	  fieldTensor_ = fieldTensor;
	  elastTensor_ = elasticityTensor;
	  couplTensor_ = couplingTensor;

    UInt couplRows, couplCols;
    if(couplTensor_ != NULL){
      couplTensor_->GetTensorSize(couplRows,couplCols);
      if(couplCols < couplRows){
        EXCEPTION("Coupling tensor has to be in non-transposed form");
      }
    }
    
	  hystCoefFunction_ = CoefFncHyst;
	
    vectorName_ = vectorName;
    transposed_ = transposed;
    onBoundary_ = onBoundary;
	}
  
	virtual ~CoefFunctionHystRHSLoad(){};

	void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
		EXCEPTION("CoefFunctionHystRHSLoad only returns vectors")
	}

	void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
		// return vector that can be put onto rhs of system; 
    // the sign of term shall be such, that it can be added with + in the
    // corresponding pdes
    // Example: in piezo systems we have (among others) -int_Omega (BN)^T*P dOmega 
    //          on rhs, so we would return -P here
    
    int timeLevelRHS;
    if(onBoundary_){
      // boundary terms usually are evluated with the current value of polarization etc
      // whereas the rhs depends on usage of deltaMat or not
      timeLevelRHS = hystCoefFunction_->GetTimeLevel("BC");
    } else {
      timeLevelRHS = hystCoefFunction_->GetTimeLevel("RHS");
    }
    
    if(onBoundary_){
       if(timeLevelRHS == -1){
//        std::cout << "Get BC wrt to the last ts" << std::endl;
      } else if (timeLevelRHS != 0){
//        std::cout << "Get BC wrt to the last "<< timeLevelRHS << "th it" << std::endl;
      } else {
//        std::cout << "Put current value of hyst to BC" << std::endl;
      }
    } else {
      if(timeLevelRHS == -1){
//        std::cout << "Get RHS Load wrt to the last ts" << std::endl;
      } else if (timeLevelRHS != 0){
//        std::cout << "Get RHS Load wrt to the last "<< timeLevelRHS << "th it" << std::endl;
      } else {
//        std::cout << "Put current value of hyst to RHS" << std::endl;
      }
    }

    if(vectorName_ == "ElecPolarization"){
//      std::cout << "ElecPolarization was requested" << std::endl;
//      if(onBoundary_){
//        std::cout << "ON BOUNDARY" << std::endl;
//      }

      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
      //outputVector = hystCoefFunction_->GetOutputOfHysteresisOperator(lpm,timeLevelRHS,false);
    } else if (vectorName_ == "PiezoLoadForMechPDE"){
//      std::cout << "PiezoLoadForMechPDE was requested" << std::endl;
      //TODO: implement correct vector function; at the moment, return simply zero
//      std::cout << couplTensor_->ToString() << std::endl;
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
//      std::cout << "outputVector: " << outputVector.ToString() << std::endl;
    } else if (vectorName_ == "PiezoLoadForElecPDE"){
//      std::cout << "PiezoLoadForElecPDE was requested" << std::endl;
      //TODO: implement correct vector function; at the moment, return simply zero
//      std::cout << couplTensor_->ToString() << std::endl;
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
//      std::cout << "outputVector: " << outputVector.ToString() << std::endl;
    } else {
      EXCEPTION("VectorName not implemented yet");
    }
	}

  virtual UInt GetVecSize() const {

    UInt numRows,numCols;
    if( (vectorName_ == "ElecPolarization") || (vectorName_ == "MagPolarization")){
      fieldTensor_->GetTensorSize(numRows,numCols);
    } else if ( (vectorName_ == "PiezoLoadForMechPDE") || (vectorName_ == "PiezoLoadForElecPDE") 
                || (vectorName_ == "MagStrictLoadForMechPDE") || (vectorName_ == "MagStrictLoadForMagPDE") ){
      couplTensor_->GetTensorSize(numRows,numCols);
    } else {
      EXCEPTION("Vector type unknown");
    }
    if(transposed_){
      return numCols;
    } else {
      return numRows;
    }
  }
      
	void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
    EXCEPTION("CoefFunctionHystRHSLoad only returns vectors")
  }

  //! Return row and columns size of tensor if coefficient function is a tensor
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    EXCEPTION("GetVecSize not used here")
  }
  
private:

  PtrCoefFct fieldTensor_;
  PtrCoefFct elastTensor_;
  PtrCoefFct couplTensor_;
  
  PtrCoefFct hystCoefFunction_;
  
  std::string vectorName_;
  bool transposed_;
  bool onBoundary_;
};


class CoefFunctionHystOutput : public CoefFunction{
  
public:
				
	CoefFunctionHystOutput(PtrCoefFct CoefFncHyst,std::string ResultName){

    if( (ResultName == "DeltaPermeability")||(ResultName == "DeltaPermittivity") ){
      dimType_ = TENSOR;
    } else {
      dimType_ = VECTOR;
    }
    
	  isAnalytic_ = false;
	  isComplex_  = false;
    dependType_ = SOLUTION;
    
	  hystCoefFunction_ = CoefFncHyst;
    resultName_ = ResultName;
	}
  
	virtual ~CoefFunctionHystOutput(){};

	void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
		EXCEPTION("CoefFunctionHystOutput only returns vectors")
	}

	void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
    // output is always the actual quantity
    int timeLevelOutput = 0;
		if(resultName_ == "ElecPolarization"){    
      outputVector = hystCoefFunction_->GetOutputOfHysteresisOperator(lpm,timeLevelOutput,false);
    } else if (resultName_ == "MagPolarization"){
      outputVector = hystCoefFunction_->GetOutputOfHysteresisOperator(lpm,timeLevelOutput,true);
    } else {
      EXCEPTION("Unknown result")
    }
	}
  
  void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
    if( (resultName_ == "DeltaPermeability")||(resultName_ == "DeltaPermittivity") ){
      //TODO: estimate permeability or permittivity around current working point
      UInt numCols,numRows;
      hystCoefFunction_->GetTensorSize(numRows,numCols);
      outputTensor.Resize(numRows,numCols);
      outputTensor.Init();
    } else {
      EXCEPTION("Unknown result")
    }
  }

  virtual UInt GetVecSize() const {
    return hystCoefFunction_->GetVecSize();
  }
  
  //! Return row and columns size of tensor if coefficient function is a tensor
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    hystCoefFunction_->GetTensorSize(numRows,numCols);
  }
  
private:
  
  PtrCoefFct hystCoefFunction_;
  std::string resultName_;
};


// ============================================================================
//  Hysteresis
// ============================================================================
//! Provide a coefficient for hysteresis modeling
//! \note This class only works for real-valued scalar data.
class CoefFunctionHyst : public CoefFunction{
public:

  //! Constructor
  CoefFunctionHyst( BaseMaterial* const material,
					shared_ptr<ElemList> actSDList,
					PtrCoefFct dependency1,
					SubTensorType tensorType,
					MaterialType matType,
					shared_ptr<FeSpace> ptFeSpace, bool performInversionTest=false);

  //! Constructor
  CoefFunctionHyst( BaseMaterial* const material,
                    shared_ptr<ElemList> actSDList,
                    PtrCoefFct dependency1,
					PtrCoefFct dependency2,
				    SubTensorType tensorType,
				    MaterialType matType,
				    shared_ptr<FeSpace> ptFeSpace, bool performInversionTest=false);

  void Init(BaseMaterial* const material,
                    shared_ptr<ElemList> actSDList,
                    PtrCoefFct dependency1,
				    SubTensorType tensorType,
				    MaterialType matType,
				    shared_ptr<FeSpace> ptFeSpace, bool performInversionTest=false);

  //! Destructor
  virtual ~CoefFunctionHyst();

//  //! Initialize with data
//  void Init( BaseMaterial* const material, shared_ptr<ElemList> actSDList);

  //! Return size of vector in case coefficient function is a vector
  virtual UInt GetVecSize() const {
    return dependCoef1_->GetVecSize();
  }

  //! Return row and columns size of tensor if coefficient function is a tensor
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      numRows = MAT_initialTensor_.GetNumRows();
      numCols = MAT_initialTensor_.GetNumCols();
  }

  void SetInputDependentFlags(UInt multiDigitInteger);

  void SetRuntimeDependentFlag(std::string flagName, UInt intState);

  void SetPreviousHystVals(bool setLastTS, bool forceMemoryLock = false);

  void GetScalar(Double& outputScalar, const LocPointMapped& lpm );

  void GetVector( Vector<Double>& outputVector,const LocPointMapped& lpm);

  void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm );

  std::string ToString() const;

  void EstimateCurrentSlope(Vector<Double> steppingDirection, Double scaling);

  Vector<Double> GetOutputOfHysteresisOperator(const LocPointMapped& lpm, int timeLevel, bool invert);

  Vector<Double> RetrieveInputToHysteresisOperator(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool invert);

  Vector<Double> RetrieveLPMSolution(LocPointMapped& actualLPM, UInt storageIdx, int timeLevel);

  Double GetOutputSaturation(){
    return MAT_ySat_;
  }
  
  PtrCoefFct GenerateMatCoefFnc(std::string tensorName,PtrCoefFct elastTensor = NULL, PtrCoefFct couplTensor = NULL ){
		PtrCoefFct ret;
    if((tensorName == "Permittivity")||(tensorName == "CouplingMechToElec")||(tensorName == "CouplingElecToMech")){
      // Note: for Piezo we need two different coupling functions
      // Reason: in case of deltaMat formulation, the coupling matrices are not 
      //         transposed of each other as the mech to elec term has deltaS/deltaE
      //         in addition to the rotated coupling operator
      PtrCoefFct eps = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType_,Global::REAL);
      
      bool transposed = false;
      if(tensorName == "CouplingElecToMech"){
        // in mech PDE we need e^T or d^T
        transposed = true;
      }
      shared_ptr<CoefFunctionHystMatTensor> c(new CoefFunctionHystMatTensor(eps, 
              elastTensor, couplTensor, hystItself_,tensorName,transposed));
      
      ret = c;
    } else if((tensorName == "Reluctivity")||(tensorName == "CouplingMechToMag")||(tensorName == "CouplingMagToMech")){
      PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY,tensorType_,Global::REAL);
      
      bool transposed = false;
      if(tensorName == "CouplingMagToMech"){
        // in mech PDE we need h^T or g^T
        transposed = true;
      }
      shared_ptr<CoefFunctionHystMatTensor> c(new CoefFunctionHystMatTensor(nu, 
              elastTensor, couplTensor, hystItself_,tensorName,transposed));
      
      ret = c;
    } else {
      EXCEPTION("tensorName not implemented yet");
    }

    return ret;
  }
  
  PtrCoefFct GenerateRHSCoefFnc(std::string vectorName,PtrCoefFct elastTensor = NULL, PtrCoefFct couplTensor = NULL, bool onBoundary = false ){
		PtrCoefFct ret;
    if( (vectorName == "ElecPolarization") || (vectorName == "PiezoLoadForMechPDE") || (vectorName == "PiezoLoadForElecPDE")){
      PtrCoefFct eps = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType_,Global::REAL);
      
      bool transposed = false;
      if(vectorName == "PiezoLoadForMechPDE"){
        // in mech PDE we need h^T or g^T
        transposed = true;
      }
      shared_ptr<CoefFunctionHystRHSLoad> c(new CoefFunctionHystRHSLoad(eps, elastTensor, couplTensor, hystItself_,vectorName,transposed,onBoundary));
      
      ret = c;
    } else if( (vectorName == "MagPolarization") || (vectorName == "MagStrictLoadForMechPDE") || (vectorName == "PMagStrictLoadForMagPDE")){
      PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY,tensorType_,Global::REAL);
      
      bool transposed = false;
      if(vectorName == "PMagStrictLoadForMechPDE"){
        // in mech PDE we need h^T or g^T
        transposed = true;
      }
      shared_ptr<CoefFunctionHystRHSLoad> c(new CoefFunctionHystRHSLoad(nu, elastTensor, couplTensor, hystItself_,vectorName,transposed,onBoundary));
      
      ret = c;
    } else {
      EXCEPTION("vectorName not implemented yet");
    }
    
    return ret;
  }
  
  PtrCoefFct GenerateOutputCoefFnc(std::string ResultName){
    
    PtrCoefFct ret;

    shared_ptr<CoefFunctionHystOutput> c(new CoefFunctionHystOutput(hystItself_,ResultName));

    ret = c;
    
    return ret;
  }
  
  
  bool deltaMatActive(){
    return deltaMatActive_;
  }
  
  int GetDeltaForm(){
    return RUN_deltaForm_;
  }
    
  int GetTimeLevel(std::string Type){
    if(Type == "Mat"){
      return timeLevel_Mat_;
    } else if (Type == "RHS"){
      return timeLevel_RHS_;
    } else if (Type == "BC"){
      return timeLevel_BC_;
    } else if (Type == "Output"){
      return timeLevel_Output_;
    } else {
      EXCEPTION("Tpye not known");
      return -2;
    }
  }

  void TestInversion(Double eps_mu);
    
private:

  void ScaleAndRotateCouplingTensor(const LocPointMapped& lpm, PtrCoefFct couplTensorCoefFnc, Matrix<Double>& rotatedCouplTensor,int timeLevel){
    // compare to Kaltenbacher "Numerical Simulation ..." 3rd Edition p. 387
	  //
	  // scaling:
	  // [e(P)] = P^i_k/P_sat[e]
	  //
	  // for vector model, an additional rotation towards P could be reasonable
    
    // get actual indices for storage array to find out if we need to reevaluated or not
    UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;

		PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx);
    
    if(rotatedCouplingTensor_requiresReeval_[storageIdx] == false){
//      std::cout << "Coupling Tensor already rotated > no reeval performed" << std::endl;
      rotatedCouplTensor = rotatedCouplingTensor_[storageIdx];
    } else {
      
      // obtain coupling tensor first
      UInt numRows, numCols;
      couplTensorCoefFnc->GetTensorSize(numRows,numCols);
      Matrix<Double> couplTensor = Matrix<Double>(numRows,numCols);
      couplTensorCoefFnc->GetTensor(couplTensor,actualLPM);
  
//      std::cout << "Retrieved coupling tensor: " << couplTensor.ToString() << std::endl;
//      std::cout << "Rotate coupling tensor" << std::endl;
      rotatedCouplTensor = Matrix<Double>(numRows,numCols);
      Matrix<Double> scaledCouplTensor = Matrix<Double>(numRows,numCols);
      
      bool invert = false;
    
      // get current polarization (elec or mag)
      Vector<Double> P = GetOutputOfHysteresisOperator(lpm, timeLevel, invert);
      
//      std::cout << "Current polarization vector " << P.ToString() << std::endl;
//      
      // calculate scaling
      assert(MAT_ySat_ != 0);
      Double scaling = P.NormL2()/MAT_ySat_;

      scaledCouplTensor = couplTensor* scaling;

//      std::cout << "Scalaed coupling tensor " << scaledCouplTensor.ToString() << std::endl;
//      
      if(P.NormL2() == 0){
//        std::cout << "Polarization is zero > perform no rotation" << std::endl;
        rotatedCouplTensor = scaledCouplTensor;
        
      } else {
        Vector<Double> dirP = Vector<Double>(P.GetSize());
        dirP = P / P.NormL2();
        
        // calculate rotation matrix
        if(numCols == 4){
          // axi case: 2x4 matrix
          WARN("Rotation for axi case not implemented yet. No rotation was peformed");
          rotatedCouplTensor = scaledCouplTensor;
          
        } else if (numCols == 6){
          // 3d plane case
          // Idea: Create rotation matrix, that maps the current direction of P onto z-axis
          //       as the z-axis is the default polarization axis in the mat file
          //       To obtained the desired behavior, we have to rotate the z-axis onto P however.
          //       This can be done by taking the transposed rotation matrix.

          Double alpha, beta, gamma;

          // 1. rotate around z-axis by angle gamma such that P lies in z-y plane
          // gamma = angle between z-y plane and dirP
          gamma = std::atan2(dirP[0],dirP[1]);

          // 2. rotate around x-axis by angle alpha such that P lies on top of z-axis
          // alpha = angle between x-y plane and z
          alpha = std::atan2(std::sqrt(dirP[1]*dirP[1]+dirP[0]*dirP[0]),dirP[2]);

          // no rotation arouind y-axis needed > beta = 0
          // WARNING: this whole procedure is only valid if coupling tensor is at least
          // transverse isotropic! Otherwise we have to figure out on which axis to rotate
          // the transverse directions
          beta = 0.0;
          Matrix<Double> R = Compute3DRotationMatrix(alpha, beta, gamma);

          // take transpose matrix to revert rotation (i.e. rotate z-axis onto dirP)
          Matrix<Double> RT;
          RT.Resize(3,3);
          R.Transpose(RT);

          assert(rotatedCouplTensor.GetNumRows() == 3);
          assert(rotatedCouplTensor.GetNumCols() == 6);

          scaledCouplTensor.PerformRotation(RT,rotatedCouplTensor);

        } else {
          // 2d plane strain or plane stress cas
          // Important remark:
          //  for 2d plane strain and stress (and somehow also for axi) the coupling tensor
          //  will be rotated by default by alpha = -90 and gamma = -90.
          //  This results in the following mapping of the coordinate axis:
          //    z > y, y > x, x > z
          //  I.e. the material is rotated such that the default polarization is in +y direction.
          //  We thus have two possibilities to align the material tensor to the 2d polarization
          //  direction:
          //  a) get original 3x6 tensor and rotate that tensor according to the 3d version above
          //      then cut out subtensor
          //  b) take cut out subtensor and perform rotation directly in 2d
          //      > Version b done in the following
          Double gamma;

          // std::atan2(dirP[0],dirP[1]) would rotate towards the y-axis
          // we want however the y-axis to rotate, so that we take -gamma instead
          gamma = -std::atan2(dirP[0],dirP[1]);

          Matrix<Double> R = Compute2DRotationMatrix(gamma);

          assert(rotatedCouplTensor.GetNumRows() == 3);
          assert(rotatedCouplTensor.GetNumCols() == 6);

          scaledCouplTensor.PerformRotation(R,rotatedCouplTensor);

        }
      }
      
      rotatedCouplingTensor_requiresReeval_[storageIdx] = false;
      rotatedCouplingTensor_[storageIdx] = rotatedCouplTensor;
      
    }
  }
  
  // helper functions for computation of rhs loads and lhs tensors
//  
//	Matrix<Double> ScaleAndRotateCouplingTensor(const LocPointMapped& lpm, bool invert, bool scaleOnly){
//	  // compare to Kaltenbacher "Numerical Simulation ..." 3rd Edition p. 387
//	  //
//	  // scaling:
//	  // [e(P)] = P^i_k/P_sat[e]
//	  //
//	  // for vector model, an additional rotation towards P could be reasonable
//
//	  // get CURRENT state of P (electric or magnetic polarization)
//	  UInt timeLevel = 0;
//
//	  Vector<Double> P = hystPol_->GetOutputOfHysteresisOperator(lpm, timeLevel, invert);
//	  Double Psat = hystPol_->GetOutputSaturation();
//
//	  // calculate scaling
//	  Double scaling = P.NormL2()/Psat;
//
//    Matrix<Double> couplTensorMat; 
//    couplTensor_->GetTensor(couplTensorMat,lpm);
//	  Matrix<Double> retMat = Matrix<Double>(couplTensorMat.GetNumRows(),couplTensorMat.GetNumCols());
//
//	  retMat = couplTensorMat* scaling;
//
//	  if(scaleOnly || (P.NormL2() == 0)){
//	    return retMat;
//	  }
//
//	  Vector<Double> dirP = Vector<Double>(P.GetSize());
//	  dirP = P / P.NormL2();
//
//	  // calculate rotation matrix
//	  if(subType_ == AXI){
//	    WARN("Rotation for axi case not implemented yet. No rotation was peformed");
//	  } else if (subType_== FULL){
//	    // Idea: Create rotation matrix, that maps the current direction of P onto z-axis
//	    //       as the z-axis is the default polarization axis in the mat file
//	    //       To obtained the desired behavior, we have to rotate the z-axis onto P however.
//	    //       This can be done by taking the transposed rotation matrix.
//
//	    Double alpha, beta, gamma;
//
//	    // 1. rotate around z-axis by angle gamma such that P lies in z-y plane
//	    // gamma = angle between z-y plane and dirP
//	    gamma = std::atan2(dirP[0],dirP[1]);
//
//	    // 2. rotate around x-axis by angle alpha such that P lies on top of z-axis
//	    // alpha = angle between x-y plane and z
//	    alpha = std::atan2(std::sqrt(dirP[1]*dirP[1]+dirP[0]*dirP[0]),dirP[2]);
//
//	    // no rotation arouind y-axis needed > beta = 0
//	    // WARNING: this whole procedure is only valid if coupling tensor is at least
//	    // transverse isotropic! Otherwise we have to figure out on which axis to rotate
//	    // the transverse directions
//	    beta = 0.0;
//	    Matrix<Double> R = Compute3DRotationMatrix(alpha, beta, gamma);
//
//	    // take transpose matrix to revert rotation (i.e. rotate z-axis onto dirP)
//	    Matrix<Double> RT;
//      RT.Resize(3,3);
//      R.Transpose(RT);
//	    Matrix<Double> rotatedState;
//
//      assert(retMat.GetNumRows() == 3);
//      assert(retMat.GetNumCols() == 6);
//
//	    rotatedState.Resize(3,6);
//
//	    retMat.PerformRotation(RT,rotatedState);
//
//	  } else {
//	    // 2d plane strain or plane stress cas
//	    // Important remark:
//	    //  for 2d plane strain and stress (and somehow also for axi) the coupling tensor
//	    //  will be rotated by default by alpha = -90 and gamma = -90.
//	    //  This results in the following mapping of the coordinate axis:
//	    //    z > y, y > x, x > z
//	    //  I.e. the material is rotated such that the default polarization is in +y direction.
//	    //  We thus have two possibilities to align the material tensor to the 2d polarization
//	    //  direction:
//	    //  a) get original 3x6 tensor and rotate that tensor according to the 3d version above
//	    //      then cut out subtensor
//	    //  b) take cut out subtensor and perform rotation directly in 2d
//	    //      > Version b done in the following
//	    Double gamma;
//
//	    // std::atan2(dirP[0],dirP[1]) would rotate towards the y-axis
//	    // we want however the y-axis to rotate, so that we take -gamma instead
//	    gamma = -std::atan2(dirP[0],dirP[1]);
//
//	    Matrix<Double> R = Compute2DRotationMatrix(gamma);
//
//      Matrix<Double> rotatedState;
//      rotatedState.Resize(2,3);
//
//      assert(retMat.GetNumRows() == 2);
//      assert(retMat.GetNumCols() == 3);
//
//      retMat.PerformRotation(R,rotatedState);
//
//	  }
//
//	  return retMat;
//	}
  
	Matrix<Double> Compute3DRotationMatrix(Double alpha, Double beta, Double gamma){
	  // Compute rotation matrix based on Kardan Angles
	  // Implementation taken from BaseMaterial class

    // Calculate rotation matrix( based on Kardan-Angles)
    // Ref.: C. Woernle, "Skript: Dynamik von Mehrkoerpersystemen,
    // Kapitel 2 "Grundlagen der Kinematik", S. 12, Univ. Rostock
    // http://iamserver.fms.uni-rostock.de/studium/mehrkoerpersysteme/unterlagen.htm

    Matrix<Double> R(3,3);
    R.Resize(3,3);
    R[0][0] =  std::cos(beta) * std::cos(gamma);
    R[0][1] = -std::cos(beta) * std::sin(gamma);
    R[0][2] =  std::sin(beta);
    R[1][0] =  std::cos(alpha)*std::sin(gamma) + std::sin(alpha)*std::sin(beta)*std::cos(gamma);
    R[1][1] =  std::cos(alpha)*std::cos(gamma) - std::sin(alpha)*std::sin(beta)*std::sin(gamma);
    R[1][2] = -std::sin(alpha)*std::cos(beta);
    R[2][0] =  std::sin(alpha)*std::sin(gamma) - std::cos(alpha)*std::sin(beta)*std::cos(gamma);
    R[2][1] =  std::sin(alpha)*std::cos(gamma) + std::cos(alpha)*std::sin(beta)*std::sin(gamma);
    R[2][2] =  std::cos(alpha)*std::cos(beta);

    return R;
	}

	Matrix<Double> Compute2DRotationMatrix(Double gamma){
    // Compute rotation matrix based on Kardan Angles
    // Here we basically rotate in the x-y-plane

    Matrix<Double> R(2,2);
    R.Resize(2,2);
    R[0][0] =  std::cos(gamma);
    R[0][1] =  -std::sin(gamma);
    R[1][0] =  std::sin(gamma);
    R[1][1] =  std::cos(gamma);

    return R;
  }
  

  void PreprocessLPM(const LocPointMapped& lpmInput, LocPointMapped& lpmOutput,
                     UInt& operatorIdx, UInt& storageIdx);

  Vector<Double> CalcOutputOfHysteresisOperator(Vector<Double> inpute, UInt operatorIdx, UInt storageIdx);

  void InitStorage();

  void CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& outputTensor, std::string evalMethod,
                         UInt storageIdx, bool intoSat, bool outofSat, bool satToSat, Vector<Double>& X_current);

  void ExtractSolutionAndInputForHystOperator(Vector<Double>& extractedSolution, Vector<Double>& extractedInput, UInt& operatorIndex,
                                                 UInt& storageIndex, const LocPointMapped& lpm, bool midpointOnly);

  UInt EvaluateHysteresisOperator( Vector<Double> inputToHystOperator, Vector<Double>& outputOfHystOperator,
                                   UInt operatorIndex, UInt storageIndex, bool overwriteMemory, Vector<Double>& currentSolution );

  bool EvaluateAtMidpointOnly();
  bool OverwriteHystMemory();

  // for usage with coefFunctionHystMat
  bool deltaMatActive_;
  int RUN_deltaForm_;
  UInt timeLevel_Mat_;
  UInt timeLevel_RHS_;
  UInt timeLevel_BC_;
  UInt timeLevel_Output_;
  PtrCoefFct hystItself_;
  
  Matrix<Double>* rotatedCouplingTensor_;
  bool* rotatedCouplingTensor_requiresReeval_;
  
  // for inversion with Levenberg Marquart and linesearch
  // for each material we store the last used linesearch stepping alpha
  // as starting value for next usage
  // starting value = -1 > first estimation done during first call to
  // computeInput_vec
  Double alphaLinesearch_;
  bool performInversionTest_;
  
  
  /*
   * ###############################################
   * ### Parameters extracted from material file ###
   * ###############################################
   */
  /*
   * Initial/small signal material tensor (permittivity / reluctivity)
   */
  Matrix<Double> MAT_initialTensor_;
  Matrix<Double> MAT_freeFieldTensor_;

  /*
   * Preisach weights and its size in form of the number of rows
   */
  Matrix<Double> MAT_PreisachWeights_;
  UInt MAT_numRows_;

  /*  elemMat.Resize( nrFncs * bOperator_->GetDimDof());
  elemMat.Init();


  // Loop over all integration points
  LocPointMapped lp;
   * Saturation values for input (x, E/H) and output (y, P/M) of
   * hysteresis operator
   */
  Double MAT_xSat_;
  Double MAT_ySat_;

  /*
   * determines whether the vector or the scalar model shall be used
   */
  CoefDimType MAT_methodType_;

  /*
   * Additional parameter for scalar Preisach model
   */
  /*
   * Direction (x,y,z) in which Polarization of material points
   */
  Directions MAT_dirP_;

  /*
   * Additional parameter for vector Preisach model
   */
  // determines wheter the vector model from 2012 (classical) or the one from#
  // 2015 (revised) is used
  bool isClassical_;
  
  /*
   * see VectorPreisach10 for more details
   */
  Double MAT_rotRes_;
  Double MAT_angResistance_;

  /*
   * if != 0, the rotation states of the vector preisach model will be clipped
   * to the provided precision (in rad)
   */
  Double MAT_angClipping_;

  /*
   *  new parameter added 03.07.2017
   *
   *  these parameter mark the minimal change (in rad and in L2-Norm)
   *  of the input to the hysteresis operator which shall be resolved;
   *  i.e. the current input will be compared against
   *  a) the old TS/iteration state
   *  b) against the previous input
   *
   *  both times, it will be checked if
   *    (currentInput - oldInput).L2Norm() < MAT_ampResolution_
   *    std::acos(currentInput.Inner(oldInput)) < MAT_angResolution_
   *
   *  if both criteria are true, the old value is taken and no reevaluation is done
   */
  Double MAT_angResolution_;
  Double MAT_ampResolution_;

  /*
   * 1,2 nested-list implementation of classical and revised vector hyst. model
   * 10,20 matrix based implementation of classical and revised vector hyst. model
   */
  UInt MAT_vecPreisachImplementationVersion_;

  /*
   * Initial input to the vector hysteresis operator;
   * and output of the vector hysteresis operator to this input
   * Note: the value specified in the mat file is the INPUT to the hyst operator
   *        not the actual polarization/magnetization state
   */
  Vector<Double> MAT_initialInput_;
  Vector<Double>* MAT_initialOutput_;

  /*
   * enables output of overlapped switching and rotation state in form of bmp
   * images
   */
  UInt MAT_bmpResolution_;
  UInt MAT_printOut_;

  /*
   * ################################################
   * ### Parameters extracted from xml input file ###
   * ################################################
   */
  /*
   * bool indicating if parameter where already set
   * Note: we cannot (maybe we can if we get a handle to the xml input during
   * the creation of the hyst operator) initialize these parameter in the constructor
   * as we need to pass a certain input from the xml flag; this should be done
   * in the StdSolveStep during the very first iteration;
   * to check if this step was done use the XMLParameterSet_ flag
   * (will be set by function SetInputDependentFlags)
   */
  bool XMLParameterSet_;
  /*
   * XML_deltaEvalVersion_ = xy
   *
   *  x > determines what values to use for deltaX and deltaY
   *  y > determines what to with deltaX and deltaY
   *
   *  XML_deltaEvalVersion_ = x0 (10,20,30,...)
   *
   *    deltaMat = addTensor +/- abs(deltaY)/abs(deltaX)
   *
   *    ->  see "DirectDivisionAbsNew" in CreateDeltaMatrix
   *
   *  XML_deltaEvalVersion_ = x1 (11,21,31,...)
   *
   *    deltaMat =
   *          addTensor +/- abs(deltaY)/abs(deltaX) (if deltaY != 0)
   *
   *          deltaMat_lastIt_[idx]/10;             (if deltaY == 0 && deltaMat_lastIt_[idx] > 10*rev_mat_fac_)
   *
   *    ->  see "DirectDivisionAbsSatStepDownNew" in CreateDeltaMatrix
   *
   *  XML_deltaEvalVersion_ = 1y (10,11)
   *
   *    deltaX = currentInput - E_B_lastIt_ or
   *           = currentInput - E_B_lastTS_ (depending on flag lastTS_)
   *
   *    deltaY = currentOutput - P_M_lastIt_ or
   *           = currentOutput - P_M_lastTS_ (depending on flag lastTS_)
   *
   *  XML_deltaEvalVersion_ = 2y (20,21)
   *
   *    a) to be done before evaluating the hysteresis operator
   *    deltaX = currentInput - E_B_lastIt_ or
   *           = currentInput - E_B_lastTS_ (depending on flag lastTS_)
   *
   *    if deltaX[i] < minimalShift  for all i
   *      -> reuse old deltaMatrix
   *    if deltaX[i] < minimalShift  for some i but not for all i
   *      -> shift input a little
   *      -> currentInput[i] += minimalShift
   *
   *    b) evaluate hystersis operator with shifted (or unshifted) input
   *
   *    deltaY = currentOutput - P_M_lastIt_ or
   *           = currentOutput - P_M_lastTS_ (depending on flag lastTS_)
   *
   *    Advantage: no issues with division by 0 during CreateDeltaMat
   *               numerical noise will be overwritten by minimalShift
   *    Disadvantage: hystoperator is not evaluated at the actual input
   *                  solution might drift over time
   *
   *  XML_deltaEvalVersion_ = 3y (30,31)
   *
   *    a) evaluate hystersis operator at currentInput and at a shiftedInput
   *      with
   *        shiftedInput = currentInput + minimalShift * (1,1,1)
   *    b)
   *      deltaX = minimalShift * (1,1,1)
   *      deltaY = shiftedOutput - currentOutput
   *
   *    Advantage: no issues wiht division by 0
   *               differential deltaMatrix computed
   *    Disadvantage: hystOperator needs to be evaluated twice
   *
   */
  UInt XML_deltaEvalVersion_;

  /*
   * kHApproxVersion_
   *
   *  used for approximating H from B in magnetic case
   *  kHApproxVersion_ = 1     -> H_new = nu0*B_new - M_old
   *  kHApproxVersion_ = 2     -> H_new = deltaMat_lastIt_ * (B_new - B_old) + H_old
   */
  UInt XML_HApproxVersion_;

  /*
   * kEvaluationDepth_ (to be set ONCE at the start of the simulation; set via
   *                    input file possible)
   *
   *  kEvaluationDepth_= 1 (standard)
   *   - each element has exactly one hysteresis operator attached
   *   - numHystOperators = numElems
   *   - each Hysteresis operator will only be evaluated at the midpoint of
   *    the corresponding element
   *   - all data structures store the values (elementSolution, current output
   *    of the Hysteresis operator, ...) at the midpoint of the corresponding
   *    element
   *
   *  kEvaluationDepth_ = 2 (extended)
   *   - each element has exactly one hysteresis operator attached
   *   - numHystOperators = numElems
   *   - for the functions GetScalar, GetVector and GetTensor (in BDB, BDU, ...
   *    integrators), the Hysteresis operator will be evaluated at the
   *    coordinates of the corresponding integration point
   *   - during the call to SetPreviousHystValues as well as during the call to
   *    GetScalar and GetVector for OUTPUT computation, the Hysteresis operator
   *    will be evaluated at the midpoint of the corresponding element
   *   - all data structures store the values at the midpoint of the element;
   *   - the Hysteresis operator will only store the state at the midpoint, i.e.
   *    all evaluation at other integrations points will be executed on temporal
   *    storage
   *
   *    -> basically, the Hysteresis operator stores the state of the material
   *      at the midpoint of the element in hope that this state does not differ
   *      too much from the actual hysteretic state at the integration points
   *
   *  kEvaluationDepth_ = 3 (full)
   *  (not yet implemented, but possible approach for future)
   *   - each INTEGRATION point + the midpoint have one hysteresis operator
   *    attached
   *   - numHystOperators = numElems * (numIntegrationsPoints + 1)
   *   - Evaluation and storage is done for each integration point
   */
  UInt XML_EvaluationDepth_;

  /*
   * activates differnet levels of runtime measurement
   * 0: no measurement
   * 1: measuerment of evaluation time of hyst operator
   * 2: detailed measurement of evaluation time of hyst operator
   */
  UInt XML_performanceMeasurement_;

  /*
   * level of text output; 0 = no output; 1 = std info; 2 = debug
   */
  UInt XML_textOutputLevel_;

  /*
   * #########################################################################
   * ### Parameters the will/might change during execution of StdSolveStep ###
   * ### these parameters are to be set from outside via SetIntegerFlag    ###
   * #########################################################################
   */

  /*
   * may be set to false to lock the rotation state of the hysteresis operator
   * changes in put will thus only affect the amplitude
   */
  bool RUN_overwriteDirection_;

  /*
   * determines whether the value from the last time step (true) or the one from
   * the last iteration (false) shall be used to compute the delta matrix
   */
  bool RUN_useLastTS_;

  /*
   * RUN_vectorToReturn_
   *
   *  used to tweak the behavior of the rhs load
   *
   *  RUN_vectorToReturn_ =
   *    1 -> evaluate Hysteresis operator and returns its value in rhs load
   *          integrator
   *    2 -> return zero vector (needed for implementations where we do not
   *          put P/M to the rhs or only for the first two iteratons)
   *    3 -> evaluate Hysteresis operator but return not its value, but the
   *          difference to the value from the last TS
   *    4 -> return only the value of the lastTS (needed for full stepping)
   */
  UInt RUN_vectorToReturn_;

  /*
   * tensorToReturn_
   *
   *  used to tweak the behavior of the coefficient function from outside
   *  needed in StdSolveStep as we need to assemble the deltaMatrix on the lhs
   *  and the freeField tensor on the rhs
   *
   *  tensorToReturn_ =
   *    1 -> return InitialTensor
   *    2 -> return FreeFieldTensor
   *    3 -> return computed deltaMatrix
   *
   */
  UInt RUN_tensorToReturn_;

  /*
   * tensorToAdd_
   *
   *  used to tweak the computation of the deltaMatrix
   *
   *  deltaMatrix =
   *    electrostatics:     tensorToAdd + deltaP/deltaE
   *    magnetics:          tensorToAdd - deltaM/deltaB
   *
   *  tensorToAdd_ =
   *    1 -> add InitialTensor
   *    2 -> add FreeFieldTensor
   *
   */
  UInt RUN_tensorToAdd_;

  /*
   * evaluationPurpose_ (to be changed during the iterative solve step)
   *
   *  flag needed in combination with kEvaluationDepth_
   *  it will determine, whether the hysteresis operator is to be evaluated
   *  at the midpoint of an element or at an integration point
   *
   *  kEvaluationDepth_   >    1 (standard)    2 (extended)          3 (full)
   *  evaluationPurpose_ v
   *    1 (assemble)          midpoint        integration point     integration point
   *    2 (store lastIT)      midpoint        midpoint              integration point + midpoint
   *    3 (store lastTS)      midpoint        midpoint              integration point + midpoint
   *    4 (output)            midpoint        midpoint              midpoint
   *
   *  (cases 2 and 3 differ in terms of memory locking; for case 2, we evaluate the
   *  hysteresis operator only on temporal storage, for case 3 we set the actual
   *  state)
   *
   */
  UInt RUN_evaluationPurpose_;

  /*
   * switch needed to activate BMP output only at the end of each time step
   * (otherwise there would be way to much bmp output if done during each iteration)
   */
  bool RUN_allowBMP_;

  /*
   * single global flag for indicating if residual could be reduced during
   * last iteration; this may be taken as an indicatior that we should adapt
   * the computatoin of the deltaMat
   * (obvious issue: the residual is a global quantity whereas the deltaMat is
   *  a local one, i.e. just because the residual did not converge, does not mean
   *  that the computation of all deltaMatrices should be changed (maybe some
   *  are already appropriate))
   */
  bool RUN_residualDecreased_;

  /*
   * ######################
   * ### Misc parameter ###
   * ######################
   */
  /*
   * For performance measurement
   */
  Timer* timer_;
  UInt totalCallingCounter_;
  UInt totalEvaluationCounter_;
  Double avgEvaluationTime_;
  Double totalEvaluationTime_;

  /*
   * Parameter passed from FE context
   */
  std::string PDEName_;

  //! global element to local element numbering
  std::map<UInt, UInt> globalElem2Local_;

  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef1_;
  PtrCoefFct dependCoef2_;

  //! type of tensor to be returned
  //! electrostatics: permittivity
  //! magnetics: reluctivity
  MaterialType matType_;

  //! type of subtensor (axi, plane)
  SubTensorType tensorType_;

  //! list of all elements
  shared_ptr<ElemList> SDList_;
  UInt numElemSD_;

  //! for standard and extended evaluation,
  //! we need one hysteresisOperator per element
  //! for full evaluation, we have N+1 hyst-Operators per element,
  //! where N is the number of integration points
  UInt numHystOperators_;

  //! for standard evaluation, we have to store values like the old
  //! input/output, the old and current deltaMatrix, etc for each element
  //! for extended and full evaluation we have to store it for each integration
  //! point + for the middlepoint
  UInt numStorageEntries_;

  //! dim of vector that is returned in GetVector
  //! basically 2 (for 2d) or 3 (for 3d)
  UInt dim_;

  /*
   * constants for permittivity and reluctivity of vacuum (eps0 and nu0)
   */
  Double rev_mat_fac_;

  Double tol_;

  //! actual hysteresis operator
  Hysteresis* hyst_;
  Hysteresis* hystTMP;

  //! material object
  BaseMaterial* material_;

  /*
   * ###############
   * ### Storage ###
   * ###############
   */
  /*
   * Storage vectors for quantities of the current iteration;
   * First letter refers to the quantity stored in the electrostatic case,
   * second letter to the quantity stored in the magnetic case, i.e.
   *
   *          electrostatics            magnetics
   * E_B_       E - ElecField             B - MagFlux         -> solution obtained directly
   *                                                              from BOperation applied to the
   *                                                              corresponding potential
   * P_M_       P - ElecPolarization      M - Magnetization   -> output of the Hysteresis operator
   *
   * E_H_       E - ElecField             H - MagField        -> input to the Hysteresis operator
   *                                        (needed as
   *                                        VectorPreisach has
   *                                        no inverse from yet)
   *
   *
   */
  Vector<Double>* E_B_;
  Vector<Double>* P_M_;
  Vector<Double>* E_H_;

  /*
   * Storage vectors for quantities of the last iteration;
   * these values will ONLY be overwritten during each iteration by calling the
   * function
   *    SetPreviousHystValues(false)
   *
   */
  Vector<Double>* E_B_lastIt_;
  Vector<Double>* P_M_lastIt_;
  Vector<Double>* E_H_lastIt_;

  /*
   * Storage vectors for quantites of the last timestep;
   * these values will ONLY be set during the FIRST iteration of each time step
   * by calling
   *    SetPreviousHystValues(true)
   */
  Vector<Double>* E_B_lastTS_;
  Vector<Double>* P_M_lastTS_;
  Vector<Double>* E_H_lastTS_;
  bool* requiresReeval_;

  /*
   * Storage matrices for current deltaMatrix and deltaMatrix of the last
   * iteration and of the last time step;
   * later one is currently only used to estimate the current magnetic field
   * H_new from the element solution B_new (we cannot directly retrieve H_new
   * from the FE system as we would have to know the correct magnetization to
   * compute H = nu0 * B - M)
   * deltaMat_lastIt_ will be set via SetPreviousHystValues
   * deltaMat_lastTS_ will be set via SetPreviousHystValues(true)
   */
  Matrix<Double>* deltaMat_;
  Matrix<Double>* deltaMat_lastIt_;
  Matrix<Double>* deltaMat_lastTS_;
  // resulting from function EstimateCurrentSlope
  Matrix<Double>* deltaMat_estimated_;

  /*
   * for deltaEvaluation version 40-50 (> cutting of input to saturation)
   * > by cutting the input to saturation, we avoid a jump in deltaMat from
   *   a higher deltaValue (e-4-e-7) to eps0/nu0 by simply reusing the last
   *   stored deltaMatrix;
   *   this has the problem, that we keep this old value of deltaMatrix as long
   *   as both the current input and the last input (that was used to compute
   *   deltaX) remain above saturation (which might be the case);
   *   if both values stay above saturation, we might actually need a slope of
   *   eps0/nu0 and the last value of deltaMat is simply to large
   * > this counter shall keep track of the double-cutting-cases; if multiple such
   *   cases appear in a row, we slowly turn the previously stored state towards
   *   eps0/nu0 (e.g. by dividing the reused deltamat by 10^cuttingApplied
   */
  Vector<UInt> cuttingApplied_;
  // during each iteration we only want to allow one step down
  // (otherwise we would step down multiple times per element one for each integration point)
  // RUN_ as it has to be reset during runtime via SetRuntimeDependentFlag
  Vector<UInt> RUN_cuttingAlreadyCounted_;


  // flags for checking if deltaMat was alreday computed during the current timestep / iteration
  // has to be reset after end of timestep / iteration
  Vector<UInt> RUN_deltaMatComputedDuringCurrentTS_;
  Vector<UInt> RUN_deltaMatComputedDuringCurrentIt_;


  // needed to get integration point
  shared_ptr<FeSpace> ptFeSpace_;
  UInt numIntegrationPoints_;
  std::map<int,UInt> locPointIndices_;



  // for bisection stuff
  Vector<Double>* dY_sol;
  Vector<Double>* X_low;
  Vector<Double>* X_up;

};


}
#endif
