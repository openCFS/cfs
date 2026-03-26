#ifndef FILE_COEFFUNCTION_HYST_HH
#define FILE_COEFFUNCTION_HYST_HH

#include "General/Environment.hh"
#include "CoefFunction.hh"
#include "Materials/BaseMaterial.hh"
#include "FeBasis/FeFunctions.hh"
#include "Utils/Timer.hh"
#include "Utils/helperStructs.hh"
//#include "Materials/Models/Preisach.hh"

namespace CoupledField {

  // forward class declaration
  class ApproxData;
  class BaseBOperator;
  class Grid;
  class FeFunctions;
  class CoefFunctionHyst;

  class CoefFunctionHystMatTensor : public CoefFunction {
  public:

    // note: coupling tensor has to be in non-transposed form (i.e. 3x6 or 2x3, 2x4)
    //CoefFunctionHystMatTensor(CoefFunctionHelper* hystHelper, std::string tensorName, bool transposed)
    CoefFunctionHystMatTensor(PtrCoefFct hystCoefFnc, std::string tensorName, bool transposed)
    :CoefFunction(){

      isAnalytic_ = false;
      isComplex_  = false;
      dependType_ = SOLUTION;

      // HystMat specific
      hystCoefFunction_ = hystCoefFnc;
      dimType_ = TENSOR;
      tensorName_ = tensorName;
      transposed_ = transposed;

      hystCoefFunction_->GetTensorSize(numRows_,numCols_);
    }

    virtual ~CoefFunctionHystMatTensor(){};

    virtual string GetName() const { return "CoefFunctionHystMatTensor"; }



    //    void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
    //      EXCEPTION("CoefFunctionHystMatTensor only returns tensors")
    //    }
    //
    //    void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
    //      EXCEPTION("CoefFunctionHystMatTensor only returns tensors")
    //    }
    //
    //    virtual UInt GetVecSize() const {
    //      EXCEPTION("GetVecSize not used here")
    //    }

    void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
      //std::cout << "Coef Function Hyst Mat - Get Tensor" << std::endl;
      //      std::cout << "GET TENSOR - 1" << std::endl;
      bool rotate = true;
      bool useAbs = false;
      std::string implementationVersion = "Division";

      /*
       * timelevel determines which value of P shall be used to evaluate the
       * scaling and rotation of the matrices;
       * the tensors returned by this functions should normally be on the current
       * time level, i.e. timelevel should be 0
       */
      hystCoefFunction_->ComputeTensor(outputTensor, lpm, tensorName_, implementationVersion, transposed_, rotate, useAbs, false );
      //      std::cout << "GET TENSOR - 2" << std::endl;
    }

    void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
      // get maximal value of diagonal mat tensor
      UInt diagSize = std::min(numRows_,numCols_);
      Matrix<Double> outputTensor = Matrix<Double>(numRows_,numCols_);
      GetTensor(outputTensor,lpm);

      Double maxVal = -1;
      Double absMaxVal = -1;
      Double curVal;
      for(UInt i = 0; i < diagSize; i++){
        curVal = outputTensor[i][i];
        if(abs(curVal)>absMaxVal){
          absMaxVal = abs(curVal);
          maxVal = curVal;
        }
      }
      outputScalar = maxVal;
    }

    //! Return row and columns size of tensor if coefficient function is a tensor
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      if( (tensorName_ == "Permittivity") || (tensorName_ == "Reluctivity")){
        hystCoefFunction_->GetTensorSize(numRows,numCols);//,tensorName_);
      } else if ( (tensorName_ == "CouplingMechToElec") || (tensorName_ == "CouplingElecToMech")
        || (tensorName_ == "CouplingMechToMag") || (tensorName_ == "ComputeMagToMech") ){
        hystCoefFunction_->GetCouplTensorSize(numRows,numCols);
      } else {
        EXCEPTION("Tensor type unknown");
      }
      if(transposed_){
        UInt tmp = numRows;
        numRows = numCols;
        numCols = tmp;
      }
    }

    std::string ToString() const {
      return "Dummy string";
    }

  private:
    UInt numRows_;
    UInt numCols_;
    PtrCoefFct hystCoefFunction_;
    std::string tensorName_;
    bool transposed_;
  };

  class CoefFunctionHystRHSLoad : public CoefFunction {

  public:

    CoefFunctionHystRHSLoad(PtrCoefFct hystCoefFnc,
      std::string vectorName, bool transposed, bool onBoundary, PtrCoefFct coefFunctionToBeIncluded):
      CoefFunctionHystRHSLoad(hystCoefFnc,
      vectorName, transposed, onBoundary){
        coefFunctionToBeIncluded_ = coefFunctionToBeIncluded;
    }

    CoefFunctionHystRHSLoad(PtrCoefFct hystCoefFnc,
      std::string vectorName, bool transposed, bool onBoundary)
    :CoefFunction(){

      isAnalytic_ = false;
      isComplex_  = false;
      dependType_ = SOLUTION;
      coefFunctionToBeIncluded_ = NULL;

      // HystRHSLoad specific
      hystCoefFunction_ = hystCoefFnc;
      dimType_ = VECTOR;
      vectorName_ = vectorName;
      transposed_ = transposed;
      onBoundary_ = onBoundary;
      usedAsRHSload_ = true;
    }

    virtual ~CoefFunctionHystRHSLoad(){};

    virtual string GetName() const { return "CoefFunctionHystRHSLoad"; }

    void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
      EXCEPTION("CoefFunctionHystRHSLoad only returns vectors")
    }

    void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
      //std::cout << "Coef Function Hyst RHS Load - Get Vector" << std::endl;
      // return vector that can be put onto rhs of system;
      // the sign of term shall be such, that it can be added with + in the
      // corresponding pdes
      // Example: in piezo systems we have (among others) -int_Omega (BN)^T*P dOmega
      //          on rhs, so we would return -P here

      int timeLevel;
      Double baseSign;

      if(onBoundary_){
        // boundary terms usually are evluated with the current value of polarization etc
        // whereas the rhs depends on usage of deltaMat or not
        timeLevel = hystCoefFunction_->GetTimeLevel("BC");
        //std::cout << "Boundary Term, timelevel = " << timeLevel << std::endl;
        // furtherrmore, the  boundary terms and the (volume) loads have opposite
        // aigns
        baseSign = -1.0;
      } else {
        timeLevel = hystCoefFunction_->GetTimeLevel("RHS");
        //std::cout << "RHS Term, timelevel = " << timeLevel << std::endl;
        baseSign = 1.0;
      }

      hystCoefFunction_->ComputeVector(outputVector,lpm, timeLevel, baseSign, vectorName_, onBoundary_, usedAsRHSload_ );

      if(coefFunctionToBeIncluded_ != NULL){
//        std::cout << "coefFunctionToBeIncluded_ found! > subtract requested tensor from the value of this function" << std::endl;
//        std::cout << "Requested vector: " << outputVector.ToString() << std::endl;
        UInt vecSize = outputVector.GetSize();
        Vector<Double> tmp = Vector<Double>(vecSize);
        coefFunctionToBeIncluded_->GetVector(tmp,lpm);
//        std::cout << "Additional vector: " << tmp.ToString() << std::endl;
        tmp.Add(-1.0,outputVector);
        outputVector = tmp;
//        std::cout << "Additional vector - requestedTensor: " << outputVector.ToString() << std::endl;
      }

    }

    virtual UInt GetVecSize() const {

      UInt numRows,numCols;
      if( (vectorName_ == "ElecPolarization") || (vectorName_ == "MagPolarization") || (vectorName_ == "MagMagnetization")){
        hystCoefFunction_->GetTensorSize(numRows,numCols);//,vectorName_);
      } else if ( (vectorName_ == "PiezoLoadForMechPDE") || (vectorName_ == "PiezoLoadForElecPDE")
        || (vectorName_ == "MagStrictLoadForMechPDE") || (vectorName_ == "MagStrictLoadForMagPDE") ){
        hystCoefFunction_->GetCouplTensorSize(numRows,numCols);
      } else {
        EXCEPTION("Vector type unknown");
      }
      if(transposed_){
        return numCols;
      } else {
        return numRows;
      }
    }

    //    void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
    //      EXCEPTION("CoefFunctionHystRHSLoad only returns vectors")
    //    }
    //

    // NEW: needed when CoefFunctionHelper::PrecomputeMaterialTensorForInverison() is called by this class
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      hystCoefFunction_->GetTensorSize(numRows,numCols);//,"None");
    }

    std::string ToString() const {
      return "Dummy string";
    }

  private:
    PtrCoefFct coefFunctionToBeIncluded_;
    PtrCoefFct hystCoefFunction_;
    std::string vectorName_;
    bool transposed_;
    bool onBoundary_;
    bool usedAsRHSload_;
  };

  class CoefFunctionHystOutput : public CoefFunction {

  public:

    CoefFunctionHystOutput(PtrCoefFct hystCoefFnc ,std::string ResultName):CoefFunction(){

      //      std::cout << "HystOutput -Generate " << ResultName << std::endl;
      if( (ResultName == "DeltaPermeability")||(ResultName == "DeltaPermittivity")
        //        ||(ResultName == "IrrStressesPiezo")||(ResultName == "IrrStrainsPiezo") ){
        ||(ResultName == "IrrStressesPiezo_TensorForm")||(ResultName == "IrrStrainsPiezo_TensorForm") ){
        dimType_ = TENSOR;
      } else {
        dimType_ = VECTOR;
      }
      //      std::cout << "dimType_ = " << dimType_ << std::endl;

      isAnalytic_ = false;
      isComplex_  = false;
      dependType_ = SOLUTION;
      onBoundary_ = false;
      usedAsRHSload_ = false;

      hystCoefFunction_ = hystCoefFnc;
      resultName_ = ResultName;
    }

    virtual ~CoefFunctionHystOutput(){};

    virtual string GetName() const { return "CoefFunctionHystOutput"; }


    void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
      EXCEPTION("CoefFunctionHystOutput returns no scalars")
    }

    void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
      //      std::cout << "HystOutput -getVector " << resultName_ << std::endl;
      int timeLevel = 0;
      Double baseSign = 1.0;
      //std::cout << "result " << resultName_ << " requested" << std::endl;
      hystCoefFunction_->ComputeVector(outputVector,lpm, timeLevel, baseSign, resultName_, onBoundary_, usedAsRHSload_ );
    }

    //    void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
    //
    //      //std::cout << "Coef Function Hyst Output - Get Vector" << std::endl;
    //      // output is always the actual quantity
    //      int timeLevelOutput = 0;
    //      if(resultName_ == "ElecPolarization"){
    //        outputVector = hystCoefFunction_->GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevelOutput);
    //      } else if (resultName_ == "MagPolarization"){
    //        outputVector = hystCoefFunction_->GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevelOutput);
    //      } else {
    //        EXCEPTION("Unknown result")
    //      }
    //      //std::cout << "Output: " << outputVector.ToString() << std::endl;
    //    }

    void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
      //      std::cout << "HystOutput -getTensor " << resultName_ << std::endl;
      if( (resultName_ == "DeltaPermeability")||(resultName_ == "DeltaPermittivity") ){
        //TODO: estimate permeability or permittivity around current working point
        UInt numCols,numRows;
        hystCoefFunction_->GetTensorSize(numRows,numCols);//,resultName_);
        outputTensor.Resize(numRows,numCols);
        outputTensor.Init();
      } else {
        hystCoefFunction_->ComputeTensor(outputTensor,lpm, resultName_, "None", false, false, false, false );
      }
    }

    virtual UInt GetVecSize() const {
      UInt vecSize = hystCoefFunction_->GetVecSize();
      if(resultName_ == "IrrStressesPiezo_VectorForm"){
        if(vecSize == 2){
          return 3;
        } else {
          return 6;
        }
      } else {
        return vecSize;
      }
    }

    //! Return row and columns size of tensor if coefficient function is a tensor
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      hystCoefFunction_->GetTensorSize(numRows,numCols);//,resultName_);
    }

    std::string ToString() const {
      return "Dummy string";
    }

  private:
    PtrCoefFct hystCoefFunction_;
    std::string resultName_;
    bool onBoundary_;
    bool usedAsRHSload_;
  };



  // ============================================================================
  //  Hysteresis
  // ============================================================================
  //! Provide a coefficient for hysteresis modeling
  //! \note This class only works for real-valued scalar data.
  class CoefFunctionHyst : public CoefFunction, public enable_shared_from_this<CoefFunctionHyst> {
  public:

    //! Constructor
    CoefFunctionHyst( BaseMaterial* const material,
      shared_ptr<ElemList> actSDList,
      PtrCoefFct dependency1,
      SubTensorType tensorType,
      MaterialType matType,
      shared_ptr<FeSpace> ptFeSpace);

    //! Constructor
    CoefFunctionHyst( BaseMaterial* const material,
    shared_ptr<ElemList> actSDList,
    PtrCoefFct dependency1,
    PtrCoefFct dependency2,
    SubTensorType tensorType,
    MaterialType matType,
    shared_ptr<FeSpace> ptFeSpace);

    virtual string GetName() const { return "CoefFunctionHyst"; }


    void Init(BaseMaterial* const material,
    shared_ptr<ElemList> actSDList,
    PtrCoefFct dependency1,
    SubTensorType tensorType,
    MaterialType matType,
    shared_ptr<FeSpace> ptFeSpace);

    /*
     * small helper to ensure same precision for directly traced hysteresis data and loaded data
     */
    inline Double restrictPrecision(Double in, UInt precision){
      std::stringstream tmp;
      tmp <<std::scientific<<std::setprecision(precision)<<in;
      return std::stod(tmp.str());
    }
    
    /*
     * taken from former hystHelper
     */
    void SetParamsOfFormerHystHelper(PtrCoefFct ptrFieldTensor,PtrCoefFct ptrElasticityTensor,
      PtrCoefFct ptrCouplingTensor){

      formerHystHelperParamsSet_ = true;
      ptrFieldTensor_ = ptrFieldTensor;
      ptrElastTensor_ = ptrElasticityTensor;
      ptrCouplTensor_ = ptrCouplingTensor;
      materialTensorsSetOnce_ = false;

      UInt couplRows, couplCols;
      if(ptrCouplTensor_ != NULL){
        ptrCouplTensor_->GetTensorSize(couplRows,couplCols);
        if(couplCols < couplRows){
          EXCEPTION("Coupling tensor has to be in non-transposed form");
        }
        isCoupled_ = true;
      } else {
        isCoupled_ = false;
      }

      strainForm_ = GetStrainForm();

      //      if(isCoupled_ == false){
      //        // no coupling tensors at all
      //        strainForm_ = -1;
      //      } else {
      //        bool useStrainForm = hystCoefFunction_->useStrainForm();
      //        if(!useStrainForm){
      //          // just take input tensors as they are, i.e. e-form or h-form
      //          strainForm_ = 0;
      //        }else{
      //          std::string PDEName = hystCoefFunction_->getPDEName();
      //          if(PDEName == "Electrostatic"){
      //            // piezo case
      //            // compute d-form as basis for coupling
      //            strainForm_ = 1;
      //          } else {
      //            // magstrict case
      //            // compute g-form
      //            strainForm_ = 2;
      //          }
      //        }
      //      }

      // !!! small material parameter are assumed to be constant over region !!!
      // i.e. for each region, just one linear material tensor of eac kind stored
      // to obtain those values, we have to evaluate the ptrCoefFunctions at an lpm
      // do this during the first call to getTensor
      tensorsInitialized_ = false;
    }

    void InitLinearTensors(const LocPointMapped& lpm){
      //std::cout << "Init Linear Tensors" << std::endl;

      if(tensorsInitialized_){
        return;
      }

      // !!! small material parameter are assumed to be constant over region !!!
      // i.e. for each region, just one linear material tensor of eac kind stored
      // to obtain those values, we have to evaluate the ptrCoefFunctions at an lpm
      // do this during the first call to getTensor
      UInt numRows, numCols;
      ptrFieldTensor_->GetTensorSize(numRows,numCols);
      Matrix<Double> epsS_nuS = Matrix<Double>(numRows,numCols);
      ptrFieldTensor_->GetTensor(epsS_nuS,lpm);

//      std::cout << "epsS_nuS: " << epsS_nuS.ToString() << std::endl;


      if(ptrCouplTensor_ == NULL){
        //        std::cout << "StrainForm = " << strainForm_ << "(due to material file)" << std::endl;
        //        std::cout << "However, no coupling was defined in .xml file." << std::endl;
        strainForm_ = -2;
        SetStrainForm(-2);
      }

      if(strainForm_ == -2){
//        std::cout << "Single field case" << std::endl;
        /*
         * single field case
         * > only fieldTensor needed
         * > remaining tensors stay unitinialized
         */
        eps_nu_base_ = epsS_nuS;

        //std::cout << "Initialized tensors: " << std::endl;
        //std::cout << "Field Tensor: " << eps_nu_base_.ToString() << std::endl;
      } else if(strainForm_ != -2){
//        std::cout << "Coupled case" << std::endl;
        /*
         * Coupled case
         */
        // get e-form / h-form
        UInt numRows1, numCols1, numRows2, numCols2;
        ptrElastTensor_->GetTensorSize(numRows1,numCols1);
        Matrix<Double> cE_B = Matrix<Double>(numRows1,numCols1);
        ptrElastTensor_->GetTensor(cE_B,lpm);

        ptrCouplTensor_->GetTensorSize(numRows2,numCols2);

        // coupling tensor should have the same number of rows as permittivity/reluctivity tensor
        // and the same number of columns as elasticity tensor
        assert(numRows2 == numRows);
        assert(numCols2 == numCols1);

        Matrix<Double> e_h = Matrix<Double>(numRows2,numCols2);
        ptrCouplTensor_->GetTensor(e_h,lpm);

        if(strainForm_ == -1){
          std::cout << "Coupled case - coupling only via irreversible parts on rhs" << std::endl;
          // just keep e-form/h-form
          eps_nu_base_ = epsS_nuS;
          couplTensor_ = e_h;
          couplTensor_.Init(); // just set coupling tensor to 0
          elastTensor_ = cE_B;
        } else if(strainForm_ == 0){
          //std::cout << "Coupled case - e/h form" << std::endl;
          // just keep e-form/h-form
          eps_nu_base_ = epsS_nuS;
          couplTensor_ = e_h;
          elastTensor_ = cE_B;
        } else if(strainForm_ == 1){
          //std::cout << "Coupled case - d form" << std::endl;
          // piezo
          // transform to d-form
          // d = cE^-1*e = sE*e
          Matrix<Double> sE = Matrix<Double>(numRows1,numCols1);
          cE_B.Invert(sE);

          Matrix<Double> d = Matrix<Double>(numRows2,numCols2);
          sE.Mult(e_h,d);

          // epsT = epsS + d e^t
          Matrix<Double> epsT = Matrix<Double>(numRows,numCols);
          Matrix<Double> e_trans = Matrix<Double>(numCols2,numRows2);
          e_h.Transpose(e_trans);
          d.Mult(e_trans,epsT);
          epsT.Add(1.0,epsS_nuS);

          eps_nu_base_ = epsT;
          couplTensor_ = d;
          // NOTE: we do not keep sE but take cE instead
          // later during the solution, we will need cE (not sE) and as cE is
          // not depending on hysteresis, we can directly use the linear tensor
          elastTensor_ = cE_B;
        } else if(strainForm_ == 2){
          //std::cout << "Coupled case - g form" << std::endl;
          // magstrict
          // transform to g-form
          // g = cE_B^-1*h = sB*h
          Matrix<Double> sB = Matrix<Double>(numRows1,numCols1);
          cE_B.Invert(sB);

          Matrix<Double> g = Matrix<Double>(numRows2,numCols2);
          sB.Mult(e_h,g);

          // nuT = nuS - g h^t
          Matrix<Double> nuT = Matrix<Double>(numRows,numCols);
          Matrix<Double> h_trans = Matrix<Double>(numCols2,numRows2);
          e_h.Transpose(h_trans);
          g.Mult(h_trans,nuT);
          nuT.Add(-1.0,epsS_nuS);
          nuT *= (-1.0);

          eps_nu_base_ = nuT;
          couplTensor_ = g;
          // NOTE: we do not keep sE but take cE instead
          // later during the solution, we will need cE (not sE) and as cE is
          // not depending on hysteresis, we can directly use the linear tensor
          elastTensor_ = cE_B;

        }
        //        std::cout << "Initialized tensors: " << std::endl;
        //        std::cout << "Field Tensor: " << eps_nu_base_.ToString() << std::endl;
        //        std::cout << "Coupling Tensor: " << couplTensor_.ToString() << std::endl;
        //        std::cout << "Elasticity Tensor: " << elastTensor_.ToString() << std::endl;
      }

      tensorsInitialized_ = true;

      bool testSettingOfFPTensor = !true;
      if(testSettingOfFPTensor){
        SetFPMaterialTensorsNEW(0);
        SetFPMaterialTensors(0);
        SetFPMaterialTensorsNEW(1);
        SetFPMaterialTensors(1);
        SetFPMaterialTensorsNEW(2);
        SetFPMaterialTensors(2);
        SetFPMaterialTensorsNEW(3);
        exit(0);
      }

    }

    void PrecomputeMaterialTensorForInverison();

    void ComputeTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm,
    std::string tensorName, std::string implementationVersion, bool transposed, bool rotate, bool useAbs,
    bool lockPrecomputationAndDeltaMat=false );

    void ComputeVector(Vector<Double>& outputVector,const LocPointMapped& lpm, int timeLevel, int baseSign, std::string vectorName, bool onBoundary, bool usedAsRHSload );

    void GetCouplTensorSize( UInt& numRows, UInt& numCols ) {
      ptrCouplTensor_->GetTensorSize(numRows,numCols);
    }

    Matrix<Double> GetBaseCouplingTensor(){
      if(tensorsInitialized_ == false){
        EXCEPTION("Coupling tensor not initialized yet");
      }
      return couplTensor_;
    }

    /*
     * New functions added April 2018
     */
    Double getWeight(ParameterPreisachWeights* weightParams, Double alpha, Double beta, Double delta){

      std::string weightType = weightParams->weightType_;
      UInt numRows = weightParams->numRows_;

      if(weightType == "Constant"){
        return weightParams->constWeight_;
      } else if((weightType == "muLorentz")||(weightType == "muLorentzExtended")){
        return LorentzianWolf(weightParams->muLorentz_A_, weightParams->muLorentz_h1_,
          weightParams->muLorentz_h2_, weightParams->muLorentz_sigma1_, 
          weightParams->muLorentz_sigma2_, alpha, beta);   
      } else if(weightType == "muDat"){
        return MuDat(weightParams->muDat_A_, weightParams->muDat_sigma1_, weightParams->muDat_h1_, weightParams->muDat_eta_, alpha, beta);
      } else if(weightType == "muDatExtended"){
        return MuDatExtended(weightParams->muDat_A_, weightParams->muDat_sigma1_, weightParams->muDat_sigma2_,
          weightParams->muDat_h1_, weightParams->muDat_h2_, weightParams->muDat_eta_, alpha, beta);
      } else if(weightType == "givenTensor"){
        // alpha = -1 shall be mapped to 0, alpha = 1.0 shall be mapped to maxNum
        UInt idxAlpha = UInt(std::round((alpha+1.0)/delta));
        UInt idxBeta = UInt(std::round((beta+1.0)/delta));
        if( idxAlpha < 0 ){
          idxAlpha = 0;
        } else if (idxAlpha >= numRows){
          idxAlpha = numRows-1;
        }
        if( idxBeta < 0 ){
          idxBeta = 0;
        } else if (idxBeta >= numRows){
          idxBeta = numRows-1;
        }

        return weightParams->weightTensor_[idxAlpha][idxBeta];
      } else {
        return 0.0;
      }
    }

    Double getWeightDerivative(ParameterPreisachWeights* weightParams, Double alpha, Double beta, Double delta, bool flipped){

      Double s, lambda;
      if(flipped){
        s = beta;
        lambda = alpha/beta;
      } else {
        s = alpha;
        lambda = beta/alpha;
      }

      std::string weightType = weightParams->weightType_;

      if(weightType == "Constant"){
        return 0.0;
      } else if((weightType == "muLorentz")||(weightType == "muLorentzExtended")){
        return dLorentzianWolf_by_ds(weightParams->muLorentz_A_, weightParams->muLorentz_h1_,
          weightParams->muLorentz_h2_, weightParams->muLorentz_sigma1_, 
          weightParams->muLorentz_sigma2_, s, lambda,flipped); 
      } else if(weightType == "muDat"){
        return dMuDat_by_ds(weightParams->muDat_A_, weightParams->muDat_sigma1_, weightParams->muDat_h1_, 
          weightParams->muDat_eta_,s,lambda,flipped);
      } else if(weightType == "muDatExtended"){
        return dMuDatExtended_by_ds(weightParams->muDat_A_, weightParams->muDat_sigma1_, weightParams->muDat_sigma2_,
          weightParams->muDat_h1_, weightParams->muDat_h2_, weightParams->muDat_eta_,s,lambda,flipped);
      } else if(weightType == "givenTensor"){
        Double s_low = s-delta/1e11;
        if( s_low < -1.0){
          s_low = -1.0;
        }
        Double s_up = s+delta/1e11;
        if( s_up > 1.0){
          s_up = 1.0;
        }

        Double deltas = s_up - s_low;
        Double wUp;
        Double wLow;
        if(flipped){
          wUp = getWeight(weightParams,lambda*s_up,s_up,delta);
          wLow = getWeight(weightParams,lambda*s_low,s_low,delta);
        } else{
          wUp = getWeight(weightParams,s_up,lambda*s_up,delta);
          wLow = getWeight(weightParams,s_low,lambda*s_low,delta);
        }
        return (wUp-wLow)/deltas;
      } else {
        return 0.0;
      }
    }

    ParameterPreisachWeights POL_weightParams_;
    ParameterPreisachWeights STRAIN_weightParams_;

    ParameterPreisachWeights POL_weightParamsForTesting_;
    ParameterPreisachWeights STRAIN_weightParamsForTesting_;
    
    ParameterPreisachOperators POL_operatorParams_;
    ParameterPreisachOperators STRAIN_operatorParams_;

    ParameterIrrStrainsAndCoupling CouplingParams_;


    bool inversionSet_;
    ParameterInversion InversionParams_;
    InitialInput POL_initial_;
    std::string STRAIN_usedHystModel_;
    int POL_weightsAlreadyAdapted_; // for Mayergoyz case
    int STRAIN_weightsAlreadyAdapted_; // for Mayergoyz case

    Double LorentzianWolf(Double A, Double h1, Double h2, Double sigma1, Double sigma2, Double alpha, Double beta){
      // Lorentzian weight function as used by F. Wolf in his PHD thesis
      // F. Wolf: Generalisiertes Preisach-Modell fuer die Simulation und Kompensation der Hysterese piezokeramischer Aktoren
      Double alphaFactor = A/(1 + ((alpha-h2)/(h2)*sigma2)*((alpha-h2)/(h2)*sigma2) );
      Double betaFactor = A/(1 + ((beta+h1)/(h1)*sigma1)*((beta+h1)/(h1)*sigma1) );
 
      return alphaFactor*betaFactor;
    }
    
    // ALREADY DONE IN XMLMaterialHandler.cc direct after reading in!
//    // move this directly to paramater reader; factors must be computed only once
//    Double LorentzianBertotti(Double C1, Double C2, Double C3, Double alpha, Double beta){
//      // Lorentzian weight function as used for Team32 and defined by Bertotti in
//      // G. Bertotti, Hysteresis in Magnetism: For Physicists, Materials Scientists, and Engineers , Academic Press, 1998
//      Double commonFactor =  C1/(C2*C3*std::sqrt(M_PI*(M_PI/2+std::atan(1/C2))));
////      Double alphaFactor = 1/(1 + ((alpha-C3)/(C2*C3))*((alpha-C3)/(C2*C3)) );
////      Double betaFactor = 1/(1 + ((-beta-C3)/(C2*C3))*((alpha-C3)/(C2*C3)) );
//// 
////      return commonFactor*commonFactor*alphaFactor*betaFactor;
////      
//      // Note: can be transferred to formulation used by Wolf in the following form:
//      /*
//       *  A = C1/(C2*C3*std::sqrt(M_PI*(M_PI/2+std::atan(1/C2))))
//       *  ((alpha-h2)/h2*sigma2)^2 == ((alpha-C3)/(C2*C3))^2
//       *  => sigma2 = 1/C2, h2 = C3
//       *  ((-beta-h1)/h1*sigma1)^2 == ((beta+C3)/()C2*C3))^2
//       *  => sigma1 = 1/C2, h1 = C3 
//       */
//      return LorentzianWolf(commonFactor,C3,C3,1/C2,1/C2,alpha,beta);
//    }
    
    Double dLorentzianWolf_by_ds(Double A, Double h1, Double h2, Double sigma1, Double sigma2, Double s, Double lambda, bool flipped){
      // Derivative of LorentzianWolf(s,lambda*s) by s
      // Used to transfer Preisach weights from scalar to vector model (Mayergoyz model, isotropic)
      // (LorentzianWolf function with alpha = s, beta = lambda*s used)
      // if flipped: alpha = lambda*s, beta = s
      if(flipped == false){
        // Case 1: alpha = s, beta = lambda*s
        // w(s,lambda*s) = A/(1 + ((s-h2)/h2*sigma2)^2) * A/(1 + ((lambda*s+h1)/h1*sigma1)^2)
        // dw(s,lambda*s)/ds = -2A/(1 + ((s-h2)/h2*sigma2)^2)^2 * (s-h2)/h2*sigma2) * sigma2/h2  * A/(1 + ((lambda*s+h1)/h1*sigma1)^2)
        //                     -2A/(1 + ((lambda*s+h1)/h1*sigma1)^2)^2 * (lambda*s+h1)/h1*sigma1) * lambda*sigma1/h1  * A/(1 + ((s-h2)/h2*sigma2)^2)
        //  = -2A^2/[ (1 + ((lambda*s+h1)/h1*sigma1)^2) * (1 + ((s-h2)/h2*sigma2)^2) ] *
        //    [ (s-h2)/h2*sigma2) * sigma2/h2 * 1/(1 + ((s-h2)/h2*sigma2)^2) + (lambda*s+h1)/h1*sigma1) * lambda*sigma1/h1 * 1/(1 + ((lambda*s+h1)/h1*sigma1)^2) ]
        //  define: tmp_s = (s-h2)/h2*sigma2; tmp_lambdas = (lambda*s+h1)/h1*sigma1
        //  = -2A^2/[ (1 + tmp_lambdas^2)*(1 + tmp_s^2) ] * 
        //    [ tmp_s*sigma2/h2*1/(1 + tmp_s^2) +  tmp_lambdas*lambda*sigma1/h1*1/(1 + tmp_lambdas^2) ]
        //  define: tmp_s_denom = 1+tmp_s^2; tmp_lambdas_denom = 1+tmp_lambdas^2
        //  = -2A^2/(tmp_s_denom * tmp_lambdas_denom) * 
        //    [ tmp_s*sigma2/h2*1/tmp_s_denom + tmp_lambdas*lambda*sigma1/h1*1/tmp_lambdas_denom ]
        Double tmp_s = (s-h2)/h2*sigma2;
        Double tmp_lambdas = (lambda*s+h1)/h1*sigma1;
        Double tmp_s_denom = 1.0+tmp_s*tmp_s;
        Double tmp_lambdas_denom = 1.0+tmp_lambdas*tmp_lambdas;
        Double factor1 = -2*A*A/(tmp_s_denom * tmp_lambdas_denom);
        Double factor2 = tmp_s*sigma2/h2*1/tmp_s_denom + tmp_lambdas*lambda*sigma1/h1*1/tmp_lambdas_denom;
        return factor1*factor2;
      } else {
        // Case 2: alpha = lambda*s, beta = s
        // w(lambda*s,s) = A/(1 + ((lambda*s-h2)/h2*sigma2)^2) * A/(1 + ((s+h1)/h1*sigma1)^2)
        // dw(lambda*s,s)/ds = -2A/(1 + ((lambda*s-h2)/h2*sigma2)^2)^2 * (lambda*s-h2)/h2*sigma2) * lambda*sigma2/h2  * A/(1 + ((s+h1)/h1*sigma1)^2)
        //                     -2A/(1 + ((s+h1)/h1*sigma1)^2)^2 * (s+h1)/h1*sigma1) * sigma1/h1  * A/(1 + ((lambda*s-h2)/h2*sigma2)^2)
        //  = -2A^2/[ (1 + ((s+h1)/h1*sigma1)^2) * (1 + ((lambda*s-h2)/h2*sigma2)^2) ] *
        //    [ (lambda*s-h2)/h2*sigma2) * lambda*sigma2/h2 * 1/(1 + ((lambda*s-h2)/h2*sigma2)^2) + (s+h1)/h1*sigma1) * sigma1/h1 * 1/(1 + ((s+h1)/h1*sigma1)^2) ]
        //  define: tmp_s = (s+h1)/h1*sigma1; tmp_lambdas = (lambda*s-h2)/h2*sigma2
        //  = -2A^2/[ (1 + tmp_lambdas^2)*(1 + tmp_s^2) ] * 
        //    [ tmp_s*sigma1/h1*1/(1 + tmp_s^2) +  tmp_lambdas*lambda*sigma2/h2*1/(1 + tmp_lambdas^2) ]
        //  define: tmp_s_denom = 1+tmp_s^2; tmp_lambdas_denom = 1+tmp_lambdas^2
        //  = -2A^2/(tmp_s_denom * tmp_lambdas_denom) * 
        //    [ tmp_s*sigma1/h1*1/tmp_s_denom + tmp_lambdas*lambda*sigma2/h2*1/tmp_lambdas_denom ]
        Double tmp_s = (s+h1)/h1*sigma1;
        Double tmp_lambdas = (lambda*s-h2)/h2*sigma2;
        Double tmp_s_denom = 1.0+tmp_s*tmp_s;
        Double tmp_lambdas_denom = 1.0+tmp_lambdas*tmp_lambdas;
        Double factor1 = -2*A*A/(tmp_s_denom * tmp_lambdas_denom);
        Double factor2 = tmp_s*sigma1/h1*1/tmp_s_denom + tmp_lambdas*lambda*sigma2/h2*1/tmp_lambdas_denom;
        return factor1*factor2;
      }
    }
    

    Double MuDat(Double A, Double sigma, Double h, Double eta, Double alpha, Double beta){
      // MuDat function for evaluating Preisach Weights
      // weights directly usable for
      //  > scalar Preisach model
      //  > vector Preisach model (Sutor)
      // after transformation, weights usable for
      //  > vector Preisach model (Mayergoyz) ISOTROPIC materials
      //
      // Source: "A Preisach-based hysteresis model for magnetic and ferroelectric hysteresis" - A. Sutor 2010
      // note: at this point we assume that A,sigma,h and eta are meant for p,e in range [-1,1]
      //       the script for determining these parameters (as well as the above mentioned
      //        artcle) assuems p,e in range [-0.5,0.5] > parameters have to be transfereed to correct range
      //                prior to calling this function; this can either be done by user or by CFS if flag in mat.xml is set
      return A/(1 + std::pow( std::pow((alpha+beta)*sigma,2) + std::pow((alpha-beta-h)*sigma,2),eta));
    }

    Double MuDatExtended(Double A, Double sigma1, Double sigma2, Double h1, Double h2, Double eta, Double alpha, Double beta){
      // extended MuDat function for evaluating Preisach Weights
      // weights directly usable for
      //  > scalar Preisach model
      //  > vector Preisach model (Sutor)
      // after transformation, weights usable for
      //  > vector Preisach model (Mayergoyz) ISOTROPIC materials
      //
      // Source: "Generalisiertes Preisach Modell für die Simulation und Kompensation der Hysterese piezokeramischer Aktoren" - F. Wolf 2014
      // note: at this point we assume that A,sigma,h and eta are meant for p,e in range [-1,1]
      //       the script for determining these parameters (as well as the above mentioned
      //        artcle) assuems p,e in range [-0.5,0.5] > parameters have to be transfereed to correct range
      //                prior to calling this function; this can either be done by user or by CFS if flag in mat.xml is set
      return A/(1 + std::pow( std::pow((alpha+beta+h1)*sigma1,2) + std::pow((alpha-beta-h2)*sigma2,2),eta));
    }

    Double dMuDat_by_ds(Double A, Double sigma, Double h, Double eta, Double s, Double lambda, bool flipped){
      // Derivative of MuDat(s,lambda*s) by s
      // Used to transfer Preisach weights from scalar to vector model (Mayergoyz model, isotropic)
      // (MuDat function with alpha = s, beta = lambda*s used)
      // if flipped: alpha = lambda*s, beta = s
      Double tmp; // tmp = (alpha-beta)/s
      if(flipped){
        tmp = lambda-1;
      } else {
        tmp = 1-lambda;
      }
      Double tmp2 = lambda+1; // tmp2 = (alpha+beta)/s

      Double denominator = std::pow(1 + std::pow( std::pow(tmp2*s*sigma,2) + std::pow((tmp*s-h)*sigma,2),eta),2);
      Double nominator1 = -A*eta;
      Double nominator2 = std::pow( std::pow(tmp2*s*sigma,2) + std::pow((tmp*s-h)*sigma,2),eta-1);
      Double nominator3 = 2*s*std::pow(tmp2*sigma,2) + 2*(tmp*s-h)*sigma*sigma*tmp;

      return nominator1*nominator2*nominator3/denominator;
    }

    Double dMuDatExtended_by_ds(Double A, Double sigma1, Double sigma2, Double h1, Double h2, Double eta, Double s, Double lambda, bool flipped){
      // Derivative of MuDat(s,lambda*s) by s
      // Used to transfer Preisach weights from scalar to vector model (Mayergoyz model, isotropic)
      // (MuDat function with alpha = s, beta = lambda*s used)
      // if flipped: alpha = lambda*s, beta = s
      Double tmp; // tmp = (alpha-beta)/s
      if(flipped){
        tmp = lambda-1;
      } else {
        tmp = 1-lambda;
      }
      Double tmp2 = lambda+1; // tmp2 = (alpha+beta)/s

      Double denominator = std::pow(1 + std::pow( std::pow((tmp2*s+h1)*sigma1,2) + std::pow((tmp*s-h2)*sigma2,2),eta),2);
      Double nominator1 = -A*eta;
      Double nominator2 = std::pow( std::pow((tmp2*s+h1)*sigma1,2) + std::pow((tmp*s-h2)*sigma2,2),eta-1);
      Double nominator3 = 2*(tmp2*s+h1)*sigma1*sigma1*tmp2 + 2*(tmp*s-h2)*sigma2*sigma2*tmp;

      return nominator1*nominator2*nominator3/denominator;
    }

    /*
     * Functions needed for numerical integration with Gauss-Tschebyscheff
     */
    inline Double x_of_y(Double y,Double xMin,Double xMax){
      return ( (1-y)*xMin/2 + (1+y)*xMax/2 );
    }

    inline Double dx_by_dy(Double xMin, Double xMax){
      return ( -xMin/2 + xMax/2 );
    }

    Matrix<Double> getTschebyscheffPointsAndWeights(UInt numPoints){
      Matrix<Double> pointsAndWeights = Matrix<Double>(numPoints,2);

      for(UInt i = 0; i < numPoints; i++){
        pointsAndWeights[i][0] = cos( (2.0*i + 1.0)/(2.0*(numPoints))*M_PI);
        // include weighting term sqrt(1-x_i^2) into stepping weight (pi/n)
        pointsAndWeights[i][1] = M_PI/(numPoints)*sin( (2.0*i + 1.0)/(2.0*(numPoints))*M_PI);
      }
      return pointsAndWeights;
    }

    Matrix<Double> evaluatePreisachWeights(ParameterPreisachWeights* weightParams){
      Double alpha,beta;
      Double delta = 2.0/weightParams->numRows_;
      Double dimHalf = weightParams->numRows_/2.0;

      if(weightParams->weightType_ == "givenTensor"){
        std::cout << "++ Use given tensor of size " << weightParams->numRows_ << " x " << weightParams->numRows_ << " as Preisach weights" << std::endl;
        return weightParams->weightTensor_;
      }

      std::cout << "++ Evaluate analytic Preisach weighting function " << weightParams->weightType_ << " for " << weightParams->numRows_ << " x " << weightParams->numRows_ << " cells" << std::endl;
      Matrix<Double> weights = Matrix<Double>(weightParams->numRows_,weightParams->numRows_);

      for(UInt i = 0; i < weightParams->numRows_; i++){
        // note: we evaluated the weights always at the element center > add deltaS/2
        alpha = (i - dimHalf)*delta + delta/2;

        for(UInt k = 0; k <= i; k++){
          beta = (k - dimHalf)*delta + delta/2;

          // note: getWeight checks for POL_weightParams_.weightType_; depending on that value
          // we evaluate muDat, muDatExt, return a const or access the already given weights (in which case
          // there is no need to call this function)
          weights[i][k] = getWeight(weightParams,alpha,beta,delta);
        }
      }
      return weights;
    }

    Matrix<Double> transformPreisachWeightsForIsotropicVectorCase(ParameterPreisachWeights* paramSet){
      // Transform Preisach weights from scalar case to vector case (isotropic, Mayergoyz model)

//      if(dim_ != 2){
//        EXCEPTION("Currently only 2D case of Preisach Weight transformation supported");
//        // 3d case starts from different formula > see "Mathematical Models of Hysteresis and their Application" - Mayergoyz 2003
//      }

      /*
       * Note: this function is independent on the kind of weights!
       * it uses the general getWeight and getWeightDerivative Functions which return the weights
       * depending on the used weight type, e.g., muDat weightTensor or NEW Lorentzian weights
       */
      UInt numRows = paramSet->numRows_;

      Matrix<Double> vectorWeights;
      Matrix<Double> vectorWeights2D;
      Matrix<Double> vectorWeights3D;

      Double dimHalf = numRows/2.0;
      UInt dimHalfInt = UInt(dimHalf);
      Double alpha,beta,s,lambda;
      Double delta = 2.0/numRows;
      bool checkWeights2D = !true;
      bool checkWeights3D = !true;
      bool stopAfterTest = !true;

      if((dim_ == 2)||checkWeights2D) {
        // 2D Isotropic
        // Source: "Analysis of Isotropic Materials with Vector Hysteresis" - O. Bottauscio 1998
        //
        // p_vect(alpha,lambda*alpha) = \int_0^alpha (3*s*s*p_scalar(alpha,beta) + s*s*s*dp_scalar(alpha,beta)/ds  )/pi*alpha*alpha*sqrt(alpha*alpha - s*s) ds
        // p_vect(0,0) = 3/4*p_scal(0,0)
        //
        // some notes: - ONLY FOR 2D
        //             - p_vect = weights for scalar model inside the Mayergoyz vector model
        //             - p_scal = weights for scalar model
        //             - lambda = beta/alpha
        //             - if p_scal(0,0) = const > integral converges towards 3/4 (found out from numerical integration)
        //             - due to denominator of pi*alpha*alpha*sqrt(alpha*alpha - s*s), perform numerical integration with
        //                Gauss-Tschebyscheff (denominator cancels out with Tschebyscheff weights; much better accuracy as with std Gauss, tested!)
        //             - formula generally valid for alpha in [alpaMin,alphaMax], beta in [-alpha,alpha]
        //                > but: for alpha = 0; lambda = beta/alpha will not work
        //                > setup Preisach plane in two parts:
        //                      a) upper part: alpha in [0,alphaMax], beta in [-alpha,alpha]
        //                      b) lower part: beta in [betaMin,0], alpha in [-beta,beta]
        //
        vectorWeights2D = Matrix<Double>(numRows,numRows);
        UInt numIntegrationPoints = 15;
        Matrix<Double> integrationPoints = getTschebyscheffPointsAndWeights(numIntegrationPoints);

        // upper part
        for(UInt i = dimHalfInt; i < numRows; i++){
          // note: we evaluated the weights always at the element center > add deltaS/2
          alpha = (i - dimHalf)*delta + delta/2;//+ 1e-6*delta;
          for(UInt k = numRows-1-i; k <= i; k++){
            beta = (k - dimHalf)*delta + delta/2;//+ 1e-6*delta;// + delta/2;

            if(alpha == 0){
              vectorWeights2D[i][k] = 0.75*getWeight(paramSet, 0, 0, delta);
            } else {
              lambda = beta/alpha;

              vectorWeights2D[i][k] = 0.0;
              // numerically integrate from 0 to alpha
              for(UInt pp = 0; pp < numIntegrationPoints; pp++){
                Double curY = integrationPoints[pp][0];
                Double curW = integrationPoints[pp][1];
                //std::cout << "s / 3*s*s*getWeight / s*s*s*getWeightDerivative / (tmp + tmp2)/(M_PI * alpha * alpha * std::sqrt(alpha*alpha - s*s))" << std::endl;
                s = x_of_y(curY,0,alpha);
                Double tmp = 3*s*s*getWeight(paramSet, s,lambda*s,delta);

                bool flipped = false; // > s takes slot of alpha, lambda*s takes place of beta

                Double tmp2 = s*s*s*getWeightDerivative(paramSet, s,lambda*s,delta, flipped);
                Double tmp3 = (tmp + tmp2)/(M_PI * alpha * alpha * std::sqrt(alpha*alpha - s*s));

                vectorWeights2D[i][k] += dx_by_dy(std::min(0.0,alpha), std::max(0.0,alpha))*curW*tmp3;
                //vectorWeights[i][k] += dx_by_dy(0.0,alpha)*curW*tmp3;
              }
            }
          }
        }

        //      std::cout << "weights so far: " << std::endl;
        //      std::cout << vectorWeights.ToString() << std::endl;
        //
        // lower part
        //      std::cout << "dimHalfInt = " << dimHalfInt << std::endl;
        for(UInt k = 0; k < dimHalfInt; k++){
          //        std::cout << "k = " << k << std::endl;
          // note: we evaluated the weights always at the element center > add deltaS/2
          beta = (k - dimHalf)*delta + delta/2; //+ 1e-6*delta;// + delta/2;

          for(UInt i = k; i < numRows-k-1; i++){
            //          std::cout << "POL_numRows_-k-1 = " << POL_numRows_-k-1 << std::endl;
            //          std::cout << "i = " << i << std::endl;
            alpha = (i - dimHalf)*delta + delta/2; //+ 1e-6*delta;// + delta/2;

            if(beta == 0){
              vectorWeights2D[i][k] = 0.75*getWeight(paramSet,alpha, beta, delta);
            } else {
              // flipped case
              lambda = alpha/beta;
              vectorWeights2D[i][k] = 0.0;
              // numerically integrate from 0 to beta > however, beta < 0 so we have to reverse the sign
              // > done vis dx_by_dy
              for(UInt pp = 0; pp < numIntegrationPoints; pp++){
                Double curY = integrationPoints[pp][0];
                Double curW = integrationPoints[pp][1];

                s = x_of_y(curY,0,beta);
                Double tmp = 3*s*s*getWeight(paramSet,lambda*s,s,delta);

                bool flipped = true; // > s takes slot of alpha, lambda*s takes place of beta
                Double tmp2 = s*s*s*getWeightDerivative(paramSet,lambda*s,s,delta, flipped);

                Double tmp3 = (tmp + tmp2)/(M_PI * beta * beta * std::sqrt(beta*beta - s*s));
                vectorWeights2D[i][k] += dx_by_dy(std::min(0.0,beta), std::max(0.0,beta))*curW*tmp3;
                //vectorWeights[i][k] += dx_by_dy(0.0,beta)*curW*tmp3;
              }
            }
          }
        }

        if(checkWeights2D){
          std::cout << "=== Check isotropic weights for Mayergoyz Vector Model - 2D" << std::endl;

          // according to "Analysis of Isotropic Materials with Vector Hysteresis" - O. Bottauscio 1998
          // int_-pi/2^pi/2 cos(phi)^3 vectorWeights(alpha*cos(phi),beta*cos(phi) dphi = scalarWeight(alpha,beta)
          // this can be checked
          // due to discretization of vectorWeights, we will not get a perfect match (except for constant weights)
          // however, we come pretty close if Preisach resolution is fine enough!
          UInt numAngles = 101;
          Double deltaPhi = M_PI/numAngles;
          Double phi;
          Double alpha,beta,c;
          UInt numTests = 100;
          for(UInt k = 0; k < numTests; k++){
            alpha = 0 + k*0.01;
            if(k > numTests/2){
              beta = 0 - k*0.01; // some random values
            } else {
              beta = 0 + k*0.01; // some random values
            }

            Double integral = 0.0;
            Double target = getWeight(paramSet,alpha+delta/2, beta+delta/2, delta);
            for(UInt i = 0; i < numAngles; i++){
              phi = -M_PI/2+i*deltaPhi;
              c = std::cos(phi);
              int idxAlphaTMP = std::floor((alpha*c+1.0)/delta);
              int idxBetaTMP = std::floor((beta*c+1.0)/delta);
              UInt idxAlpha = UInt(std::max(0,idxAlphaTMP));
              UInt idxBeta = UInt(std::max(0,idxBetaTMP));
              integral += c*c*c*vectorWeights2D[idxAlpha][idxBeta];
            }
            //Double deltaAngle = M_PI/numDirections_;
            // Note: here we do not average out over the halfspace so not multiplication by 2/M_PI
            integral *= M_PI/numAngles;

            std::cout << "Alpha / Beta = " << alpha << " / " << beta << std::endl;
            std::cout << "absErr / relError = " << integral-target << " / " << integral/target - 1.0 << std::endl;
//            std::cout << "Scalar weight: " << target << std::endl;
//            std::cout << "int_-pi/2^pi/2 cos(phi)^3 vectorWeights(alpha*cos(phi),beta*cos(phi) dphi - scalarWeights(alpha,beta): " << integral-target << std::endl;
          }
        }
      }

      if((dim_ == 3)||checkWeights3D) {
        // 3D Isotropic
        // Source: Mathematical Models of Hysteresis and their Application" - Mayergoyz 2003
        //
        // P(alpha,lambda*alpha) = 1/(2*pi*alpha) * d/dalpha (alpha*alpha F(alpha,lambda*alpha) )
        // where P and F are p_vect and p_scal integrated over alpha,lambda*alpha
        // --> as described in Bottauscio: p_vect = -d/dalpha d/dbeta P(alpha,beta) and p_scal = -d/dalpha d/dbeta F(alpha,beta)
        //
        // after some calculus:
        // p_vect(alpha,lambda*alpha) = 1/pi (2*p_scalar(alpha,lambda*alpha) + alpha/2 dp_scalar(alpha,lambda*alpha)/dalpha)
        //
        // some notes: - to obtain the (hopefully) correct derivatives, we have to consider that alpha = beta/lambda and/or beta = lambda*alpha
        //             - 3D isotropic case actually much simpler (from mathematical point of view) than 2D isotropic case
        //             - formula generally valid for alpha in [alpaMin,alphaMax], beta in [-alpha,alpha]
        //                > but: for alpha = 0; lambda = beta/alpha will not work
        //                > setup Preisach plane in two parts:
        //                      a) upper part: alpha in [0,alphaMax], beta in [-alpha,alpha]
        //                      b) lower part: beta in [betaMin,0], alpha in [-beta,beta]
        vectorWeights3D = Matrix<Double>(numRows,numRows);

        // upper part
        for(UInt i = dimHalfInt; i < numRows; i++){
          // note: we evaluated the weights always at the element center > add deltaS/2
          alpha = (i - dimHalf)*delta + delta/2;//+ 1e-6*delta;
          for(UInt k = numRows-1-i; k <= i; k++){
            beta = (k - dimHalf)*delta + delta/2;//+ 1e-6*delta;// + delta/2;

            if(alpha == 0){
              vectorWeights3D[i][k] = 2.0/M_PI*getWeight(paramSet, 0, 0, delta);
            } else {
              // p_vect(alpha,lambda*alpha) = 1/pi (2*p_scalar(alpha,lambda*alpha) + alpha/2 dp_scalar(alpha,lambda*alpha)/dalpha)
              bool flipped = false;
              vectorWeights3D[i][k] = 2.0*getWeight(paramSet, alpha, beta, delta);
              vectorWeights3D[i][k] += alpha/2.0*getWeightDerivative(paramSet, alpha, beta, delta, flipped);
              vectorWeights3D[i][k] *= 1.0/M_PI;
            }
          }
        }

        //      std::cout << "weights so far: " << std::endl;
        //      std::cout << vectorWeights.ToString() << std::endl;
        //
        // lower part
        //      std::cout << "dimHalfInt = " << dimHalfInt << std::endl;
        for(UInt k = 0; k < dimHalfInt; k++){
          //        std::cout << "k = " << k << std::endl;
          // note: we evaluated the weights always at the element center > add deltaS/2
          beta = (k - dimHalf)*delta + delta/2; //+ 1e-6*delta;// + delta/2;

          for(UInt i = k; i < numRows-k-1; i++){
            //          std::cout << "POL_numRows_-k-1 = " << POL_numRows_-k-1 << std::endl;
            //          std::cout << "i = " << i << std::endl;
            alpha = (i - dimHalf)*delta + delta/2; //+ 1e-6*delta;// + delta/2;

            if(alpha == 0){
              vectorWeights3D[i][k] = 2.0/M_PI*getWeight(paramSet, 0, 0, delta);
            } else {
              // p_vect(alpha,lambda*alpha) = 1/pi (2*p_scalar(alpha,lambda*alpha) + alpha/2 dp_scalar(alpha,lambda*alpha)/dalpha)
              bool flipped = false;
              vectorWeights3D[i][k] = 2.0*getWeight(paramSet, alpha, beta, delta);
              vectorWeights3D[i][k] += alpha/2.0*getWeightDerivative(paramSet, alpha, beta, delta, flipped);
              vectorWeights3D[i][k] *= 1.0/M_PI;
            }
          }
        }

        if(checkWeights3D){
          std::cout << "=== Check isotropic weights for Mayergoyz Vector Model - 3D" << std::endl;
          // from Mayergoyz (p. 190; eq. 3.163)
          // int_phi=0^2pi int_theta=0^pi/2 cos(theta) P(alpha*cos(theta),beta*cos(theta))sin(theta) dtheta dphi = F(alpha,beta)
          // --> with p_vect = -d/dalpha d/dbeta P(alpha,beta) and p_scal = -d/dalpha d/dbeta F(alpha,beta)
          // int_phi=0^2pi int_theta=0^pi/2 cos(theta) -d/dalpha d/dbeta P(alpha*cos(theta),beta*cos(theta))sin(theta) dtheta dphi = p_scal(alpha,beta)
          // --> change of variables for derivative
          // int_phi=0^2pi int_theta=0^pi/2 cos^3(theta) -d/d(alpha*cos(theta)) d/d(beta*cos(theta)) P(alpha*cos(theta),beta*cos(theta))sin(theta) dtheta dphi = p_scal(alpha,beta)
          // int_phi=0^2pi int_theta=0^pi/2 cos^3(theta) p_vect(alpha*cos(theta),beta*cos(theta))sin(theta) dtheta dphi = p_scal(alpha,beta)
          // --> integrand independent of phi
          // 2*pi*int_theta=0^pi/2 cos^3(theta) sin(theta) p_vect(alpha*cos(theta),beta*cos(theta)) dtheta = p_scal(alpha,beta)
          // --> evaluate this integral for checking
          UInt numAngles = 101;
          Double deltaTheta = M_PI/numAngles;
          Double theta;
          Double alpha,beta,c;
          UInt numTests = 100;
          for(UInt k = 0; k < numTests; k++){
            alpha = 0 + k*0.01;
            if(k > numTests/2){
              beta = 0 - k*0.01; // some random values
            } else {
              beta = 0 + k*0.01; // some random values
            }

            Double integral = 0.0;
            Double target = getWeight(paramSet,alpha+delta/2, beta+delta/2, delta);
            for(UInt i = 0; i < numAngles; i++){
              theta = i*deltaTheta/2.0; // from 0 to M_PI/2
              c = std::cos(theta);
              s = std::sin(theta);
              UInt idxAlpha = UInt(std::round((alpha*c+1.0)/delta));
              UInt idxBeta = UInt(std::round((beta*c+1.0)/delta));
              integral += c*c*c*s*vectorWeights3D[idxAlpha][idxBeta];
            }
            //Double deltaAngle = M_PI/numDirections_;
            // Note: here we do not average out over the halfspace so not multiplication by 2/M_PI
            integral *= M_PI/2.0/numAngles;

            integral *= 2*M_PI; // integral over phi

            std::cout << "Alpha / Beta = " << alpha << " / " << beta << std::endl;
            std::cout << "absErr / relError = " << integral-target << " / " << integral/target - 1.0 << std::endl;
//            std::cout << "Scalar weight: " << target << std::endl;
//            std::cout << "2*pi*int_theta=0^pi/2 cos^3(theta) sin(theta) p_vect(alpha*cos(theta),beta*cos(theta)) dtheta - scalarWeights(alpha,beta): " << integral-target << std::endl;
//
//            std::cout << "Alpha = " << alpha << std::endl;
//            std::cout << "Beta = " << beta << std::endl;
//            std::cout << "Scalar weight: " << target << std::endl;
//            std::cout << "2*pi*int_theta=0^pi/2 cos^3(theta) sin(theta) p_vect(alpha*cos(theta),beta*cos(theta)) dtheta: " << integral << std::endl;
          }
        }
      }

      if(dim_ == 2){
        vectorWeights = vectorWeights2D;
      } else if(dim_ == 3){
        vectorWeights = vectorWeights3D;
      } else {
        EXCEPTION("dim_ " << dim_ << " not supported");
      }

      if( (checkWeights2D || checkWeights3D) && stopAfterTest ){
        std::cout << "End after check" << std::endl;
        exit(0);
      }

//            std::cout << "weights final: " << std::endl;
//            std::cout << vectorWeights.ToString() << std::endl;
      return vectorWeights;
    }

    void ReadInMaterial(BaseMaterial * const material);

    void ReadAndSetParamsForHystOperator(BaseMaterial * const material, bool setForStrains);

    void ReadAndSetWeights(BaseMaterial * const material, bool setForStrains, int forcedPreisachResolutionForTests = -1);

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
      numRows = eps_nu_base_.GetNumRows();
      numCols = eps_nu_base_.GetNumCols();
    }

    void SetFlag(std::string flagName, Integer intState);
    void SetDoubleFlag(std::string flagName, Double doubleState);

    Hysteresis* getTemporalHystOperator(bool forStrains, bool forceScalarDirection = false, UInt dimOfHystOperator = 2, int forcedPreisachResolutionForTests = -1, bool enforceContinuousEvaluation=false);

    void TraceHystOperator(UInt baseSteps, Double& maxSlope, Double& minSlope, Double& negCoercivity, Double& maxPolarization, bool dedicatedOperatorForStrains = false);
    void TraceHystOperatorVector(UInt baseSteps, Double& maxSlope, Double& minSlope, Double& negCoercivity, Double& maxPolarization, bool dedicatedOperatorForStrains = false);
 
    void SetFPMaterialTensorsNEW(Integer intState){
      EXCEPTION("SetFPMaterialTensorsNEW is Deprecated and should not be used");
      /*
       * NEW version after reading
       * "Analysis of the Convergence of the Fixed-Point Method Used for Solving Nonlinear Rotational Magnetic Field Problems"
       * - Dlala et al; 2008
       *
       * Global Fixpoint:
       *
       * nu_FP = 0.5 * ( dH/dB|min + dH/dB|max )
       *
       * where dH/dB|min / dH/dB|max are the smallest / maximal partial derivatives dH_i/dB_i over the whole possible
       * range of B, i.e. for the actual global case we have to find min and max slope in between negative and positve
       * saturation (or even beyond if polarization grows beyond saturation which can be the case if Mayergoyz is used
       * or if the anhysteretic paramter are set)
       *
       * > if polarization goes to real saturation, i.e. dB/dH grows only with nu_0 in saturation, we have to satisfy
       *    nu_0/2 < nu_FP <= nu_0
       *   in order to ensure the condition for contraction
       *    abs( 1/nu_FP dH/dB - 1 ) < 1
       * > if polarization does not saturate (due to anhysteretic parts or due to the excitation of the system staying
       *    away from saturation), nu_FP can be smaller than nu_0/2 (in fact I had some very nice convergence for
       *    Mayergoyz with anhysteretic parts for nu_0/8); however, it is risky to simply set a smaller nu_FP without
       *    knowing where the route goes
       * > in electrostatics, we have to use eps_FP which I assume should act similarily; here we have to be careful
       *    as we do not know the actual dD/dE as eps_0 is the minimum (whereas nu_0 is the maximum in magnetics)
       *
       * conclusion for setting nu_FP so small as possible and to find out dD/dE|max in electrostatics:
       *  > perform a complete run over one hysteresis operator going from N times negative saturation to N times positve
       *  saturation and determine the maximal and minimal dP/dH (dP/dE); from these values, determine maximal and minimal slopes
       *  > for this estimate, we have to actually traverse the major loop which means that the memory of the hyst-operator
       *  has to be set, i.e. perform this on a temporal hyst operator to not destroy the storage of an actual operator
       *  > perform this complete run only once!
       *
       * Lcoal Fixpoint
       *
       * nu_FP = C * ( dH/dB^n|min + dH/dB^n|max )
       *
       * where dH/dB^n is the slope at the last computed solution (last time step)
       * > factor C is a convergence factor that optimally should be 1 but as the slope is only estimated, C > 1 required
       * > the larger C, the larger the area of convergence but the slower the convergence
       * > the Authors of the mentioned paper suggest to start with a small C and perform a full (period) of simulation
       *  (they simulate periodic systems); in case of no-convergence C is increased
       *  --> this procedure seems not usable in general but one can show, that a linesearch-stepping with alpha basically
       *      results in C_effective = C/alpha, i.e. if alpha is smaller than 1, the convergence range is increased (as
       *      one would suspect from performing a linesearch)
       * > here we start at C = 2 (which seems to work quite well in practice according to the authors) and hope that
       *    linesearch will do the rest
       */
      Matrix<Double> FP_Tensor = Matrix<Double>(dim_,dim_);
      Matrix<Double> baseTensor = Matrix<Double>(dim_,dim_);
      Matrix<Double> tmp = Matrix<Double>(dim_,dim_);
      Matrix<Double> diagForm;

      UInt FPStateBAK = GetFPMaterialState();
      // temporarily disable FPMaterial state, such that GetTensor does not return the fixpoint tensor but
      // eps_S / nu_S (can be computed from eps_T or nu_T)
      SetFPMaterialState(0);

      // eps or mu might depend on solution, i.e. we have to trace each local point separately (unfortunately!)
      std::map<UInt, LocPointMapped >::iterator LPMit;
      std::map<UInt, LocPointMapped >::iterator LPMitStart = allLPMmap_.begin();
      std::map<UInt, LocPointMapped >::iterator LPMitEnd = allLPMmap_.end();

      std::string requiredTensorName;
      if(PDEName_ == "Electromagnetics"){
        requiredTensorName = "Reluctivity";
      } else {
        requiredTensorName = "Permittivity";
      }

      if(!hystModelTraced_){
        UInt baseSteps = 100;
        // get some additinal info if we are already traverse loop
        Double negCoercivity = 0.0;
        Double maxPolarization = 0.0;

        TraceHystOperator(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
        hystModelTraced_ = true;
      }

      bool testoutput = false;
      UInt operatorIdx, storageIdx, idxElem;
      LocPointMapped actualLPM;
      bool onBoundary;

      for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
        /*
         * check for boundary elements; on these surface elements we cannot eval the material relation and this
         * should not be required either
         */
        onBoundary = PreprocessLPM(LPMit->second, actualLPM, operatorIdx, storageIdx, idxElem);
        if(onBoundary){
          continue;
        }

        // last flag important -> no delta mat may be returned, so set it to true!
        // due to setting FPMaterialState to 0 we should get eps_S or nu_S here!
        ComputeTensor(baseTensor, LPMit->second, requiredTensorName, "DEFAULT", false, false, false, true );

        if(testoutput){
          std::cout << "baseTensor " << baseTensor.ToString() << std::endl;
        }

        if((intState == 0)||(intState == 3)){
          // > FP Tensor basically the current material tensor
          // > for H-Version (intState == 3)
          FP_Tensor = baseTensor;
        } else if(intState == 1){
          if(!globalFPset_){
            // GLOBAL FP Method
            // > use information of maximal and minimal slope of hyst operator
            // Electrostatics / Piezo:
            // FP_Tensor = C*0.5*( max(dD/dE) + min(dD/dE) )
            //  max(dD/dE) = max(diag(baseTensor)) + maxSlopeGlobal_
            //  min(dD/dE) = min(diag(baseTensor)) + minSlopeGlobal_
            // Magnetics / Magstrict:
            // FP_Tensor = C*0.5*( max(dH/dB) + min(dH/dB) )
            //  max(dH/dB) = 1/min(dB/dH) = 1/(min(diag(inv(baseTensor))) + minSlopeGlobal_)
            //  min(dH/dB) = 1/max(dB/dH) = 1/(max(diag(inv(baseTensor))) + maxSlopeGlobal_)

            // DEFAULT: 1
            //          Double globalFPFactor_C = 0.0625*2; //1.0;
//            std::cout << "Using globalFPFactor_C = " << globalFPFactor_C << std::endl;

            if(requiredTensorName == "Reluctivity"){
              baseTensor.Invert(tmp);

              if(testoutput){
                std::cout << "baseTensor inverted " << tmp.ToString() << std::endl;
              }

              tmp.GetDiagInMatrix(diagForm);
            } else {
              baseTensor.GetDiagInMatrix(diagForm);
            }

            Double slope1,slope2,FPSlope;
            slope1 = minSlopeGlobal_ + diagForm.GetMin();
            slope2 = maxSlopeGlobal_ + diagForm.GetMax();


            if(testoutput){
              std::cout << "minSlopeGlobal_ " << minSlopeGlobal_ << std::endl;
              std::cout << "maxSlopeGlobal_ " << maxSlopeGlobal_ << std::endl;
              std::cout << "diagForm.GetMin() " << diagForm.GetMin() << std::endl;
              std::cout << "diagForm.GetMax() " << diagForm.GetMax() << std::endl;
              std::cout << "slope1 " << slope1 << std::endl;
              std::cout << "slope2 " << slope2 << std::endl;
            }

            if(requiredTensorName == "Reluctivity"){
              slope1 = 1.0/slope1;
              slope2 = 1.0/slope2;
//              std::cout << "slope1 " << slope1 << std::endl;
//              std::cout << "slope2 " << slope2 << std::endl;
            }
            FPSlope = globalFPFactor_C*0.5*(slope1 + slope2);

            FP_Tensor.Init();
            for(UInt i = 0; i < dim_; i++){
              FP_Tensor[i][i] = FPSlope;
            }
            globalFPset_ = true;
          }
        } else {
          // LOCAL FP Method
          // > use information of maximal and minimal slope at current point
          // Electrostatics / Piezo:
          // FP_Tensor = C*0.5*( max(dD/dE) + min(dD/dE) )
          //  max(dD/dE) = max(diag(currentJacobian))
          //  min(dD/dE) = min(diag(currentJacobian))
          // Magnetics / Magstrict:
          // FP_Tensor = C*0.5*( max(dH/dB) + min(dH/dB) )
          //  max(dH/dB) = max(diag(currentJacobian))
          //  min(dH/dB) = min(diag(currentJacobian))
          // > note that currentJacobian is the Jacobian of the material relation, i.e. in case of magnetics it will
          //    will already return dH/dB
          // > the higher the factor C the slower convergence but the larger the convergence radius
          //    if we start with a low value here, we might compensate by linesearch (with maximal factor 1)
          //    if we start with a high value, we cannot compensate so easy

          // DEFAULT > 1
//          Double localFPFactor_C = 2.0;//8.0;//2.0;
//          std::cout << "Using localFPFactor_C = " << localFPFactor_C << std::endl;

          // new; set timelevel old to 3 forces jacobian; do not use storage values here but force new evaluation
          // (Note: the value would not be saved either, but this is not imporatant as the deltaMat storage is not required here)
          if(requiredTensorName == "Reluctivity"){
            baseTensor.Invert(tmp);
            // we have to pass mu here
            FP_Tensor = GetMaterialRelation(LPMit->second, 0, 3, tmp, "Reluctivity",true);
          } else {
            // we have to pass eps here
            FP_Tensor = GetMaterialRelation(LPMit->second, 0, 3, baseTensor, "Permittivity",true);
          }

          if(testoutput){
            std::cout << "Resulting Jacobian: " << FP_Tensor.ToString() << std::endl;
          }

          if(intState == 3){
//            std::cout << "Local FP B - take full Jacobian" << std::endl;
            // take full Jacobian as FP_Tensor
            for(UInt i = 0; i < dim_; i++){
              for(UInt j = 0; j < dim_; j++){
                // just for testing
                FP_Tensor[i][j] *= localFPFactor_C;
              }
            }
          } else {
            // compute FP Tensor as stated in literature
//            std::cout << "Local FP B - take min and max diagonal entries of Jacobian" << std::endl;
            FP_Tensor.GetDiagInMatrix(diagForm);

            if(testoutput){
              std::cout << "diagForm.GetMin() " << diagForm.GetMin() << std::endl;
              std::cout << "diagForm.GetMax() " << diagForm.GetMax() << std::endl;
            }

            Double FPSlope = localFPFactor_C*0.5*(diagForm.GetMax() + diagForm.GetMin());
            FP_Tensor.Init();
            for(UInt i = 0; i < dim_; i++){
              FP_Tensor[i][i] = FPSlope;
            }
          }
        }
        FPMaterialTensor_[LPMit->first] = FP_Tensor;

        if(testoutput){
          std::cout << "SetFPMaterialTensorsNEW" << std::endl;
          std::cout << "Calculated FPMaterialTensor according to input-flag " << intState << std::endl;
          std::cout << "Resulting FPMaterialTensor: " << FP_Tensor.ToString() << std::endl;
          testoutput = false;
        }

      }

      // restore FP state
      SetFPMaterialState(FPStateBAK);
    }

    void SetFPMaterialTensorsTMP(Integer intState){
      /*
       * NEW version after reading
       * "Analysis of the Convergence of the Fixed-Point Method Used for Solving Nonlinear Rotational Magnetic Field Problems"
       * - Dlala et al; 2008
       *
       * Global Fixpoint:
       *
       * nu_FP = 0.5 * ( dH/dB|min + dH/dB|max )
       *
       * where dH/dB|min / dH/dB|max are the smallest / maximal partial derivatives dH_i/dB_i over the whole possible
       * range of B, i.e. for the actual global case we have to find min and max slope in between negative and positve
       * saturation (or even beyond if polarization grows beyond saturation which can be the case if Mayergoyz is used
       * or if the anhysteretic paramter are set)
       *
       * > if polarization goes to real saturation, i.e. dB/dH grows only nu_0 = nu_max in saturation, we have to satisfy
       *    nu_0/2 < nu_FP <= nu_0
       *   in order to ensure the condition for contraction
       *    abs( 1/nu_FP dH/dB - 1 ) < 1
       * > if polarization does not saturate (due to anhysteretic parts or due to the excitation of the system staying
       *    away from saturation), nu_FP can be smaller than nu_0/2 (in fact I had some very nice convergence for
       *    Mayergoyz with anhysteretic parts for nu_0/8); however, it is risky to simply set a smaller nu_FP without
       *    knowing where the route goes
       * > in electrostatics, we have to use eps_FP which I assume should act similarily; here we have to be careful
       *    as we do not know the actual dD/dE as eps_0 is the minimum (whereas nu_0 is the maximum in magnetics)
       *
       * conclusion for setting nu_FP so small as possible and to find out dD/dE|max in electrostatics:
       *  > perform a complete run over one hysteresis operator going from N times negative saturation to N times positve
       *  saturation and determine the maximal and minimal dP/dH (dP/dE); from these values, determine maximal and minimal slopes
       *  > for this estimate, we have to actually traverse the major loop which means that the memory of the hyst-operator
       *  has to be set, i.e. perform this on a temporal hyst operator to not destroy the storage of an actual operator
       *  > perform this complete run only once!
       *
       * Lcoal Fixpoint
       *
       * nu_FP = C * ( dH/dB^n|min + dH/dB^n|max )
       *
       * where dH/dB^n is the slope at the last computed solution (last time step)
       * > factor C is a convergence factor that optimally should be 1 but as the slope is only estimated, C > 1 required
       * > the larger C, the larger the area of convergence but the slower the convergence
       * > the Authors of the mentioned paper suggest to start with a small C and perform a full (period) of simulation
       *  (they simulate periodic systems); in case of no-convergence C is increased
       *  --> this procedure seems not usable in general but one can show, that a linesearch-stepping with alpha basically
       *      results in C_effective = C/alpha, i.e. if alpha is smaller than 1, the convergence range is increased (as
       *      one would suspect from performing a linesearch)
       * > here we start at C = 2 (which seems to work quite well in practice according to the authors) and hope that
       *    linesearch will do the rest
       *
       * 21.02.2019
       * Additions and modifications:
       * a) Contraction / Convergence factors (called C above) are read from input.xml and passed to this function
       *    from solveStepHyst
       * b) Intstate = 0: set FP Tensor to default material tensor nu_0/nu_s > should not be required at all; check this case
       *      and stop
       *    Intstate = 1: global FP method (B); trace hyst operator if not already done; estimate maximal and minimal slope
       *      from the traced data; set FP Tensor for all LPM to nu_FP*Identity
       *    Intstate = 2: local FP method (B); at each LPM, evaluate Jacobian; from Jacobian retrieve min and max slopes
       *      and compute nu_FP from it; set FP_Tensor to nu_FP*Identity; min and max slope are only taken from DIAGONAL
       *      entries
       *    Intstate = 21: NEW; as state 2 but min and max entries are taken from WHOLE Jacobian
       *    Intstate = 3: NEW; global FP method (H); trace hyst operator if not already done; estimate maximal and minimal
       *      slope from the traced data; set FP Tensor for all LPM to nu_0/nu_s; use maximal and minimal slope in function
       *      retrieveInputToHystOperator
       *    Intstate = 4: NEW; local FP (H); at each LPM, evaluate Jacobian; from Jacobian retrieve min and max slopes;
       *      set FP_Tensor to default material tensor; use slopes later in retrieveInputToHystOperator; take entries only
       *      from DIAGONAL of Jacobian
       *    Intstate = 41: NEW; local FP (H); at each LPM, evaluate Jacobian; from Jacobian retrieve min and max slopes;
       *      set FP_Tensor to default material tensor; use slopes later in retrieveInputToHystOperator; take entries only
       *      from WHOLE Jacobian
       *
       *
       */
      Matrix<Double> Identity = Matrix<Double>(dim_,dim_);
      for(UInt i = 0; i < dim_; i++){
        for(UInt j = 0; j < dim_; j++){
          if(i == j){
            Identity[i][j] = 1.0;
          } else {
            Identity[i][j] = 0.0;
          }
        }
      }
      Matrix<Double> diagForm = Matrix<Double>(dim_,dim_);
      Matrix<Double> FP_Tensor = Matrix<Double>(dim_,dim_);
      Matrix<Double> Jacobian = Matrix<Double>(dim_,dim_);
      Matrix<Double> JacobianAbs = Matrix<Double>(dim_,dim_);
      Matrix<Double> baseTensor = Matrix<Double>(dim_,dim_);
      Matrix<Double> baseTensorInverse = Matrix<Double>(dim_,dim_);
      Matrix<Double> baseTensorUsed = Matrix<Double>(dim_,dim_);

      bool testoutput = false;

      if(testoutput){
        std::cout << "SetFPMaterialTensorsTMP with flag " << intState << " called" << std::endl;
      }

      if( (intState == 1) || (intState == 3) ){
        // global B or global H both require a traced hyst operator
        if(!hystModelTraced_){
          UInt baseSteps = 100;
          // get some additinal info if we are already traverse loop
          Double negCoercivity = 0.0;
          Double maxPolarization = 0.0;
          TraceHystOperator(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
          hystModelTraced_ = true;
        }
      }

      UInt FPStateBAK = GetFPMaterialState();
      // temporarily disable FPMaterial state, such that GetTensor does not return the fixpoint tensor but
      // eps_S / nu_S (can be computed from eps_T or nu_T)
      SetFPMaterialState(0);

      bool Hversion = false;
      if((intState == 3)||(intState == 4)||(intState == 41)){
        Hversion = true;
      }
      
      bool outputOnce = false;
      
      // eps or mu might depend on solution, i.e. we have to trace each local point separately (unfortunately!)
      std::map<UInt, LocPointMapped >::iterator LPMit;
      std::map<UInt, LocPointMapped >::iterator LPMitStart = allLPMmap_.begin();
      std::map<UInt, LocPointMapped >::iterator LPMitEnd = allLPMmap_.end();

      std::string requiredTensorName;
      std::string requiredTensorNameJacobian;
      bool invertTrace = false;
      bool invertBaseTensor = false;
      bool invertSlope = false;
      bool evalJacobian = false;

      if((intState == 2)||(intState == 21)||
          (intState == 4)||(intState == 41)){
        evalJacobian = true;
      }

      if(PDEName_ == "Electromagnetics"){
        requiredTensorName = "Reluctivity";
        if(Hversion){
          // we take dP/dH from estimation of hyst operator and add mu_0 or nu_S^-1 to get dB/dH
          invertTrace = false; // we take dP/dH directly
          invertBaseTensor = true; // we have to invert nu_0/nu_S to get mu
          invertSlope = false; // take dB/dH
          requiredTensorNameJacobian = "Permeability";
        } else {
          // we take dP/dH from estimation of hyst operator and add mu_0 or nu_S^-1 to get dB/dH
          // then we invert the sum
          invertTrace = false; // take dP/dH first
          invertBaseTensor = true; // we have to invert nu_0/nu_S to get mu
          invertSlope = true; // get dH/dB
          requiredTensorNameJacobian = "Reluctivity";
        }
      } else {
        if(Hversion){
          EXCEPTION("HVersion not implemented for electrostatics/piezos");
        }
        requiredTensorName = "Permittivity";
        invertTrace = false; // we take dP/dE directly
        invertBaseTensor = false; // we take eps directly
        invertSlope = false; // take dD/dE
      }

      UInt operatorIdx, storageIdx, idxElem;
      LocPointMapped actualLPM;
      bool onBoundary;
      
      std::map<UInt, Matrix<Double> > jacobianEvaluatedAtElement;
      std::map<UInt, Matrix<Double> >::iterator jacobianEvaluatedAtElement_it;
      
      // trace all LPM
      for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
        /*
         * check for boundary elements; on these surface elements we cannot eval the material relation and this
         * should not be required either
         */
        onBoundary = PreprocessLPM(LPMit->second, actualLPM, operatorIdx, storageIdx, idxElem);
        if(onBoundary){
          continue;
        }

        // last flag important -> no delta mat may be returned, so set it to true!
        // due to setting FPMaterialState to 0 we should get eps_S or nu_S here!
        ComputeTensor(baseTensor, LPMit->second, requiredTensorName, "DEFAULT", false, false, false, true );

        if(invertBaseTensor){
          baseTensor.Invert(baseTensorInverse);
          baseTensorUsed = baseTensorInverse;
        } else {
          baseTensorUsed = baseTensor;
        }

        if(testoutput){
          std::cout << "baseTensor " << baseTensor.ToString() << std::endl;
          std::cout << "baseTensorInverse " << baseTensorInverse.ToString() << std::endl;
          std::cout << "baseTensorUsed " << baseTensorUsed.ToString() << std::endl;
        }

        if(evalJacobian){
          if(evalJacAtMidpointOnly_){
//            std::cout << "evalJacAtMidpointOnly" << std::endl;
            // check if Jacobian was already evaluated for the current element; if so > reload Jacobian
            jacobianEvaluatedAtElement_it = jacobianEvaluatedAtElement.find(idxElem);
            if(jacobianEvaluatedAtElement_it != jacobianEvaluatedAtElement.end()){
              Jacobian = jacobianEvaluatedAtElement_it->second;
//              std::cout << "Jacobian found for element "<<idxElem<<": "<<Jacobian.ToString() << std::endl;
            } else {
              // enforce midpoint, update LPM
              int RUN_evaluationPurpose_BAK = RUN_evaluationPurpose_;
              // by setting the evaluationPurpose globally, we ensure that all functions that preprocess
              // the LPM and which are not manually forced to the midpoint will nevertheless only return the midpoint
              // due to manual enforcing, however, this step may be superfluous
              RUN_evaluationPurpose_ = 4; // output; only at midpoint
              PreprocessLPM(LPMit->second, actualLPM, operatorIdx, storageIdx, idxElem);
              
              // evaluate and store Jacobian at midpoint
              // note: here we are not required to pass a flag forceMidpoint=true as GetMaterialRelation checks the
              // flag evalJacAtMidpointOnly_ by itself
              Jacobian = GetMaterialRelation(LPMit->second, 0, 3, baseTensorUsed, requiredTensorNameJacobian,true);
              jacobianEvaluatedAtElement[idxElem] = Jacobian;
//              std::cout << "Computed and stored Jacobian for element "<<idxElem<<": "<<Jacobian.ToString() << std::endl;
              
              // revert LPM
              RUN_evaluationPurpose_ = RUN_evaluationPurpose_BAK;
              PreprocessLPM(LPMit->second, actualLPM, operatorIdx, storageIdx, idxElem);
            }
          } else {
            // new; set parameter timelevel_old to 3 forces jacobian; do not use storage values here but force new evaluation
            // (Note: the value would not be saved either, but this is not imporatant as the deltaMat storage is not required here)
            Jacobian = GetMaterialRelation(LPMit->second, 0, 3, baseTensorUsed, requiredTensorNameJacobian,true);
          }
          if(testoutput){
            std::cout << "Resulting Jacobian: " << Jacobian.ToString() << std::endl;
          }
        }

//        std::cout << "SET FP TENSOR WITH FLAG " << intState << std::endl;
        
        if(intState == 0){
          WARN("SetFPTensor with flag 0 should not be used");
          // baseTensor is the normal tensor for the system (eps for electrostatics/piezo; nu for magnetics/magnetostriction)
          FP_Tensor = baseTensor;
        } else if( (intState == 1) || (intState == 3) ){
          /*
           * Global H or B Version
           * > get base tensor
           * > get dP/dH from traced hyst operator
           * > compute nu_FP or eps_FP
           * B-Version: set FP Tensor to nu_FP*Identity (or eps_FP*Identity)
           * H-Version: set FP Tensor to baseTensor (nu), set muFPH_ to mu_FP
           */
          if(testoutput){
            std::cout << std::setprecision(16) << "Using globalFPFactor_C = " << globalFPFactor_C << std::endl;
          }

          Double slope1,slope2,FPSlope,FPSlopeComp;
          slope1 = minSlopeGlobal_;
          slope2 = maxSlopeGlobal_;
          if(invertTrace){
            Double tmp;
            // switch position as min and max change
            tmp = 1.0/slope1;
            slope1 = 1.0/slope2;
            slope2 = tmp;
          }

          baseTensorUsed.GetDiagInMatrix(diagForm);
          if(testoutput){
            std::cout << "diagForm: " << diagForm.ToString() << std::endl;
          }

          slope1 += diagForm.GetMin();
          slope2 += diagForm.GetMax();

          if(invertSlope){
            slope1 = 1.0/slope1;
            slope2 = 1.0/slope2;
          }

          // should also work for H-version with globalFPFactor_C = 1 but apparently it does not > maybe slope estimate is wrong?!
          UInt restrictedPrecision = 9;
          Double safetyFactor = 1.00001;
          FPSlopeComp = safetyFactor*globalFPFactor_C*0.5*(slope1 + slope2);
          FPSlope = restrictPrecision(FPSlopeComp, restrictedPrecision); // restrict precision to hopefully increase stability
          
          if(outputOnce){
            std::cout << std::scientific << std::setprecision(12) << "###Global globalFPFactor_C and additional safetyFactor: " << globalFPFactor_C << " / " << safetyFactor << std::endl;
            std::cout << std::scientific << std::setprecision(12) << "###Global FP Slope (computed/truncated)### " << FPSlopeComp << " / " << FPSlope << std::endl;
            outputOnce = false;
          }
          FP_Tensor.Init();

          // nu_max/2 < nu_FP <= nu_max
          // > maybe should be checked
//          std::cout << "baseTensorUsed (mu/nu/eps) = " << baseTensorUsed.ToString() << std::endl;
//          std::cout << "FPSlope = " << FPSlope << std::endl;

          if(intState == 1){
            // B-version
            FP_Tensor.Add(FPSlope,Identity);
          } else {
            // H-version
            // > use baseTensor for system
            FP_Tensor.Add(1.0,baseTensor);
            // set muFP to FPSlope
            muFPH_[LPMit->first] = FPSlope;

            if(testoutput){
              std::cout << "GLOBAL FP H - Set FP-Tensor to baseTensor; set muFPH to " << FPSlope << std::endl;
            }

          }

        } else if((intState == 2)||(intState == 21) || (intState == 4)||(intState == 41)){
          /*
           * Local H or B Version
           * > get dD/dE, dH/dB or dB/dH from Jacobian (yes, Jacobian contains already the correct base tensor)
           * > compute nu_FP or eps_FP
           * B-Version: set FP Tensor to nu_FP*Identity (or eps_FP*Identity)
           * H-Version: set FP Tensor to baseTensor (nu), set muFPH_ to mu_FP
           */
          if(testoutput){
            std::cout << "Using localFPFactor_C = " << localFPFactor_C << std::endl;
          }

          Double slope1,slope2,FPSlope;
          if( (intState == 2) || (intState == 4) ){
            Jacobian.GetDiagInMatrix(diagForm);
            // compute FP factor only from diagonal entries
            slope1 = diagForm.GetMin();
            slope2 = diagForm.GetMax();
          } else {
            Jacobian.GetAbsValues(JacobianAbs);
            // use min/max of whole Jacobian; build abs values first as side diagonal elements mostly are negative
            slope1 = JacobianAbs.GetMin();
            slope2 = JacobianAbs.GetMax();
          }

          if(testoutput){
            std::cout << "idx="<<LPMit->first<<"-slope1=" << slope1 <<"-slope2="<<slope2<< std::endl;
            std::cout << "Jacobian = " << Jacobian.ToString() << std::endl;
          }

          // not required as Jacobian is already the correct tensor (nu or mu)
//          if(invertSlope){
//            if(testoutput){
//              std::cout << "Invert slope" << std::endl;
//            }
//            slope1 = 1.0/slope1;
//            slope2 = 1.0/slope2;
//          }

          FPSlope = localFPFactor_C*0.5*(slope1 + slope2);
          FP_Tensor.Init();

          if( (intState == 2) || (intState == 21) ){
            // B-version
            FP_Tensor.Add(FPSlope,Identity);
          } else {
            // H-version
            // > use baseTensor for system
            FP_Tensor.Add(1.0,baseTensor);
            // set muFP to FPSlope
            muFPH_[LPMit->first] = FPSlope;
            if(testoutput){
              std::cout << "LOCAL FP H - Set FP-Tensor to baseTensor; set muFPH to " << FPSlope << std::endl;
            }

          }

        } else {
          EXCEPTION("FP Tensor unknown");
        }

        FPMaterialTensor_[LPMit->first] = FP_Tensor;

        if(testoutput){
          std::cout << "SetFPMaterialTensors" << std::endl;
          std::cout << "Calculated FPMaterialTensor according to input-flag " << intState << std::endl;
          std::cout << "Resulting FPMaterialTensor: " << FP_Tensor.ToString() << std::endl;
          testoutput = false;
        }
      }

      // restore FP state
      SetFPMaterialState(FPStateBAK);
    }

    void SetFPMaterialTensors(Integer intState){
      EXCEPTION("old version of SetFPMaterialTensors is deprecated");
//
//      Matrix<Double> FP_Tensor = Matrix<Double>(dim_,dim_);
//      if(intState == 0){
//        FP_Tensor = rev_mat_fac_MATRIX_;
//      } else if(intState == 1){
//        // global fp approach;
//        // magnetics (B-version): take (nu_min + nu_max)/2 \approx nu_max/2 = nu_0/2 = rev_mat_fac_/2
//        // electrostatics (E-version): take (eps_min + eps_max)/2 but the problem is, that we do not know eps_max
//        // TODO: find way to compute estimate global eps_FP
//        if(material_->GetMaterialDatabaseName() == "Electrostatic"){
//          EXCEPTION("Global FP approach for electrostatics not implemented yet.");
//        } else {
//          FP_Tensor.Init();
//          Double factor = 0.0625;
//          std::cout << "Used scaling factor for GLOBAL FP: " << factor << std::endl;
//          FP_Tensor.Add(factor,rev_mat_fac_MATRIX_);
//        }
//      } else if(intState == 2){
//        // approximate nu/eps from Jacobian at each storage location
//      } else {
//        EXCEPTION("FPMaterial flag has to be 0,1 or 2.");
//      }
//
//      std::map<UInt, LocPointMapped >::iterator LPMit;
//      std::map<UInt, LocPointMapped >::iterator LPMitStart = allLPMmap_.begin();
//      std::map<UInt, LocPointMapped >::iterator LPMitEnd = allLPMmap_.end();
//
//      bool testoutput = true;
//      if(intState < 2){
//        // iterate over all LPM and set FP_Tensor accordingly
//        for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
//          FPMaterialTensor_[LPMit->first] = FP_Tensor;
//
//          if(testoutput){
//            std::cout << "SetFPMaterialTensors" << std::endl;
//            std::cout << "Calculated FPMaterialTensor according to input-flag " << intState << std::endl;
//            std::cout << "Resulting FPMaterialTensor: " << FP_Tensor.ToString() << std::endl;
//            testoutput = false;
//          }
//
////          if(LPMit->first == 1){
////            std::cout << "Global estimate for FP tensor" << std::endl;
////            std::cout << "FPMaterialTensor_["<<LPMit->first<<"] = " << FPMaterialTensor_[LPMit->first].ToString() << std::endl;
////          }
//
//        }
//      } else {
//        Matrix<Double> tmpTensor = Matrix<Double>(dim_,dim_);
//        // according to the source
//        // "Locally Convergent Fixed-Point Method for Solving Time-Stepping Nonlinear Field Problems" - E. Dlala et al. 2007
//        // the slope estimation has to be scaled by a convergenceFactor that is strictly larger than 1
//        // the actual value has/can be retrieved during iteration (similar to linesearch); here just set 2 for testing
//        Double convergenceFactor = 8.0;
//        std::cout << "Used convergence factor for LOCAL FP: " << convergenceFactor << std::endl;
//        Double  sign = 1.0;
//
//        Matrix<Double> rev_mat_fac_MATRIX_inverted = Matrix<Double>(dim_,dim_);
//        rev_mat_fac_MATRIX_.Invert(rev_mat_fac_MATRIX_inverted);
//
//        // iterate over all LPM, estimate Jacobian and set the value
//        for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
//          // note: estimateSlope_Jacobian returns dJ/dB in magnetics and dP/dE in electrostatics
//          // in magnetics, we need dH/dB and in electrostatics dD/dE
//          // > magnetics:
//          //    compute nu_0*(I - dJ/dB) > ist doch mist! besser invert(mu + dB/dH)
//          // > electrostatics:
//          //    compute eps_ + dP/dE
//          //
//          // additional note / TODO:
//          // for local FP we do not actually use the Jacobian but only its trace (somehow) > change and test
//
//          if(material_->GetMaterialDatabaseName() == "Electrostatic"){
//            tmpTensor = EstimateSlope_Jacobian(LPMit->second, false, false, "DEFAULT");
//            FP_Tensor = tmpTensor;
//            FP_Tensor.Add(1.0,rev_mat_fac_MATRIX_);
//          } else {
//            // new; set timelevel old to 3 forces jacobian; do not use storage values here but force new evaluation
//            // (Note: the value would not be saved either, but this is not imporatant as the deltaMat storage is not required here)
//            FP_Tensor = GetMaterialRelation(LPMit->second, 0, 3, rev_mat_fac_MATRIX_inverted, "Reluctivity",true);
//            sign = 1.0;
//            // old
////            tmpTensor.Mult(rev_mat_fac_MATRIX_,FP_Tensor);
////            FP_Tensor.Add(-1.0,rev_mat_fac_MATRIX_);
////            sign = -1.0;
//          }
//          FPMaterialTensor_[LPMit->first].Init();
//          FPMaterialTensor_[LPMit->first].Add(sign*convergenceFactor,FP_Tensor);
//
//          if(testoutput){
//            std::cout << "SetFPMaterialTensors" << std::endl;
//            std::cout << "Calculated FPMaterialTensor according to input-flag " << intState << std::endl;
//            std::cout << "Resulting FPMaterialTensor: " << FPMaterialTensor_[LPMit->first].ToString() << std::endl;
//            testoutput = false;
//          }
//
////          if(LPMit->first == 1){
////            std::cout << "Local estimate for FP tensor" << std::endl;
////            std::cout << "FPMaterialTensor_["<<LPMit->first<<"] = " << FPMaterialTensor_[LPMit->first].ToString() << std::endl;
////          }
//        }
//      }
    }

    UInt GetFPMaterialState(){
      return FPMaterialTensorSet_;
    }

     void SetFPMaterialState(UInt newState){
      FPMaterialTensorSet_ = newState;
    }

    UInt GetFPMaterialStateRHS(){
      return FPMaterialRHSSet_;
    }

    Matrix<Double> GetFPMaterialTensor(const LocPointMapped& OriginalLPM){
      UInt operatorIdx, storageIdx, idxElem;
      LocPointMapped actualLPM;
      PreprocessLPM(OriginalLPM, actualLPM, operatorIdx, storageIdx, idxElem);
      return FPMaterialTensor_[storageIdx];
    }

    Vector<Double> GetFPCorrectionVector(const LocPointMapped& OriginalLPM, Integer timeLevel){
      Vector<Double> returnVec = Vector<Double>(dim_);
      if(FPMaterialRHSSet_ == 0){
        // rev_mat_fac_MATRIX_ _ rev_mat_fac_MATRIX_ = zeroMat
        returnVec.Init();
      } else {
        EXCEPTION("DEPRECATED; USE ONLY STANDARD RESIDUAL ON RHS; NO CORRECTION REQUIRED");
//        UInt operatorIdx, storageIdx;
//        LocPointMapped actualLPM;
//
//        bool onBoundary = PreprocessLPM(OriginalLPM, actualLPM, operatorIdx, storageIdx);
//        Vector<Double> E_B = RetrieveLPMSolution(actualLPM, storageIdx, timeLevel, onBoundary);
//        // note: correctionMatrix = eps_FP - eps_0 > add
//        Matrix<Double> correctionTensor = Matrix<Double>(dim_,dim_);
//
//        // TODO: muss geänder werden im Falle von Piezo/Magstrict
//        // here we use not eps_0 but eps_S or eps_T
//        // eps_S would be the easier case as it is assumed constant
//        // eps_T however depends on P, so we should use the current value of eps_T
//        // only last flag relevant, locking the deltaMatrix computation
//        if(PDEName_ == "Electromagnetics"){
//          ComputeTensor(correctionTensor, OriginalLPM, "Reluctivity", "DEFAULT", false, false, false, true );
//        } else {
//          ComputeTensor(correctionTensor, OriginalLPM, "Permittivty", "DEFAULT", false, false, false, true );
//        }
////        correctionTensor = rev_mat_fac_MATRIX_;
//        correctionTensor.Add(-1.0,FPMaterialTensor_[storageIdx]);
//
//        correctionTensor.Mult(E_B,returnVec);
      }
      return returnVec;
    }

    void PrecomputePolarization(bool overwriteMemory);
    void RestorePrevHystVals(UInt storageLocation);
    void SetPreviousHystVals(UInt storageLocation);
    void SetPreviousHystValsOLD(bool setLastTS, bool forceMemoryLock = false);

    void EvaluateRotDirectionForVectorExtension();

    void GetScalar(Double& outputScalar, const LocPointMapped& lpm );

    void GetVector( Vector<Double>& outputVector,const LocPointMapped& lpm);

    void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm );

    void checkSaturationStateAllElements(Double& lastTSSatAvg, Double& lastItSatAvg, Double& curItSatAvg,
    Double& oppositeDirAsTSAvg, Double& oppositeDirAsItAvg){

      // go over all elements and check for saturation;
      // Note: here we return Double values in order to give a percentage over the number of saturated elements
      UInt lastTSSat, lastItSat, curItSat, oppositeDirAsTS, oppositeDirAsIt;
      UInt lastTSSatAvgInt, lastItSatAvgInt, curItSatAvgInt, oppositeDirAsTSAvgInt, oppositeDirAsItAvgInt;
      lastTSSatAvgInt = 0;
      lastItSatAvgInt = 0;
      curItSatAvgInt = 0;
      oppositeDirAsTSAvgInt = 0;
      oppositeDirAsItAvgInt = 0;
      for(UInt i = 0; i < numElemSD_; i++){
        checkSaturationState(i, lastTSSat, lastItSat, curItSat, oppositeDirAsTS, oppositeDirAsIt);
        lastTSSatAvgInt += lastTSSat;
        lastItSatAvgInt += lastItSat;
        curItSatAvgInt += curItSat;
        oppositeDirAsTSAvgInt += oppositeDirAsTS;
        oppositeDirAsItAvgInt += oppositeDirAsIt;
      }
      lastTSSatAvg = ((Double)lastTSSatAvgInt)/((Double)numElemSD_);
      lastItSatAvg = ((Double)lastItSatAvgInt)/((Double)numElemSD_);
      curItSatAvg = ((Double)curItSatAvgInt)/((Double)numElemSD_);
      oppositeDirAsTSAvg = ((Double)oppositeDirAsTSAvgInt)/((Double)numElemSD_);
      oppositeDirAsItAvg = ((Double)oppositeDirAsItAvgInt)/((Double)numElemSD_);
    }

    void checkSaturationState(UInt operatorIdx, UInt& lastTSSat, UInt& lastItSat, UInt& curItSat,
    UInt& oppositeDirAsTS, UInt& oppositeDirAsIt){

      /*
       * Little helper function checking if element solution was/is in saturation
       * -> checks stored array values at storageIdx = operatorIdx+0 (i.e. midpoint of element)
       * -> return integer instead of bools to allow simple addition (e.g. when summing over all
       *      elements)
       */
      // for magnetics, we check for B, i.e. we have to compare with ySat + mu*xSat (ySat is saturation for P!)
      // for electrostatics, we check for E, i.e. we have to compare with xSat
      // > both are stored in the same array, however
      Double satValue;
      if(PDEName_ == "Electromagnetics"){
        // will not  be exactly correct when taking the norm but gives beter estimate than skipping it
        satValue = POL_operatorParams_.outputSat_ + eps_nu_base_.NormL2()*POL_operatorParams_.inputSat_;
      } else {
        satValue = POL_operatorParams_.inputSat_;
      }

      if(E_B_lastTS_[operatorIdx].NormL2() >= satValue){
        lastTSSat = 1;
      } else {
        lastTSSat = 0;
      }
      if(E_B_lastIt_[operatorIdx].NormL2() >= satValue){
        lastItSat = 1;
      } else {
        lastItSat = 0;
      }
      if(E_B_[operatorIdx].NormL2() >= satValue){
        curItSat = 1;
      } else {
        curItSat = 0;
      }

      Double dotProduct;
      E_B_lastTS_[operatorIdx].Inner(E_B_[operatorIdx],dotProduct);
      if(dotProduct < 0){
        // direction changed at least be 90 degree
        // > possible jump from pos into neg saturation
        oppositeDirAsTS = 1;
      } else {
        oppositeDirAsTS = 0;
      }

      E_B_lastIt_[operatorIdx].Inner(E_B_[operatorIdx],dotProduct);
      if(dotProduct < 0){
        // direction changed at least be 90 degree
        // > possible jump from pos into neg saturation
        oppositeDirAsIt = 1;
      } else {
        oppositeDirAsIt = 0;
      }

    }

    std::string ToString() const;

    void ActiveOneShotSlopeEstimation(Double steppingLength, Double scaling){
      ES_steppingLength_ = steppingLength;
      ES_scaling_ = scaling;
      if(storageInitialized_){
        for(UInt i = 0; i < numStorageEntries_; i++){
          takeEstimatedSlope_[i] = true;
        }
      }
    }

    Matrix<Double> EstimateSlope(const LocPointMapped& Originallpm, bool includeStrains, bool useAbs, std::string estimationType, std::string implementationVersion,
    Double steppingLength=1e-10, Double scaling=1.0, bool forceMidpoint=false);

    Matrix<Double> EstimateSlope_DeltaMat(const LocPointMapped& Originallpm, bool includeStrains, bool useAbs, std::string implementationVersion, bool forceMidpoint);

    Matrix<Double> EstimateSlope_Jacobian(const LocPointMapped& Originallpm, bool includeStrains, bool useAbs, std::string implementationVersion, bool forceMidpoint);

    Matrix<Double> CalcDeltaMat(Vector<Double> E_B_diff, Vector<Double> P_J_diff, bool useAbs,std::string implementationVersion, Double cuttingTol = -1.0);
    Matrix<Double> CalcDeltaMatStrains(Vector<Double> E_B_diff, Vector<Double> S_diff, bool useAbs,std::string implementationVersion, Double cuttingTol = -1.0);

    Matrix<Double> GetDeltaMat(const LocPointMapped& Originallpm, int timelevel1, int timelevel2, bool useStrains, bool useAbs, std::string implementationVersion);

//    bool PrecomputeSlopes(const LocPointMapped& Originallpm, std::string implementationVersion,
//          bool strainVersion, bool forceReevaluation);

    Vector<Double> GetIrreversibleStrains(const LocPointMapped& Originallpm, int timeLevel, bool forceMidpoint=false);
    Matrix<Double> ConvertFromVoigtToTensor(Vector<Double> Si_voigt);
    Matrix<Double> GetIrreversibleStrainTensor(const LocPointMapped& Originallpm, int timeLevel, bool forceMidpoint=false);
    Vector<Double> GetIrreversibleStrains(UInt storageIdx, int timeLevel);

    Vector<Double> GetPrecomputedInputToHysteresisOperator(const LocPointMapped& Originallpm, int timeLevel, bool forceMidpoint=false);

    Vector<Double> GetPrecomputedOutputOfHysteresisOperator(const LocPointMapped& Originallpm, int timeLevel, bool forStrain, bool forceMidpoint=false);

    void PrecomputeScaledAndRotatedCouplingTensor(Matrix<Double> baseTensor, bool rotate);

    Vector<Double> ComputeIrreversibleStrains(Vector<Double> P,Vector<Double> E_H,Vector<Double> Pold);

    Vector<Double> RetrieveInputToHysteresisOperator(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool onBoundary, bool overwriteMemory=false);
    Vector<Double> RetrieveInputToHysteresisOperatorOLD(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool onBoundary);

    Vector<Double> RetrieveLPMSolution(LocPointMapped& actualLPM, UInt storageIdx, int timeLevel, bool onBoundary);

    Double GetOutputSaturation(){
      return POL_operatorParams_.outputSat_;
    }

    PtrCoefFct GenerateMatCoefFnc(std::string tensorName){
      //      std::cout << "Generate mat coef function " << tensorName << std::endl;
      PtrCoefFct ret;
      if((tensorName == "Permittivity")||(tensorName == "CouplingMechToElec")||(tensorName == "CouplingElecToMech")){
        // Note: for Piezo we need two different coupling functions
        // Reason: in case of deltaMat formulation, the coupling matrices are not
        //         transposed of each other as the mech to elec term has deltaS/deltaE
        //         in addition to the rotated coupling operator
        PtrCoefFct eps = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY_TENSOR,tensorType_,Global::REAL);

        if(formerHystHelperParamsSet_ == false){
          SetParamsOfFormerHystHelper(eps,elastTensorFct_,couplTensorFct_);
          formerHystHelperParamsSet_ = true;
        }

        bool transposed = false;
        if(tensorName == "CouplingElecToMech"){
          // in mech PDE we need e^T or d^T
          transposed = true;
        }
        shared_ptr<CoefFunctionHystMatTensor> c(new CoefFunctionHystMatTensor(this->shared_from_this(),tensorName,transposed));

        ret = c;
      } else if((tensorName == "Reluctivity")||(tensorName == "LinearReluctivity")||
          (tensorName == "CouplingMechToMag")||(tensorName == "CouplingMagToMech")){
        PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY_TENSOR,tensorType_,Global::REAL);

        if(formerHystHelperParamsSet_ == false){
//          std::cout << "Create hystHelper" << std::endl;
//          std::cout << "- where? GenerateMatCoefFnc" << std::endl;

          SetParamsOfFormerHystHelper(nu,elastTensorFct_,couplTensorFct_);
          formerHystHelperParamsSet_ = true;
        }

        bool transposed = false;
        if(tensorName == "CouplingMagToMech"){
          // in mech PDE we need h^T or g^T
          transposed = true;
        }
        shared_ptr<CoefFunctionHystMatTensor> c(new CoefFunctionHystMatTensor(this->shared_from_this(),tensorName,transposed));

        ret = c;
      } else {
        EXCEPTION("tensorName not implemented yet");
      }

      return ret;
    }

    PtrCoefFct GenerateRHSCoefFnc(std::string vectorName, PtrCoefFct coefFunctionToBeIncluded){
      // currently only used for magnetic pde where we want to modify the rhs fluxExcitation
      // > might be required in a similar way for elecPDE
      // Background: if flux is prescribed, this flux should be decomposed into linear and hysteretic part
      // without this modification, the flux will be added above the hysteretic part
      PtrCoefFct ret;
      bool onBoundary = false;
      bool transposed = false;

      if(vectorName == "MagPolarization"){
        PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY_TENSOR,tensorType_,Global::REAL);

        if(formerHystHelperParamsSet_ == false){
          SetParamsOfFormerHystHelper(nu,elastTensorFct_,couplTensorFct_);
          formerHystHelperParamsSet_ = true;
        }

        shared_ptr<CoefFunctionHystRHSLoad> c(new CoefFunctionHystRHSLoad(this->shared_from_this(),vectorName,transposed,onBoundary,coefFunctionToBeIncluded));
        ret = c;
      }
      return ret;
    }


    PtrCoefFct GenerateRHSCoefFnc(std::string vectorName, bool onBoundary = false ){
      //      std::cout << "Generate rhs coef function " << vectorName << std::endl;
      PtrCoefFct ret;
      if( (vectorName == "ElecPolarization") || (vectorName == "IrrStressForMechPDE") || (vectorName == "CoupledIrrStrainForElecPDE")){
        PtrCoefFct eps = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY_TENSOR,tensorType_,Global::REAL);

        if(formerHystHelperParamsSet_ == false){
          SetParamsOfFormerHystHelper(eps,elastTensorFct_,couplTensorFct_);
          formerHystHelperParamsSet_ = true;
        }

        bool transposed = false;
        // NOTE: for mechPDE we do not need the coupling tensor when we apply the strains > no transpose needed
        // note further, that only the irr-strains are put to both rhs; polarization only in elec pde
        //        if(vectorName == "IrrStrainForMechPDE"){
        //          // in mech PDE we need h^T or g^T
        //          transposed = true;
        //        }
        shared_ptr<CoefFunctionHystRHSLoad> c(new CoefFunctionHystRHSLoad(this->shared_from_this(),vectorName,transposed,onBoundary));

        ret = c;
      } else if( (vectorName == "MagMagnetization") || (vectorName == "MagStrictLoadForMechPDE") || (vectorName == "MagStrictLoadForMagPDE")){
        PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY_TENSOR,tensorType_,Global::REAL);

        if(formerHystHelperParamsSet_ == false){
//          std::cout << "Create hystHelper" << std::endl;
//          std::cout << "- where? CoefFunctionHystRHSLoad" << std::endl;

          SetParamsOfFormerHystHelper(nu,elastTensorFct_,couplTensorFct_);
          formerHystHelperParamsSet_ = true;
        }

        bool transposed = false;
        if(vectorName == "MagStrictLoadForMechPDE"){
          // in mech PDE we need h^T or g^T
          transposed = true;
        }
        shared_ptr<CoefFunctionHystRHSLoad> c(new CoefFunctionHystRHSLoad(this->shared_from_this(),vectorName,transposed,onBoundary));

        ret = c;
      } else {
        //std::cout << "RHSCoefFunction for " << vectorName << " requested." << std::endl;
        EXCEPTION("vectorName not implemented yet");
      }

      return ret;
    }

    PtrCoefFct GenerateOutputCoefFnc(std::string ResultName){
      //      std::cout << "Generate output coef function " << ResultName << std::endl;
      PtrCoefFct ret;

      if( (ResultName == "ElecPolarization") || (ResultName == "IrrStressesPiezo_TensorForm")
        || (ResultName == "IrrStressesPiezo_VectorForm")
        || (ResultName == "IrrStrainsPiezo_TensorForm") || (ResultName == "CoupledIrrStrainsPiezo") ){
        PtrCoefFct eps = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY_TENSOR,tensorType_,Global::REAL);

        if(formerHystHelperParamsSet_ == false){
          SetParamsOfFormerHystHelper(eps,elastTensorFct_,couplTensorFct_);
          formerHystHelperParamsSet_ = true;
        }

        shared_ptr<CoefFunctionHystOutput> c(new CoefFunctionHystOutput(this->shared_from_this(),ResultName));
        ret = c;
      } else if ( (ResultName == "MagPolarization") || (ResultName == "MagMagnetization") || (ResultName == "MagFieldIntensityHyst") ){
        PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY_TENSOR,tensorType_,Global::REAL);

        if(formerHystHelperParamsSet_ == false){
//          std::cout << "Create hystHelper" << std::endl;
//          std::cout << "- where? GenerateOutputCoefFnc" << std::endl;

          SetParamsOfFormerHystHelper(nu,elastTensorFct_,couplTensorFct_);
          formerHystHelperParamsSet_ = true;
        }

        shared_ptr<CoefFunctionHystOutput> c(new CoefFunctionHystOutput(this->shared_from_this(),ResultName));
        ret = c;

      } else {
        EXCEPTION("Result not implemented");
      }
      return ret;
    }

    std::string getPDEName(){
      return PDEName_;
    }

    bool useStrainForm(){
      return CouplingParams_.useStrainForm_;
    }

    int GetStrainForm(){
      return CouplingParams_.strainForm_;
    }

    void SetStrainForm(int strainForm){
      CouplingParams_.strainForm_ = strainForm;
    }

    bool deltaMatActive(){
      return deltaMatActive_;
    }

    int GetDeltaForm(){
      return RUN_deltaForm_;
    }

    int GetTimeLevel(std::string Type){
      if(forceCurrentTS_){
        return 0;
      }
      // keep old names until all new ones are implemented
      if(Type == "Material"){
        return timeLevel_Mat_;
      } else if (Type == "DeltaMat"){
        return timeLevel_deltaMatPol_;
      } else if (Type == "DeltaMatPol"){
        return timeLevel_deltaMatPol_;
      } else if (Type == "DeltaMatStrain"){
        return timeLevel_deltaMatStrain_;
      } else if (Type == "DeltaMatCoupling"){
        return timeLevel_deltaMatCoupling_;
      } else if (Type == "RHS"){
        return timeLevel_rhsPol_;
      } else if (Type == "RHSPol"){
        return timeLevel_rhsPol_;
      } else if (Type == "RHSStrain"){
        return timeLevel_rhsStrain_;
      } else if (Type == "RHSCoupling"){
        return timeLevel_rhsCoupling_;
      } else if (Type == "BC"){
        return timeLevel_boundaryTerms_;
      } else if (Type == "Output"){
        return timeLevel_results_;
      } else {
        EXCEPTION("Type not known");
        return -2;
      }
    }
    //
    //    Integer timeLevel_Mat_;
    //    Integer timeLevel_rhsPol_;
    //    Integer timeLevel_rhsStrain_; // piezo / magstrict only
    //    Integer timeLevel_rhsCoupling_; // piezo / magstrict only
    //    Integer timeLevel_deltaMatPol_;
    //    Integer timeLevel_deltaMatStrain_;
    //    Integer timeLevel_deltaMatCoupling_;
    //    Integer timeLevel_results_;
    //    Integer timeLevel_boundaryTerms_;

    void TestJacobianApproximations();
    void TestInversion(PtrParamNode testNode, PtrParamNode infoNode);

    void CreatePeriodicTestSignal(std::string name, Double amplitudeScaling, Double numPeriods, UInt stepsPerPeriods, Vector<Double>& xVals, Vector<Double>& yVals, Double initialAmplitude = std::numeric_limits<double>::infinity());

    void CreateNonPeriodicTestSignal(std::string name, Double amplitudeScaling, UInt numberOfSteps, Vector<Double>& xVals, Vector<Double>& yVals);

    void TestHystOperatorWithSignal(std::string name, Vector<Double> xVals, Vector<Double> yVals,
    Vector<Double> zVals, UInt dimHystOperator, int forcedPreisachResolution,
    bool testInversion, bool printStatistics, bool writeResultsToFile, bool measurePerformance,
    std::string commonPerfFile, bool test1D, bool outputIrrStrains, std::string nameTagForPerfFile,
    std::string additionalTag1,std::string additionalTag2,std::string additionalTag3);

    void WriteSignalToFile(std::string name, Vector<Double> xVals, Vector<Double> yVals, Vector<Double> zVals);
    void WriteTestStepToInfoXML(const UInt testNum,
          const Double Px, const Double Py, const Double Pz,
          const Double Pabs, const Double Pangle1, const Double Pangle2);
    
    void SetElastAndCouplTensor(PtrCoefFct elastTensor, PtrCoefFct couplTensor){
      //std::cout << "Elast and Coupl tensor were passed by coupled pde!" << std::endl;
      elastTensorFct_ = elastTensor;
      couplTensorFct_ = couplTensor;
      COUPLED_inXMLFile_ = true;
    }

    bool couplingTensorSet(){
      if (couplTensorFct_ == NULL){
        return false;
      } else {
        return true;
      }
    }

    /*
     * Helper function that allows to add additional SD lists to the hyst operator
     * > needed to add storage space for surface elements
     * > has to be executed before storage was initialized!
     */
    void AddAdditionalSDList(shared_ptr<EntityList> actSDList, RegionIdType volReg, bool isSurface);

    /*
     * Helper function; creates mapping between operators and their corresponding
     * storage entries;
     * for standard evaluation:
     *    operatorIdx = storageIdx; evaluation only at midpoint
     *  > each entry of map has form
     *      N, <N, true>
     *  for extended evaluation:
     *    each operator corresponds to an element; operator gets evaluated at
     *    each integration point inside the element + at midpoint
     *  > entries of map
     *      M, <n1, false>
     *      M, <n2, true>
     *      M, <n3, false>
     *      ...
     *  for full evaluation:
     *    each integration point + midpoint get own hyst operator
     *    operatorIdx = storageIdx
     *  > entries of map
     *      m1, <m1, false>
     *      m2, <m2, true>
     *      m3, <m3, false>
     *      ...
     */
    void MapStorageIndexToOperator(UInt operatorIdx, UInt storageIdx, bool isMidpoint){

      if(operatorToStorage_.find(operatorIdx) == operatorToStorage_.end()){
        std::map<UInt, bool> innerMap;
        innerMap[storageIdx] = isMidpoint;
        operatorToStorage_[operatorIdx] = innerMap;
      } else {
        operatorToStorage_[operatorIdx][storageIdx] = isMidpoint;
      }
    }

    bool anyMatrixForLocalInversionRequiresComputation(){
      for(UInt i = 0; i < numStorageEntries_; i++){
        if(!eps_mu_local_Set_[i]){
          return true;
        }
      }
      return false;
    }

    void getMatrixForLocalInversion(const LocPointMapped& Originallpm, Matrix<Double>& matrixForInversion, Matrix<Double>& matrixForInversionInverted){
      UInt operatorIdx, storageIdx, idxElem;
      LocPointMapped actualLPM;

      PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx, idxElem);
      matrixForInversion = eps_mu_local_[storageIdx];
      matrixForInversionInverted = epsInv_nu_local_[storageIdx];
    }

    void getLPMMaps(std::map<UInt, LocPointMapped >& allLPM, std::map<UInt, LocPointMapped >& midpointLPM){
      allLPM = allLPMmap_;
      midpointLPM = midpointLPMmap_;
    }

    void setMatrixForLocalInversion(Matrix<Double> matrixForInversion, Matrix<Double> matrixForInversionInverted, UInt storageIdx, bool reuse){
      if(reuse){
        eps_mu_local_Set_[storageIdx] = true;
      } else {
        // eps_mu_local_Set_ shall be set only once in each iteration
        // > use flag which gets reset once during setprevioushystvals function
        if(eps_mu_local_Set_[storageIdx] == false){
          eps_mu_local_[storageIdx] = matrixForInversion;
          epsInv_nu_local_[storageIdx] = matrixForInversionInverted;
          eps_mu_local_Set_[storageIdx] = true;
        }
      }
    }

    void ClipDirection(Vector<Double>& targetVector);

    void GetScaledAndRotatedCouplingTensor(const LocPointMapped& lpm, Matrix<Double>& couplTensor, Matrix<Double>& rotatedCouplTensor,int timeLevel,
    bool rotate = true){
      // compare to Kaltenbacher "Numerical Simulation ..." 3rd Edition p. 387
      //
      // scaling:
      // [e(P)] = P^i_k/P_sat[e]
      //
      // for vector model, an additional rotation towards P could be reasonable

      // get actual indices for storage array to find out if we need to reevaluate or not
      UInt operatorIdx, storageIdx, idxElem;
      LocPointMapped actualLPM;

      PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx, idxElem);

      rotatedCouplTensor = rotatedCouplingTensor_[storageIdx];
      //
      //      if( (rotatedCouplingTensor_requiresReeval_[storageIdx] == false) && (timeLevel == lastUsedTimeLevelForRotation_[storageIdx]) ){
      //        //      std::cout << "Coupling Tensor already rotated > no reeval performed" << std::endl;
      //        rotatedCouplTensor = rotatedCouplingTensor_[storageIdx];
      //      } else {
      //
      //        // obtain coupling tensor first
      //        UInt numRows, numCols;
      //        numCols = couplTensor.GetNumCols();
      //        numRows = couplTensor.GetNumRows();
      //        //couplTensorCoefFnc->GetTensorSize(numRows,numCols);
      //        //Matrix<Double> couplTensor = Matrix<Double>(numRows,numCols);
      //        //couplTensorCoefFnc->GetTensor(couplTensor,actualLPM);
      //
      //        //      std::cout << "Retrieved coupling tensor: " << couplTensor.ToString() << std::endl;
      //        //      std::cout << "Rotate coupling tensor" << std::endl;
      //        rotatedCouplTensor.Resize(numRows,numCols);
      //        Matrix<Double> scaledCouplTensor = Matrix<Double>(numRows,numCols);
      //
      //        // get current polarization (elec or mag)
      //        Vector<Double> P = GetPrecomputedOutputOfHysteresisOperator(lpm, timeLevel);
      //
      //        //      std::cout << "Current polarization vector " << P.ToString() << std::endl;
      //        //
      //        // calculate scaling
      //        assert(POL_operatorParams_.outputSat_ != 0);
      //        Double scaling = P.NormL2()/POL_operatorParams_.outputSat_;
      //
      //        scaledCouplTensor = couplTensor* scaling;
      //
      //        //      std::cout << "Scalaed coupling tensor " << scaledCouplTensor.ToString() << std::endl;
      //        //
      //        if(rotate == true){
      //          if(P.NormL2() == 0){
      //            //        std::cout << "Polarization is zero > perform no rotation" << std::endl;
      //            rotatedCouplTensor = scaledCouplTensor;
      //
      //          } else {
      //            Vector<Double> dirP = Vector<Double>(P.GetSize());
      //            dirP = P / P.NormL2();
      //
      //            // calculate rotation matrix
      //            if(numCols == 4){
      //              // axi case: 2x4 matrix
      //              WARN("Rotation for axi case not implemented yet. No rotation was peformed");
      //              rotatedCouplTensor = scaledCouplTensor;
      //
      //            } else if (numCols == 6){
      //              // 3d plane case
      //              // Idea: Create rotation matrix, that maps the current direction of P onto z-axis
      //              //       as the z-axis is the default polarization axis in the mat file
      //              //       To obtained the desired behavior, we have to rotate the z-axis onto P however.
      //              //       This can be done by taking the transposed rotation matrix.
      //
      //              Double alpha, beta, gamma;
      //
      //              // 1. rotate around z-axis by angle gamma such that P lies in z-y plane
      //              // gamma = angle between z-y plane and dirP
      //              gamma = std::atan2(dirP[0],dirP[1]);
      //
      //              // 2. rotate around x-axis by angle alpha such that P lies on top of z-axis
      //              // alpha = angle between x-y plane and z
      //              alpha = std::atan2(std::sqrt(dirP[1]*dirP[1]+dirP[0]*dirP[0]),dirP[2]);
      //
      //              // no rotation arouind y-axis needed > beta = 0
      //              // WARNING: this whole procedure is only valid if coupling tensor is at least
      //              // transverse isotropic! Otherwise we have to figure out on which axis to rotate
      //              // the transverse directions
      //              beta = 0.0;
      //              Matrix<Double> R = Compute3DRotationMatrix(alpha, beta, gamma);
      //
      //              // take transpose matrix to revert rotation (i.e. rotate z-axis onto dirP)
      //              Matrix<Double> RT;
      //              RT.Resize(3,3);
      //              R.Transpose(RT);
      //
      //              assert(rotatedCouplTensor.GetNumRows() == 3);
      //              assert(rotatedCouplTensor.GetNumCols() == 6);
      //
      //              scaledCouplTensor.PerformRotation(RT,rotatedCouplTensor);
      //
      //            } else {
      //              // 2d plane strain or plane stress case
      //              // Important remark:
      //              //  for 2d plane strain and stress (and somehow also for axi) the coupling tensor
      //              //  will be rotated by default by alpha = -90 and gamma = -90.
      //              //  This results in the following mapping of the coordinate axis:
      //              //    z > y, y > x, x > z
      //              //  I.e. the material is rotated such that the default polarization is in +y direction.
      //              //  We thus have two possibilities to align the material tensor to the 2d polarization
      //              //  direction:
      //              //  a) get original 3x6 tensor and rotate that tensor according to the 3d version above
      //              //      then cut out subtensor
      //              //  b) take cut out subtensor and perform rotation directly in 2d
      //              //      > Version b done in the following
      //              Double gamma;
      //
      //              // std::atan2(dirP[0],dirP[1]) would rotate towards the y-axis
      //              // we want however the y-axis to rotate, so that we take -gamma instead
      //              gamma = -std::atan2(dirP[0],dirP[1]);
      //
      //              Matrix<Double> R = Compute2DRotationMatrix(gamma);
      //
      //              assert(rotatedCouplTensor.GetNumRows() == 2);
      //              assert(rotatedCouplTensor.GetNumCols() == 3);
      //
      //              scaledCouplTensor.PerformRotation(R,rotatedCouplTensor);
      //
      //            }
      //          }
      //        } else {
      //          // no rotation > take scale tensor
      //          rotatedCouplTensor = scaledCouplTensor;
      //        }
      //
      //        rotatedCouplingTensor_requiresReeval_[storageIdx] = false;
      //        rotatedCouplingTensor_[storageIdx] = rotatedCouplTensor;
      //        lastUsedTimeLevelForRotation_[storageIdx] = timeLevel;
      //      }
    }

    /*
     * static member for tracking traced hyst data over multiple instances
     */
    static std::map< std::string, TracedData > tracedOperatorData_;


  private:
    void DetermineChunkSizeAndThreads(UInt size, UInt& chunksize, UInt& numThreads);


    // helper functions for computation of rhs loads and lhs tensors
    //
    //	Matrix<Double> GetScaledAndRotatedCouplingTensor(const LocPointMapped& lpm, bool invert, bool scaleOnly){
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
    //	  Vector<Double> P = hystPol_->GetPrecomputedOutputOfHysteresisOperator(lpm, timeLevel, invert);
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
      // Ref.: C. Woernle, Mehrkörpersysteme "Eine Einführung in die Kinematik und Dynamik von Systemen starrer Körper" 
      // Springer 2016; Kapitel 3.6

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


    void EvaluateHystOperatorsInt(Integer intFlag);
    void EvaluateHystOperators(bool setMatForInversion, bool resetSolToZeroFirst, bool overwriteMemory, bool setRotStateForVecExtension);

    void PrecomputeIrrStrains();

    bool PreprocessLPM(const LocPointMapped& lpmInput, LocPointMapped& lpmOutput,
    UInt& operatorIdx, UInt& storageIdx, UInt& idxElem, bool forceMidpoint = false);

    Vector<Double> CalcOutputOfHysteresisOperator(Vector<Double> inpute, UInt operatorIdx, UInt storageIdx,
    bool forceMemoryLock = false, bool forceMemoryWrite = false, bool useOperatorForStrain = false);

    void ForceRemanence();
    void InitStorage();
    void InitSpecificStorage(int storageLocation);
    void InitLPMMaps();

    Matrix<Double> GetMaterialRelation(const LocPointMapped& Originallpm, int timelevel_new, int timelevel_old,
          Matrix<Double> eps_mu, std::string tensorName, bool prohibitUseOfStorage);

    Matrix<Double> DeltaMatApproximationOfMaterialRelation(Vector<Double>& E_H_diff, Vector<Double>& P_J_diff,
          Matrix<Double>& eps_mu, bool needsInversion);

    Matrix<Double> JacobianApproximationOfMaterialRelation(Vector<Double>& E_H, Matrix<Double>& eps_mu,
          UInt operatorIdx, UInt version, Double steppingSign = 1.0);

    void CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& outputTensor, std::string evalMethod,
    UInt storageIdx, bool intoSat, bool outofSat, bool satToSat, Vector<Double>& X_current);

//    void ExtractSolutionAndInputForHystOperator(Vector<Double>& extractedSolution, Vector<Double>& extractedInput, UInt& operatorIndex,
//    UInt& storageIndex, const LocPointMapped& lpm, bool midpointOnly);

    bool EvaluateAtMidpointOnly();
    bool OverwriteHystMemory();

    // for usage with coefFunctionHystMat
    bool deltaMatActive_;
    // if set to true, timelevel for evaluation will be at the current timestep
    // independent of the actual value of timeLevel_xxx_
    bool forceCurrentTS_;
    int RUN_deltaForm_;

    // NEW timelevel to allow for more flexibility
    Integer timeLevel_Mat_;
    Integer timeLevel_rhsPol_;
    Integer timeLevel_rhsStrain_; // piezo / magstrict only
    Integer timeLevel_rhsCoupling_; // piezo / magstrict only
    Integer timeLevel_deltaMatPol_;
    Integer timeLevel_deltaMatStrain_;
    Integer timeLevel_deltaMatCoupling_;
    Integer timeLevel_results_;
    Integer timeLevel_boundaryTerms_;

    IntegOrder IntegOrder_;
		IntScheme::IntegMethod IntegMethod_;

    // for coupled pdes, we need to know about the mechanical material as well
    // as the coupling tensor;
    // these parameters have to be set from inside the coupled pde (e.g. PiezoCoupling.cc)
    // and have to be set prior to the integrator definition in the single pdes;
    // best place to set these parameter: inside the define integrator function of the coupled
    // pde (as this one will be executed prior to the integrator definition of the single pdes)
    PtrCoefFct elastTensorFct_;
    PtrCoefFct couplTensorFct_;

    Matrix<Double>* rotatedCouplingTensor_;
    bool* rotatedCouplingTensor_requiresReeval_;
    Vector<int> lastUsedTimeLevelForRotation_;

    // for inversion
    bool performInversionTest_;

    /*
     * for slope estimation
     * (can be set in solvestep)
     */
    bool* takeEstimatedSlope_;
    bool ES_includeStrains_;
    bool ES_useAbs_;
    std::string ES_implementationVersion_;
    Double ES_steppingLength_;
    Double ES_scaling_;

    /*
     * ###############################################
     * ### Parameters extracted from material file ###
     * ###############################################
     */
    bool COUPLED_inXMLFile_;
    bool COUPLED_inMatFile_;
    bool COUPLED_ownOperatorForStrains_;

    /*
     * Preisach weights and its size in form of the number of rows
     */
    //    Matrix<Double> POL_weightParams_.weightTensor_;
    Matrix<Double> POL_PreisachWeightsStrains_;
    UInt POL_numRows_;
    UInt POL_numRowsStrain_;

    /*
     * Saturation values for input (x, E/H) and output (y, P/M) of
     * hysteresis operator
     */
    //    Double POL_operatorParams_.inputSat_;
    //    Double POL_operatorParams_.outputSat_;
    Double COUPLED_sSat_; // new strain in saturation

    /*
     * determines whether the vector or the scalar model shall be used
     */
    //    CoefDimType POL_operatorParams_.methodType_;

    /*
     * Additional parameter for scalar Preisach model
     */
    /*
     * Direction (x,y,z) in which Polarization of material points
     * NEW: instead of specifying only x,y,z we allow to specify a direction in the mat.xml file
     */
    //    Vector<Double> POL_operatorParams_.fixDirection_;
    bool POL_useExtension_;
    bool POL_setWithFlux_;

    /*
     * Additional parameter for vector Preisach model
     */
    // determines wheter the vector model from 2012 (classical) or the one from#
    // 2015 (revised) is used
    //    bool POL_operatorParams_.isClassical_;

    /*
     * see VectorPreisach10 for more details
     */
    //    Double POL_operatorParams_.angularDistance_;

    /*
     * if != 0, the rotation states of the vector preisach model will be clipped
     * to the provided precision (in rad)
     */
    //    Double POL_operatorParams_.angularClipping_;

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
     *    (currentInput - oldInput).L2Norm() < POL_operatorParams_.amplitudeResolution_
     *    std::acos(currentInput.Inner(oldInput)) < POL_operatorParams_.angularResolution_
     *
     *  if both criteria are true, the old value is taken and no reevaluation is done
     */
    //    Double POL_operatorParams_.angularResolution_;
    //    Double POL_operatorParams_.amplitudeResolution_;
    //
    /*
     * 1,2 nested-list implementation of classical and revised vector hyst. model
     * 10,20 matrix based implementation of classical and revised vector hyst. model
     */
    //    UInt POL_operatorParams_.evalVersion_;
    //
    /*
     * evaluation of Hyst operator for inputs larger than Saturation
     * > for Mayergoyz hyst model no bisection possible
     */
    //    bool POL_operatorParams_.fieldsAlignedAboveSat_;
    // mayergoyz vector without clipping produces values above PSaturated_
    // when input grows above XSaturated_ > set flag to false
    // otherwise, flag = true
    //    bool POL_operatorParams_.hystOutputRestrictedToSat_;
    //    int POL_operatorParams_.numDirections_;
    //    int POL_operatorParams_.outputClipping_;

    /*
     * Initial input to the vector hysteresis operator;
     * and output of the vector hysteresis operator to this input
     * Note: the value specified in the mat file is the INPUT to the hyst operator
     *        not the actual polarization/magnetization state
     */
    Vector<Double> POL_initialInput_;
    bool POL_prescribeInitialOutput_;

    /*
     * enables output of overlapped switching and rotation state in form of bmp
     * images
     */
    UInt POL_bmpResolution_;
    //    UInt POL_operatorParams_.printOut_;

    bool storageInitialized_;
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
     * EvaluationDepth_ (to be set ONCE at the start of the simulation; set via
     *                    input file possible)
     *
     *  EvaluationDepth_= 1 (standard)
     *   - each element has exactly one hysteresis operator attached
     *   - numHystOperators = numElems
     *   - each Hysteresis operator will only be evaluated at the midpoint of
     *    the corresponding element
     *   - all data structures store the values (elementSolution, current output
     *    of the Hysteresis operator, ...) at the midpoint of the corresponding
     *    element
     *
     *  EvaluationDepth_ = 2 (extended)
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
     *  EvaluationDepth_ = 3 (full)
     *  (not yet implemented, but possible approach for future)
     *   - each INTEGRATION point + the midpoint have one hysteresis operator
     *    attached
     *   - numHystOperators = numElems * (numIntegrationsPoints + 1)
     *   - Evaluation and storage is done for each integration point
     */
    Integer XML_EvaluationDepth_;

    /*
     * activates differnet levels of runtime measurement
     * 0: no measurement
     * 1: measuerment of evaluation time of hyst operator
     * 2: detailed measurement of evaluation time of hyst operator
     * 3: measure performance ONLY during testInversion
     */
    Integer XML_performanceMeasurement_;

    /*
     * #########################################################################
     * ### Parameters the will/might change during execution of StdSolveStep ###
     * ### these parameters are to be set from outside via SetIntegerFlag    ###
     * #########################################################################
     */

    /*
     * may be set to false to lock the rotation state of the hysteresis operator
     * changes in put will thus only affect the amplitude
     * > does not work if false; set it to true and leave it true
     * REMOVED 1.6.2018
     */
    //    bool RUN_overwriteDirection_;

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
    Integer RUN_evaluationPurpose_;

    /*
     * switch needed to activate BMP output only at the end of each time step
     * (otherwise there would be way to much bmp output if done during each iteration)
     */
    bool RUN_allowBMP_;

    /*
     * ######################
     * ### Misc parameter ###
     * ######################
     */
    
    /*
     * For testing of hyst operator in general
     */
    std::string targetDirForResultFiles_;
    bool writeResultsToInfoXML_;
    PtrParamNode infoNodeLoc_;
      
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

    // volume regions
    std::set<RegionIdType> volRegions_;

    //! global element to local element numbering
    std::map<UInt, UInt> globalElem2Local_;
    std::map<UInt, bool> globalElemOnSurf_;

    //! mapping from storageIndex to LPM (original, needs to be preprocessed later)
    std::map<UInt, LocPointMapped > allLPMmap_;
    std::map<UInt, LocPointMapped > midpointLPMmap_;

    //! mappting from operatorIndex to list of all corresponding storageIndices
    //! for later pre-evaluation of hyst-operators; operators can be evaluated in
    //! parallel, single storage entries not (if they use the same operator)
    //! use map for storageIndices such that we can additionally store the flag
    //! isCenter which tells if storageIndex corresponds to element midpoint
    std::map<UInt, std::map<UInt, bool> > operatorToStorage_;

    //! Coefficient function which this one depends one
    PtrCoefFct dependCoef1_;
    PtrCoefFct dependCoef1Surf_;

    //! type of tensor to be returned
    //! electrostatics: permittivity
    //! magnetics: reluctivity
    MaterialType matType_;

    //! type of subtensor (axi, plane)
    SubTensorType tensorType_;

    //! list of all elements
    std::list<shared_ptr<EntityList> > SDList_List_;
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

    //! actual hysteresis operator
    Hysteresis* hyst_;
    Hysteresis* hystTMP;

    Hysteresis* hystStrain_;
    Hysteresis* hystStrainTMP;

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
     * P_J_       P - ElecPolarization      J - mag Polarization -> output of the Hysteresis operator
     *
     * E_H_       E - ElecField             H - MagField        -> input to the Hysteresis operator
     *                                        (needed as
     *                                        VectorPreisach has
     *                                        no inverse from yet)
     *
     *
     */
    Vector<Double>* E_B_;
    Vector<Double>* P_J_;
    Vector<Double>* P_J_forStrains_;
    Vector<Double>* E_H_;
    Vector<Double>* Si_;

    /*
     * Storage vectors for quantities of the last iteration;
     * these values will ONLY be overwritten during each iteration by calling the
     * function
     *    SetPreviousHystValues(0)
     *
     */
    Vector<Double>* E_B_lastIt_;
    Vector<Double>* P_J_lastIt_;
    Vector<Double>* E_H_lastIt_;
    Vector<Double>* Si_lastIt_;

    /*
     * Storage vectors for quantities of optional estimator step;
     * Idea: perform estimator step at beginning of each timestep to check
     *  if current boundary conditions/rhs load would lead to a reversal of
     *  hyst operator; this helps to find the correct stepping direction for
     *  approximation of Jacobian with help of finite differences
     *    SetPreviousHystValues(2)
     *
     * NOTE: check which of these are actually necessay; most probably only E_B
     *
     */
    bool estimatorSet_;
    Vector<Double>* E_B_estimatorStep_;
    Vector<Double>* P_J_estimatorStep_;
    Vector<Double>* E_H_estimatorStep_;
    Vector<Double>* Si_estimatorStep_;

    Vector<Double>* E_B_BAK_;
    Vector<Double>* P_J_BAK_;
    Vector<Double>* E_H_BAK_;
    Vector<Double>* Si_BAK_;


    /*
     * Storage vectors for quantites of the last timestep;
     * these values will ONLY be set during the FIRST iteration of each time step
     * by calling
     *    SetPreviousHystValues(1)
     */
    Vector<Double>* E_B_lastTS_;
    Vector<Double>* P_J_lastTS_;
    Vector<Double>* P_J_forStrains_lastTS_;
    Vector<Double>* E_H_lastTS_;
    Vector<Double>* Si_lastTS_;
    bool* requiresReeval_;
    bool* Si_requiresReeval_;


    /*
     * New idea: if convergence using deltaMat is bad, we might
     * bring in some even older iterates to the computation
     * e.g. deltaMat = deltaMatCenteral = (P_J_ - P_J_secondtolastIt_ )/(E_B_ - E_B_secondtolastIt_ )
     * which would be "second order" around P_J_lastIt_
     * or more general
     * gamma*(P_J_ - P_J_lastIt_ )/(E_B_ - E_B_lastIt_ )+(1-gamma)(P_J_lastIt_ - P_J_secondtolastIt_ )/(E_B_lastIt_ - E_B_secondtolastIt_ )
     *
     *  = alpha*deltaMat_ + (1-alpha)*deltaMatPrev_
     * > only create storage for deltaMatPrev_
     * > set previous values using SetPreviousHystValues
     */
    Matrix<Double>* deltaMat_;
    Matrix<Double>* deltaMatPrev_;
    bool* deltaMat_requiresReeval_;

    Matrix<Double>* deltaMatStrain_;
    Matrix<Double>* deltaMatStrainPrev_;
    bool* deltaMatStrain_requiresReeval_;

    Double includeOldDelta_;

    Matrix<Double>* eps_mu_local_; // used for local inversion of preisach operator > eps/mu; might depend on P itself
    Matrix<Double>* epsInv_nu_local_; // used for local inversion of preisach operator > eps/mu; might depend on P itself
    bool* eps_mu_local_Set_;

    // needed to get integration point
    shared_ptr<FeSpace> ptFeSpace_;
    UInt numIntegrationPoints_;
    std::map<int,UInt> locPointIndices_;

    // flags indicating if we need inversion (for e.g. for magnetics)
    // and if we have an inverse model for that
    bool needsInversion_;
    //    bool POL_operatorParams_.hasInverseModel_;
    // For inversion via Levenberg Marquardt
    //    UInt InversionParams_.maxNumIts;
    //    Double InversionParams_.tolH;
    //    Double InversionParams_.tolB;
    //    Double InversionParams_.jacRes;
    //    bool InversionParams_.useTikhonov;
    //    Double InversionParams_.alphaLSStart;
    //    Double InversionParams_.alphaLSMin;
    //    Double InversionParams_.alphaLSMax;

    // when setting the system into remanence for the case of magnetics
    // we do this by forcing the input to the hystoperator to become 0;
    // for electrostatics, we can simply set the solution vector to 0 but
    // for magnetics, this would cause B = 0 and thus bring the system to
    // coercive points; nevertheless we want to reset solVec to 0; to avoid
    // evaluating the hyst operator (after bringing it to remanence by hand)
    // we simply lock its evaluation (via setFlag)
    bool hystOperatorLocked_;


    /*
     * new: initialize storage on demand
     */
    bool lastTSStorageInitialized_;
    bool lastItStorageInitialized_;
    bool backupStorageInitialized_;
    bool estimatorStorageInitialized_;
    bool strainRequired_;

    bool DEBUG_infoOutput_;
    bool DEBUG_debugOutput_;
    bool DEBUG_printCallingOrder_;

    /*
     * FIXPOINT ITERATION (for magnetics)
     * I) B-Version
     *  B_k+1 = B_k + nu_FP^-1 (H_J - H(B_k))
     *        = B_k + nu_FP^-1 (H_J - nu_0 B_k + nu_0 P(H_k))
     *  P(H_k) = Preisach Operator
     *  rot(H_J) = J_i
     *
     *  FE Procedure:
     *    1a. [FULL STEP] > NOT USED
     *       Solve rot(nu_FP B_k+1) = J + rot( (nu_FP - nu_0) B_k ) + rot( nu_0 P(H_k) )
     *        or
     *    1b. [INCREMENTAL STEP] > USED
     *       Solve rot(nu_FP DeltaB) = J + rot( (nu_FP - nu_0) B_k ) + rot( nu_0 P(H_k) ) - rot( nu_FP B_k )
     *                                = J - rot( nu_0 B_k ) + rot( nu_0 P(H_k) )
     *                                = residual(B_k)
     *       B_k+1 = B_k + eta_LS DeltaB
     *    2. Solve via local inversion (i.e. on element level)
     *       H_k+1 = nu_0 (B_k+1 - P(H_k+1))
     *
     * I.1) Global version
     *  nu_FP = globalFPFactor_C * 0.5 * (nu_max_global + nu_min_global)
     *    > get evaluated only once at beginning of transient simulation
     *    > nu_max and nu_min must be the global maxima and minima; estimation done by tracing hyst-operator
     *    > globalFPFactor_C usually 1.0 (only introduced for sake of completeness)
     *
     * I.2) Local version
     *  nu_FP = localFPFactor_C * 0.5 * (nu_max_currentTS + nu_min_currentTS)
     *    > get evaluated at the beginning of each time-step
     *    > nu_max and nu_min must be the maxima and minima at the last known solution (previous time-step); estimation
     *        from computing the Jacobian at the previous time-step solution
     *    > localFPFactor_C > 1; usually 1-2; the larger localFPFactor_C the larger the range in which solution can con-
     *        verge but the slower the convergence rate
     *
     * Required data structure:
     *  Global version:
     *    globalFPFactor_C > defined via input.xml
     *    nu_FP_global > computed from global minima and maxima of nu > retrieved from tracing hyst operator once
     *    nu_0 or nu_S for computation of residual
     *
     *  Local version:
     *    localFPFactor_C > defined via input.xml
     *    nu_FP_local[numStorageIdx] > computed from jacobian at last time-step solution; different for each storage idx
     *    nu_0 or nu_S for computation of residual
     */


    /* II) H-Version
     *  H_k+1 = H_k + mu_FP^-1 (B_A - B(H_k))
     *        = H_k + mu_FP^-1 (B_A - mu_0 H_k - P(H_k))
     *  P(H_k) = Preisach Operator
     *  B_A = Solution of FE-System
     *
     *  FE Procedure:
     *    1a. [FULL STEP] > NOT USED
     *       Solve rot(mu_0^-1 B_k+1) = J + rot( mu_0^-1 P(H_k) )
     *        or
     *    1b. [INCREMENTAL STEP] > USED
     *       Solve rot(mu_0^-1 DeltaB) = J - rot( mu_0^-1 B_k ) + rot( mu_0^-1 P(H_k) )
     *                                 = residual(B_k)
     *       B_k+1 = B_k + eta_LS DeltaB
     *    2. Update H_k+1 using FP method
     *       H_k+1 = H_k + mu_FP^-1 (B_k+1 - B_k)
     *             = H_k + mu_FP^-1 (B_k+1 - mu_0 H_k - P_k)
     *             = H_k - (mu_0/mu_FP)H_k + mu_0/mu_FP (B_k+1 - P_k)/mu_0
     *             = (1-omega) H_k + omega (B_k+1 - P_k)/mu_0 > derzeitige implementierung
     *
     * II.1) Global version (see e.g. Kurz, Fetzer, Lehner)
     *  0 < omega < 2mu_0/mu_max
     *  with mu_max = (deltaB/deltaH)_max
     *    --> mu_FP = mu_0/omega > mu_max/2
     *    --> mu_FP = C*mu_max/2 with C > 1
     *    or with mu_min >= mu_0 > 0
     *        mu_FP = globalFPFactor_C*0.5*(mu_max + mu_min) with globalFPFactor_C >= 1
     *
     * II.2) Local version > similar as in B-Version
     *    > mu_FP = localFPFactor_C*0.5*(mu_max + mu_min)
     *    > mu_max and mu_min are retrieved from Jacobian at beginning of time-step
     *
     * Required data structure:
     *  Global version:
     *    globalFPFactor_C > defined via input.xml
     *    mu_FP_global > computed from global minima and maxima of mu > retrieved from tracing hyst operator once
     *    mu_0 or mu_S for computation of residual and system
     *
     *  Local version:
     *    localFPFactor_C > defined via input.xml
     *    mu_FP_local[numStorageIdx] > computed from jacobian at last time-step solution; different for each storage idx
     *    mu_0 or mu_S for computation of residual and system
     */




    /*
     * new for fixpoint iterations:
     * the original material tensor nu_0/eps_0 will not lead to (good) convergence(
     * (see e.g. "Locally Convergent Fixed-Point Method for Solving Time-Stepping Nonlinear Field Problems" - E. Dlala et al. 2007)
     * instead use nu_FP/eps_FP which leads to better convergence
     * the value of nu_FP can be constant for all timesteps or can be adapted in each timestep (e.g. to the current Jacobian)
     * note: if nu != nu_0 but nu_FP the material law
     * H = nu*B - nu_0*J_Preisach
     * is changed, too! to ensure, that we obtain the same material relation independent on nu_FP we have to adapt J
     * i.e.
     * H = nu_0*B - nu_FP*B + nu_FP*B - nu_0*J_Preisach
     * H = nu_FP*B - nu_0*( J_Preisach - (1 - nu_FP/nu_0)*B )
     *   = nu_FP*B - nu_0*J_FP
     *  > the term (1 - nu_FP/nu_0)*B gets subtracted from rhs term if flag is set
     */
    // TODO: should be synchronized with coupled case where we generally have nuS or nuT instead of nu_0
    //        > in this case, the FPMaterialTensor should not be rev_mat_fac_MATRIX but adapted to nuS / nuT
    Matrix<Double>* FPMaterialTensor_;
    UInt FPMaterialTensorSet_; // 0 > not active; 1 > global (const for all TS); 2 > local (changing in each TS)
    UInt FPMaterialRHSSet_; // correction for rhs; we basically have two options: a) use system matrix with FP tensor
    // to subtract old solution from rhs > this requires correction of hyst value
    // b) use system matrix with standard material tensor (instead of FP one) > no correction of hyst value needed
    //  (as basically the rhs remains the same for incremental stepping, even if system matrix uses FP tensor)
    Matrix<Double> rev_mat_fac_MATRIX_;

    bool useFPH_;
    bool skipStorage_; // new parameter introduced 18.4.2020; see function PrecomputePolarization for more info
    fixpointFlag selectedFPApproach_;

    // convergence factors for global and local fixpoint methods
    // > can be set via input.xml
    Double globalFPFactor_C;
    Double localFPFactor_C;
    
    // parameters to adapt tracing
    Double trace_JacResolution_;
    bool trace_forceCentral_;
    bool trace_forceRetracing_;

    // information from hyst operator
    bool hystModelTraced_;
    bool globalFPset_;
    Double maxSlopeGlobal_;
    Double minSlopeGlobal_;
    Vector<Double> maxSlopeLocal_;
    Vector<Double> minSlopeLocal_;
    // for H-version; evaluate muFP_ in function set FPMaterialTensor
    //  and use it in RetrieveInputToHystOperator
    // (for B-version include nuFP_ directly in FPMaterialTensor_)
    Vector<Double> muFPH_;
    // 24.6.2020: introduced new flag
    // if set to true, the slope for the localized fixpoint methods will only be computed from the Jacobian at the
    // elements center; if false, the Jacobian will be evaluated at each integration point and each integration point
    // will gets its own FP-tensor
    bool evalJacAtMidpointOnly_;

    /*
     * inherited from HystHelper
     */

    PtrCoefFct ptrFieldTensor_;
    PtrCoefFct ptrElastTensor_;
    PtrCoefFct ptrCouplTensor_;
    bool isCoupled_;

    // base tensors for computation of (eventually) hyst-dependent material tensors
    // e.g. nu(P) = P/P_sat*eps_nu_base_
    Matrix<Double> eps_nu_base_ ;
    Matrix<Double> eps_mu_base_ ;
    Matrix<Double> elastTensor_;
    Matrix<Double> couplTensor_;
    
    Matrix<Double> eps_nu_baseFULL_ ;
    Matrix<Double> eps_mu_baseFULL_ ;
    Matrix<Double> elastTensorFULL_;
    Matrix<Double> couplTensorFULL_;

    bool tensorsInitialized_;
    bool materialTensorsSetOnce_;

    /*
     * -1 : not coupled at all
     *  0 : coupled e-form/h-form (piezo/magstrict)
     *  1 : coupled d-form (piezo)
     *  2 : coupled g-form (magstrict)
     */
    int strainForm_;
    bool formerHystHelperParamsSet_;


  };


}
#endif
