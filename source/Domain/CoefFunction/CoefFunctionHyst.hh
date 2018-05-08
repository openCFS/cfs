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
  
  class CoefFunctionHelper : public CoefFunction{
  public:
    
    // note: coupling tensor has to be in non-transposed form (i.e. 3x6 or 2x3, 2x4)
    CoefFunctionHelper(PtrCoefFct ptrFieldTensor,PtrCoefFct ptrElasticityTensor,
      PtrCoefFct ptrCouplingTensor, PtrCoefFct CoefFncHyst)
    :CoefFunction(){
      
      isAnalytic_ = false;
      isComplex_  = false;
      dependType_ = SOLUTION;
      
      ptrFieldTensor_ = ptrFieldTensor;
      ptrElastTensor_ = ptrElasticityTensor;
      ptrCouplTensor_ = ptrCouplingTensor;
      hystCoefFunction_ = CoefFncHyst;
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
      
      if(isCoupled_ == false){
        // no coupling tensors at all
        strainForm_ = -1;
      } else {
        bool useStrainForm = hystCoefFunction_->useStrainForm();
        if(!useStrainForm){
          // just take input tensors as they are, i.e. e-form or h-form
          strainForm_ = 0;
        }else{
          std::string PDEName = hystCoefFunction_->getPDEName();
          if(PDEName == "Electrostatic"){
            // piezo case
            // compute d-form as basis for coupling
            strainForm_ = 1;
          } else {
            // magstrict case
            // compute g-form 
            strainForm_ = 2;
          }
        }
      }
      
      // !!! small material parameter are assumed to be constant over region !!!
      // i.e. for each region, just one linear material tensor of eac kind stored
      // to obtain those values, we have to evaluate the ptrCoefFunctions at an lpm
      // do this during the first call to getTensor
      tensorsInitialized_ = false;
    }
    
    virtual ~CoefFunctionHelper(){};
    
    void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
      EXCEPTION("CoefFunctionHelper returns no scalars")
    }
    
    void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
      EXCEPTION("CoefFunctionHelper returns no vectors")
    }
    
    void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
      EXCEPTION("CoefFunctionHelper returns no tensors")
    }
    
    UInt GetVecSize() const {
      return hystCoefFunction_->GetVecSize();
    }
    
    void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      hystCoefFunction_->GetTensorSize(numRows, numCols);
    }
    
    void GetCouplTensorSize( UInt& numRows, UInt& numCols ) const {
      ptrCouplTensor_->GetTensorSize(numRows,numCols);
    }
    
    int GetTimeLevel(std::string name){
      return hystCoefFunction_->GetTimeLevel(name);
    }
    
    std::string ToString() const {
      return "Dummy string";
    }
        
    void PrecomputeMaterialTensorForInverison();
    
    void ComputeTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm, int timeLevel_to_diff, 
      std::string tensorName, std::string implementationVersion, bool transposed, bool rotate, bool useAbs, bool lockPrecomputationAndDeltaMat=false );
        
    void ComputeVector(Vector<Double>& outputVector,const LocPointMapped& lpm, int timeLevel, int baseSign, std::string vectorName, bool onBoundary );
    
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
      
      if(strainForm_ == -1){
        //std::cout << "Single field case" << std::endl;
        /*
         * single field case
         * > only fieldTensor needed
         * > remaining tensors stay unitinialized
         */
        fieldTensor_ = epsS_nuS;
        
        //std::cout << "Initialized tensors: " << std::endl;
        //std::cout << "Field Tensor: " << fieldTensor_.ToString() << std::endl;
      } else if(strainForm_ != -1){
        //std::cout << "Coupled case" << std::endl;
        /*
         * Coupled case
         */
        // get e-form / h-form
        UInt numRows1, numCols1, numRows2, numCols2;
        ptrElastTensor_->GetTensorSize(numRows1,numCols1);
        Matrix<Double> cE_B = Matrix<Double>(numRows1,numCols1);
        ptrElastTensor_->GetTensor(cE_B,lpm);
        
        ptrCouplTensor_->GetTensorSize(numRows2,numCols2);
        Matrix<Double> e_h = Matrix<Double>(numRows2,numCols2);
        ptrCouplTensor_->GetTensor(e_h,lpm);
        
        if(strainForm_ == 0){
          //std::cout << "Coupled case - e/h form" << std::endl;
          // just keep e-form/h-form
          fieldTensor_ = epsS_nuS;
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
          
          fieldTensor_ = epsT;
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
          
          fieldTensor_ = nuT;
          couplTensor_ = g;
          // NOTE: we do not keep sE but take cE instead
          // later during the solution, we will need cE (not sE) and as cE is
          // not depending on hysteresis, we can directly use the linear tensor
          elastTensor_ = cE_B;
          
        }
        //        std::cout << "Initialized tensors: " << std::endl;
        //        std::cout << "Field Tensor: " << fieldTensor_.ToString() << std::endl;
        //        std::cout << "Coupling Tensor: " << couplTensor_.ToString() << std::endl;
        //        std::cout << "Elasticity Tensor: " << elastTensor_.ToString() << std::endl;
      }
      
      tensorsInitialized_ = true;
    }
    
    PtrCoefFct ptrFieldTensor_;
    PtrCoefFct ptrElastTensor_;
    PtrCoefFct ptrCouplTensor_;
    bool isCoupled_;
    
    PtrCoefFct hystCoefFunction_;
    
    Matrix<Double> fieldTensor_ ;
    Matrix<Double> elastTensor_;
    Matrix<Double> couplTensor_;
    
    bool tensorsInitialized_;
    bool materialTensorsSetOnce_;
    
    /*
     * -1 : not coupled at all
     *  0 : coupled e-form/h-form (piezo/magstrict)
     *  1 : coupled d-form (piezo)
     *  2 : coupled g-form (magstrict)
     */
    int strainForm_;
  };
  
  class CoefFunctionHystMatTensor : public CoefFunction {
  public:
		
    // note: coupling tensor has to be in non-transposed form (i.e. 3x6 or 2x3, 2x4)
    CoefFunctionHystMatTensor(CoefFunctionHelper* hystHelper, std::string tensorName, bool transposed)
    :CoefFunction(){
          
      isAnalytic_ = false;
      isComplex_  = false;
      dependType_ = SOLUTION;
        
      // HystMat specific
      hystHelper_ = hystHelper;
      dimType_ = TENSOR;      
      tensorName_ = tensorName;
      transposed_ = transposed;
    }
    
    virtual ~CoefFunctionHystMatTensor(){};
    
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
      int timeLevel_to_diff = hystHelper_->GetTimeLevel("DeltaMat"); 
      
      hystHelper_->ComputeTensor(outputTensor, lpm, timeLevel_to_diff, tensorName_, implementationVersion, transposed_, rotate, useAbs );
//      std::cout << "GET TENSOR - 2" << std::endl;
    }
    
    //! Return row and columns size of tensor if coefficient function is a tensor
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      if( (tensorName_ == "Permittivity") || (tensorName_ == "Reluctivity")){
        hystHelper_->GetTensorSize(numRows,numCols);
      } else if ( (tensorName_ == "CouplingMechToElec") || (tensorName_ == "CouplingElecToMech") 
        || (tensorName_ == "CouplingMechToMag") || (tensorName_ == "ComputeMagToMech") ){
        hystHelper_->GetCouplTensorSize(numRows,numCols);
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
    
    CoefFunctionHelper* hystHelper_;
    std::string tensorName_;
    bool transposed_;
  };
  
  class CoefFunctionHystRHSLoad : public CoefFunction {
    
  public:
    
    CoefFunctionHystRHSLoad(CoefFunctionHelper* hystHelper,
      std::string vectorName, bool transposed, bool onBoundary)
    :CoefFunction(){
      
      isAnalytic_ = false;
      isComplex_  = false;
      dependType_ = SOLUTION;
        
      // HystRHSLoad specific
      hystHelper_ = hystHelper;
      dimType_ = VECTOR;            
      vectorName_ = vectorName;
      transposed_ = transposed;
      onBoundary_ = onBoundary;
    }
    
    virtual ~CoefFunctionHystRHSLoad(){};
    
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
        timeLevel = hystHelper_->GetTimeLevel("BC");
        //std::cout << "Boundary Term, timelevel = " << timeLevel << std::endl;
        // furtherrmore, the  boundary terms and the (volume) loads have opposite
        // aigns
        baseSign = -1.0;
      } else {
        timeLevel = hystHelper_->GetTimeLevel("RHS");
        //std::cout << "RHS Term, timelevel = " << timeLevel << std::endl;
        baseSign = 1.0;
      }
      
      hystHelper_->ComputeVector(outputVector,lpm, timeLevel, baseSign, vectorName_, onBoundary_ );
    }
    
    virtual UInt GetVecSize() const {
      
      UInt numRows,numCols;
      if( (vectorName_ == "ElecPolarization") || (vectorName_ == "MagPolarization") || (vectorName_ == "MagMagnetization")){
        hystHelper_->GetTensorSize(numRows,numCols);
      } else if ( (vectorName_ == "PiezoLoadForMechPDE") || (vectorName_ == "PiezoLoadForElecPDE") 
        || (vectorName_ == "MagStrictLoadForMechPDE") || (vectorName_ == "MagStrictLoadForMagPDE") ){
        hystHelper_->GetCouplTensorSize(numRows,numCols);
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
      hystHelper_->GetTensorSize(numRows,numCols);
    }
    
    std::string ToString() const {
      return "Dummy string";
    }
    
  private:
    CoefFunctionHelper* hystHelper_;
    std::string vectorName_;
    bool transposed_;
    bool onBoundary_;
  };
  
  
  //class CoefFunctionHystMatTensor : public CoefFunction{
  //  public:
  //		
  //    // note: coupling tensor has to be in non-transposed form (i.e. 3x6 or 2x3, 2x4)
  //    CoefFunctionHystMatTensor(PtrCoefFct ptrFieldTensor,PtrCoefFct ptrElasticityTensor,
  //      PtrCoefFct ptrCouplingTensor, PtrCoefFct CoefFncHyst,
  //      std::string tensorName, bool transposed)
  //    :CoefFunction(){
  //      
  //      dimType_ = TENSOR;
  //      isAnalytic_ = false;
  //      isComplex_  = false;
  //      dependType_ = SOLUTION;
  //      
  //      ptrFieldTensor_ = ptrFieldTensor;
  //      ptrElastTensor_ = ptrElasticityTensor;
  //      ptrCouplTensor_ = ptrCouplingTensor;
  //      hystCoefFunction_ = CoefFncHyst;
  //      
  //      UInt couplRows, couplCols;
  //      if(ptrCouplTensor_ != NULL){
  //        ptrCouplTensor_->GetTensorSize(couplRows,couplCols);
  //        if(couplCols < couplRows){
  //          EXCEPTION("Coupling tensor has to be in non-transposed form");
  //        }
  //        isCoupled_ = true;
  //      } else {
  //        isCoupled_ = false;
  //      }
  //      
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
  //      
  //      tensorName_ = tensorName;
  //      transposed_ = transposed;
  //      
  //      // !!! small material parameter are assumed to be constant over region !!!
  //      // i.e. for each region, just one linear material tensor of eac kind stored
  //      // to obtain those values, we have to evaluate the ptrCoefFunctions at an lpm
  //      // do this during the first call to getTensor
  //      tensorsInitialized_ = false;
  //    }
  //    
  //    virtual ~CoefFunctionHystMatTensor(){};
  //    
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
  //    
  //    void InitLinearTensors(const LocPointMapped& lpm){
  //      std::cout << "Init Linear Tensors" << std::endl;
  //      
  //      // !!! small material parameter are assumed to be constant over region !!!
  //      // i.e. for each region, just one linear material tensor of eac kind stored
  //      // to obtain those values, we have to evaluate the ptrCoefFunctions at an lpm
  //      // do this during the first call to getTensor
  //      UInt numRows, numCols;
  //      ptrFieldTensor_->GetTensorSize(numRows,numCols);
  //      Matrix<Double> epsS_nuS = Matrix<Double>(numRows,numCols);
  //      ptrFieldTensor_->GetTensor(epsS_nuS,lpm);
  //      
  //      if(strainForm_ == -1){
  //        std::cout << "Single field case" << std::endl;
  //        /*
  //         * single field case
  //         * > only fieldTensor needed
  //         * > remaining tensors stay unitinialized
  //         */
  //        fieldTensor_ = epsS_nuS;
  //        
  //        std::cout << "Initialized tensors: " << std::endl;
  //        std::cout << "Field Tensor: " << fieldTensor_.ToString() << std::endl;
  //      } else if(strainForm_ != -1){
  //        std::cout << "Coupled case" << std::endl;
  //        /*
  //         * Coupled case
  //         */
  //        // get e-form / h-form
  //        UInt numRows1, numCols1, numRows2, numCols2;
  //        ptrElastTensor_->GetTensorSize(numRows1,numCols1);
  //        Matrix<Double> cE_B = Matrix<Double>(numRows1,numCols1);
  //        ptrElastTensor_->GetTensor(cE_B,lpm);
  //        
  //        ptrCouplTensor_->GetTensorSize(numRows2,numCols2);
  //        Matrix<Double> e_h = Matrix<Double>(numRows2,numCols2);
  //        ptrCouplTensor_->GetTensor(e_h,lpm);
  //        
  //        if(strainForm_ == 0){
  //          std::cout << "Coupled case - e/h form" << std::endl;
  //          // just keep e-form/h-form
  //          fieldTensor_ = epsS_nuS;
  //          couplTensor_ = e_h;
  //          elastTensor_ = cE_B;
  //        } else if(strainForm_ == 1){
  //          std::cout << "Coupled case - d form" << std::endl;
  //          // piezo
  //          // transform to d-form
  //          // d = cE^-1*e = sE*e
  //          Matrix<Double> sE = Matrix<Double>(numRows1,numCols1);
  //          cE_B.Invert(sE);
  //          
  //          Matrix<Double> d = Matrix<Double>(numRows2,numCols2);
  //          sE.Mult(e_h,d);
  //          
  //          // epsT = epsS + d e^t
  //          Matrix<Double> epsT = Matrix<Double>(numRows,numCols);
  //          Matrix<Double> e_trans = Matrix<Double>(numCols2,numRows2);
  //          e_h.Transpose(e_trans);
  //          d.Mult(e_trans,epsT);
  //          epsT.Add(1.0,epsS_nuS);
  //          
  //          fieldTensor_ = epsT;
  //          couplTensor_ = d;
  //          // NOTE: we do not keep sE but take cE instead
  //          // later during the solution, we will need cE (not sE) and as cE is
  //          // not depending on hysteresis, we can directly use the linear tensor
  //          elastTensor_ = cE_B;
  //        } else if(strainForm_ == 2){
  //          std::cout << "Coupled case - g form" << std::endl;
  //          // magstrict
  //          // transform to g-form
  //          // g = cE_B^-1*h = sB*h
  //          Matrix<Double> sB = Matrix<Double>(numRows1,numCols1);
  //          cE_B.Invert(sB);
  //          
  //          Matrix<Double> g = Matrix<Double>(numRows2,numCols2);
  //          sB.Mult(e_h,g);
  //          
  //          // nuT = nuS - g h^t
  //          Matrix<Double> nuT = Matrix<Double>(numRows,numCols);
  //          Matrix<Double> h_trans = Matrix<Double>(numCols2,numRows2);
  //          e_h.Transpose(h_trans);
  //          g.Mult(h_trans,nuT);
  //          nuT.Add(-1.0,epsS_nuS);
  //          nuT *= (-1.0);
  //          
  //          fieldTensor_ = nuT;
  //          couplTensor_ = g;
  //          // NOTE: we do not keep sE but take cE instead
  //          // later during the solution, we will need cE (not sE) and as cE is
  //          // not depending on hysteresis, we can directly use the linear tensor
  //          elastTensor_ = cE_B;
  //          
  //        }
  //        std::cout << "Initialized tensors: " << std::endl;
  //        std::cout << "Field Tensor: " << fieldTensor_.ToString() << std::endl;
  //        std::cout << "Coupling Tensor: " << couplTensor_.ToString() << std::endl;
  //        std::cout << "Elasticity Tensor: " << elastTensor_.ToString() << std::endl;
  //      }
  // 
  //      tensorsInitialized_ = true;
  //    }
  //    
  //    void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
  //      std::cout << "Coef Function Hyst Mat - Get Tensor" << std::endl;
  //      UInt numCols, numRows;
  //      // GetTensorSize will return the size of the actual tensor that shall
  //      // be returned, i.e. if it is the transposed of the couplTensor_
  //      // numCols = numRows(coupTensor_)
  //      GetTensorSize(numRows,numCols);
  //      Matrix<Double> tmp;
  //      
  //      outputTensor.Resize(numRows,numCols);
  //      outputTensor.Init();
  //      
  //      if(tensorsInitialized_ == false){
  //        InitLinearTensors(lpm);
  //      }
  //      
  //      if(transposed_){
  //        tmp = Matrix<Double>(numCols,numRows);
  //      } else {
  //        tmp = Matrix<Double>(numRows,numCols);
  //      }
  //      tmp.Init();
  //      
  //      int deltaForm_ = hystCoefFunction_->GetDeltaForm();
  //      bool deltaFormActive_ = hystCoefFunction_->deltaMatActive();
  //      bool rotate = true;
  //      bool useAbs = false;
  //      std::string implementationVersion = "Division";
  //      
  //      /*
  //       * timelevel determines which value of P shall be used to evaluate the
  //       * scaling and rotation of the matrices; 
  //       * the tensors returned by this functions should normally be on the current
  //       * time level, i.e. timelevel should be 0
  //       */
  //      int timelevel_for_diff = hystCoefFunction_->GetTimeLevel("Mat"); 
  //      int timelevel_cur = 0;
  //      
  //      if(tensorName_ == "Permittivity"){
  //        std::cout << "Get Permittivity" << std::endl;
  //        /*
  //         * The following cases are possible:
  //         *  I. pure electrostatics:
  //         *     a) no-deltaMat > return epsS
  //         *     b) deltamat    > return epsS + dP/dE
  //         *  II. coupled piezo-electric case
  //         *     a) e-form as basis
  //         *        return epsS + dP/dE - e*dS/dE
  //         *     b) d-form as basis
  //         *        return epsT - d [c] d^T + dP/dE - d [c] dS/dE
  //         */
  //
  //        // in case of piezoelectricity, we have to check for e-form or d-form
  //        // in case of e-form, fieldTensor is seens (and set to be) epsS which needs
  //        // no further treatment
  //        Matrix<Double> e_scaled;
  //        Matrix<Double> tmp2;
  //        if(strainForm_ == 1){
  //          std::cout << "Get Permittivity - Compute epsS from epsT" << std::endl;
  //          // d-form shall be used as basis; fieldTensor such represents epsT
  //          // to get the required epsS we have to compute
  //          //  epsS = epsT - d(P)[c]d(P)^T
  //          // where d(P) is the scaled and rotated coupling tensor in d-form
  //          Matrix<Double> d_scaled, d_scaled_transposed;
  //          
  //          hystCoefFunction_->ScaleAndRotateCouplingTensor(lpm,couplTensor_,d_scaled,timelevel_cur,rotate);
  //
  //          // e = d*c
  //          d_scaled.Mult(elastTensor_,e_scaled);
  //          d_scaled.Transpose(d_scaled_transposed);
  //          // tmp2 = e*d^T = dcd^T
  //          e_scaled.Mult(d_scaled_transposed,tmp2);
  //
  //          // tmp = epsT
  //          tmp = fieldTensor_;
  //
  //          // tmp = epsT - d*c*d^T
  //          tmp.Add(-1.0,tmp2);
  //          
  //        } else {
  //          std::cout << "Get Permittivity - Take epsS directly" << std::endl;
  //          // uncoupled case or e-form
  //          // just take epsS directly
  //          tmp = fieldTensor_;
  //        }
  //        
  //        // check if deltaMat shall be added
  //        // here we have two flags:
  //        //  deltaForm_ is used to indicate if we are using a delta formulation in general
  //        //  deltaFormActive_ indicates if we actually need this deltaMatrix for the current
  //        //    evaluation (the issue is, that we need a deltaMatrix on the lhs and a non-deltaMatrix
  //        //    on the rhs
  //        if(deltaFormActive_ && (deltaForm_ != 0) ) {
  //          std::cout << "Get Permittivity - Compute DeltaMatrix" << std::endl;
  //          //        std::cout << "Compute DeltaMatrix" << std::endl;
  //          bool useStrains = false;
  //          
  //          // deltaMat will be computed using the current value (timelevel_cur = 0)
  //          // and the value at timelevel (-1 > last ts; +1 > last iteration)
  //          Matrix<Double> deltaMat = hystCoefFunction_->GetDeltaMat(lpm, timelevel_cur, timelevel_for_diff, useStrains, useAbs,implementationVersion);
  //          
  //          tmp.Add(1.0,deltaMat);
  //          std::cout << "DeltaMat: " << deltaMat.ToString() << std::endl;
  //          if(strainForm_ != -1){
  //            std::cout << "Get Permittivity - Add dS/dE" << std::endl;
  //            // coupled case
  //            // here we have to add -e*dS/dE in addition
  //            useStrains = true;
  //            Matrix<Double> deltaMat_strains = hystCoefFunction_->GetDeltaMat(lpm, timelevel_cur, timelevel_for_diff, useStrains, useAbs,implementationVersion);
  //            e_scaled.Mult(deltaMat_strains,tmp2);
  //            
  //            tmp.Add(-1.0,tmp2);
  //            std::cout << "DeltaMatStrain: " << deltaMat_strains.ToString() << std::endl;
  //          }
  //          
  //        } else if (deltaForm_ != 0){
  //          //        std::cout << "deltaFormSet_ == true, but deltaMatNotActive!" << std::endl;
  //        } else {
  //          //        std::cout << "deltaFormSet_ == false" << std::endl;
  //        }
  //        
  //      } else if (tensorName_ == "CouplingMechToElec"){
  //        //      std::cout << "ComputeMechToElec" << std::endl;
  //        Matrix<Double> rotatedCouplTensor;
  //     
  //        if(strainForm_ == 1){
  //          // use d-form as basis, i.e. the followings steps have to be applied
  //          // 1. scale and rotate d
  //          hystCoefFunction_->ScaleAndRotateCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
  //          // 2. compute e from d
  //          // e = d*c
  //          rotatedCouplTensor.Mult(elastTensor_,tmp);
  //        } else {
  //          // use e-form as basis, i.e. the following step has to be done
  //          // 1. acale and rotate e
  //          hystCoefFunction_->ScaleAndRotateCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
  //          tmp = rotatedCouplTensor;
  //        }
  //        
  //      } else if (tensorName_ == "CouplingElecToMech"){ 
  //        // this is basically the same tensor as mechToElec except for the case of
  //        // deltaFormulation; in the later case, we have to add c*dS/dE
  //        
  //        //      std::cout << "ComputeElecToMech" << std::endl;
  //        Matrix<Double> rotatedCouplTensor;
  //     
  //        if(strainForm_ == 1){
  //          // use d-form as basis, i.e. the followings steps have to be applied
  //          // 1. scale and rotate d
  //          hystCoefFunction_->ScaleAndRotateCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
  //          // 2. compute e from d
  //          // e = d*c
  //          rotatedCouplTensor.Mult(elastTensor_,tmp);
  //        } else {
  //          // use e-form as basis, i.e. the following step has to be done
  //          // 1. acale and rotate e
  //          hystCoefFunction_->ScaleAndRotateCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
  //          tmp = rotatedCouplTensor;
  //        }
  //        
  //        if(deltaFormActive_ && (deltaForm_ != 0) ) {
  //          //        std::cout << "Compute DeltaMatrix" << std::endl;
  //          bool useStrains = true;
  //          
  //          // deltaMat will be computed using the current value (timelevel_cur = 0)
  //          // and the value at timelevel (-1 > last ts; +1 > last iteration)
  //          Matrix<Double> deltaMat = hystCoefFunction_->GetDeltaMat(lpm, timelevel_cur, timelevel_for_diff, useStrains, useAbs,implementationVersion);
  //          Matrix<Double> tmp2;
  //          // tmp2 = dS/dE*c
  //          deltaMat.Mult(elastTensor_,tmp2);
  //          
  //          // tmp = e_scaled_rotated + dS/dE*c
  //          tmp.Add(1.0,tmp2);
  //        }        
  //      } else {
  //        std::cout << "Tensor " << tensorName_ << " requested, but not yet implemented " << std::endl;
  //      }
  //
  //      if(transposed_){
  //        std::cout << "Perform transpose" << std::endl;
  //        tmp.Transpose(outputTensor);
  //        std::cout << "Transposed tensor: " << tmp.ToString() << std::endl;
  //      } else {
  //        outputTensor = tmp;
  //      }
  //      std::cout << "Return the following tensor: " << outputTensor.ToString() << std::endl;
  //    }
  //    
  //    //! Return row and columns size of tensor if coefficient function is a tensor
  //    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
  //      if( (tensorName_ == "Permittivity") || (tensorName_ == "Reluctivity")){
  //        ptrFieldTensor_->GetTensorSize(numRows,numCols);
  //      } else if ( (tensorName_ == "CouplingMechToElec") || (tensorName_ == "CouplingElecToMech") 
  //        || (tensorName_ == "CouplingMechToMag") || (tensorName_ == "ComputeMagToMech") ){
  //        ptrCouplTensor_->GetTensorSize(numRows,numCols);
  //      } else {
  //        EXCEPTION("Tensor type unknown");
  //      }
  //      if(transposed_){
  //        UInt tmp = numRows;
  //        numRows = numCols;
  //        numCols = tmp;
  //      }
  //    }
  //    
  //  private:
  //    
  //    PtrCoefFct ptrFieldTensor_;
  //    PtrCoefFct ptrElastTensor_;
  //    PtrCoefFct ptrCouplTensor_;
  //    bool isCoupled_;
  //    
  //    PtrCoefFct hystCoefFunction_;
  //    
  //    std::string tensorName_;
  //    bool transposed_;
  //    
  //    Matrix<Double> fieldTensor_ ;
  //    Matrix<Double> elastTensor_;
  //    Matrix<Double> couplTensor_;
  //    
  //    bool tensorsInitialized_;
  //    
  //    /*
  //     * -1 : not coupled at all
  //     *  0 : coupled e-form/h-form (piezo/magstrict)
  //     *  1 : coupled d-form (piezo)
  //     *  2 : coupled g-form (magstrict)
  //     */
  //    int strainForm_;
  //  };
  
  //  class CoefFunctionHystRHSLoad : public CoefFunction{
  //    
  //  public:
  //    
  //    CoefFunctionHystRHSLoad(PtrCoefFct fieldTensor,PtrCoefFct elasticityTensor,
  //      PtrCoefFct couplingTensor, PtrCoefFct CoefFncHyst,
  //      std::string vectorName, bool transposed, bool onBoundary)
  //    :CoefFunction(){
  //      
  //      dimType_ = VECTOR;
  //      isAnalytic_ = false;
  //      isComplex_  = false;
  //      dependType_ = SOLUTION;
  //      
  //      fieldTensor_ = fieldTensor;
  //      elastTensor_ = elasticityTensor;
  //      couplTensor_ = couplingTensor;
  //      
  //      UInt couplRows, couplCols;
  //      if(couplTensor_ != NULL){
  //        couplTensor_->GetTensorSize(couplRows,couplCols);
  //        if(couplCols < couplRows){
  //          EXCEPTION("Coupling tensor has to be in non-transposed form");
  //        }
  //      }
  //      
  //      hystCoefFunction_ = CoefFncHyst;
  //      
  //      vectorName_ = vectorName;
  //      transposed_ = transposed;
  //      onBoundary_ = onBoundary;
  //    }
  //    
  //    virtual ~CoefFunctionHystRHSLoad(){};
  //    
  //    void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
  //      EXCEPTION("CoefFunctionHystRHSLoad only returns vectors")
  //    }
  //    
  //    void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){
  //      std::cout << "Coef Function Hyst RHS Load - Get Vector" << std::endl;
  //      // return vector that can be put onto rhs of system; 
  //      // the sign of term shall be such, that it can be added with + in the
  //      // corresponding pdes
  //      // Example: in piezo systems we have (among others) -int_Omega (BN)^T*P dOmega 
  //      //          on rhs, so we would return -P here
  //      
  //      int timeLevel;
  //      Double specificSign;
  //      Double baseSign;
  //      
  //      if(onBoundary_){
  //        // boundary terms usually are evluated with the current value of polarization etc
  //        // whereas the rhs depends on usage of deltaMat or not
  //        timeLevel = hystCoefFunction_->GetTimeLevel("BC");
  //        // furtherrmore, the  boundary terms and the (volume) loads have opposite
  //        // aigns
  //        baseSign = -1.0;
  //      } else {
  //        timeLevel = hystCoefFunction_->GetTimeLevel("RHS");
  //        baseSign = 1.0;
  //      }
  //      std::cout << "Timelevel: " << timeLevel << " (0 = current, 1 = last it, -1 = last ts)" << std::endl;
  //      //      if(onBoundary_){
  //      //        if(timeLevelRHS == -1){
  //      //          //        std::cout << "Get BC wrt to the last ts" << std::endl;
  //      //        } else if (timeLevelRHS != 0){
  //      //          //        std::cout << "Get BC wrt to the last "<< timeLevelRHS << "th it" << std::endl;
  //      //        } else {
  //      //          //        std::cout << "Put current value of hyst to BC" << std::endl;
  //      //        }
  //      //      } else {
  //      //        if(timeLevelRHS == -1){
  //      //          //        std::cout << "Get RHS Load wrt to the last ts" << std::endl;
  //      //        } else if (timeLevelRHS != 0){
  //      //          //        std::cout << "Get RHS Load wrt to the last "<< timeLevelRHS << "th it" << std::endl;
  //      //        } else {
  //      //          //        std::cout << "Put current value of hyst to RHS" << std::endl;
  //      //        }
  //      //      }
  //      
  //      if(vectorName_ == "ElecPolarization"){
  //        std::cout << "ElecPolarization was requested" << std::endl;
  //        if(onBoundary_){
  //          std::cout << "ON BOUNDARY" << std::endl;
  //        }
  //        // rhs (electrostatics, w = testfunction): 
  //        // + int_Volume nabla(w)^T P dOmega
  //        // - int_Surface w P^T n dGamma
  //        // > specificSign = +
  //        specificSign = 1.0;
  //        
  //        outputVector = Vector<Double>(GetVecSize());
  //        outputVector.Init();
  //        outputVector = hystCoefFunction_->GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevel);
  //        
  //        outputVector.ScalarMult(specificSign*baseSign);
  //        
  //      } else if (vectorName_ == "PiezoLoadForMechPDE"){
  //        //      std::cout << "PiezoLoadForMechPDE was requested" << std::endl;
  //        //TODO: implement correct vector function; at the moment, return simply zero
  //        //      std::cout << couplTensor_->ToString() << std::endl;
  //        outputVector = Vector<Double>(GetVecSize());
  //        outputVector.Init();
  //        //      std::cout << "outputVector: " << outputVector.ToString() << std::endl;
  //      } else if (vectorName_ == "PiezoLoadForElecPDE"){
  //        //      std::cout << "PiezoLoadForElecPDE was requested" << std::endl;
  //        //TODO: implement correct vector function; at the moment, return simply zero
  //        //      std::cout << couplTensor_->ToString() << std::endl;
  //        outputVector = Vector<Double>(GetVecSize());
  //        outputVector.Init();
  //        //      std::cout << "outputVector: " << outputVector.ToString() << std::endl;
  //      } else if(vectorName_ == "MagPolarization"){
  //        std::cout << "MagPolarization was requested" << std::endl;
  //        // rhs (magnetics, w = testfunction): 
  //        // + int_Volume rot(w)^T M dOmega
  //        // > specificSign = +
  //        specificSign = 1.0;
  //        
  //        outputVector = Vector<Double>(GetVecSize());
  //        outputVector.Init();
  //        outputVector = hystCoefFunction_->GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevel);
  //        
  //        outputVector.ScalarMult(specificSign*baseSign);
  //      }
  //      else {
  //        EXCEPTION("VectorName not implemented yet");
  //      }
  //      
  //      std::cout << "Coef Function Hyst RHS Load - return: " << std::endl;
  //      std::cout << outputVector.ToString() << std::endl;
  //    }
  //    
  //    virtual UInt GetVecSize() const {
  //      
  //      UInt numRows,numCols;
  //      if( (vectorName_ == "ElecPolarization") || (vectorName_ == "MagPolarization")){
  //        fieldTensor_->GetTensorSize(numRows,numCols);
  //      } else if ( (vectorName_ == "PiezoLoadForMechPDE") || (vectorName_ == "PiezoLoadForElecPDE") 
  //        || (vectorName_ == "MagStrictLoadForMechPDE") || (vectorName_ == "MagStrictLoadForMagPDE") ){
  //        couplTensor_->GetTensorSize(numRows,numCols);
  //      } else {
  //        EXCEPTION("Vector type unknown");
  //      }
  //      if(transposed_){
  //        return numCols;
  //      } else {
  //        return numRows;
  //      }
  //    }
  //    
  //    void GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm ){
  //      EXCEPTION("CoefFunctionHystRHSLoad only returns vectors")
  //    }
  //    
  //    //! Return row and columns size of tensor if coefficient function is a tensor
  //    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
  //      EXCEPTION("GetVecSize not used here")
  //    }
  //    
  //  private:
  //    
  //    PtrCoefFct fieldTensor_;
  //    PtrCoefFct elastTensor_;
  //    PtrCoefFct couplTensor_;
  //    
  //    PtrCoefFct hystCoefFunction_;
  //    
  //    std::string vectorName_;
  //    bool transposed_;
  //    bool onBoundary_;
  //  };
  //  
  
  class CoefFunctionHystOutput : public CoefFunction {
    
  public:
    
    CoefFunctionHystOutput(CoefFunctionHelper* hystHelper ,std::string ResultName):CoefFunction(){
      
      if( (ResultName == "DeltaPermeability")||(ResultName == "DeltaPermittivity") ){
        dimType_ = TENSOR;
      } else {
        dimType_ = VECTOR;
      }
      
      isAnalytic_ = false;
      isComplex_  = false;
      dependType_ = SOLUTION;
      onBoundary_ = false;
      
      hystHelper_ = hystHelper;
      resultName_ = ResultName;
    }
    
    virtual ~CoefFunctionHystOutput(){};
    
    void GetScalar(Double& outputScalar, const LocPointMapped& lpm ){
      EXCEPTION("CoefFunctionHystOutput returns no scalars")
    }
    
    void GetVector(Vector<Double>& outputVector,const LocPointMapped& lpm){

      int timeLevel = 0;
      Double baseSign = 1.0;
      //std::cout << "result " << resultName_ << " requested" << std::endl;
      hystHelper_->ComputeVector(outputVector,lpm, timeLevel, baseSign, resultName_, onBoundary_ );
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
      if( (resultName_ == "DeltaPermeability")||(resultName_ == "DeltaPermittivity") ){
        //TODO: estimate permeability or permittivity around current working point
        UInt numCols,numRows;
        hystHelper_->GetTensorSize(numRows,numCols);
        outputTensor.Resize(numRows,numCols);
        outputTensor.Init();
      } else {
        EXCEPTION("Unknown result")
      }
    }
    
    virtual UInt GetVecSize() const {
      return hystHelper_->GetVecSize();
    }
    
    //! Return row and columns size of tensor if coefficient function is a tensor
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      hystHelper_->GetTensorSize(numRows,numCols);
    }
    
    std::string ToString() const {
      return "Dummy string";
    }
    
  private:
    
    CoefFunctionHelper* hystHelper_;
    std::string resultName_;
    bool onBoundary_;
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
      shared_ptr<FeSpace> ptFeSpace);
    
    //! Constructor
    CoefFunctionHyst( BaseMaterial* const material,
    shared_ptr<ElemList> actSDList,
    PtrCoefFct dependency1,
    PtrCoefFct dependency2,
    SubTensorType tensorType,
    MaterialType matType,
    shared_ptr<FeSpace> ptFeSpace);
    
    void Init(BaseMaterial* const material,
    shared_ptr<ElemList> actSDList,
    PtrCoefFct dependency1,
    SubTensorType tensorType,
    MaterialType matType,
    shared_ptr<FeSpace> ptFeSpace);
    
    /*
     * New functions added April 2018
     */
    Double getWeight(Double alpha, Double beta, Double delta){
      
      if(weightType_ == "Constant"){
        return constWeight_;
      } else if(weightType_ == "muDat"){
        return MuDat(muDat_A_, muDat_sigma1_, muDat_h1_, muDat_eta_, alpha, beta);
      } else if(weightType_ == "muDatExtended"){
        return MuDatExtended(muDat_A_, muDat_sigma1_, muDat_sigma2_, muDat_h1_, muDat_h2_, muDat_eta_, alpha, beta);
      } else if(weightType_ == "givenTensor"){
        // alpha = -1 shall be mapped to 0, alpha = 1.0 shall be mapped to maxNum
        UInt idxAlpha = UInt(std::round((alpha+1.0)/delta));
        UInt idxBeta = UInt(std::round((beta+1.0)/delta));
        if( idxAlpha < 0 ){
          idxAlpha = 0;
        } else if (idxAlpha >= MAT_numRows_){
          idxAlpha = MAT_numRows_-1;
        }
        if( idxBeta < 0 ){
          idxBeta = 0;
        } else if (idxBeta >= MAT_numRows_){
          idxBeta = MAT_numRows_-1;
        }
        
        return weightTensor_[idxAlpha][idxBeta];
      } else {
        return 0.0;
      }
    }
    
    Double getWeightDerivative(Double alpha, Double beta, Double delta, bool flipped){
      
      Double s, lambda;
      if(flipped){
        s = alpha;
        lambda = beta/alpha;
      } else {
        s = beta;
        lambda = alpha/beta;
      }
      
      if(weightType_ == "Constant"){
        return 0.0;
      } else if(weightType_ == "muDat"){
        return dMuDat_by_ds(muDat_A_, muDat_sigma1_, muDat_h1_, muDat_eta_,s,lambda,flipped);
      } else if(weightType_ == "muDatExtended"){
        return dMuDatExtended_by_ds(muDat_A_, muDat_sigma1_, muDat_sigma2_, muDat_h1_, muDat_h2_, muDat_eta_,s,lambda,flipped);
      } else if(weightType_ == "givenTensor"){
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
            wUp = getWeight(lambda*s_up,s_up,delta);
            wLow = getWeight(lambda*s_low,s_low,delta);
        } else{
            wUp = getWeight(s_up,lambda*s_up,delta);
            wLow = getWeight(s_low,lambda*s_low,delta);
        }
        return (wUp-wLow)/deltas;
      } else {
        return 0.0;
      }
    }
    
    std::string weightType_;
    Double constWeight_;
    Double muDat_A_;
    Double muDat_sigma1_;
    Double muDat_sigma2_;
    Double muDat_h1_;
    Double muDat_h2_;
    Double muDat_eta_;
    Double MAT_anhysteretic_a_;
    Double MAT_anhysteretic_b_;
    Double MAT_anhysteretic_c_;
    Matrix<Double> weightTensor_;
    std::string usedHystModel_;
    bool anhystOnly_;
    
    Double MuDat(Double A, Double sigma, Double h, Double eta, Double alpha, Double beta){
      // MuDat function for evaluating Preisach Weights
      // weights directly usable for 
      //  > scalar Preisach model
      //  > vector Preisach model (Sutor)
      // after transformation, weights usable for
      //  > vector Preisach model (Mayergoyz) ISOTROPIC materials
      //
      // Source: "A Preisach-based hysteresis model for magnetic and ferroelectric hysteresis" - A. Sutor 2010
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
    
    Matrix<Double> evaluatePreisachWeights(){
      Double alpha,beta;
      Double delta = 2.0/MAT_numRows_;
      Double dimHalf = MAT_numRows_/2.0;
      
      if(weightType_ == "givenTensor"){
        return weightTensor_;
      }
      
      Matrix<Double> weights = Matrix<Double>(MAT_numRows_,MAT_numRows_);
      
      for(UInt i = 0; i < MAT_numRows_; i++){
        // note: we evaluated the weights always at the element center > add deltaS/2
        alpha = (i - dimHalf)*delta + delta/2;
        
        for(UInt k = 0; k <= i; k++){
          beta = (k - dimHalf)*delta + delta/2;
          
          // note: getWeight checks for weightType_; depending on that value
          // we evaluate muDat, muDatExt, return a const or access the already given weights (in which case
          // there is no need to call this function)
          weights[i][k] = getWeight(alpha,beta,delta);
        }
      }
      return weights;
    }
    
    
    Matrix<Double> transformPreisachWeightsForIsotropicVectorCase(){
      // Transform Preisach weights from scalar case to vector case (isotropic, Mayergoyz model)
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
      if(dim_ != 2){
        EXCEPTION("Currently only 2D case of Preisach Weight transformation supported");
        // 3d case starts from different formula > see "Mathematical Models of Hysteresis and their Application" - Mayergoyz 2003
      }
      Matrix<Double> vectorWeights = Matrix<Double>(MAT_numRows_,MAT_numRows_);
      
      Double dimHalf = MAT_numRows_/2.0;
      UInt dimHalfInt = UInt(dimHalf);
      Double alpha,beta,s,lambda;
      Double delta = 2.0/MAT_numRows_;
      
      UInt numIntegrationPoints = 15;
      Matrix<Double> integrationPoints = getTschebyscheffPointsAndWeights(numIntegrationPoints);

      // upper part
      for(UInt i = dimHalfInt; i < MAT_numRows_; i++){
        // note: we evaluated the weights always at the element center > add deltaS/2
        alpha = (i - dimHalf)*delta + 1e-6*delta; 
        for(UInt k = MAT_numRows_-1-i; k <= i; k++){
          beta = (k - dimHalf)*delta + 1e-6*delta;// + delta/2;
          
          if(alpha == 0){
            vectorWeights[i][k] = 0.75*getWeight(alpha, beta, delta);
          } else {
            lambda = beta/alpha;
            vectorWeights[i][k] = 0.0;
            // numerically integrate from 0 to alpha
            for(UInt pp = 0; pp < numIntegrationPoints; pp++){
              Double curY = integrationPoints[pp][0];
              Double curW = integrationPoints[pp][1];
              
              s = x_of_y(curY,0,alpha);
              Double tmp = 3*s*s*getWeight(s,lambda*s,delta);
              
              bool flipped = false; // > s takes slot of alpha, lambda*s takes place of beta
              Double tmp2 = s*s*s*getWeightDerivative(s,lambda*s,delta, flipped);
              Double tmp3 = (tmp + tmp2)/(M_PI * alpha * alpha * std::sqrt(alpha*alpha - s*s));

              vectorWeights[i][k] += dx_by_dy(std::min(0.0,alpha), std::max(0.0,alpha))*curW*tmp3;
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
        beta = (k - dimHalf)*delta + 1e-6*delta;// + delta/2;
        
        for(UInt i = k; i < MAT_numRows_-k-1; i++){
//          std::cout << "MAT_numRows_-k-1 = " << MAT_numRows_-k-1 << std::endl;
//          std::cout << "i = " << i << std::endl;
          alpha = (i - dimHalf)*delta + 1e-6*delta;// + delta/2;
          
          if(beta == 0){
            vectorWeights[i][k] = 0.75*getWeight(alpha, beta, delta);
          } else {
            // flipped case
            lambda = alpha/beta;
            vectorWeights[i][k] = 0.0;
            // numerically integrate from 0 to beta > however, beta < 0 so we have to reverse the sign
            // > done vis dx_by_dy
            for(UInt pp = 0; pp < numIntegrationPoints; pp++){
              Double curY = integrationPoints[pp][0];
              Double curW = integrationPoints[pp][1];
              
              s = x_of_y(curY,0,beta);
              Double tmp = 3*s*s*getWeight(lambda*s,s,delta);
              
              bool flipped = true; // > s takes slot of alpha, lambda*s takes place of beta
              Double tmp2 = s*s*s*getWeightDerivative(lambda*s,s,delta, flipped);
              
              Double tmp3 = (tmp + tmp2)/(M_PI * beta * beta * std::sqrt(beta*beta - s*s));
              vectorWeights[i][k] += dx_by_dy(std::min(0.0,beta), std::max(0.0,beta))*curW*tmp3;
              //vectorWeights[i][k] += dx_by_dy(0.0,beta)*curW*tmp3;
            }
          }
        } 
      }
      
//      std::cout << "weights final: " << std::endl;
//      std::cout << vectorWeights.ToString() << std::endl;
      
      bool checkWeights = !true;
      if(checkWeights){
        // according to "Analysis of Isotropic Materials with Vector Hysteresis" - O. Bottauscio 1998
        // int_-pi/2^pi/2 cos(phi)^3 vectorWeights(alpha*cos(phi),beta*cos(phi) dphi != scalarWeight(alpha,beta)
        // this can be checked
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
          Double target = getWeight(alpha, beta, delta);
          for(UInt i = 0; i < numAngles; i++){
            phi = -M_PI/2+i*deltaPhi;
            c = std::cos(phi);
            UInt idxAlpha = UInt(std::round((alpha*c+1.0)/delta));
            UInt idxBeta = UInt(std::round((beta*c+1.0)/delta));
            integral += c*c*c*vectorWeights[idxAlpha][idxBeta];
          }
          //Double deltaAngle = M_PI/numDirections_;
          // Note: here we do not average out over the halfspace so not multiplication by 2/M_PI
          integral *= M_PI/numAngles;

          std::cout << "Alpha = " << alpha << std::endl;
          std::cout << "Beta = " << beta << std::endl;
          std::cout << "Scalar weight: " << target << std::endl;
          std::cout << "int_-pi/2^pi/2 cos(phi)^3 vectorWeights(alpha*cos(phi),beta*cos(phi) dphi: " << integral << std::endl;
          }
        
        EXCEPTION("End after check");
      }

      return vectorWeights;
    }

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
      numRows = MAT_eps_mu_SmallSignal_.GetNumRows();
      numCols = MAT_eps_mu_SmallSignal_.GetNumCols();
    }
    
    void SetFlag(std::string flagName, Integer intState);
    
    void SetCurrentHystVals(bool overwriteMemory);
    void SetPreviousHystVals(bool setLastTS);
    void SetPreviousHystValsOLD(bool setLastTS, bool forceMemoryLock = false);
    
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
        satValue = MAT_ySat_ + MAT_eps_mu_SmallSignal_.NormL2()*MAT_xSat_;
      } else {
        satValue = MAT_xSat_;
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

    Matrix<Double> EstimateSlope(const LocPointMapped& Originallpm, bool includeStrains, bool useAbs, std::string implementationVersion,
        Double steppingLength=1e-10, Double scaling=1.0);
    
    Matrix<Double> CalcDeltaMat(Vector<Double> E_B_diff, Vector<Double> P_J_diff, bool useAbs,std::string implementationVersion, Double cuttingTol = -1.0);
    Matrix<Double> CalcDeltaMatStrains(Vector<Double> E_B_diff, Vector<Double> S_diff, bool useAbs,std::string implementationVersion, Double cuttingTol = -1.0);
    
    Matrix<Double> GetDeltaMat(const LocPointMapped& Originallpm, int timelevel1, int timelevel2, bool useStrains, bool useAbs, std::string implementationVersion );
    
    Vector<Double> GetIrreversibleStrains(const LocPointMapped& Originallpm, int timeLevel);
    
    Vector<Double> GetPrecomputedInputToHysteresisOperator(const LocPointMapped& Originallpm, int timeLevel);

    Vector<Double> GetPrecomputedOutputOfHysteresisOperator(const LocPointMapped& Originallpm, int timeLevel);
    
    Vector<Double> RetrieveInputToHysteresisOperator(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool onBoundary);
    
    Vector<Double> RetrieveLPMSolution(LocPointMapped& actualLPM, UInt storageIdx, int timeLevel, bool onBoundary);
    
    Double GetOutputSaturation(){
      return MAT_ySat_;
    }
    
    PtrCoefFct GenerateMatCoefFnc(std::string tensorName){
      
      PtrCoefFct ret;
      if((tensorName == "Permittivity")||(tensorName == "CouplingMechToElec")||(tensorName == "CouplingElecToMech")){
        // Note: for Piezo we need two different coupling functions
        // Reason: in case of deltaMat formulation, the coupling matrices are not 
        //         transposed of each other as the mech to elec term has deltaS/deltaE
        //         in addition to the rotated coupling operator
        PtrCoefFct eps = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType_,Global::REAL);
        
        if(hystHelper_ == NULL){
          hystHelper_ = new CoefFunctionHelper(eps,elastTensorFct_,
          couplTensorFct_, hystItself_);
        }
        
        bool transposed = false;
        if(tensorName == "CouplingElecToMech"){
          // in mech PDE we need e^T or d^T
          transposed = true;
        }
        shared_ptr<CoefFunctionHystMatTensor> c(new CoefFunctionHystMatTensor(hystHelper_,tensorName,transposed));
        
        ret = c;
      } else if((tensorName == "Reluctivity")||(tensorName == "CouplingMechToMag")||(tensorName == "CouplingMagToMech")){
        PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY,tensorType_,Global::REAL);
        
        if(hystHelper_ == NULL){
          hystHelper_ = new CoefFunctionHelper(nu,elastTensorFct_,
          couplTensorFct_, hystItself_);
        }
        
        bool transposed = false;
        if(tensorName == "CouplingMagToMech"){
          // in mech PDE we need h^T or g^T
          transposed = true;
        }
        shared_ptr<CoefFunctionHystMatTensor> c(new CoefFunctionHystMatTensor(hystHelper_,tensorName,transposed));
        
        ret = c;
      } else {
        EXCEPTION("tensorName not implemented yet");
      }
      
      return ret;
    }
    
    PtrCoefFct GenerateRHSCoefFnc(std::string vectorName, bool onBoundary = false ){
      PtrCoefFct ret;
      if( (vectorName == "ElecPolarization") || (vectorName == "PiezoLoadForMechPDE") || (vectorName == "PiezoLoadForElecPDE")){
        PtrCoefFct eps = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType_,Global::REAL);
        
        if(hystHelper_ == NULL){
          hystHelper_ = new CoefFunctionHelper(eps,elastTensorFct_,
          couplTensorFct_, hystItself_);
        }
        
        bool transposed = false;
        if(vectorName == "PiezoLoadForMechPDE"){
          // in mech PDE we need h^T or g^T
          transposed = true;
        }
        shared_ptr<CoefFunctionHystRHSLoad> c(new CoefFunctionHystRHSLoad(hystHelper_,vectorName,transposed,onBoundary));
        
        ret = c;
      } else if( (vectorName == "MagMagnetization") || (vectorName == "MagStrictLoadForMechPDE") || (vectorName == "PMagStrictLoadForMagPDE")){
        PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY,tensorType_,Global::REAL);
        
        if(hystHelper_ == NULL){
          hystHelper_ = new CoefFunctionHelper(nu,elastTensorFct_,
          couplTensorFct_, hystItself_);
        }
        
        bool transposed = false;
        if(vectorName == "MagStrictLoadForMechPDE"){
          // in mech PDE we need h^T or g^T
          transposed = true;
        }
        shared_ptr<CoefFunctionHystRHSLoad> c(new CoefFunctionHystRHSLoad(hystHelper_,vectorName,transposed,onBoundary));
        
        ret = c;
      } else {
        //std::cout << "RHSCoefFunction for " << vectorName << " requested." << std::endl;
        EXCEPTION("vectorName not implemented yet");
      }
      
      return ret;
    }
    
    PtrCoefFct GenerateOutputCoefFnc(std::string ResultName){
      
      PtrCoefFct ret;
      
      if(ResultName == "ElecPolarization") {
        PtrCoefFct eps = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType_,Global::REAL);
        
        if(hystHelper_ == NULL){
          hystHelper_ = new CoefFunctionHelper(eps,elastTensorFct_,
          couplTensorFct_, hystItself_);
        }
        
        shared_ptr<CoefFunctionHystOutput> c(new CoefFunctionHystOutput(hystHelper_,ResultName));
        ret = c;
      } else if ( (ResultName == "MagPolarization") || (ResultName == "MagMagnetization") || (ResultName == "MagFieldIntensityHyst") ){
        PtrCoefFct nu = material_->GetTensorCoefFnc( MAG_RELUCTIVITY,tensorType_,Global::REAL);
        
        if(hystHelper_ == NULL){
          hystHelper_ = new CoefFunctionHelper(nu,elastTensorFct_,
          couplTensorFct_, hystItself_);
        }
        
        shared_ptr<CoefFunctionHystOutput> c(new CoefFunctionHystOutput(hystHelper_,ResultName));
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
      return MAT_useStrainForm_;
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
      if(Type == "Material"){
        return timeLevel_Material_;
      } else if (Type == "DeltaMat"){
        return timeLevel_DeltaMat_;
      } else if (Type == "RHS"){
        return timeLevel_RHS_;
      } else if (Type == "BC"){
        return timeLevel_BC_;
      } else if (Type == "Output"){
        return timeLevel_Output_;
      } else {
        EXCEPTION("Type not known");
        return -2;
      }
    }
    
    void TestInversion(UInt testNumber, bool stopAfterTests, bool printStatistics, bool writeResultsToFiles);
    
    void SetElastAndCouplTensor(PtrCoefFct elastTensor, PtrCoefFct couplTensor){
      //std::cout << "Elast and Coupl tensor were passed by coupled pde!" << std::endl;
      elastTensorFct_ = elastTensor;
      couplTensorFct_ = couplTensor;
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
    void AddAdditionalSDList(shared_ptr<EntityList> actSDList, bool isSurface);
    
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
        if(!matrixForInversionWasSet_[i]){
          return true;
        }
      }
      return false;
    }
    
    void getMatrixForLocalInversion(const LocPointMapped& Originallpm, Matrix<Double>& matrixForInversion, Matrix<Double>& matrixForInversionInverted){
      UInt operatorIdx, storageIdx;
      LocPointMapped actualLPM;
      
      PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
      matrixForInversion = matrixForInversion_[storageIdx];
      matrixForInversionInverted = matrixForInversionInverted_[storageIdx];
    }
    
    void getLPMMaps(std::map<UInt, LocPointMapped >& allLPM, std::map<UInt, LocPointMapped >& midpointLPM){
      allLPM = allLPMmap_;
      midpointLPM = midpointLPMmap_;
    }
    
    void setMatrixForLocalInversion(Matrix<Double> matrixForInversion, Matrix<Double> matrixForInversionInverted, UInt storageIdx, bool reuse){
      if(reuse){
        matrixForInversionWasSet_[storageIdx] = true;
      } else {
        // matrixForInversionWasSet_ shall be set only once in each iteration
        // > use flag which gets reset once during setprevioushystvals function
        if(matrixForInversionWasSet_[storageIdx] == false){
          matrixForInversion_[storageIdx] = matrixForInversion;
          matrixForInversionInverted_[storageIdx] = matrixForInversionInverted;
          matrixForInversionWasSet_[storageIdx] = true;
        } 
      }
    }
     
    void ClipDirection(Vector<Double>& targetVector);
    
    void ScaleAndRotateCouplingTensor(const LocPointMapped& lpm, Matrix<Double>& couplTensor, Matrix<Double>& rotatedCouplTensor,int timeLevel,
      bool rotate = true){
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
        numCols = couplTensor.GetNumCols();
        numRows = couplTensor.GetNumRows();
        //couplTensorCoefFnc->GetTensorSize(numRows,numCols);
        //Matrix<Double> couplTensor = Matrix<Double>(numRows,numCols);
        //couplTensorCoefFnc->GetTensor(couplTensor,actualLPM);
        
        //      std::cout << "Retrieved coupling tensor: " << couplTensor.ToString() << std::endl;
        //      std::cout << "Rotate coupling tensor" << std::endl;
        rotatedCouplTensor.Resize(numRows,numCols);
        Matrix<Double> scaledCouplTensor = Matrix<Double>(numRows,numCols);
        
        // get current polarization (elec or mag)
        Vector<Double> P = GetPrecomputedOutputOfHysteresisOperator(lpm, timeLevel);
        
        //      std::cout << "Current polarization vector " << P.ToString() << std::endl;
        //      
        // calculate scaling
        assert(MAT_ySat_ != 0);
        Double scaling = P.NormL2()/MAT_ySat_;
        
        scaledCouplTensor = couplTensor* scaling;
        
        //      std::cout << "Scalaed coupling tensor " << scaledCouplTensor.ToString() << std::endl;
        //      
        if(rotate == true){
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
              // 2d plane strain or plane stress case
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
              
              assert(rotatedCouplTensor.GetNumRows() == 2);
              assert(rotatedCouplTensor.GetNumCols() == 3);
              
              scaledCouplTensor.PerformRotation(R,rotatedCouplTensor);
              
            }
          }
        } else {
          // no rotation > take scale tensor
          rotatedCouplTensor = scaledCouplTensor;
        }
        
        rotatedCouplingTensor_requiresReeval_[storageIdx] = false;
        rotatedCouplingTensor_[storageIdx] = rotatedCouplTensor;
        
      }
    }
    
    CoefFunctionHelper* hystHelper_;
  
  private:
        
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
    

    void EvaluateHystOperatorsInt(Integer intFlag);
    void EvaluateHystOperators(bool setMatForInversion, bool resetSolToZeroFirst, bool overwriteMemory);
    
    bool PreprocessLPM(const LocPointMapped& lpmInput, LocPointMapped& lpmOutput,
    UInt& operatorIdx, UInt& storageIdx, bool forceMidpoint = false);
    
    Vector<Double> CalcOutputOfHysteresisOperator(Vector<Double> inpute, UInt operatorIdx, UInt storageIdx,
    bool forceMemoryLock = false, bool forceMemoryWrite = false);
    
    void ForceRemanence();
    void InitStorage();
    void InitLPMMaps();
    
    void CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& outputTensor, std::string evalMethod,
    UInt storageIdx, bool intoSat, bool outofSat, bool satToSat, Vector<Double>& X_current);
    
    void ExtractSolutionAndInputForHystOperator(Vector<Double>& extractedSolution, Vector<Double>& extractedInput, UInt& operatorIndex,
    UInt& storageIndex, const LocPointMapped& lpm, bool midpointOnly);
        
    bool EvaluateAtMidpointOnly();
    bool OverwriteHystMemory();
    
    // for usage with coefFunctionHystMat
    bool deltaMatActive_;
    // if set to true, timelevel for evaluation will be at the current timestep
    // independent of the actual value of timeLevel_xxx_
    bool forceCurrentTS_;
    int RUN_deltaForm_;
    Integer timeLevel_Material_;
    Integer timeLevel_DeltaMat_;
    Integer timeLevel_RHS_;
    Integer timeLevel_BC_;
    Integer timeLevel_Output_;
    PtrCoefFct hystItself_;
    
    
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
    /*
     * Initial/small signal material tensor (permittivity / reluctivity)
     */
    //Matrix<Double> MAT_initialTensor_; // small signal tensor of permittivity / reluctivity > no longer used
    Matrix<Double> MAT_eps_mu_SmallSignal_; // small signal tensor of permittivity / permeability from mat file
    //Matrix<Double> MAT_freeFieldTensor_; // eps0, nu0 > no longer used
    
    /*
     * Preisach weights and its size in form of the number of rows
     */
    Matrix<Double> MAT_PreisachWeights_;
    UInt MAT_numRows_;
    
    /* 
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
     * NEW: instead of specifying only x,y,z we allow to specify a direction in the mat.xml file
     */
    Vector<Double> MAT_dirP_;
    bool MAT_useExtension_;
    
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
     */
    bool RUN_overwriteDirection_;
    
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
    
    //! material object
    BaseMaterial* material_;
    
    bool MAT_useStrainForm_;
    int MAT_dim_beta_;
    Matrix<Double> MAT_betaCoefs_;
    
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
    Vector<Double>* E_H_;
    Vector<Double>* Si_;
    
    /*
     * Storage vectors for quantities of the last iteration;
     * these values will ONLY be overwritten during each iteration by calling the
     * function
     *    SetPreviousHystValues(false)
     *
     */
    Vector<Double>* E_B_lastIt_;
    Vector<Double>* P_J_lastIt_;
    Vector<Double>* E_H_lastIt_;
    Vector<Double>* Si_lastIt_;
        
    /*
     * Storage vectors for quantites of the last timestep;
     * these values will ONLY be set during the FIRST iteration of each time step
     * by calling
     *    SetPreviousHystValues(true)
     */
    Vector<Double>* E_B_lastTS_;
    Vector<Double>* P_J_lastTS_;
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
    
    Matrix<Double>* matrixForInversion_; // used for local inversion of preisach operator > eps/mu; might depend on P itself
    Matrix<Double>* matrixForInversionInverted_; // used for local inversion of preisach operator > eps/mu; might depend on P itself
    bool* matrixForInversionWasSet_;
    
    // needed to get integration point
    shared_ptr<FeSpace> ptFeSpace_;
    UInt numIntegrationPoints_;
    std::map<int,UInt> locPointIndices_;
    
    // flags indicating if we need inversion (for e.g. for magnetics)
    // and if we have an inverse model for that
    bool needsInversion_;
    bool hasInverseModel_;
    // For inversion via Levenberg Marquardt
    UInt INV_maxIter_;
    Double INV_resTolH_;
    Double INV_resTolB_;
    Double INV_jacobiResolution_;
    bool INV_useTikhonov_;
    Double INV_alphaLSStart_;
    Double INV_alphaLSMin_;
    Double INV_alphaLSMax_;
    
    // when setting the system into remanence for the case of magnetics
    // we do this by forcing the input to the hystoperator to become 0;
    // for electrostatics, we can simply set the solution vector to 0 but
    // for magnetics, this would cause B = 0 and thus bring the system to
    // coercive points; nevertheless we want to reset solVec to 0; to avoid
    // evaluating the hyst operator (after bringing it to remanence by hand)
    // we simply lock its evaluation (via setFlag)
    bool hystOperatorLocked_;
    
  };
  
}
#endif
