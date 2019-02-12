#include "CoefFunctionHyst.hh"

// classes for function / spline approximation
#include "Materials/Models/Preisach.hh"
#include "Materials/Models/SimplePreisachInv.hh"
#include "Materials/Models/VectorPreisachv10.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include <cmath>
#include "Utils/Timer.hh"

namespace CoupledField {
	DECLARE_LOG(coeffcthyst)
	DEFINE_LOG(coeffcthyst, "coeffcthyst")
  
	CoefFunctionHyst::CoefFunctionHyst(BaseMaterial * const material,
          shared_ptr<ElemList> actSDList,
          PtrCoefFct dependency1,
          SubTensorType tensorType, MaterialType matType, shared_ptr<FeSpace> ptFeSpace) 
  : CoefFunction() {
    
		dependCoef1Surf_ = NULL;
		Init(material, actSDList, dependency1, tensorType, matType, ptFeSpace);
	}
  
	CoefFunctionHyst::CoefFunctionHyst(BaseMaterial * const material,
          shared_ptr<ElemList> actSDList,
          PtrCoefFct dependency1, PtrCoefFct dependCoef1Surf,
          SubTensorType tensorType, MaterialType matType, shared_ptr<FeSpace> ptFeSpace) 
  : CoefFunction() {
    
		WARN("Currently we support only single-input hysteresis operators! The second dependency is only allowed if it is the surface version of the original dependency!")
		dependCoef1Surf_ = dependCoef1Surf;
		Init(material, actSDList, dependency1, tensorType, matType, ptFeSpace);
	}
  
	void CoefFunctionHyst::Init(BaseMaterial * const material,
          shared_ptr<ElemList> actSDList,
          PtrCoefFct dependency1,
          SubTensorType tensorType,
          MaterialType matType,
          shared_ptr<FeSpace> ptFeSpace) {
    
		// this type of coefficient is nonlinear (i.e. solution dependent)
		dimType_ = VECTOR;
		dependType_ = SOLUTION;
		isAnalytic_ = false;
		isComplex_ = false;
		dependCoef1_ = dependency1;
    storageInitialized_ = false;    
		tensorType_ = tensorType;
    
    if(tensorType_ == AXI){
      EXCEPTION("Coef functions hyst does only support plane 2d and full 3d models (at the moment)")
    }
    
		matType_ = matType;
		material_ = material;
		tol_ = 1e-14;
		ptFeSpace_ = ptFeSpace;
    
		/*
     * for performance measurement
     */
		timer_ = new Timer();
		totalCallingCounter_ = 0;
		totalEvaluationCounter_ = 0;
		avgEvaluationTime_ = 0.0;
		totalEvaluationTime_ = 0.0;
    
    /*
     * for linesearch in VectorPreisach
     */
    alphaLinesearch_ = -1.0;
    
		// for interaction with coefFunctionHystMat
		deltaMatActive_ = true;
		hystItself_ = PtrCoefFct(this);
    
    // per default the coef functtions for coupling are NULL
    // for coupled pdes, these functions have to be set BEFORE both coupled
    // and single field integrators get defined!
    elastTensorFct_ = NULL;
    couplTensorFct_ = NULL;
    
    /*
     * Default case: deltaMat towards last TS
     * > material tensors at current value (0)
     * > deltaMat towards last ts (-1)
     * > rhs hyst on last ts (-1)
     * > bc on current value (0)
     * > output on current value (0)
     */
    RUN_deltaForm_ = -1;
		timeLevel_Material_ = 0; // i.e. time step for which to evaluate material tensors (eps, nu, ...)
		timeLevel_DeltaMat_ = -1; // timestep to which delta is computed
		timeLevel_RHS_ = -1; // rhs term has to be last ts value as delta is towards that value
		timeLevel_BC_ = 0; // BC term has to use the last known, i.e. current value
		timeLevel_Output_ = 0; // output should always be the current value
    forceCurrentTS_ = false;
    
		/*
     * set initial values for runtime dependent parameter
     */
		// before anything can be computed, a system has to be assembled first -> 1
		RUN_evaluationPurpose_ = 1;

		// do not output switching state; that is very costly and should only be done
		// for debugging or special figure computation
		RUN_allowBMP_ = false;
    
		// is set to false, the hyst operator will keep its rotation state
		// if the initial state is unset, this initial rotation state will be 0 0
		// unless the flag is set to true, it will stay 0 0 and so the output of the
		// hyst operator will be 0 0, too
		// > at least during the first inputs, this flag should be true so that the
		//   rotation list gets initialized properly
		// > tests showed that this flag should always be true
		RUN_overwriteDirection_ = true;
		
		/*
     * set temporary values for xml dependent parameter
     */
		// currently, we have no direct handle of the xml input during the creation
		// of this object; the initialization of the XML dependent parameter is therefore
		// postponed to the first call of StdSolveStep (prior to that we do not need the
		// paramters either way)
		XMLParameterSet_ = false;
    
		XML_EvaluationDepth_ = 0;
		XML_performanceMeasurement_ = 0;
    XML_testInversion_ = 0;
    
		/*
     * set values for material dependent parameter
     */
		material->GetScalar(MAT_xSat_, X_SATURATION, Global::REAL);
		material->GetScalar(MAT_ySat_, Y_SATURATION, Global::REAL);
		material->GetTensor(MAT_PreisachWeights_, PREISACH_WEIGHTS, Global::REAL);
		MAT_numRows_ = MAT_PreisachWeights_.GetNumRows();
    
		/*
     * get elements and integration points
     */
		// set map: global to local element number
		EntityIterator it = actSDList->GetIterator();
		UInt iel = 0;
		UInt globalElNr;
		for (it.Begin(); !it.IsEnd(); it++, iel++) {
			globalElNr = it.GetElem()->elemNum;
			globalElem2Local_[globalElNr] = iel;
			globalElemOnSurf_[globalElNr] = false;
		}
    
		//store subdomain list of elements
		SDList_List_.push_back(actSDList);
		numElemSD_ = actSDList->GetSize();
    
		// pick out the first element (even though any of them would do as they share the
		// same material) and extract midpoint
		it.Begin();
		const Elem * el = it.GetElem();
		LocPoint lp = Elem::shapes[el->type].midPointCoord;
		LocPointMapped lpm;
		shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
		lpm.Set(lp, esm, 0.0);
    
    
		// get number of integration points
		IntegOrder order;
		IntScheme::IntegMethod method;
		StdVector<LocPoint> intPoints;
		StdVector<Double> weights;
    
		ptFeSpace_->GetFe(it, method, order);
		// store method and order for later on so that we can retrieve the
		// integration points without getting the fe
		// > this is useful for surface elements where we do not necessarily know the
		//   neighboring volume element (GetFe will fail in this case)
		IntegMethod_ = method;
		IntegOrder_ = order;
		
		ptFeSpace_->GetIntScheme()->GetIntPoints(Elem::GetShapeType(el->type), method, order,
            intPoints, weights);
    
		
		numIntegrationPoints_ = intPoints.GetSize();
    
		// create a map that gives to each integration point + to midpoint an index
		// we cannot have a map of loc points directly (not sortable) but each point
		// has also a number which we can use for sorting; midpoint gets index -1
		locPointIndices_.insert(std::pair<int, UInt>(-1, 0));
    
		for (UInt i = 0; i < numIntegrationPoints_; i++) {
			//std::cout << intPoints[i].number << std::endl;
			locPointIndices_.insert(std::pair<int, UInt>(intPoints[i].number, i + 1));
		}
    
		//std::cout << "NumIntegrationPoints: " << numIntegrationPoints_ << std::endl;
		//std::cout << "Integration points:" << std::endl;
		//std::cout << intPoints.ToString() << std::endl;
    
		//std::cout << "Mapping: " << std::endl;
		//std::map<int,UInt>::iterator mapit;
		//for(mapit=locPointIndices_.begin(); mapit!=locPointIndices_.end(); mapit++){
		//  std::cout << mapit->first << "; " << mapit->second << std::endl;
		//}
    
		// dim_ is the dim of the output retrieved by GetVector
		// not tha same as Preisach_Dim that determines whether we use scalar or vector model
		dim_ = dependCoef1_->GetVecSize();
    
		PtrCoefFct matCoef = material_->GetTensorCoefFnc(matType_, tensorType_,
            Global::REAL, false);
    
		/*
     * MAT_initialTensor_ is the small signal tensor from the mat file
     */
		matCoef->GetTensor(MAT_initialTensor_, lpm);
    
		// to calculate differential material properties, we need to know e0 / nu0
		if (material_->GetMaterialDatabaseName() == "Electrostatic") {
			rev_mat_fac_ = 8.854187817e-12; //eps0
			PDEName_ = "Electrostatic";
      needsInversion_ = false;
      
      // small signal tensor = permittivity
      MAT_smallSignalTensor_ = MAT_initialTensor_;
      
		} else if (material_->GetMaterialDatabaseName() == "Electromagnetics") {
			rev_mat_fac_ = 795774.7155; //nu0
			PDEName_ = "Electromagnetics";
      needsInversion_ = true;
      
      PtrCoefFct permeability = material_->GetTensorCoefFnc(MAG_PERMEABILITY,tensorType_,
              Global::REAL, false);
      
      MAT_smallSignalTensor_ = Matrix<Double>(dim_, dim_);
      permeability->GetTensor(MAT_smallSignalTensor_, lpm);
      
		} else {
			EXCEPTION("Currently only Electrostatics and Electromagnetics are supported");
		}
    
		/*
     * MAT_freeFieldTensor_ is the tensor to be returned, when
     * RUN_tensorToReturn_ = 2
     * and the tensor to be added to deltaMat if
     * RUN_tensorToAdd_ = 2
     * (only ONE needed per PDE)
     */
		MAT_freeFieldTensor_ = Matrix<Double>(dim_, dim_);
		MAT_freeFieldTensor_.Init();
		for (UInt i = 0; i < dim_; i++) {
			MAT_freeFieldTensor_[i][i] = rev_mat_fac_;
		}
    
    // for piezoelectric / magnetostrictive coupling
    double strainForm;
		material_->GetScalar(strainForm, HYST_STRAIN_FORM);
    
    if(((int) strainForm) != 0){
      MAT_useStrainForm_ = true;
    } else {
      MAT_useStrainForm_ = false;
    }
    
    
    //    std::cout << "StrainForm: " << strainForm << std::endl;
    //    std::cout << "MAT_useStrainForm_: " << MAT_useStrainForm_ << std::endl;
    //    
    //    
	}
  
	CoefFunctionHyst::~CoefFunctionHyst() {
    
		delete timer_;
    
		delete[] E_B_;
		delete[] P_J_;
		delete[] E_H_;
    
		delete[] E_B_lastIt_;
		delete[] P_J_lastIt_;
		delete[] E_H_lastIt_;
    
		delete[] E_B_lastTS_;
		delete[] P_J_lastTS_;
		delete[] E_H_lastTS_;
    
		delete[] deltaMat_;
		delete[] deltaMat_lastIt_;
		delete[] deltaMat_lastTS_;
		delete[] deltaMat_estimated_;
    
		delete[] dY_sol;
		delete[] X_low;
		delete[] X_up;
    
		delete[] MAT_initialOutput_;
	}
  
	void CoefFunctionHyst::InitStorage() {
		/*
     * this function sets up the actual storage vectors and the the hysteretic
     * objects
     * it has to be called AFTER the XML dependent parameter are known
     * as the size of the storage depends on the evaluation depth (standard,
     * extended, full) which is retrieved from the non-linear parameters (see
     * stdSolveStep::ReadNonLinData)
     */
    if(storageInitialized_){
      return;
    }
    
    
		// 1. determine the number of required storages
		if (XML_EvaluationDepth_ == 1) {
			// standard evaluation
			// > one hyst operator per element
			// > only one storage for each element
			numHystOperators_ = numElemSD_;
			numStorageEntries_ = numElemSD_;
		} else if (XML_EvaluationDepth_ == 2) {
			// extended evaluation
			// > one hyst operator per element
			// > each integration point (+ midpoint) needs one storage
			numHystOperators_ = numElemSD_;
			numStorageEntries_ = numElemSD_ * (numIntegrationPoints_ + 1);
		} else if (XML_EvaluationDepth_ == 3) {
			// full evaluation
			// > each integration point (+ midpoint) get a hyst operator
			// > each integration point (+ midpoint) needs one storage
			numHystOperators_ = numElemSD_ * (numIntegrationPoints_ + 1);
			numStorageEntries_ = numElemSD_ * (numIntegrationPoints_ + 1);
		} else {
			EXCEPTION("Evaluation depth < 1 or > 3 not allowed")
		}
    
		//std::cout << "XML_EvaluationDepth_: " << XML_EvaluationDepth_ << std::endl;
    //		std::cout << "numIntegrationPoints_: " << numIntegrationPoints_ << std::endl;
    //		std::cout << "numStorageEntries_: " << numStorageEntries_ << std::endl;
    //		std::cout << "numHystOperators_: " << numHystOperators_ << std::endl;
    
    /*
     * Input for piezoelectric / magnetostrictive setups
     */   

		double tmp;
		material_->GetScalar(tmp, DIM_BETA_COEFS);
		MAT_dim_beta_ = (int) tmp;
    material_->GetTensor(MAT_betaCoefs_, HYST_BETA_COEFS, Global::REAL);
    
		/*
     * setup hysteresis operator
     * > in total numHystOperators_
     */
		std::string dimTypeStr;
		material_->GetScalar(dimTypeStr, PREISACH_DIM);
    
		// use variable MAT_methodType_ to distinguish between scalar and vector model
		// do not confuse this with dimType_ !
		if (dimTypeStr == "SCALAR") {
			MAT_methodType_ = SCALAR;
		} else if (dimTypeStr == "VECTOR") {
			MAT_methodType_ = VECTOR;
		}
    
		if (MAT_methodType_ == SCALAR) {
			//get direction
			std::string str;
			material_->GetScalar(str, P_DIRECTION);
			Directions dir;
			String2Enum(str, dir);
			MAT_dirP_ = dir;
      bool isVirgin = true;
      
      hyst_ = new Preisach(numHystOperators_, MAT_xSat_, MAT_ySat_, MAT_PreisachWeights_, isVirgin);
      hasInverseModel_ = false;
      
			/*
       * currently we do not support initial input for Scalar model
       */
			MAT_initialInput_ = Vector<Double>(dim_);
			MAT_initialInput_.Init();
			MAT_initialOutput_ = new Vector<Double>[numStorageEntries_];
			for (UInt k = 0; k < numStorageEntries_; k++) {
				MAT_initialOutput_[k] = Vector<Double>(dim_);
				MAT_initialOutput_[k].Init();
			}
      
		} else if (MAT_methodType_ == VECTOR) {
      hasInverseModel_ = false;
			material_->GetScalar(MAT_rotRes_, ROT_RESISTANCE, Global::REAL);
			material_->GetScalar(MAT_angResistance_, ANG_DISTANCE, Global::REAL);
			material_->GetScalar(MAT_angResolution_, ANG_RESOLUTION, Global::REAL);
			material_->GetScalar(MAT_angClipping_, ANG_CLIPPING, Global::REAL);
			material_->GetScalar(MAT_ampResolution_, AMP_RESOLUTION, Global::REAL);
      
			double printOut;
			double bmpResolution;
			material_->GetScalar(printOut, PRINT_PREISACH);
			material_->GetScalar(bmpResolution, PRINT_PREISACH_RESOLUTION);
      
			/*
       * if printOut > 0 -> activate output of overlaid switching and rotation state; output every printOut timestep
       * bmpResolution -> number of pixels (std = 1000)
       */
			MAT_printOut_ = (UInt) printOut;
			MAT_bmpResolution_ = (UInt) bmpResolution;
      
			double eval;
			material_->GetScalar(eval, EVAL_VERSION);
      
			MAT_vecPreisachImplementationVersion_ = (UInt) eval;
      
			// this flag is not used currently;
			// if we have an initial state, isVirgin should be false, but as already set,
			// it is not used at the momement
			bool isVirgin = true;
      
			if (MAT_vecPreisachImplementationVersion_ == 1) {
				isClassical_ = true; // original vector preisach model -> sutor2012
        
				hyst_ = new VectorPreisachv10_ListApproach(numHystOperators_, MAT_xSat_, MAT_ySat_,
                MAT_PreisachWeights_, MAT_rotRes_, dim_, isVirgin,
                isClassical_, MAT_angResistance_, MAT_angClipping_);
			} else if (MAT_vecPreisachImplementationVersion_ == 2) {
				isClassical_ = false; // revised vector preisach model -> sutor2015
        
				hyst_ = new VectorPreisachv10_ListApproach(numHystOperators_, MAT_xSat_, MAT_ySat_,
                MAT_PreisachWeights_, MAT_rotRes_, dim_, isVirgin,
                isClassical_, MAT_angResistance_, MAT_angClipping_);
			} else if (MAT_vecPreisachImplementationVersion_ == 10) {
				isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
        
				hyst_ = new VectorPreisachv10_MatrixApproach(numHystOperators_, MAT_xSat_, MAT_ySat_,
                MAT_PreisachWeights_, MAT_rotRes_, dim_, isVirgin,
                isClassical_, MAT_angResistance_, MAT_angClipping_);
			} else if (MAT_vecPreisachImplementationVersion_ == 20) {
				isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
        
				hyst_ = new VectorPreisachv10_MatrixApproach(numHystOperators_, MAT_xSat_, MAT_ySat_,
                MAT_PreisachWeights_, MAT_rotRes_, dim_, isVirgin,
                isClassical_, MAT_angResistance_, MAT_angClipping_);
			} else {
				EXCEPTION("MAT_vecPreisachImplementationVersion_ has to be one of the following: \n "
                "1: classical vector model (sutor2012) \n"
                "2: revised vector model (sutor2015) [DEFAULT] \n"
                "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
                "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
			}
      
			// initial input that shall be feeded to the vectorPreisach operator
			MAT_initialInput_ = Vector<Double>(dim_);
      
			material_->GetScalar(MAT_initialInput_[0], INITIAL_STATE_X, Global::REAL);
			material_->GetScalar(MAT_initialInput_[1], INITIAL_STATE_Y, Global::REAL);
			if (dim_ == 3) {
				material_->GetScalar(MAT_initialInput_[2], INITIAL_STATE_Z, Global::REAL);
			}
      
			if (MAT_initialInput_.NormL2() > 1e-16) {
				WARN("Currently the treatment of initial states is not working properly! \n"
                "Depending on the selected evaluation approach, the initial state of \n"
                "the hysteresis operator will act as an excitation from the first iteration on or not. \n"
                "Furthermore, the initial state has to fit to the boundary condition \n"
                "(i.e. fluxParallel with initial state standing perpendicular on the boundary will not work).\n"
                "The only save usage is in context of the debugging fixPoint iteration \n"
                "(i.e. output of hysteresis operator does not couple back).")
			}
			MAT_initialOutput_ = new Vector<Double>[numStorageEntries_];
      
			if (MAT_initialInput_.NormL2() > 1e-16) {
				for (UInt k = 0; k < numHystOperators_; k++) {
					/*
           * although all hyst operators should produce the same output, we have to
           * call the evaluation for each one of them to get its internal memory set
           */
					Vector<Double> tmp;
					tmp = hyst_->computeValue_vec(MAT_initialInput_, k, true, true);
          
					if (numHystOperators_ != numStorageEntries_) {
						// for standard and extended evaluation, we use one hyst operator for multiple
						// intregration points but for extended evalutation, each integration point
						// has its own storage
						for (UInt i = 0; i < numIntegrationPoints_ + 1; i++) {
							MAT_initialOutput_[k * (numIntegrationPoints_ + 1) + i] = tmp;
						}
					} else {
						// for standard evaluation, we have one hyst operator and one storage per element
						// (numHystOperators = numStorageEntries) and for full evaluation
						// we have one hyst operator and one storage for each integration point (+midpoint)
						// (numHystOperators = numStorageEntries)
						MAT_initialOutput_[k] = tmp;
					}
				}
			} else {
				for (UInt k = 0; k < numStorageEntries_; k++) {
					MAT_initialOutput_[k] = Vector<Double>(dim_);
					MAT_initialOutput_[k].Init();
				}
			}
		}
    
		if (MAT_ampResolution_ == 0) {
			// for scalar case and in case thatt no resolution was set in mat file
			MAT_ampResolution_ = 1e-17;
		}
    
    
		/*
     * finally initialize storage
     * NEW: use same storage for vector and scalar model
     * (even though the scalar model could use smaller storage, too)
     * NEW: create storage for each hystOperator or for each integration point
     * depending on evaluation depth
     */
		Vector<Double> zeroVec = Vector<Double>(dim_);
		zeroVec.Init();
    
		E_B_ = new Vector<Double>[numStorageEntries_];
		P_J_ = new Vector<Double>[numStorageEntries_];
		E_H_ = new Vector<Double>[numStorageEntries_];
    
		E_B_lastIt_ = new Vector<Double>[numStorageEntries_];
		P_J_lastIt_ = new Vector<Double>[numStorageEntries_];
		E_H_lastIt_ = new Vector<Double>[numStorageEntries_];
    
		E_B_lastTS_ = new Vector<Double>[numStorageEntries_];
		P_J_lastTS_ = new Vector<Double>[numStorageEntries_];
		E_H_lastTS_ = new Vector<Double>[numStorageEntries_];
    
		deltaMat_ = new Matrix<Double>[numStorageEntries_];
		deltaMat_lastIt_ = new Matrix<Double>[numStorageEntries_];
		deltaMat_lastTS_ = new Matrix<Double>[numStorageEntries_];
		deltaMat_estimated_ = new Matrix<Double>[numStorageEntries_];
    
		requiresReeval_ = new bool[numStorageEntries_];
    Si_requiresReeval_ = new bool[numStorageEntries_];
    deltaMat_requiresReeval_ = new bool[numStorageEntries_];
    deltaMatStrain_requiresReeval_ = new bool[numStorageEntries_];
    
		rotatedCouplingTensor_ = new Matrix<Double>[numStorageEntries_];
		rotatedCouplingTensor_requiresReeval_ = new bool[numStorageEntries_];
		//TODO: coupling tensor einlesen und rotatedCouplingTensor entsprechend initialisieren
    
		Vector<Double> slope = Vector<Double>(dim_);
		slope.Init();
    
		dY_sol = new Vector<Double>[numStorageEntries_];
		X_low = new Vector<Double>[numStorageEntries_];
		X_up = new Vector<Double>[numStorageEntries_];
    
		for (UInt k = 0; k < numStorageEntries_; k++) {
			dY_sol[k] = zeroVec;
			X_low[k] = zeroVec;
			X_up[k] = zeroVec;
			requiresReeval_[k] = true;
      
			/*
       * General note:
       *  We assume that our material was in virgin state, i.e.
       *    all _lastIt_ and _lastTS_ vectors are zeroVec
       *    and the _lastIt_ and _lastTS_ matrices are MAT_InitialTensor
       *  The possible initial input provided from the matrial file shall
       *    be treated as if it was an evaluation before the first actual iteration
       *    i.e. the arrays for the current values (E_H_, E_B_ and P_J_) store
       *    the initial input and output
       */
      
			// E_B_ stores (as the name says) E for electrostatics and B for magnetics
			// E acts as both, the element solution and the input to the hysteresis operator
			//  in case of electrostatics
			// B acts as solution in magnetics, but not as input
			//  i.e. we cannot set E_B_ to MAT_initialInput_ in that case
			//  as MAT_initialInput_ = H
			if (PDEName_ == "Electromagnetics") {
				// here we have to estimate the initial solution first from H (MAT_initialInput_)
				// and M (MAT_initialOutput_
				E_B_[k] = MAT_initialInput_ + MAT_initialOutput_[k]; //H+M
				E_B_[k] = E_B_[k] / rev_mat_fac_; //B = (H+M)/nu
			} else {
				E_B_[k] = MAT_initialInput_;
			}
			E_B_lastIt_[k] = zeroVec;
			E_B_lastTS_[k] = zeroVec;
      
			// E_H_ stores E for electrostatics and H for magnetics
			// -> both are the input to the hyst operator, so simply initialize these
			//    arrays with MAT_initialInput_
			E_H_[k] = MAT_initialInput_;
			E_H_lastIt_[k] = zeroVec;
			E_H_lastTS_[k] = zeroVec;
      
			// P_J_ stores P for electrostatics and M for magnetics
			// -> both are the output of the hyst operator, so simply initialize those
			//    arrays with MAT_initialOutput_[k]
			P_J_[k] = MAT_initialOutput_[k];
			//std::cout << "k: " << k << std::endl;
			//std::cout << "P_J_[k]: " << P_J_[k].ToString() << std::endl;
			P_J_lastIt_[k] = zeroVec;
			P_J_lastTS_[k] = zeroVec;
      
			deltaMat_lastIt_[k] = MAT_initialTensor_;
			deltaMat_lastTS_[k] = MAT_initialTensor_;
      
			deltaMat_estimated_[k] = Matrix<Double>(dim_, dim_);
      
			if (MAT_initialInput_.NormL2() == 0) {
				// no initial input > use MAT_initialTensor_
				deltaMat_[k] = MAT_initialTensor_;
			} else {
        // TODO: compute first deltaMat or start without it
				// create a first deltaMatrix    
				//CreateDeltaMatrix(MAT_initialInput_, MAT_initialOutput_[k], deltaMat_[k], evalMethodName, k, intoSat, outofSat, satToSat, MAT_initialInput_);
			}
		}
 
    storageInitialized_ = true;
	}
  
	
	void CoefFunctionHyst::AddAdditionalSDList(shared_ptr<EntityList> actSDList, bool isSurface){
		if(storageInitialized_ == true){
			EXCEPTION("Storage was already initialized! Cannot add further SDLists");
		}
		//std::cout << "Add additional SD list to HystOperator" << std::endl;
		SDList_List_.push_back(actSDList);
		
		EntityIterator it = actSDList->GetIterator();
		
		UInt globalElNr;
		for (it.Begin(); !it.IsEnd(); it++) {
			globalElNr = it.GetElem()->elemNum;
			//std::cout << "Add element with #"<<globalElNr << " to list" << std::endl;
			globalElem2Local_[globalElNr] = numElemSD_;
			globalElemOnSurf_[globalElNr] = isSurface;
			numElemSD_++;
		}
		
		// add integration points to locPointIndices_
		
		// pick out the first element (even though any of them would do as they share the
		// same material) and extract midpoint
		it.Begin();
		const Elem * el = it.GetElem();
		
		
		//std::cout << "Elem Type: " << el->type << std::endl;
		LocPoint lp = Elem::shapes[el->type].midPointCoord;
		
		
		
		//std::cout << "Midpoint coords: " << Elem::shapes[el->type].midPointCoord.ToString() << std::endl;
		LocPointMapped lpm;
		shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
		//std::cout << "0" << std::endl;
		lpm.Set(lp, esm, 0.0);
		//std::cout << "1" << std::endl;
		// get number of integration points
		
		StdVector<LocPoint> intPoints;
		StdVector<Double> weights;
		// here we use the information about the integration method and order that we retrieved in the
		// constructor
		ptFeSpace_->GetIntScheme()->GetIntPoints(Elem::GetShapeType(el->type), IntegMethod_, IntegOrder_,
            intPoints, weights);
    
		UInt numIntegrationPointsLocal = intPoints.GetSize();
		if(numIntegrationPointsLocal > numIntegrationPoints_){
			// should not be the case as the number of integration points on surf elements usually is smaller than the
			// number of integration points on the volume
			//std::cout << "On the newly added element, more integration points are defined as on the already stored ones." << std::endl;
			numIntegrationPoints_ = numIntegrationPointsLocal;
		}
		
		// create a map that gives to each integration point + to midpoint an index
		// we cannot have a map of loc points directly (not sortable) but each point
		// has also a number which we can use for sorting; midpoint gets index -1
		// midpoint already inserted
		
		//locPointIndices_.insert(std::pair<int, UInt>(-1, 0));
		
		for (UInt i = 0; i < numIntegrationPointsLocal; i++) {
			//std::cout << intPoints[i].number << std::endl;
			//std::cout << "Add integration point pair (" << intPoints[i].number << ", " << i+1 << ")" << std::endl;
			locPointIndices_.insert(std::pair<int, UInt>(intPoints[i].number, i + 1));
		}
	}
	
	bool CoefFunctionHyst::PreprocessLPM(const LocPointMapped& lpmInput,
          LocPointMapped& lpmOutput, UInt& operatorIdx, UInt& storageIdx, bool forceMidpoint) {
		/*
     *	Input:
     *		LocPointMapped lpmInput coming from numerical integration
     *
     *	Output:
     *		LocPointMapped lpmOutput specifying the position at which the hysteresis operator
     *		shall be evaluated
     *
     *		UInt operatorIdx reffering to the hysteresis operator that has to be evaluated
     *
     *		UInt storageIdx specifying the position in the storage arrays
     */
		const Elem * el = lpmInput.ptEl;
    //std::cout << "Preprocess LPM" << std::endl;
    //std::cout << "Input LPM on surface? " << lpmInput.isSurface << std::endl;
    //std::cout << "Corresponding element: " << el->ToString() << std::endl;
    //std::cout << "Coordinates of LP: " << lpmInput.lp.coord.ToString() << std::endl;
    
		/*
     * standard evaluation: midpointOnly ALWAYS used
     * extended evaluation: midpointOnly = true for output computation and
     *                                      during SetPreviousHystValues
     * full evaluation: midpointOnly = true only for output computation
     *
     * > use function EvaluateAtMidpointOnly() to determine the flag
     * 
     * for setprevioushystvalues, we might want to evaluate BOTH at integration points
     * and the midpoint; however the function EvaluateAtMidpointOnly would only give us
     * the actualy integrations points by then; by the additional input parameter
     * forceMidpoint, we can still retrieve the midpoint, though
     */
    if(forceMidpoint == false){
      forceMidpoint = EvaluateAtMidpointOnly();
    }
    
		if (forceMidpoint) {
      //std::cout << "Evaluate only at midpoint" << std::endl;
			/*
       * get solution at midpoint of element (elemSolution)
       */
			LocPoint lp = Elem::shapes[el->type].midPointCoord;
			shared_ptr<ElemShapeMap> esm = lpmInput.shapeMap;
      
			lpmOutput.Set(lp, esm, 0.0);
      
		} else {
      //std::cout << "Evaluate at integration point" << std::endl;
			/*
       * get solution at actual integration point
       */
			lpmOutput = lpmInput;
		}
    
		UInt idxElem = globalElem2Local_[el->elemNum];
		bool onBoundary = globalElemOnSurf_[el->elemNum];
		UInt idxPoint = locPointIndices_[lpmOutput.lp.number];
    //std::cout << "Global Indices: " << std::endl;
    //std::cout << "ElementI Number: " << el->elemNum << std::endl;
    //std::cout << "LPM Number: " << lpmOutput.lp.number << std::endl;		
		
		//std::cout << "Mapped indices: " << std::endl;
    //std::cout << "ElementIndex: " << idxElem << std::endl;
    //std::cout << "PointIndex: " << idxPoint << std::endl;
    
		// storage index = index where to store results, input, ...
		// operator index = index of hyst operator to evaluate
		if (XML_EvaluationDepth_ == 1) {
			// standard evaluation
			// storageIndex = operatorIndex = elementIndex
			storageIdx = idxElem;
			operatorIdx = idxElem;
		} else if (XML_EvaluationDepth_ == 2) {
			// extended evaluation
			// operatorIndex = elementIndex (only one operator per element)
			// storageIndex = combinedIndex = elementIndex*(numIntegrationPoints+1)+pointIndex
			operatorIdx = idxElem;
			storageIdx = idxElem * (numIntegrationPoints_ + 1) + idxPoint;
		} else {
			// full evaluation
			// one operator and one storage per integration point
			// operatorIndex = storageIndex = combinedIndex
			operatorIdx = idxElem * (numIntegrationPoints_ + 1) + idxPoint;
			storageIdx = operatorIdx;
		}
		
		return onBoundary;
	}
  
  Matrix<Double> CoefFunctionHyst::GetDeltaMat(const LocPointMapped& Originallpm, int timelevel_new, int timelevel_old, bool useStrains, bool useAbs,
          std::string implementationVersion){
    
    
    LOG_DBG(coeffcthyst) << "GetDeltaMat";
    LOG_DBG(coeffcthyst) << "Compute delta between timelevel " << timelevel_new << " and time level " << timelevel_old << "(-2 = deactivated; -1 = lastTS; 0 = current IT; 1 = lastIT)";
    
    if(tensorType_ == AXI){
      EXCEPTION("GetDeltaMat only implemented for 2d plane and full 3d setups");
    }
    
    Matrix<Double> deltaMat;
    UInt numCols = dim_;
    UInt numRows;
    
    Double absTol = 1e-16;
    Double relTol = 1e-16;
    
    UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;
    
		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
	  // deltaMat can only be computed on volume elements, not on boundaries
		//bool onBoundary = false;
		if(onBoundary == true){
			EXCEPTION("DeltaMat not defined on boundary");
		}
		
    if(useStrains){
      
      if(deltaMatStrain_requiresReeval_[storageIdx] == false){
        return deltaMatStrain_[storageIdx];
      }
      
      /*
       * Compute matrix dM, such that
       * dS = dM*dE (or dB)
       * 
       * with dS = Strains(timeLevel1)-Strains(timeLevel2)
       *      dE = E(timeLevel1)-E(timeLevel2)
       */
      
      if(dim_ == 3){
        // 6x3 tensor
        numRows = 6;
      } else {
        // 3x2 tensor
        numRows = 3;
      }
      deltaMat = Matrix<Double>(numRows,numCols);
      
      Vector<Double> E_B_new = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_new, onBoundary);
      Vector<Double> E_B_old = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_old, onBoundary);
      
      Vector<Double> E_B_diff = E_B_new;
      E_B_diff -= E_B_old;
      Double E_B_diff_norm = E_B_diff.NormL2();
      
      if(E_B_diff_norm <= absTol){
        // variation of solution is very small
        // return zero matrix
        return deltaMat;
      }
      Vector<Double> S_new = GetIrreversibleStrains(Originallpm,timelevel_new);
      Vector<Double> S_old = GetIrreversibleStrains(Originallpm,timelevel_old);
      Vector<Double> S_diff = S_new;
      S_diff -= S_old;
      
      if(implementationVersion == "Division"){
        if(dim_ == 2){
          
          if(abs(E_B_diff[0])/E_B_diff_norm > relTol){
            // division can be used
            deltaMat[0][0] = S_diff[0]/E_B_diff[0];
          } else if(abs(E_B_diff[1])/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = S_diff[0]/E_B_diff[1];
            deltaMat[1][0] = S_diff[0]/E_B_diff[1];
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }
          
          if(abs(E_B_diff[1])/E_B_diff_norm > relTol){
            deltaMat[1][1] = S_diff[1]/E_B_diff[1];
          } else if(abs(E_B_diff[0])/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = S_diff[1]/E_B_diff[0];
            deltaMat[1][0] = S_diff[1]/E_B_diff[0];
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }
          Double sum,diff;
          sum = E_B_diff[1]+E_B_diff[2];
          diff = E_B_diff[1]-E_B_diff[2];
          
          if(abs(sum)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[2][0] = S_diff[1]/sum;
            deltaMat[2][1] = S_diff[1]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[2][0] = S_diff[1]/diff;
            deltaMat[2][1] = -S_diff[1]/diff;
          } 
          
        } else {
          Double sum,diff;
          sum = E_B_diff[1]+E_B_diff[2];
          diff = E_B_diff[1]-E_B_diff[2];
          
          if(abs(E_B_diff[0])/E_B_diff_norm > relTol){
            // division can be used
            deltaMat[0][0] = S_diff[0]/E_B_diff[0];
          } else if(abs(sum)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = S_diff[0]/sum;
            deltaMat[1][0] = S_diff[0]/sum;
            deltaMat[0][2] = S_diff[0]/sum;
            deltaMat[2][0] = S_diff[0]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = S_diff[0]/diff;
            deltaMat[1][0] = S_diff[0]/diff;
            deltaMat[0][2] = -S_diff[0]/diff;
            deltaMat[2][0] = -S_diff[0]/diff;
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }
          
          sum = E_B_diff[0]+E_B_diff[2];
          diff = E_B_diff[0]-E_B_diff[2];
          
          if(abs(E_B_diff[1])/E_B_diff_norm > relTol){
            // division can be used
            deltaMat[1][1] = S_diff[1]/E_B_diff[1];
          } else if(abs(sum)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = S_diff[1]/sum;
            deltaMat[1][0] = S_diff[1]/sum;
            deltaMat[1][2] = S_diff[1]/sum;
            deltaMat[2][1] = S_diff[1]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = S_diff[1]/diff;
            deltaMat[1][0] = S_diff[1]/diff;
            deltaMat[1][2] = -S_diff[1]/diff;
            deltaMat[2][1] = -S_diff[1]/diff;
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }
          
          sum = E_B_diff[0]+E_B_diff[1];
          diff = E_B_diff[0]-E_B_diff[1];
          
          if(abs(E_B_diff[2])/E_B_diff_norm > relTol){
            // division can be used
            deltaMat[2][2] = S_diff[2]/E_B_diff[2];
          } else if(abs(sum)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][2] = S_diff[2]/sum;
            deltaMat[2][0] = S_diff[2]/sum;
            deltaMat[1][2] = S_diff[2]/sum;
            deltaMat[2][1] = S_diff[2]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][2] = S_diff[2]/diff;
            deltaMat[2][0] = S_diff[2]/diff;
            deltaMat[1][2] = -S_diff[2]/diff;
            deltaMat[2][1] = -S_diff[2]/diff;
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }     
          
          sum = E_B_diff[1]+E_B_diff[2];
          diff = E_B_diff[1]-E_B_diff[2];
          
          if(abs(sum)/E_B_diff_norm > relTol) {
            // ds_zy/(de_y + de_z)
            deltaMat[3][1] = S_diff[3]/sum;
            deltaMat[3][2] = S_diff[3]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // ds_zy/(de_y + de_z)
            deltaMat[3][1] = S_diff[3]/diff;
            deltaMat[3][2] = -S_diff[3]/diff;
          } else if(abs(E_B_diff[0])/E_B_diff_norm > relTol) {
            deltaMat[3][0] = S_diff[3]/E_B_diff[0];
          }
          
          sum = E_B_diff[0]+E_B_diff[2];
          diff = E_B_diff[0]-E_B_diff[2];
          
          if(abs(sum)/E_B_diff_norm > relTol) {
            // ds_zy/(de_y + de_z)
            deltaMat[4][0] = S_diff[4]/sum;
            deltaMat[4][2] = S_diff[4]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // ds_zy/(de_y + de_z)
            deltaMat[4][0] = S_diff[4]/diff;
            deltaMat[4][2] = -S_diff[4]/diff;
          } else if(abs(E_B_diff[1])/E_B_diff_norm > relTol) {
            deltaMat[4][1] = S_diff[4]/E_B_diff[1];
          }
          
          sum = E_B_diff[0]+E_B_diff[1];
          diff = E_B_diff[0]-E_B_diff[1];
          
          if(abs(sum)/E_B_diff_norm > relTol) {
            // ds_zy/(de_y + de_z)
            deltaMat[5][0] = S_diff[5]/sum;
            deltaMat[5][1] = S_diff[5]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // ds_zy/(de_y + de_z)
            deltaMat[5][0] = S_diff[5]/diff;
            deltaMat[5][1] = -S_diff[5]/diff;
          } else if(abs(E_B_diff[2])/E_B_diff_norm > relTol) {
            deltaMat[5][2] = S_diff[5]/E_B_diff[2];
          } 
          
        }
      } else {
        EXCEPTION("Delta mat version not implemented yet");
      }
      
      deltaMatStrain_[storageIdx] = deltaMat;
      deltaMatStrain_requiresReeval_[storageIdx] = false;
      
    } else {
      
      if(deltaMat_requiresReeval_[storageIdx] == false){
        LOG_DBG(coeffcthyst) << "Reuse old deltaMat";
        //std::cout << "Reuse old deltaMat" << std::endl;
        return deltaMat_[storageIdx];
      }
      
      /*
       * Compute matrix dM, such that
       * dP = dM*dE (or dB)
       * 
       * with dP = Polarization(timeLevel1)-Polarizatoin(timeLevel2)
       *      dE = E(timeLevel1)-E(timeLevel2)
       */
      numRows = dim_;
      deltaMat = Matrix<Double>(numRows,numCols);
      
      Vector<Double> E_B_new = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_new, onBoundary);
      Vector<Double> E_B_old = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_old, onBoundary);
      
      LOG_DBG(coeffcthyst) << "Old solution of system (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString();
      LOG_DBG(coeffcthyst) << "Current solution of system (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString();
      //std::cout << "Old solution (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString() << std::endl;
      //std::cout << "Current solution (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString() << std::endl;
      
      Vector<Double> E_B_diff = E_B_new;
      E_B_diff -= E_B_old;
      Double E_B_diff_norm = E_B_diff.NormL2();
      
      //std::cout << "Difference Vector: " << E_B_diff.ToString() << std::endl;
      
      if(E_B_diff_norm <= absTol){
        // variation of solution is very small
        // return zero matrix
        return deltaMat;
      }
      
      Vector<Double> P_J_new = GetOutputOfHysteresisOperator(Originallpm,timelevel_new);
      Vector<Double> P_J_old = GetOutputOfHysteresisOperator(Originallpm,timelevel_old);
      
      LOG_DBG(coeffcthyst) << "Old output of hyst operator: " << P_J_old.ToString();
      LOG_DBG(coeffcthyst) << "Current output of hyst operator: " << P_J_new.ToString();
      
      //std::cout << "Old output of hystoperator (solution at timelevel = " << timelevel_old << "): " << P_J_old.ToString() << std::endl;
      //std::cout << "Current output of hystoperator (solution at timelevel = " << timelevel_new << "): " << P_J_new.ToString() << std::endl;
      
      Vector<Double> P_J_diff = P_J_new;
      P_J_diff -= P_J_old;
      //std::cout << "Difference Vector: " << P_J_diff.ToString() << std::endl;
      if(implementationVersion == "Division"){
        if(dim_ == 2){
          
          if(abs(E_B_diff[0])/E_B_diff_norm > relTol){
            // division can be used
            deltaMat[0][0] = P_J_diff[0]/E_B_diff[0];
          } else if(abs(E_B_diff[1])/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = P_J_diff[0]/E_B_diff[1];
            deltaMat[1][0] = P_J_diff[0]/E_B_diff[1];
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }
          
          if(abs(E_B_diff[1])/E_B_diff_norm > relTol){
            deltaMat[1][1] = P_J_diff[1]/E_B_diff[1];
          } else if(abs(E_B_diff[0])/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = P_J_diff[1]/E_B_diff[0];
            deltaMat[1][0] = P_J_diff[1]/E_B_diff[0];
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }
          
        } else {
          Double sum,diff;
          sum = E_B_diff[1]+E_B_diff[2];
          diff = E_B_diff[1]-E_B_diff[2];
          
          if(abs(E_B_diff[0])/E_B_diff_norm > relTol){
            // division can be used
            deltaMat[0][0] = P_J_diff[0]/E_B_diff[0];
          } else if(abs(sum)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = P_J_diff[0]/sum;
            deltaMat[1][0] = P_J_diff[0]/sum;
            deltaMat[0][2] = P_J_diff[0]/sum;
            deltaMat[2][0] = P_J_diff[0]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = P_J_diff[0]/diff;
            deltaMat[1][0] = P_J_diff[0]/diff;
            deltaMat[0][2] = -P_J_diff[0]/diff;
            deltaMat[2][0] = -P_J_diff[0]/diff;
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }
          
          sum = E_B_diff[0]+E_B_diff[2];
          diff = E_B_diff[0]-E_B_diff[2];
          
          if(abs(E_B_diff[1])/E_B_diff_norm > relTol){
            // division can be used
            deltaMat[1][1] = P_J_diff[1]/E_B_diff[1];
          } else if(abs(sum)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = P_J_diff[1]/sum;
            deltaMat[1][0] = P_J_diff[1]/sum;
            deltaMat[1][2] = P_J_diff[1]/sum;
            deltaMat[2][1] = P_J_diff[1]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][1] = P_J_diff[1]/diff;
            deltaMat[1][0] = P_J_diff[1]/diff;
            deltaMat[1][2] = -P_J_diff[1]/diff;
            deltaMat[2][1] = -P_J_diff[1]/diff;
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }
          
          sum = E_B_diff[0]+E_B_diff[1];
          diff = E_B_diff[0]-E_B_diff[1];
          
          if(abs(E_B_diff[2])/E_B_diff_norm > relTol){
            // division can be used
            deltaMat[2][2] = P_J_diff[2]/E_B_diff[2];
          } else if(abs(sum)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][2] = P_J_diff[2]/sum;
            deltaMat[2][0] = P_J_diff[2]/sum;
            deltaMat[1][2] = P_J_diff[2]/sum;
            deltaMat[2][1] = P_J_diff[2]/sum;
          } else if(abs(diff)/E_B_diff_norm > relTol) {
            // distribute entry to non-diagonal entries
            deltaMat[0][2] = P_J_diff[2]/diff;
            deltaMat[2][0] = P_J_diff[2]/diff;
            deltaMat[1][2] = -P_J_diff[2]/diff;
            deltaMat[2][1] = -P_J_diff[2]/diff;
          } 
          //          else {
          //            // all entries are too small (relatively > return zero matrx
          //            // leave entry zero
          //          }                   
        }
      } else {
        EXCEPTION("Delta mat version not implemented yet");
      }
      
      deltaMat_[storageIdx] = deltaMat;
      deltaMat_requiresReeval_[storageIdx] = false;
      
    }
    LOG_DBG(coeffcthyst) << "Compute deltaMat for storageIDX " << storageIdx << ": " << deltaMat.ToString();
    return deltaMat;
    
  }
  
  // compute irreversible strains Si
  // Note: Return Vector in VoigtNotation!
  Vector<Double> CoefFunctionHyst::GetIrreversibleStrains(const LocPointMapped& Originallpm, int timeLevel) {
    
    if(tensorType_ == AXI){
      EXCEPTION("ComputeIrreversibleStrains only implemented for 2d plane and full 3d setups");
    }
    /*
     * Model irreversible strains by polynomial of polarization/output of hyst operator
     * ( see Numerical Simulation of Mechatronic Sensors and Actuators 3rd Edition p.386
     * 
     * Si = 3/2*(beta0 + beta1*|P| + beta2*|P|^2 + ... + betan*|P|^n)*(dirP*dirP^T - 1/3[I])
     * 
     */
    
    // Initial step: 
    // check timeLevel
    // only if timeLevel == 0 (current value) and reevaluation is required
    // we evaluate Si; otherwise reuse value
    UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;
    
		PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
    
		if (timeLevel == -1) {
			// return value from last time step
			return Si_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
			// return value from last iteration
			return Si_lastIt_[storageIdx];
		} else {
			// get current state; here we have to check if a reevaluation is needed
			// note: requiresReeval_ has one entry for each storage, not for each operator 
			if (false == Si_requiresReeval_[storageIdx]) {
        //				std::cout << "no reeval needed as flag is false!" << std::endl;
				// return current value
				return Si_[storageIdx];
			}
		}
    
    // computation required
    // setup storage
    Matrix<Double> Si_tensor = Matrix<Double>(dim_,dim_);
    Si_tensor.Init();
    
    UInt sizeVoigt;
    if(dim_ == 3){
      sizeVoigt = 6;
    } else {
      // rememeber: axi not supported
      sizeVoigt = 3;
    }
    
    Vector<Double> Si_voigt = Vector<Double>(sizeVoigt);
    Si_voigt.Init(0.0);
    
    if(MAT_dim_beta_ <= 0){
      WARN("No beta coefficients were defined. Cannot approximate Si; return empty vector");
      return Si_voigt;
    }
    
    Matrix<Double> negeye = Matrix<Double>(dim_,dim_);
    for(UInt i = 0; i < dim_; i++){
      negeye[i][i] = -1.0;
    }
    
    // get polarizaton / output of hyst operator first
    Vector<Double> P = GetOutputOfHysteresisOperator(Originallpm, timeLevel);
    
    // Si = 3/2*(beta0 + beta1*|P| + beta2*|P|^2 + ... + betan*|P|^n)*(dirP*dirP^T - 1/3[I])
    Double normP = P.NormL2();
    
    if(P.NormL2() == 0){
      // only constant offset 
      // as we do not have dirP neither, just return
      // -3/2*beta0*1/3[I] = -1/2*beta0*[I]
      Si_tensor.Add(0.5*MAT_betaCoefs_[0][0],negeye);
    } else {
      // get dirP
      Vector<Double>dirP = P;
      dirP = dirP/normP;
      Matrix<Double>dyadic = Matrix<Double>(dim_,dim_);
      dyadic.DyadicMult(dirP);
      dyadic.Add(1.0/3.0,negeye);
      
      // use Horner scheme
      // start with beta n
      Double poly = MAT_betaCoefs_[0][MAT_dim_beta_-1];
      for(UInt i = MAT_dim_beta_-2; i >= 0; i--){
        poly = poly*normP + MAT_betaCoefs_[0][i];
      }
      Si_tensor.Add(1.5*poly,dyadic);
    }
    
    // transform matrix to voigt notation
    if(dim_ == 2){
      Si_voigt[0] = Si_tensor[0][0];
      Si_voigt[1] = Si_tensor[1][1];
      Si_voigt[2] = Si_tensor[0][1];
    } else {
      Si_voigt[0] = Si_tensor[0][0];
      Si_voigt[1] = Si_tensor[1][1];
      Si_voigt[2] = Si_tensor[2][2];
      Si_voigt[3] = Si_tensor[1][2];
      Si_voigt[4] = Si_tensor[0][2];
      Si_voigt[5] = Si_tensor[0][1];
    }
    
    // flag gets reset at end of each iteration (so the stored value can beu
    // used during one iteration by several terms/function calls)
    Si_requiresReeval_[storageIdx] = false;
    Si_[storageIdx] = Si_voigt;
    return Si_voigt;
    
  }
  
  
  
	Vector<Double> CoefFunctionHyst::GetOutputOfHysteresisOperator(const LocPointMapped& Originallpm, int timeLevel) {
    //std::cout << "Get Output of Hyst Operator" << std::endl;
    //std::cout << "Timelevel: " << timeLevel << " (0 = current, -1 = lastTS, +1 = lastIt)" << std::endl;
		/*
     * Input:
     *  LocPointMapped specifying the integration point where the state of the hyst operator
     *  shall be obtained
     *
     *  UInt timeLevel:
     *		-1: get value from last ts
     *    0: get current value
     *    1: get value from last iteration
     *
     * Flags:
     *  bool invert_
     *    true: use inverse hyst operator if possible; otherwise perform iterative inversion
     *    false: use standard hyst operator
     *
     * Return value:
     *  Vector defining the state of the hyst operator
     */
    
		assert(XMLParameterSet_);
		totalCallingCounter_++;
    
		/*
     * 1. preprocess lpm to get information about the actual lpm where we have to
     *  evaluate the hyst operator, the index of the operator and the index of the storage array
     */
		UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;
    
		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
    
		/*std::cout << "timeLevel " << timeLevel << std::endl;
     std::cout << "storageIdx: " << storageIdx << std::endl;
     std::cout << "size of arrays: " << P_J_lastTS_->GetSize() << std::endl;
     */
		/*
     * 2. check time level and need for evaluation
     */
		if(timeLevel == -2){
			Vector<Double> zeroVec = Vector<Double>(dim_);
			zeroVec.Init();
			return zeroVec;
		} else if (timeLevel == -1) {
			// return value from last time step
      //std::cout << "Get value from last TS: " << P_J_lastTS_[storageIdx].ToString() << std::endl;
			return P_J_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
			// return value from last iteration
			return P_J_lastIt_[storageIdx];
		} else {
			// get current state; here we have to check if a reevaluation is needed
			// note: requiresReeval_ has one entry for each storage, not for each operator
			if (false == requiresReeval_[storageIdx]) {
				//std::cout << "no reeval needed as flag is false!" << std::endl;
				// return current value
				return P_J_[storageIdx];
			}
		}
    
		/*
     * 3. Evaluate hysteresis operator
     */
		// first get input
		Vector<Double> input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx, onBoundary);
    
		// then compute output
		return CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx);
    
	}
  
	Vector<Double> CoefFunctionHyst::RetrieveLPMSolution(LocPointMapped& actualLPM, UInt storageIdx, int timeLevel, bool onBoundary) {
    
		if (timeLevel == -1) {
			// get solution from last ts
			return E_B_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
			// get solution from last iteration
			return E_B_lastIt_[storageIdx];
		} else {
			// get current input from system
			Vector<Double> LPMSolution;
			if(onBoundary){
				if(dependCoef1Surf_ != NULL){
					LPMSolution = Vector<Double>(dependCoef1Surf_->GetVecSize());
					dependCoef1Surf_->GetVector(LPMSolution, actualLPM);
				} else {
					EXCEPTION("LPM solution on boundary requested, but coefFunction on boundary not set.")
				}
				//std::cout << "Retrieved input ON boundary: " << LPMSolution.ToString() << std::endl;
				
				// for fieldParallel boundary conditions, En should be zero, Pn not
				// (En should be zero due to NeumannBC but this does not seem to work, so 
				// we subtract might have to subtract normal part)
				//bool tangentialOnly = true;
				
			} else {
				LPMSolution = Vector<Double>(dependCoef1_->GetVecSize());
				dependCoef1_->GetVector(LPMSolution, actualLPM);
			}
			return LPMSolution;
		}
	}
  
	Vector<Double> CoefFunctionHyst::RetrieveInputToHysteresisOperator(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool onBoundary) {
    //std::cout << "RetrieveInputToHysteresisOperator" << std::endl;
		// get current solution first
		// as we want the current input only, set timeLevel to 0
		UInt timeLevel = 0;
		Vector<Double> curLPMSolution = RetrieveLPMSolution(actualLPM, storageIdx, timeLevel, onBoundary);
    
		if (needsInversion_ && (hasInverseModel_ == false) ) {
      //std::cout << "Inversion needed" << std::endl;
      // we have an inverse field (e.g. magnetics) but no inverse model (i.e. VectorPreisach)
      //
			// check if solution did change since last time > if not, hope that input / output pair of hyst
			// operator can be reused
			Vector<Double> diff = curLPMSolution;
			Vector<Double> toDiff = E_B_[storageIdx];
			diff -= toDiff;
      
      //std::cout << "StorageIDX: " << storageIdx << std::endl;
      //std::cout << "current Solution (input): " << curLPMSolution.ToString() << std::endl;
      //std::cout << "old Solution (loaded): " << E_B_[storageIdx] << std::endl;
      //std::cout << "difference: " << diff.ToString() << std::endl;
      //std::cout << "diff.NormL2(): " << diff.NormL2() << std::endl;
      Vector<Double> retrievedInput;
			if (diff.NormL2() < MAT_ampResolution_) {
        //std::cout << "(Nearly) no difference in amplitude to previous element solution" << std::endl;
        retrievedInput = E_H_[storageIdx];
				//}
			} else {
        //bool overwriteMemory = OverwriteHystMemory();
        // for input computation never overwrite storage
        // we basically search in the current state for a fitting solution
        // then, during the step calcoutput we evaluate the hyst operator with
        // the retrieved input and here we actually overwrite the memory, if necessary
        bool overwriteMemory = false;
        
        if(MAT_methodType_ == SCALAR){
          retrievedInput = Vector<Double>(dim_);
          retrievedInput[MAT_dirP_] = hyst_->computeInputAndUpdate(curLPMSolution[MAT_dirP_], MAT_smallSignalTensor_[MAT_dirP_][MAT_dirP_], operatorIdx, overwriteMemory);
        } else {
          // Idea: use Levenberg Marquart to estimate  inversion of hyst operator
          retrievedInput = hyst_->computeInput_vec(curLPMSolution, operatorIdx, MAT_smallSignalTensor_,
                  alphaLinesearch_, overwriteMemory, RUN_overwriteDirection_ );
        }
        //std::cout << "Computed input: " << retrievedInput.ToString() << std::endl;
        // we should not store the value here but first after evaluating the hyst operator in calcOutput function
        // otherwise we will never get to the point where P_J_ gets overwritten
        // and thus we get no real stepping at all
        //E_H_[storageIdx] = retrievedInput;
			}
      // store current solution for later checks
      E_B_[storageIdx] = curLPMSolution;
      
      //std::cout << "Returned input: " << retrievedInput << std::endl;
      return retrievedInput;
      
		} else {
      // store current solution for later checks
      E_B_[storageIdx] = curLPMSolution;
      //std::cout << "Returned input: " << curLPMSolution << std::endl;
			// for electrostatic / piezoelectric material, we can use the element solution (E) as input
			return curLPMSolution;
		} 
	}
  
	Vector<Double> CoefFunctionHyst::CalcOutputOfHysteresisOperator(Vector<Double> inputToHystOperator, UInt operatorIdx, UInt storageIdx, 
          bool forceMemoryLock, bool forceMemoryWrite) {
    
    //std::cout << "CalcOutputOfHysteresisOperator" << std::endl;
    //std::cout << "InputToHystOperator: " << inputToHystOperator.ToString() << std::endl;
    
		Vector<Double> outputOfHystOperator;
    
		/*
     *  overwrite the hysteresis memory
     *          or
     *  work on temporal storage
     */
    // during standard evaluation, the helperfunction OverwriteHystMemory()
    // will return the correct flag; however, for setting the previous hyst value
    // it might be that we are forced to write or lock the memory
    bool overwriteMemory;
    if(forceMemoryLock == true){
      overwriteMemory = false;
    } else if(forceMemoryWrite == true){
      overwriteMemory = true;
    } else {
      overwriteMemory = OverwriteHystMemory();
    }
    
		/*
     * Check if a reevaluation really is required by comparing the actual input to
     * the last used input
     */
		Vector<Double> diff = inputToHystOperator;
		Vector<Double> toDiff = E_H_[storageIdx];
		diff -= toDiff;
    
		if ((diff.NormL2() < MAT_ampResolution_)&&(overwriteMemory == false)) {
      //	if (XML_textOutputLevel_ == 2) {
      //std::cout << "(Nearly) no difference in amplitude to previous input" << std::endl;
      //std::cout << "Input: " << inputToHystOperator.ToString() << std::endl;
      //std::cout << "Old: " << toDiff.ToString() << std::endl;
      //std::cout << "Apmplitude resolution: " << MAT_ampResolution_ << std::endl;
      //	}
      //std::cout << "return P_J_[storageIdx]: " << P_J_[storageIdx].ToString() << std::endl;
			return P_J_[storageIdx];
		}
    //std::cout << "Evaluate Hysteresis operator" << std::endl;
		/*
     * Call hystersis operator
     */
		totalEvaluationCounter_++;
    
		if (MAT_methodType_ == SCALAR) {
      
			outputOfHystOperator = Vector<Double>(dim_);
			outputOfHystOperator.Init();
      
			if (XML_performanceMeasurement_) {
				timer_->Start();
			}
      
      //std::cout << "Used input for SCALAR model: " << inputToHystOperator[MAT_dirP_] << std::endl;
      
			outputOfHystOperator[MAT_dirP_] = hyst_->computeValueAndUpdate(inputToHystOperator[MAT_dirP_], operatorIdx, overwriteMemory);
      
			if (XML_performanceMeasurement_) {
				timer_->Stop();
				//        std::cout << "Scalar Preisach operator" << std::endl;
				//        std::cout << "Total time (" << totalEvaluationCounter_ << "calls): " << timer_->GetCPUTime() << std::endl;
				//        std::cout << "Avg Time: " << timer_->GetCPUTime()/totalEvaluationCounter_ << std::endl;
				totalEvaluationTime_ = timer_->GetCPUTime();
				avgEvaluationTime_ = timer_->GetCPUTime() / totalEvaluationCounter_;
			}
      
		} else {
			if (XML_performanceMeasurement_) {
				timer_->Start();
			}
      //std::cout << "Comput hyst output using Vector model " << std::endl;
			//std::cout << "OperatorIdx: " << operatorIdx << std::endl;
			outputOfHystOperator = hyst_->computeValue_vec(inputToHystOperator, operatorIdx, overwriteMemory, RUN_overwriteDirection_);
      //std::cout << "Output: " << outputOfHystOperator.ToString(2) << std::endl;
			if (XML_performanceMeasurement_) {
				timer_->Stop();
				//        std::cout << "Vector Preisach operator" << std::endl;
				//        std::cout << "Total time (" << totalEvaluationCounter_ << "calls): " << timer_->GetCPUTime() << std::endl;
				//        std::cout << "Avg Time: " << timer_->GetCPUTime()/totalEvaluationCounter_ << std::endl;
				totalEvaluationTime_ = timer_->GetCPUTime();
				avgEvaluationTime_ = timer_->GetCPUTime() / totalEvaluationCounter_;
			}
      
			/*!
       * Print out of Hysteresis state
       * ONLY FOR DEBUGGING REASONS!
       * USE CAREFULLY! MASSIVE IMAGE OUTPUT!
       * only for vector case
       */
			static UInt cnt = 0;
			static UInt firstIdx = 1;
      
			if ((MAT_printOut_ > 0) && (RUN_allowBMP_ == true)) {
				if ((cnt % MAT_printOut_ == 0)&&(operatorIdx == firstIdx)) {
					//std::cout << "Outputting bmp" << std::endl;
					std::stringstream filenamebuf;
					filenamebuf << "Switch_Elem" << firstIdx << "_Step" << std::setfill('0') << std::setw(5) << cnt << "_v" << MAT_vecPreisachImplementationVersion_ << "_numRows" << MAT_numRows_ << ".bmp";
					hyst_->switchingStateToBmp(MAT_bmpResolution_, filenamebuf.str(), operatorIdx, true);
				}
        
				if (operatorIdx == firstIdx) {
					cnt++;
					/*
           * disable output until reset
           * -> otherwise we would get two images for each timestep if P and D are computed
           */
					RUN_allowBMP_ = false;
				}
			}
		}
    
    //std::cout << "Store input/output combination for later usage" << std::endl;
		E_H_[storageIdx] = inputToHystOperator;
    //std::cout << "P_J_[storageIdx] before setting: " << P_J_[storageIdx].ToString() << std::endl;
		P_J_[storageIdx] = outputOfHystOperator;
    //std::cout << "P_J_[storageIdx] after setting: " << P_J_[storageIdx].ToString() << std::endl;
    
    //std::cout << "Computed output: " << outputOfHystOperator.ToString() << std::endl;
		// this flag has to be reset by hand!
		requiresReeval_[storageIdx] = false;
    
		return outputOfHystOperator;
    
	}
  
	void CoefFunctionHyst::SetPreviousHystVals(bool setLastTS, bool forceMemoryLock) {
    
//		if(XML_textOutputLevel_ == 2){
//      std::cout << "++ SetPreviousHystVals ++" << std::endl;
//      std::cout << "lastTS? " << setLastTS << std::endl;
//		}
    
		/*
     * Function to backup input/output pair as well as current deltaMatrix;
     * This function is called at the beginning of each iteration BEFORE the
     * system gets solved. Therewith, the backup always store the values of
     * the last iteration. These values are needed for the computation of the
     * new deltaMatrix. Evaluation of the hysteresis operator is done onto a temp.
     * copy, i.e. the actual memory is not changed).
     *
     * If setLastTS_ is true, the current state is instead stored to
     * the arrays ..._lastTS_ and the evaluation is performed on the actual storage
     * (i.e. the hysteresis operator is set). This call should only be performed
     * after the end of a timestep.
     *
     * Depending on the XML_EvaluationDepth_ this function does the following steps:
     * depth = 1 or 2 (standard or extended):
     *  iterate over all elements
     *    extract element index
     *    extract solution and input to hysteresis operator at midpoint of element
     *    evaluate hysteresis operator with the extracted input
     *    store extracted solution, extracted input and evaluated output of hystoperator
     *     and store current state of deltaMatrix (no reevaluation)
     *
     * depth = 3 (full evalution)
     *  iterate over all elements
     *    extract element index
     *    get integration points for element
     *
     *    iterate over all integration points
     *      compute corresponding index (unique combination of element index and integration point index)
     *      extract solution and input to hystoperator at integration point
     *      evaluate hyst operator with the determined index
     *      store extracted solution, input and output and store current state of
     *        deltaMattrix
     *
     *    additionally: get coordinates of element midpoint (is normally no integration point)
     *    extract solution and input at midpoint
     *    evaluate hyst operator with index of midpoint (again a unique combination of element index + point index)
     *    store extracted solution, ...
     *
     * Note: the additional evaluation at the midpoint is needed for the output computation as we
     *  evaluate the hyst operator at the midpoint; therefore we need to know the state at the midpoint
     *
     */
    
		Vector<Double> solution = Vector<Double>(dim_);
		Vector<Double> input = Vector<Double>(dim_);
		Vector<Double> output = Vector<Double>(dim_);
		UInt storageIdx;
		UInt operatorIdx;
    
		/*
     * 0.1 backup current evaluationPurpose
     */
		UInt evalPurposeBackup = RUN_evaluationPurpose_;
    
		/*
     * 0.2 set evaluation purpose to
     *  2 if setLastTS is false
     *  3 if setLastTS is true
     *
     * by setting the evaluation purpose, we automatically will get the
     * right values for midpointOnly and OverwriteHystMemory from the
     * corresponding functions
     */
		if (setLastTS) {
			// overwriteMemory = true
			// for extended evaluation (where we also evaluate at each integration point)
			// we only overwrite memory for the midpoint (where the hystoperator is located)
			RUN_evaluationPurpose_ = 3;
		} else {
			// overwriteMemory = false
			RUN_evaluationPurpose_ = 2;
		}
    
		bool overwriteMemoryIntPoints = false;
		if (XML_EvaluationDepth_ == 3) {
      // only for full integraation (where we have a hyst operator per integration
      // point) we overwrite the memory at each of these points
			overwriteMemoryIntPoints = true;
		}
		bool overwriteMemoryMidPoint = OverwriteHystMemory();
		bool midpointOnly = EvaluateAtMidpointOnly();
    
		if (forceMemoryLock) {
			overwriteMemoryIntPoints = false;
			overwriteMemoryMidPoint = false;
		}
    
		bool onBoundary;
		
		//std::cout << "EvaluationDepth: " << XML_EvaluationDepth_ << std::endl;
		//std::cout << "RUN_evaluationPurpose_: " << RUN_evaluationPurpose_ << std::endl;
		//std::cout << "midpointOnly: " << midpointOnly << std::endl;
		//std::cout << "overwriteMemoryMidPoint: " << overwriteMemoryMidPoint << std::endl;
    UInt outerListCNT = 0;
		
		
		/*
     * 1. iterate over all elements
     */
    for (std::list<shared_ptr<EntityList> >::iterator listIt=SDList_List_.begin(); listIt != SDList_List_.end(); ++listIt) {
      shared_ptr<EntityList> SDList_ = *listIt;
			//std::cout << "outerListCNT: " << outerListCNT << std::endl;
			outerListCNT++;
			
      EntityIterator it = SDList_->GetIterator();
      LocPoint lp;
      LocPointMapped lpm;
      const Elem * el;
      shared_ptr<ElemShapeMap> esm;
      Vector<Double> zeroVec = Vector<Double>(dim_);
      zeroVec.Init();
      LocPointMapped actualLPM;
			
      for (it.Begin(); !it.IsEnd(); it++) {
				
        /*
         * 1.0 get element and shape map
         */
        el = it.GetElem();
        esm = it.GetGrid()->GetElemShapeMap(el, true);
				//std::cout << "ElemType: " << el->type << std::endl;
        /*
         * 1.1a full evaluation > get integration points and iterate over them
         *                        this is done in addition to the evalution at the midpoint
         */
        if (midpointOnly == false) {
          // Get integration points
          IntegOrder order;
          IntScheme::IntegMethod method;
          StdVector<LocPoint> intPoints;
          StdVector<Double> weights;
					
          ptFeSpace_->GetFe(it, method, order);
					
          ptFeSpace_->GetIntScheme()->GetIntPoints(Elem::GetShapeType(el->type), method, order,
                  intPoints, weights);
					
          // Loop over all integration points
          const UInt numIntPts = intPoints.GetSize();
          for (UInt i = 0; i < numIntPts; i++) {
						
            // Calculate for each integration point the LocPointMapped
            lpm.Set(intPoints[i], esm, weights[i]);
            
            //   preprocess LPM to get correct operator and storage indices
            onBoundary = PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx);
						//std::cout << "Integration point NR: " << i << std::endl;
						//std::cout << "On boundary? " << onBoundary << std::endl;
						
            // get input
            input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx,onBoundary);
            // new: retrieveinput will store the actual solution directly to E_B_
            solution = E_B_[storageIdx];
						
            // then compute output
            output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, overwriteMemoryIntPoints);
						
            // OLD
            //          					// Calculate for each integration point the LocPointMapped
            //          					lpm.Set(intPoints[i], esm, weights[i]);
            //          
            //          					// extract solution and input for hystoperator
            //          					ExtractSolutionAndInputForHystOperator(solution, input, operatorIdx, storageIdx, lpm, false);
            //          
            //          					// evaluate at corresponding operator index
            //				EvaluateHysteresisOperator(input, output, operatorIdx, storageIdx, overwriteMemoryIntPoints, solution);
						
            // store at storage index
            if (setLastTS) {
              //std::cout << "Set values for integration point " << storageIdx << std::endl;
              E_B_lastTS_[storageIdx] = solution;
              //std::cout << "P_J_lastTS_[storageIdx] before setting: " << P_J_lastTS_[storageIdx] << std::endl;
              P_J_lastTS_[storageIdx] = output;
              //std::cout << "P_J_lastTS_[storageIdx] after setting: " << P_J_lastTS_[storageIdx] << std::endl;
              E_H_lastTS_[storageIdx] = input;
              deltaMat_lastTS_[storageIdx] = deltaMat_[storageIdx];
              dY_sol[storageIdx] = zeroVec;
              X_low[storageIdx] = zeroVec;
              X_up[storageIdx] = zeroVec;
            } else {
              E_B_lastIt_[storageIdx] = solution;
              P_J_lastIt_[storageIdx] = output;
              E_H_lastIt_[storageIdx] = input;
              deltaMat_lastIt_[storageIdx] = deltaMat_[storageIdx];
            }
          }
        }
        /*
         * 1.1b evaluation at midpoint > done for standard, extended and full evalutation
         */
        lp = Elem::shapes[el->type].midPointCoord;
        lpm.Set(lp, esm, 0.0);
        
				//std::cout << "Midpoint " << std::endl;
				//std::cout << "On boundary? " << onBoundary << std::endl;
        //      // preprocess LPM to get correct operator and storage indices
        bool forceMidpoint = true;
        onBoundary = PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx,forceMidpoint);
				
        // get input
        input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx, onBoundary);
				
        // new: retrieveinput will store the actual solution directly to E_B_
        solution = E_B_[storageIdx];
				
        // then compute output
        output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, overwriteMemoryMidPoint);
				
        // OLD
        //      			// extract solution and input for hystoperator
        //      			ExtractSolutionAndInputForHystOperator(solution, input, operatorIdx, storageIdx, lpm, true);
        //      
        //      			// evaluate at corresponding operator index
        //      			EvaluateHysteresisOperator(input, output, operatorIdx, storageIdx, overwriteMemoryMidPoint, solution);
				
        // store at storage index
        if (setLastTS) {
          E_B_lastTS_[storageIdx] = solution;
          P_J_lastTS_[storageIdx] = output;
          E_H_lastTS_[storageIdx] = input;
          deltaMat_lastTS_[storageIdx] = deltaMat_[storageIdx];
          dY_sol[storageIdx] = zeroVec;
          X_low[storageIdx] = zeroVec;
          X_up[storageIdx] = zeroVec;
        } else {
          E_B_lastIt_[storageIdx] = solution;
          P_J_lastIt_[storageIdx] = output;
          E_H_lastIt_[storageIdx] = input;
          deltaMat_lastIt_[storageIdx] = deltaMat_[storageIdx];
        }
      } //actsdlist
    } // sdlistlist
    //	if (setLastTS) {
    //EstimateCurrentSlope(Vector<Double>(dim_), 1e-4);
    //EstimateCurrentSlope(Vector<Double>(dim_), 1e-8);
    //EstimateCurrentSlope(Vector<Double>(dim_), 1e-10);
    //	}
    
    
		/*
     * 2.0 restore old state of evaluationPurpose
     */
		RUN_evaluationPurpose_ = evalPurposeBackup;
    
	}
  
	void CoefFunctionHyst::GetScalar(Double& outputScalar, const LocPointMapped& lpm) {
		EXCEPTION("GetScalar not implemented for coefFncHyst");
	}
  
	void CoefFunctionHyst::GetVector(Vector<Double>& outputVector, const LocPointMapped& lpm) {
    
		// return current state of hyst operator
		UInt timeLevel = 0;
		outputVector = GetOutputOfHysteresisOperator(lpm, timeLevel);
	}
  
	void CoefFunctionHyst::GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm) {
    EXCEPTION("GetTensor not implemented for coefFncHyst");
	}
  
  void CoefFunctionHyst::SetFlag(std::string flagName, Integer intState) {
    
    /*
     * Flags that (may) change during runtime
     */
    if(flagName == "evaluationPurpose"){
      RUN_evaluationPurpose_ = intState;
    }
    else if(flagName == "allowBMP"){
      RUN_allowBMP_ = bool(intState);
    } 
    else if(flagName == "forceCurrentTS"){
      forceCurrentTS_ = bool(intState);
    } 
    else if(flagName == "SetTimeLevelRHSHyst"){
      timeLevel_RHS_ = intState;
    }
    else if(flagName == "SetTimeLevelMaterial"){
      timeLevel_Material_ = intState;
    }
    else if(flagName == "SetTimeLevelDeltaMat"){
      timeLevel_DeltaMat_ = intState;
      if(timeLevel_DeltaMat_ == -2){
        deltaMatActive_ = false;
      } else {
        deltaMatActive_ = true;
      }
    }
    else if(flagName == "resetReeval"){
      //std::cout << "Reset reeval flag" << std::endl;
      for (UInt i = 0; i < numStorageEntries_; i++) {
        requiresReeval_[i] = true;
        Si_requiresReeval_[i] = true;
        deltaMat_requiresReeval_[i] = true;
        deltaMatStrain_requiresReeval_[i] = true;
        rotatedCouplingTensor_requiresReeval_[i] = true;
      }
    }
    /*
     * Flags that are to be set from xml file
     */
    else if(flagName == "deltaForm"){
      // ==0 > no deltaMat 
      // !=0 > deltaMat gets computed
      RUN_deltaForm_ = intState;
    }
    else if(flagName == "evaluationDepth"){
      /*
       * 1 -> standard evaluation
       * 2 -> extended evaluation
       * 3 -> full evaluation
       *
       * -> details: see .hh file
       */
      if (intState < 1) {
        XML_EvaluationDepth_ = 1;
      } else if (intState > 3) {
        XML_EvaluationDepth_ = 3;
      } else {
        XML_EvaluationDepth_ = intState;
      }

      // the evaluation depth is crucial for the number of storage entries
      // > init only after evalDepth is set!
      InitStorage();
      XMLParameterSet_ = true;
    }
    else if(flagName == "measurePerformance"){
      XML_performanceMeasurement_ = intState; 
      hyst_->setFlags(XML_performanceMeasurement_);
    }    
    else if(flagName == "testInversion"){
      XML_testInversion_ = intState; 
      if(storageInitialized_ == true){
        TestInversion(MAT_smallSignalTensor_);
      } else {
        WARN("Could not perform inversion tests as storage was not iniitalized yet.");
      }
    } 
    else {
      std::stringstream except;
      except << "Flag " << flagName << " unknown to CoefFncHyst.";
      EXCEPTION(except.str())
    } 
  }

	bool CoefFunctionHyst::EvaluateAtMidpointOnly() {
		/*
     *  XML_EvaluationDepth_   >    1 (standard)    2 (extended)                          3 (full)
     *  RUN_evaluationPurpose_ v
     *    1 (assemble)          midpoint            integration point                     integration point
     *    2 (store lastIt)      midpoint            integration point + midpoint          integration point + midpoint
     *    3 (store lastTS)      midpoint            integration point + midpoint          integration point + midpoint
     *    4 (output)            midpoint            midpoint                              midpoint
     *
     */
		if (XML_EvaluationDepth_ == 1) {
			return true;
		} else if (XML_EvaluationDepth_ == 2) {
			// extended
			// NEW: even for extended evaluation, we have to store values at each integration point
			// (to get correct deltaMatrices)
			// but in difference to full evaluation, we lock the hysteresis operator for all points except the midpoint!
			if (RUN_evaluationPurpose_ == 4) {
				return true;
			} else {
				return false;
			}
		} else {
			// full
			if (RUN_evaluationPurpose_ <= 3) {
				return false;
			} else {
				return true;
			}
		}
	}
  
	bool CoefFunctionHyst::OverwriteHystMemory() {
		/*
     *  RUN_evaluationPurpose_
     *    1 (assemble)                false
     *    2 (store lastIT)            false
     *    3 (store lastTS)            true
     *    4 (output)                  false
     *
     *  Motivation:
     *    During the nonlinear solution process, we oftentimes overestimate
     *    the input to the hysteresis operator and have to step back.
     *    These stepping (especially during linesearch) lead to permanent changes
     *    of the hysteresis memory which are (normally) not wanted.
     *    Therefore, the hysteresis memory should be locked, i.e. all changes
     *    are done to a temporal copy of the current state.
     *    This leads to a nonlinear but nonhysteretic stepping.
     *
     *    At the end of each timestep (where we assume to have found a valid solution)
     *    we finally unlock the hysteresis memory and evaluate the hysteresis operator
     *    once more at the last computed solution to bring the operator to the correct
     *    state.
     *
     *    Conclusion:
     *      during assemble and storage (of last iteration values) and during
     *      output computation, we lock the hysteresis memory
     *      (output computation would not be the problem as we would simply
     *      apply the same input as during the storage of the lastTS value)
     *
     *      for storing the lastTS values), we unlock it
     *
     *
     */
		if (RUN_evaluationPurpose_ == 3) {
			return true;
		} else {
			return false;
		}
    
	}
  
  void CoefFunctionHyst::TestInversion(Matrix<Double> eps_mu){
    
    if(XML_testInversion_ == 0){
      return;
    }
    
    // Tests inversion of VECTORPreisach Operator and returns 
    //  estimation for alpha (linesearch stepping paramater) that can be used
    //  for later computation
    // > can be used in actual simulations to setup an appropriate linesearch factor alpha
    //   to accelerate the actual computations later on!
    
    Vector<Double>* xIn;
    bool vector;
    bool isVirgin = true;
    
    if (MAT_methodType_ == SCALAR) {
      vector = false;
      
      // create one temoporary ScaparPreisach operator with the same
      // paramter as the later used ones, except, we just need one entry      
      hystTMP = new Preisach(1, MAT_xSat_, MAT_ySat_, MAT_PreisachWeights_, isVirgin);
      
    } else {
      vector = true;
      
      // create one temporary VectorPreisach operator with the same
      // parameter as the later used ones, except, we just need it to store
      // entries for one single element / point
      if (MAT_vecPreisachImplementationVersion_ == 1) {
        isClassical_ = true; // original vector preisach model -> sutor2012
        
        hystTMP = new VectorPreisachv10_ListApproach(1, MAT_xSat_, MAT_ySat_,
                MAT_PreisachWeights_, MAT_rotRes_, dim_, isVirgin,
                isClassical_, MAT_angResistance_, MAT_angClipping_);
      } else if (MAT_vecPreisachImplementationVersion_ == 2) {
        isClassical_ = false; // revised vector preisach model -> sutor2015
        
        hystTMP = new VectorPreisachv10_ListApproach(1, MAT_xSat_, MAT_ySat_,
                MAT_PreisachWeights_, MAT_rotRes_, dim_, isVirgin,
                isClassical_, MAT_angResistance_, MAT_angClipping_);
      } else if (MAT_vecPreisachImplementationVersion_ == 10) {
        isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
        
        hystTMP = new VectorPreisachv10_MatrixApproach(1, MAT_xSat_, MAT_ySat_,
                MAT_PreisachWeights_, MAT_rotRes_, dim_, isVirgin,
                isClassical_, MAT_angResistance_, MAT_angClipping_);
      } else if (MAT_vecPreisachImplementationVersion_ == 20) {
        isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
        
        hystTMP = new VectorPreisachv10_MatrixApproach(1, MAT_xSat_, MAT_ySat_,
                MAT_PreisachWeights_, MAT_rotRes_, dim_, isVirgin,
                isClassical_, MAT_angResistance_, MAT_angClipping_);
      } else {
        EXCEPTION("MAT_vecPreisachImplementationVersion_ has to be one of the following: \n "
                "1: classical vector model (sutor2012) \n"
                "2: revised vector model (sutor2015) [DEFAULT] \n"
                "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
                "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
      }
    }
    
    UInt numCases = 2;
    UInt totalSteps;
    for(UInt ca = 0; ca < numCases; ca++){
      if(ca == 0){
        /*
         * Testcase 1 --- (self-designed) benchmark signal
         * 
         * 1d-Signal shape:
         *  exp-increase, cos decrease, linear increase, 1/x decrease, log increase, x^3 increase, hold state
         * 
         * 2d/3d-Signal shape:
         *  Multiply 1d-Signal with chaning direction vector
         * 
         */
        
        UInt steps1,steps2,steps3,steps4,steps5,steps6,steps7;
        steps1=50;
        steps2=50;
        steps3=50;
        steps4=50;
        steps5=50;
        steps6=50;
        steps7=10;
        totalSteps=steps1+steps2+steps3+steps4+steps5+steps6+steps7;
        
        xIn = new Vector<Double>[totalSteps]; 
        Vector<Double> xScal = Vector<Double>(totalSteps);
        xScal.Init();
        Double delta;
        
        // 1. expontential increase from [10^-15 to 10^0[
        delta = 15.0/steps1;
        
        for(UInt i = 0; i < steps1; i++){
          xScal[i] = std::pow(10,-15.0+i*delta);
        }
        // 2. cos decrease towards -0.5[
        delta = -std::acos(-0.5)/steps2;
        
        for(UInt i = 0; i < steps2; i++){
          xScal[steps1+i] = std::cos(i*delta);
        }
        // 3. linear increase to 0.4
        delta = 0.9/steps3;
        
        for(UInt i = 0; i < steps3; i++){
          xScal[steps1+steps2+i] = -0.5+i*delta;
        }
        // 4. 1/x decrease towards -0.3
        Double Xstart = 1.0/(0.7 + 1.0/steps4);
        delta = (steps4 - Xstart)/steps4;
        
        for(UInt i = 0; i < steps4; i++){
          xScal[steps1+steps2+steps3+i] = 1.0/(Xstart+i*delta)-1.0/steps4-0.3;
        }
        // 5. logarithmic increase towards 0.2
        Xstart = std::exp(-3);
        Double Xend = std::exp(2);
        delta = (Xend-Xstart)/steps5;
        
        for(UInt i = 0; i < steps5; i++){
          xScal[steps1+steps2+steps3+steps4+i] = 0.1*std::log(Xstart+i*delta);
        }
        // 6. x^3 increase to 0.8
        Double denom = std::pow(steps6-1,3);
        for(UInt i = 0; i < steps6; i++){
          xScal[steps1+steps2+steps3+steps4+steps5+i] = 0.6*std::pow(i,3)/denom + 0.2;
        }
        
        // 7. hold signal at value 0.8
        for(UInt i = 0; i < steps7; i++){
          xScal[steps1+steps2+steps3+steps4+steps5+steps6+i] = 0.8;
        }
        
        //        std::ofstream testOutput;
        //        testOutput.open("GeneratedInput");
        //        for(UInt i = 0; i < totalSteps; i++){
        //          testOutput << i << " " << xScal[i] << std::endl;
        //        }
        //        testOutput.close();
        Vector<Double> dir = Vector<Double>(dim_);
        dir[0] = 1.0;
        dir[1] = 0.0;
        for(UInt i = 0; i < totalSteps; i++){
          xIn[i] = Vector<Double>(dim_);
          xIn[i].Init();
          dir[0] = (Double) (totalSteps-i)/totalSteps;
          dir[1] = (Double) i/totalSteps;
          xIn[i].Add(MAT_xSat_*xScal[i],dir);
        }
        
      } else {
        totalSteps = 200;
        xIn = new Vector<Double>[totalSteps]; 
        /*
         * Testcase: Put material into full saturation in y-direction; from that state on
         * apply a decreasing triangular signal in x-direction; for initial period hold
         * signla in y-direction; during first real period decrease to 1/2 of value;
         * then to 1/4 and finally to 0
         */
        xIn[0] = Vector<Double>(dim_);
        xIn[0].Init(0.0);
        xIn[0][1] = -MAT_xSat_;
        
        UInt numStepsFirstIncrease = (UInt) ( (Double)totalSteps/8.0 ) +1;
        Double incr = MAT_xSat_/((Double) numStepsFirstIncrease-1);
        
        for(UInt i = 1; i < numStepsFirstIncrease; i++){
          xIn[i] = Vector<Double>(dim_);
          xIn[i].Init(0.0);
          xIn[i][0] = xIn[i-1][0] + incr; 
          if(vector){
            xIn[i][1] = -MAT_xSat_;
          }
        }
        
        UInt cnt = numStepsFirstIncrease;
        
        UInt remainingSteps = totalSteps - numStepsFirstIncrease;
        UInt stepsPerPeriod = 20;
        UInt numPeriods = remainingSteps/stepsPerPeriod;
        remainingSteps = remainingSteps%stepsPerPeriod;
        
        incr = 2*MAT_xSat_/((Double)(stepsPerPeriod-1));
        Double sign = -1.0;
        Double decrFactor = 1.0 - 1.0/((Double) numPeriods);
        
        for(UInt j = 0; j < numPeriods; j++){
          for(UInt i = 0; i < stepsPerPeriod; i++){
            xIn[cnt] = Vector<Double>(dim_);
            xIn[cnt].Init(0.0);
            xIn[cnt][0] = xIn[cnt-1][0] + sign*incr;
            
            if(vector){
              if(j==0){
                xIn[cnt][1] = -MAT_xSat_/2.0;
              } else if(j==1){
                xIn[cnt][1] = -MAT_xSat_/4.0;
              }
            }
            
            cnt++;
          }
          sign = -1.0*sign;
          incr = incr*decrFactor;
        }
        for(UInt i = 0; i < remainingSteps; i++){
          xIn[cnt] = Vector<Double>(dim_);
          xIn[cnt].Init(0.0);
          xIn[cnt][0] = xIn[cnt-1][0] + sign*incr;
          cnt++;
        }
      }
      
      std::stringstream stringstream1;
      stringstream1 << "xInput_" << ca;
      std::string filenameInput = stringstream1.str();
      std::stringstream stringstream2;
      stringstream2 << "xRetrieved_" << ca;
      std::string filenameOutput = stringstream2.str();
      std::stringstream stringstream3;
      stringstream3 << "AllResults_" << ca;
      std::string filenameOutputFull = stringstream3.str();
      
      std::ofstream testOutputSig;
      testOutputSig.open(filenameInput);
      for(UInt i = 0; i < totalSteps; i++){
        testOutputSig << i << " " << xIn[i][0] << " " << xIn[i][1] << std::endl;
      }
      testOutputSig.close();
      
      
      Vector<Double>* xRetrieved = new Vector<Double>[totalSteps];
      Vector<Double>* yRetrieved = new Vector<Double>[totalSteps];
      Vector<Double>* yError = new Vector<Double>[totalSteps];
      Vector<Double>* xError = new Vector<Double>[totalSteps];
      Vector<Double>* yIn = new Vector<Double>[totalSteps]; 
      
      std::ofstream testOutput;
      std::ofstream testOutput2;
      testOutput2.open(filenameOutput);
      testOutput.open(filenameOutputFull);
      testOutput << "## Results of TestInversion" << std::endl; 
      testOutput << "# Number of tests: " << totalSteps << std::endl;
      testOutput << "# eps_mu: " << eps_mu.ToString() << std::endl;
      testOutput << "# xSat: " << MAT_xSat_ << std::endl;
      testOutput << "# ySat: " << MAT_ySat_ << std::endl;
      testOutput << "# Starting value of inversion " << xIn[0][0] << " " << xIn[0][1] << std::endl;
      testOutput << "# " << std::endl;
      testOutput << "# Testnumber    x_sol[0]    x_sol[1]    x_retrieved[0]    x_retrieved[1]    x_error[0]    x_error[1]    y_target[0]    y_target[1]    y_retrieved[0]    y_retrieved[1]    y_error[0]    y_error[1]" << std::endl;
      
      // this was true all the time but I don't get why
      // when we evaluate the hyst operator we just want to know which value
      // of y we have to solve for but we do not want to set the material already
      // if we would do so, the mateiral would already have the solution inside
      bool overwriteMemory = false;
      bool overwriteDirection = true;
      yIn[0] = Vector<Double>(dim_);
      yIn[0].Init(0.0);
      // Get polarization P
      if(vector){
        yIn[0] = hystTMP->computeValue_vec(xIn[0], 0, overwriteMemory, overwriteDirection);
      } else {
        yIn[0][0] = hystTMP->computeValueAndUpdate(xIn[0][0], 0, overwriteMemory);
      }
      
      // D = eps*E + P
      for(UInt j = 0; j < dim_; j++){
        yIn[0][j] += eps_mu[j][j]*xIn[0][j];
      }
      
      xRetrieved[0] = Vector<Double>(dim_);
      xRetrieved[0].Init(0.0);
      
      overwriteMemory = false;
      if(vector){
        xRetrieved[0] = hystTMP->computeInput_vec(yIn[0], 0, eps_mu, alphaLinesearch_, overwriteMemory, overwriteDirection);
      } else {
        xRetrieved[0][0] = hystTMP->computeInputAndUpdate(yIn[0][0], eps_mu[0][0], 0, overwriteMemory);
      }
      // evaluate Preisach OP to get the actual output; DO NOT OVERWRITE MEMORY
      yRetrieved[0] = Vector<Double>(dim_);
      yRetrieved[0].Init(0.0);
      
      xError[0] = Vector<Double>(dim_);
      xError[0].Init(0.0);
      
      xError[0] = xRetrieved[0];
      xError[0] -= xIn[0];
      
      yError[0] = Vector<Double>(dim_);
      yError[0].Init(0.0);
      
      overwriteMemory = true;
      if(vector){
        yRetrieved[0] = hystTMP->computeValue_vec(xRetrieved[0], 0, overwriteMemory, overwriteDirection);
      } else {
        yRetrieved[0][0] = hystTMP->computeValueAndUpdate(xRetrieved[0][0], 0, overwriteMemory);
      }
      
      yError[0] = yRetrieved[0];
      yError[0] -= yIn[0];
      
      testOutput << std::setprecision(9) <<  "0    "<< xIn[0][0] <<"    "<< xIn[0][1] <<"    "<< xRetrieved[0][0] <<"    "<< xRetrieved[0][1]<<"    "<<xError[0][0]<<"    "<<xError[0][1]<<"    "<<yIn[0][0]<<"    "<<yIn[0][1]<<"    "<<yRetrieved[0][0]<<"    "<<yRetrieved[0][1]<<"    "<<yError[0][0]<<"    "<<yError[0][1]<< std::endl;
      testOutput2 << std::setprecision(9) <<  "0   "<< xRetrieved[0][0] <<"    "<< xRetrieved[0][1]<<std::endl;
      //    Double inc = 2*MAT_xSat_/( (Double) numTests);
      
      for(UInt i = 1; i < totalSteps; i++){
        std::cout << "##### TEST NR " << i << "#####" << std::endl;
        //      xIn[i] = Vector<Double>(dim_);
        //      xIn[i][0] = xIn[i-1][0] + inc;
        //      xIn[i][1] = xIn[i-1][0] + inc/2;
        overwriteMemory = false;
        
        if(vector){
          yIn[i] = hystTMP->computeValue_vec(xIn[i], 0, overwriteMemory, overwriteDirection);
        } else {
          yIn[i] = Vector<Double>(dim_);
          yIn[i].Init(0.0);
          yIn[i][0] = hystTMP->computeValueAndUpdate(xIn[i][0], 0, overwriteMemory);
        }
        
        for(UInt j = 0; j < dim_; j++){
          yIn[i][j] += eps_mu[j][j]*xIn[i][j];
        }
        
        overwriteMemory = false;
        
        if(vector){
          xRetrieved[i] = hystTMP->computeInput_vec(yIn[i], 0, eps_mu, alphaLinesearch_, overwriteMemory, overwriteDirection);
        } else {
          xRetrieved[i] = Vector<Double>(dim_);
          xRetrieved[i].Init(0.0);
          xRetrieved[i][0] = hystTMP->computeInputAndUpdate(yIn[i][0], eps_mu[0][0], 0, overwriteMemory);
        }
        
        yRetrieved[i] = Vector<Double>(dim_);
        yRetrieved[i].Init(0.0);
        
        yError[i] = Vector<Double>(dim_);
        yError[i].Init(0.0);
        
        overwriteMemory = true;
        if(vector){
          yRetrieved[i] = hystTMP->computeValue_vec(xRetrieved[i], 0, overwriteMemory, overwriteDirection);
        } else {
          yRetrieved[i] = Vector<Double>(dim_);
          yRetrieved[i].Init(0.0);
          yRetrieved[i][0] = hystTMP->computeValueAndUpdate(xRetrieved[i][0], 0, overwriteMemory);
        }
        
        for(UInt j = 0; j < dim_; j++){
          yRetrieved[i][j] += eps_mu[j][j]*xRetrieved[i][j];
        }
        yError[i] = yRetrieved[i];
        yError[i] -= yIn[i];
        
        xError[i] = Vector<Double>(dim_);
        xError[i].Init(0.0);
        
        xError[i] = xRetrieved[i];
        xError[i] -= xIn[i];
        
        testOutput << i<< std::setprecision(9) << "    "<<xIn[i][0]<<"    "<<xIn[i][1]<<"    "<<xRetrieved[i][0]<<"    "<<xRetrieved[i][1]<<"    "<<xError[i][0]<<"    "<<xError[i][1]<<"    "<<yIn[i][0]<<"    "<<yIn[i][1]<<"    "<<yRetrieved[i][0]<<"    "<<yRetrieved[i][1]<<"    "<<yError[i][0]<<"    "<<yError[i][1]<< std::endl;
        testOutput2 << i << std::setprecision(9) << "   "<< xRetrieved[i][0] <<"    "<< xRetrieved[i][1]<<std::endl;
        
        std::cout << "##### TARGET Y-VECTOR #####" << std::endl;
        std::cout << yIn[i].ToString() << std::endl;
        std::cout << "##### RETRIEVED Y-VECTOR #####" << std::endl;
        std::cout << yRetrieved[i].ToString() << std::endl;
        std::cout << "##### CORRECT X-VECTOR #####" << std::endl;
        std::cout << xIn[i].ToString() << std::endl;
        std::cout << "##### RETRIEVED X-VECTOR #####" << std::endl;
        std::cout << xRetrieved[i].ToString() << std::endl;
        std::cout << "##### ERROR VECTOR wrt X #####" << std::endl;
        std::cout << xError[i].ToString() << std::endl;
        std::cout << "##### ERROR VECTOR wrt Y #####" << std::endl;
        std::cout << yError[i].ToString() << std::endl;
      }
      testOutput.close();
    }
    EXCEPTION("STOP");
  }
  
  
  std::string CoefFunctionHyst::ToString() const {
    
    static UInt initialOutput = 0;
    std::ostringstream oss;
    
    if (initialOutput == 0) {
      oss << "+++ CoefFunctionHyst +++\n";
      oss << "++ General information from FE context: \n"
              << "  coefFunctionHyst created for pde " << PDEName_ << "\n"
              << "  coefFunctionHyst depends on " << dependCoef1_->ToString() << "\n";
      oss << "\n";
      
      oss << "++ Parameter extracted from input: \n"
              << "+ From material file: \n"
              << "  initial/smallfield tensor " << MAT_initialTensor_.ToString() << "\n"
              << "  xSaturation " << MAT_xSat_ << "\n"
              << "  ySaturation " << MAT_ySat_ << "\n";
      
      if (MAT_methodType_ == SCALAR) {
        oss << "  dirP " << MAT_dirP_ << "\n";
      } else {
        oss << "  rotRes " << MAT_rotRes_ << "\n"
                << "  angResistance " << MAT_angResistance_ << "\n"
                << "  angClipping " << MAT_angClipping_ << "\n"
                << "  angResolution " << MAT_angResolution_ << "\n"
                << "  amplitudeResolution " << MAT_ampResolution_ << "\n"
                << "  evalVersion " << MAT_vecPreisachImplementationVersion_ << "\n";
        oss << "  initial input to hyst operator " << MAT_initialInput_.ToString() << "\n";
        oss << "  initial output of hyst operator (exemplary for element 0) " << MAT_initialOutput_[0].ToString() << "\n";
        oss << "  printOut " << MAT_printOut_ << "\n";
        if (MAT_printOut_) {
          oss << "  bmpResolution " << MAT_bmpResolution_ << "\n";
        }
      }
      oss << "\n";
      
      oss << "+ From xml file: \n"
              << "  evaluationDepth " << XML_EvaluationDepth_ << "\n"
              << "  performanceMeasurement " << XML_performanceMeasurement_ << "\n"
              << "\n";
    }
    
    if (initialOutput > 0) {
        oss << "++ Current state of runtime parameter: \n"
                << "  evaluationPurpose " << RUN_evaluationPurpose_ << "\n"
                << "  allowBMP " << RUN_allowBMP_ << "\n";
      
      if (XML_performanceMeasurement_ >= 1) {
        oss << "++ Performance measurements:\n "
                << "  Total number of calls to EvaluateHysteresisOperator: " << totalCallingCounter_ << "\n"
                << "  Total number of actual evaluations: " << totalEvaluationCounter_ << "\n"
                << "  Total evaluation time for Hysteresis operator: " << totalEvaluationTime_ << "\n"
                << "  Average evaluation time for Hysteresis operator: " << avgEvaluationTime_ << "\n";
        if (XML_performanceMeasurement_ >= 2) {
          oss << hyst_->runtimeToString();
        }
      }
    }
    initialOutput++;
    
    return oss.str();
  }
  
}// namespace
