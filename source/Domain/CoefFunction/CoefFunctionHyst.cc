#include "CoefFunctionHyst.hh"

// classes for function / spline approximation
#include "Materials/Models/Preisach.hh"
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
		SubTensorType tensorType, MaterialType matType, shared_ptr<FeSpace> ptFeSpace,
    bool performInversionTest) : CoefFunction() {

		Init(material, actSDList, dependency1, tensorType, matType, ptFeSpace, performInversionTest);
	}

	CoefFunctionHyst::CoefFunctionHyst(BaseMaterial * const material,
		shared_ptr<ElemList> actSDList,
		PtrCoefFct dependency1, PtrCoefFct dependency2,
		SubTensorType tensorType, MaterialType matType, shared_ptr<FeSpace> ptFeSpace,
    bool performInversionTest) : CoefFunction() {

		WARN("Currently we support only single-input hysteresis operators! The second dependency will be ignored!")
		dependCoef2_ = dependency2;
		Init(material, actSDList, dependency1, tensorType, matType, ptFeSpace, performInversionTest);
	}

	void CoefFunctionHyst::Init(BaseMaterial * const material,
		shared_ptr<ElemList> actSDList,
		PtrCoefFct dependency1,
		SubTensorType tensorType,
		MaterialType matType,
		shared_ptr<FeSpace> ptFeSpace, bool performInversionTest) {

		// this type of coefficient is nonlinear (i.e. solution dependent)
		dimType_ = VECTOR;
		dependType_ = SOLUTION;
		isAnalytic_ = false;
		isComplex_ = false;
		dependCoef1_ = dependency1;
		tensorType_ = tensorType;
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
    performInversionTest_ = performInversionTest;

		// for interaction with coefFunctionHystMat
		deltaMatActive_ = true;
		hystItself_ = PtrCoefFct(this);

		/*
		 * RUN_deltaForm_ = -1 > delta = current - last ts
		 *                       rhs = last ts
		 * RUN_deltaForm_ = 0  > no delta mat
		 *                       rhs = current
		 * RUN_deltaForm_ = x  > delta between currernt and last xth iteration value
		 *                       rhs = last xth iteration value
		 */
		RUN_deltaForm_ = -1; // i.e. compute delta between last known value and last ts value
		timeLevel_Mat_ = -1; // i.e. delta wrt to last ts
		timeLevel_RHS_ = -1; // rhs term has to be last ts value as delta is towards that value
		timeLevel_BC_ = 0; // BC term has to use the last known, i.e. current value
		timeLevel_Output_ = 0; // output should always be the current value






		/*
		 * set initial values for runtime dependent parameter
		 */
		// before anything can be computed, a system has to be assembled first -> 1
		RUN_evaluationPurpose_ = 1;

		// we assume that we want to get the output of the hystoperator when we call
		// this function -> 1
		RUN_vectorToReturn_ = 1;

		// we have no solution yet, so we cannot determine a deltaMatrix
		// return initial tensor -> 1
		RUN_tensorToReturn_ = 1;

		// although it has no relevance unless RUN_tensorToReturn_ == 3
		// suppose we use the initial tensor as kick-off addition to the deltaMatrix
		RUN_tensorToAdd_ = 1;

		// use last ts value for first comutation of delta matrix
		RUN_useLastTS_ = true;

		// do not output switching state; that is very costly and should only be done
		// for debugging or special figure computation
		RUN_allowBMP_ = false;

		/*
		 * set temporary values for xml dependent parameter
		 */
		// currently, we have no direct handle of the xml input during the creation
		// of this object; the initialization of the XML dependent parameter is therefore
		// postponed to the first call of StdSolveStep (prior to that we do not need the
		// paramters either way)
		XMLParameterSet_ = false;

		XML_EvaluationDepth_ = 0;
		XML_deltaEvalVersion_ = 0;
		XML_HApproxVersion_ = 0;
		XML_performanceMeasurement_ = 0;
		XML_textOutputLevel_ = 0;

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
		}

		//store subdomain list of elements
		SDList_ = actSDList;
		numElemSD_ = SDList_->GetSize();

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
		 * MAT_initialTensor_ is the tensor to be returned, if
		 * RUN_tensorToReturn_ = 1
		 * and the tensor to be added to deltaMat if
		 * RUN_tensorToAdd_ = 1
		 * (actually it is the small signal tensor from the mat file)
		 */
		matCoef->GetTensor(MAT_initialTensor_, lpm);

		// to calculate differential material properties, we need to know e0 / nu0
		if (material_->GetMaterialDatabaseName() == "Electrostatic") {
			rev_mat_fac_ = 8.854187817e-12; //eps0
			PDEName_ = "Electrostatic";
		} else if (material_->GetMaterialDatabaseName() == "Electromagnetics") {
			rev_mat_fac_ = 795774.7155; //nu0
			PDEName_ = "Electromagnetics";
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

	}

	CoefFunctionHyst::~CoefFunctionHyst() {

		delete timer_;

		delete[] E_B_;
		delete[] P_M_;
		delete[] E_H_;

		delete[] E_B_lastIt_;
		delete[] P_M_lastIt_;
		delete[] E_H_lastIt_;

		delete[] E_B_lastTS_;
		delete[] P_M_lastTS_;
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

			material_->GetScalar(MAT_rotRes_, ROT_RESISTANCE, Global::REAL);
			material_->GetScalar(MAT_angResistance_, ANG_DISTANCE, Global::REAL);
			material_->GetScalar(MAT_angResolution_, ANG_RESOLUTION, Global::REAL);
			material_->GetScalar(MAT_angClipping_, ANG_CLIPPING, Global::REAL);
			material_->GetScalar(MAT_ampResolution_, AMP_RESOLUTION, Global::REAL);

			int printOut;
			int bmpResolution;
			material_->GetScalar(printOut, PRINT_PREISACH);
			material_->GetScalar(bmpResolution, PRINT_PREISACH_RESOLUTION);

			/*
			 * if printOut > 0 -> activate output of overlaid switching and rotation state; output every printOut timestep
			 * bmpResolution -> number of pixels (std = 1000)
			 */
			MAT_printOut_ = (UInt) printOut;
			MAT_bmpResolution_ = (UInt) bmpResolution;

			int eval;
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

		cuttingApplied_ = Vector<UInt>(numStorageEntries_);
		cuttingApplied_.Init(0);

		RUN_residualDecreased_ = true;

		RUN_cuttingAlreadyCounted_ = Vector<UInt>(numStorageEntries_);
		RUN_cuttingAlreadyCounted_.Init(0);

		RUN_deltaMatComputedDuringCurrentTS_ = Vector<UInt>(numStorageEntries_);
		RUN_deltaMatComputedDuringCurrentTS_.Init(0);

		RUN_deltaMatComputedDuringCurrentIt_ = Vector<UInt>(numStorageEntries_);
		RUN_deltaMatComputedDuringCurrentIt_.Init(0);

		E_B_ = new Vector<Double>[numStorageEntries_];
		P_M_ = new Vector<Double>[numStorageEntries_];
		E_H_ = new Vector<Double>[numStorageEntries_];

		E_B_lastIt_ = new Vector<Double>[numStorageEntries_];
		P_M_lastIt_ = new Vector<Double>[numStorageEntries_];
		E_H_lastIt_ = new Vector<Double>[numStorageEntries_];

		E_B_lastTS_ = new Vector<Double>[numStorageEntries_];
		P_M_lastTS_ = new Vector<Double>[numStorageEntries_];
		E_H_lastTS_ = new Vector<Double>[numStorageEntries_];

		deltaMat_ = new Matrix<Double>[numStorageEntries_];
		deltaMat_lastIt_ = new Matrix<Double>[numStorageEntries_];
		deltaMat_lastTS_ = new Matrix<Double>[numStorageEntries_];
		deltaMat_estimated_ = new Matrix<Double>[numStorageEntries_];

		requiresReeval_ = new bool[numStorageEntries_];


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
			 *    i.e. the arrays for the current values (E_H_, E_B_ and P_M_) store
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

			// P_M_ stores P for electrostatics and M for magnetics
			// -> both are the output of the hyst operator, so simply initialize those
			//    arrays with MAT_initialOutput_[k]
			P_M_[k] = MAT_initialOutput_[k];
			//std::cout << "k: " << k << std::endl;
			//std::cout << "P_M_[k]: " << P_M_[k].ToString() << std::endl;
			P_M_lastIt_[k] = zeroVec;
			P_M_lastTS_[k] = zeroVec;

			deltaMat_lastIt_[k] = MAT_initialTensor_;
			deltaMat_lastTS_[k] = MAT_initialTensor_;

			deltaMat_estimated_[k] = Matrix<Double>(dim_, dim_);

			if (MAT_initialInput_.NormL2() == 0) {
				// no initial input > use MAT_initialTensor_
				deltaMat_[k] = MAT_initialTensor_;
			} else {
				// create a first deltaMatrix
				std::string evalMethodName;
				UInt lastDigit = XML_deltaEvalVersion_ % 10;

				if (lastDigit == 0) {
					evalMethodName = "DirectDivisionAbsNew";
				} else if (lastDigit == 1) {
					evalMethodName = "DirectDivisionAbsSatStepDownNew";
				} else if (lastDigit == 2) {
					evalMethodName = "DirectDivisionAbsCut";
				} else if (lastDigit == 3) {
					evalMethodName = "GeometricMean";
				} else if (lastDigit == 4) {
					evalMethodName = "VeryNewApproach";
				} else {
					EXCEPTION("Evaluation approach not implemented yet");
				}

				bool outofSat = false;
				bool intoSat = false;
				bool satToSat = false;
				if (abs(MAT_initialOutput_[k].NormL2() - MAT_ySat_) >= tol_) {
					// hyst operator is in saturation
					intoSat = true;
				}
				CreateDeltaMatrix(MAT_initialInput_, MAT_initialOutput_[k], deltaMat_[k], evalMethodName, k, intoSat, outofSat, satToSat, MAT_initialInput_);
			}
		}

    if(performInversionTest_){
      std::cout << "Test inversion of vector hyst operator" << std::endl;
      TestInversion(rev_mat_fac_);
    }
		//EstimateCurrentSlope(slope,1e-6);

	}

	void CoefFunctionHyst::PreprocessLPM(const LocPointMapped& lpmInput,
		LocPointMapped& lpmOutput, UInt& operatorIdx, UInt& storageIdx) {
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

		/*
		 * standard evaluation: midpointOnly ALWAYS used
		 * extended evaluation: midpointOnly = true for output computation and
		 *                                      during SetPreviousHystValues
		 * full evaluation: midpointOnly = true only for output computation
		 *
		 * > use function EvaluateAtMidpointOnly() to determine the flag
		 */
		bool forceMidpoint = EvaluateAtMidpointOnly();

		if (forceMidpoint) {
			if (XML_textOutputLevel_ == 2) {
				std::cout << "Evaluate only at midpoint" << std::endl;
			}
			/*
			 * get solution at midpoint of element (elemSolution)
			 */
			LocPoint lp = Elem::shapes[el->type].midPointCoord;
			shared_ptr<ElemShapeMap> esm = lpmInput.shapeMap;

			lpmOutput.Set(lp, esm, 0.0);

		} else {
			if (XML_textOutputLevel_ == 2) {
				std::cout << "Evaluate at integration point" << std::endl;
			}
			/*
			 * get solution at actual integration point
			 */
			lpmOutput = lpmInput;
		}

		UInt idxElem = globalElem2Local_[el->elemNum];
		UInt idxPoint = locPointIndices_[lpmOutput.lp.number];

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
	}

	Vector<Double> CoefFunctionHyst::GetOutputOfHysteresisOperator(const LocPointMapped& Originallpm, int timeLevel, bool invert) {

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

		if (XML_textOutputLevel_ == 2) {
		}

		assert(XMLParameterSet_);
		totalCallingCounter_++;

		/*
		 * 1. preprocess lpm to get information about the actual lpm where we have to
		 *  evaluate the hyst operator, the index of the operator and the index of the storage array
		 */
		UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;

		PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);

		/*std::cout << "timeLevel " << timeLevel << std::endl;
		std::cout << "storageIdx: " << storageIdx << std::endl;
		std::cout << "blub" << std::endl;
		std::cout << "size of arrays: " << P_M_lastTS_->GetSize() << std::endl;
		 */
		/*
		 * 2. check time level and need for evaluation
		 */
		if (timeLevel == -1) {
			// return value from last time step
			return P_M_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
			// return value from last iteration
			return P_M_lastIt_[storageIdx];
		} else {
			// get current state; here we have to check if a reevaluation is needed
			// note: requiresReeval_ has one entry for each storage, not for each operator
			if (false == requiresReeval_[storageIdx]) {
//				std::cout << "no reeval needed as flag is false!" << std::endl;
				// return current value
				return P_M_[storageIdx];
			}
		}

		/*
		 * 3. Evaluate hysteresis operator
		 */
		// first get input
		Vector<Double> input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx, invert);

		// then compute output
		return CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx);

	}

	Vector<Double> CoefFunctionHyst::RetrieveLPMSolution(LocPointMapped& actualLPM, UInt storageIdx, int timeLevel) {

		if (timeLevel == -1) {
			// get solution from last ts
			return E_B_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
			// get solution from last iteration
			return E_B_lastIt_[storageIdx];
		} else {
			// get current input from system
			Vector<Double> LPMSolution = Vector<Double>(dependCoef1_->GetVecSize());
			dependCoef1_->GetVector(LPMSolution, actualLPM);
			return LPMSolution;
		}
	}

	Vector<Double> CoefFunctionHyst::RetrieveInputToHysteresisOperator(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool invert) {

		// get current solution first
		// as we want the current input only, set timeLevel to 0
		UInt timeLevel = 0;
		Vector<Double> curLPMSolution = RetrieveLPMSolution(actualLPM, storageIdx, timeLevel);

		if (invert) {
			// check if solution did change since last time > if not, hope that input / output pair of hyst
			// operator can be reused
			Vector<Double> diff = curLPMSolution;
			Vector<Double> toDiff = E_B_[storageIdx];
			diff -= toDiff;

			if (diff.NormL2() < MAT_ampResolution_) {
				if (XML_textOutputLevel_ == 2) {
					std::cout << "(Nearly) no difference in amplitude to previous element solution" << std::endl;
				}
				// reuse previous input
				return E_H_[storageIdx];
			} else {
				EXCEPTION("Implemenent me");
			}

			// Idea: use Levenberg Marquart to estimate  inversion of hyst operator//
		} else {
			// for electrostatic / piezoelectric material, we can use the element solution (E) as input
			return curLPMSolution;
		}
	}

	Vector<Double> CoefFunctionHyst::CalcOutputOfHysteresisOperator(Vector<Double> inputToHystOperator, UInt operatorIdx, UInt storageIdx) {

		if (XML_textOutputLevel_ == 2) {
		}

		Vector<Double> outputOfHystOperator;

		/*
		 *  overwrite the hysteresis memory
		 *          or
		 *  work on temporal storage
		 */
		bool overwriteMemory = OverwriteHystMemory();

		/*
		 * Check if a reevaluation really is required by comparing the actual input to
		 * the last used input
		 */
		Vector<Double> diff = inputToHystOperator;
		Vector<Double> toDiff = E_H_[storageIdx];
		diff -= toDiff;

		if ((diff.NormL2() < MAT_ampResolution_)&&(overwriteMemory == false)) {
			if (XML_textOutputLevel_ == 2) {
				std::cout << "(Nearly) no difference in amplitude to previous input" << std::endl;
				std::cout << "Input: " << inputToHystOperator.ToString() << std::endl;
				std::cout << "Old: " << toDiff.ToString() << std::endl;
			}
			return P_M_[storageIdx];
		}

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

			outputOfHystOperator = hyst_->computeValue_vec(inputToHystOperator, operatorIdx, overwriteMemory, RUN_overwriteDirection_);

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
					std::cout << "Outputting bmp" << std::endl;
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

		E_H_[storageIdx] = inputToHystOperator;
		E_B_[storageIdx] = inputToHystOperator;
		P_M_[storageIdx] = outputOfHystOperator;

		// this flag has to be reset by hand!
		requiresReeval_[storageIdx] = false;

		if (XML_textOutputLevel_ == 2) {

		}

		return outputOfHystOperator;

	}

	void CoefFunctionHyst::ExtractSolutionAndInputForHystOperator(Vector<Double>& extractedSolution, Vector<Double>& extractedInput, UInt& operatorIndex, UInt& storageIndex,
		const LocPointMapped& lpm, bool midpointOnly) {

		if (XML_textOutputLevel_ == 2) {
			std::cout << "++ ExtractSolutionAndInputForHystOperatorNew ++" << std::endl;
		}

		/*
		 * This function will extract the current solution from the FE-system.
		 * If midpointOnly is true, the solution will be extracted at the midpoint
		 * of the element corresponding to lpm
		 *  -> elementSolution
		 * If midpointOnly is false, the solution will be extracted at lpm directly
		 *  -> solution at integration point
		 *
		 * The gathered solution (E for electrostatics, B for magnetics)
		 *  will be stored to vector E_B_
		 *
		 * Additionally, this function will prepare the input for the hysteresis operator
		 * For magnetics, H will be approximated from B and will be stored to E_H_
		 * For electrostatics, E will additionally be stored to E_H_
		 *
		 * Return values:
		 *  operatorIndex = index of hyst operator that has to be evaluated for this integration point
		 *  storageIndex = index in storage vectors where input, output, deltaMat has to be stored
		 *
		 *  extractedSolution = E in case of electrostatics, B in case of magnetics
		 *  extractedInput = E in case of electrostatics, H in case of magnetics
		 *
		 */

		/*
		 * First part - get actual solution
		 * Electrostatics: elecFieldIntensity
		 * Magnetics: magFluxDensity (not magFieldIntensity)
		 */
		Vector<Double> sol = Vector<Double>(dependCoef1_->GetVecSize());
		const Elem * el = lpm.ptEl;

		/*
		 * standard evaluation: midpointOnly ALWAYS used
		 * extended evaluation: midpointOnly = true for output computation and
		 *                                      during SetPreviousHystValues
		 * full evaluation: midpointOnly = true only for output computation
		 */
		if (midpointOnly) {
			if (XML_textOutputLevel_ == 2) {
				std::cout << "Evaluate only at midpoint" << std::endl;
			}
			/*
			 * get solution at midpoint of element (elemSolution)
			 */
			LocPoint lp = Elem::shapes[el->type].midPointCoord;
			LocPointMapped copylpm;
			shared_ptr<ElemShapeMap> esm = lpm.shapeMap;

			copylpm.Set(lp, esm, 0.0);
			dependCoef1_->GetVector(sol, copylpm);

		} else {
			if (XML_textOutputLevel_ == 2) {
				std::cout << "Evaluate at integration point" << std::endl;
			}
			/*
			 * get solution at actual integration point
			 */
			dependCoef1_->GetVector(sol, lpm);
		}
		extractedSolution = sol;

		UInt idxElem = globalElem2Local_[el->elemNum];
		UInt idxPoint = locPointIndices_[lpm.lp.number];

		// storage index = index where to store results, input, ...
		// operator index = index of hyst operator to evaluate
		if (XML_EvaluationDepth_ == 1) {
			// standard evaluation
			// storageIndex = operatorIndex = elementIndex
			storageIndex = idxElem;
			operatorIndex = idxElem;
		} else if (XML_EvaluationDepth_ == 2) {
			// extended evaluation
			// operatorIndex = elementIndex (only one operator per element)
			// storageIndex = combinedIndex = elementIndex*(numIntegrationPoints+1)+pointIndex
			operatorIndex = idxElem;
			storageIndex = idxElem * (numIntegrationPoints_ + 1) + idxPoint;
		} else {
			// full evaluation
			// one operator and one storage per integration point
			// operatorIndex = storageIndex = combinedIndex
			operatorIndex = idxElem * (numIntegrationPoints_ + 1) + idxPoint;
			storageIndex = operatorIndex;
		}

		//
		//  std::cout << "LP: " << lpm.lp.number << ": " << lpm.lp.coord.ToString() << std::endl;
		//  std::cout << "Retrieved index: " << std::endl;
		//  std::cout << "Elem index: " << idx << std::endl;
		//  std::cout << "Point index original and after mapping: " << lpm.lp.number << " / " << locPointIndices_[lpm.lp.number] << std::endl;
		//  std::cout << "Combined index (elemIndex*(numIntegrationPoints+1)+pointIndex): " << idx*(numIntegrationPoints_+1)+locPointIndices_[lpm.lp.number] << std::endl;


		/*
		 * Second part - get input for hysteresis operator
		 * Electrostatics: elecFieldIntensity
		 * Magnetics: magFluxDensity (not magFieldIntensity)
		 */
		Vector<Double> input;

		if (PDEName_ == "Electromagnetics") {
			/*
			 * two approaches to get H from B
			 * 1. H_new = nu0*B_new - M_old
			 *
			 * 2. H_new = deltaH_old/deltaB_old * deltaB_new + H_old
			 *
			 *      with
			 *        deltaH_old/deltaB_old = deltaMatrix from last iteration
			 *        deltaB_new = B_new - B_old
			 *
			 */
			if ((XML_HApproxVersion_ == 1) || (XML_HApproxVersion_ == 3)) {
				input = sol*rev_mat_fac_; //H = nu0*B
				if (RUN_useLastTS_) {
					input -= P_M_lastTS_[storageIndex]; //subtract old M
				} else {
					input -= P_M_lastIt_[storageIndex]; //subtract old M
				}
			} else {
				Vector<Double> deltaB = sol;
				if (RUN_useLastTS_) {
					deltaB -= E_B_lastTS_[storageIndex];
					input = deltaMat_lastTS_[storageIndex] * deltaB; //deltaH_new = deltaH_old/deltaB_old * deltaB_new
					input += E_H_lastTS_[storageIndex]; //H = deltaH_new + H_old
				} else {
					deltaB -= E_B_lastIt_[storageIndex];
					input = deltaMat_lastIt_[storageIndex] * deltaB; //H = deltaH_old/deltaB_old * deltaB_new
					input += E_H_lastIt_[storageIndex]; //H = deltaH_new + H_old
				}
			}
			extractedInput = input;

		} else {
			/*
			 * electrostatics: return elemSolution directly
			 */
			extractedInput = sol;
		}

		//E_H_[storageIndex] = extractedInput;
		E_B_[storageIndex] = extractedSolution;

		if (XML_textOutputLevel_ == 2) {
			std::cout << "extractedSolution: " << extractedSolution.ToString() << std::endl;
			std::cout << "extractedInput: " << extractedInput.ToString() << std::endl;
			std::cout << "Index for storage: " << storageIndex << std::endl;
			std::cout << "Index of operator: " << operatorIndex << std::endl;
		}

	}

	UInt CoefFunctionHyst::EvaluateHysteresisOperator(Vector<Double> inputToHystOperator,
		Vector<Double>& outputOfHystOperator, UInt operatorIndex, UInt storageIndex,
		bool overwriteMemory, Vector<Double>& currentSolution) {

		if (XML_textOutputLevel_ == 2) {
			std::cout << "++ EvaluateHysteresisOperator ++" << std::endl;
			std::cout << "Index for storage: " << storageIndex << std::endl;
			std::cout << "Index of operator: " << operatorIndex << std::endl;
			//std::cout << "Overwrite memroy? " << overwriteMemory << std::endl;
			//std::cout << "StoreVectors? " << storeVectors << std::endl;
		}

		assert(XMLParameterSet_);
		totalCallingCounter_++;

		/*
		 * Check if the current input to the hyst operator with index operatorIndex
		 * is the same as the one that was stored in vector E_H_ at index storageIndex
		 *
		 * for standard evaluation:
		 *  a) we always evaluate at the midpoint of an element
		 *  b) we only store entries for lastTS, lastIt, currentIt at the midpoint
		 *  c) operatorIndex and storageIndex are equal to the index of the element
		 *
		 * for extended evaluation:
		 *  a) we evaluate the same hystoperator at different points (with locked
		 *      memory of course)
		 *  b) we store entries for each integration / evaluation point individually
		 *  c) operatorIndex = index of element; storageIndex = combined element + point index
		 *
		 * for full evaluation:
		 *  a) we evaluate different hystoperators for each point
		 *  b) we store entries for each point
		 *  c) operatorIndex = storageIndex = combined element + point index
		 */
		Vector<Double> diff = inputToHystOperator;
		Vector<Double> toDiff = E_H_[storageIndex];
		diff -= toDiff;

		//std::cout << "diff: " << diff.ToString() << std::endl;

		if ((diff.NormL2() < MAT_ampResolution_)&&(overwriteMemory == false)) {
			if (XML_textOutputLevel_ == 2) {
				std::cout << "(Nearly) no difference in amplitude to previous input" << std::endl;
				std::cout << "Input: " << inputToHystOperator.ToString() << std::endl;
				std::cout << "Old: " << toDiff.ToString() << std::endl;
			}
			outputOfHystOperator = P_M_[storageIndex];
			//std::cout << "outputOfHystOperator: " << outputOfHystOperator.ToString() << std::endl;
			return 1;
		}
		// Note: we do not have to compare to the values of the lastIteration or lastTS
		// Reason: at the beginning of each iteration, E_H_ and E_H_lastIt_ should contain
		//         the same values; at the beginning of each timestep, E_H_ and E_H_lastTS_
		//         should contain the same values; both are the consequence of using
		//         this function in SetPreviousHystValues; by doing so, the evaluated value
		//         input will be written to E_H_ first (done in this function) and directly
		//         afterwards to E_H_lastIt (or E_H_lastTS) in the SetPreviousHystValues fnc.

		totalEvaluationCounter_++;

		if (MAT_methodType_ == SCALAR) {

			outputOfHystOperator = Vector<Double>(dim_);
			outputOfHystOperator.Init();

			if (XML_performanceMeasurement_) {
				timer_->Start();
			}

			outputOfHystOperator[MAT_dirP_] = hyst_->computeValueAndUpdate(inputToHystOperator[MAT_dirP_], operatorIndex, overwriteMemory);

			if (XML_performanceMeasurement_) {
				timer_->Stop();
				//        std::cout << "Scalar Preisach opeator" << std::endl;
				//        std::cout << "Total time (" << totalEvaluationCounter_ << "calls): " << timer_->GetCPUTime() << std::endl;
				//        std::cout << "Avg Time: " << timer_->GetCPUTime()/totalEvaluationCounter_ << std::endl;
				totalEvaluationTime_ = timer_->GetCPUTime();
				avgEvaluationTime_ = timer_->GetCPUTime() / totalEvaluationCounter_;
			}

		} else {
			if (XML_performanceMeasurement_) {
				timer_->Start();
			}

			/*
			 * currenty not used; idea behind this flag was to lock the direction of the
			 * rotational operator during the evaluation of the VectorPreisach model so
			 * that only the switching state is updated;
			 * did not work so well but maybe it will have some usage in the future
			 */
			RUN_overwriteDirection_ = true;

			if (XML_HApproxVersion_ == 3) {
				std::cout << "use iterative refinement" << std::endl;
				std::cout << "initialInput " << std::endl;
				std::cout << inputToHystOperator.ToString() << std::endl;
				UInt numRefinementSteps = 6;
				Double tol = 1e-3;
				Vector<Double> diff;
				Vector<Double> tmpInput = inputToHystOperator;
				for (UInt i = 0; i < numRefinementSteps; i++) {

					std::cout << "input H_" << i << std::endl;
					std::cout << tmpInput.ToString() << std::endl;

					diff = tmpInput;

					// compute current M
					// do not overwrite memory here
					outputOfHystOperator = hyst_->computeValue_vec(tmpInput, operatorIndex, false, RUN_overwriteDirection_);

					std::cout << "hyst output M_" << i << std::endl;
					std::cout << outputOfHystOperator.ToString() << std::endl;

					// check if H changes due to M
					tmpInput = currentSolution*rev_mat_fac_; // nu0*B
					tmpInput -= outputOfHystOperator; // -M

					std::cout << "H_" << i + 1 << " = nu0*B - M_" << i << " = " << std::endl;
					std::cout << tmpInput.ToString() << std::endl;

					diff -= tmpInput;

					Double relDiff = diff.NormL2() / tmpInput.NormL2();
					std::cout << "|H_" << i << " - H_" << i + 1 << "| / |H_" << i + 1 << "| = " << relDiff << std::endl;

					if (relDiff < tol) {
						break;
					}
				}

				std::cout << "Initial input: " << std::endl;
				std::cout << inputToHystOperator.ToString() << std::endl;

				inputToHystOperator = tmpInput;

				std::cout << "Refined input: " << std::endl;
				std::cout << inputToHystOperator.ToString() << std::endl;

				if (overwriteMemory) {
					// if we actually want to set the hystoperator itself, we have to do this only with the actual refined input
					outputOfHystOperator = hyst_->computeValue_vec(inputToHystOperator, operatorIndex, overwriteMemory, RUN_overwriteDirection_);
				}

			} else {
				outputOfHystOperator = hyst_->computeValue_vec(inputToHystOperator, operatorIndex, overwriteMemory, RUN_overwriteDirection_);
			}



			if (XML_performanceMeasurement_) {
				timer_->Stop();
				//        std::cout << "Vector Preisach opeator" << std::endl;
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
				if ((cnt % MAT_printOut_ == 0)&&(operatorIndex == firstIdx)) {
					std::cout << "Outputting bmp" << std::endl;
					std::stringstream filenamebuf;
					filenamebuf << "Switch_Elem" << firstIdx << "_Step" << std::setfill('0') << std::setw(5) << cnt << "_v" << MAT_vecPreisachImplementationVersion_ << "_numRows" << MAT_numRows_ << ".bmp";
					hyst_->switchingStateToBmp(MAT_bmpResolution_, filenamebuf.str(), operatorIndex, true);
				}

				if (operatorIndex == firstIdx) {
					cnt++;
					/*
					 * disable output until reset
					 * -> otherwise we would get two images for each timestep if P and D are computed
					 */
					RUN_allowBMP_ = false;
				}
			}
		}

		E_H_[storageIndex] = inputToHystOperator;
		P_M_[storageIndex] = outputOfHystOperator;

		//if(XML_textOutputLevel_ == 2){
		// std::cout << "Evaluated output: " << outputOfHystOperator << std::endl;
		//}

		return 0;
	}

	void CoefFunctionHyst::EstimateCurrentSlope(Vector<Double> steppingDirection, Double scaling) {

		//std::cout << "Estimate actual slope around current point" << std::endl;
		//std::cout << "scaling: " << scaling << std::endl;

		// Problem: the first iteration of the non-linear solve step is crucial; if the slope is too
		//          too small or too large, the solution update might be too small or large that the
		//          following computation of the deltaMatrix will result in values which lead to even
		//          worse results or even oscillating solutions;
		// Idea: to get a better start into the solution process, estimate the slope around the current
		//       point X,Y by evaluating Y(X-scaling*steppingDirection) and Y(X+scaling*steppingDirection)
		//       then create a deltaMatrix using these to Y and 2*scaling*steppingDirection as nominator and
		//       denominator, respectively

		/*
		 * backup current evaluationPurpose
		 */
		UInt evalPurposeBackup = RUN_evaluationPurpose_;

		RUN_evaluationPurpose_ = 1; // assemble, hyst operator locked, evaluation at midpoint only or at each int. point

		bool overwriteMemory = OverwriteHystMemory();
		bool midpointOnly = EvaluateAtMidpointOnly();

		/*
		 * 1. iterate over all elements
		 */
		EntityIterator it = SDList_->GetIterator();
		LocPoint lp;
		LocPointMapped lpm;
		const Elem * el;
		shared_ptr<ElemShapeMap> esm;
		Vector<Double> zeroVec = Vector<Double>(dim_);
		zeroVec.Init();
		UInt idxElem;
		UInt idxPoint;
		UInt operatorIndex;
		UInt storageIndex;

		Vector<Double> Y_up = Vector<Double>(dim_);
		Vector<Double> Y_down = Vector<Double>(dim_);
		Vector<Double> X_up = Vector<Double>(dim_);
		Vector<Double> X_down = Vector<Double>(dim_);
		Vector<Double> dY;
		Vector<Double> dX = Vector<Double>(dim_);
		Vector<Double> steppingDirectionUsed;

		for (it.Begin(); !it.IsEnd(); it++) {

			/*
			 * 1.0 get element and shape map
			 */
			el = it.GetElem();
			esm = it.GetGrid()->GetElemShapeMap(el, true);
			idxElem = globalElem2Local_[el->elemNum];

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

					idxPoint = locPointIndices_[lpm.lp.number];

					if (XML_EvaluationDepth_ == 2) {
						// extended evaluation
						// operatorIndex = elementIndex (only one operator per element)
						// storageIndex = combinedIndex = elementIndex*(numIntegrationPoints+1)+pointIndex
						operatorIndex = idxElem;
						storageIndex = idxElem * (numIntegrationPoints_ + 1) + idxPoint;
					} else if (XML_EvaluationDepth_ == 3) {
						// full evaluation
						// one operator and one storage per integration point
						// operatorIndex = storageIndex = combinedIndex
						operatorIndex = idxElem * (numIntegrationPoints_ + 1) + idxPoint;
						storageIndex = operatorIndex;
					} else {
						EXCEPTION("EvaluationDepth_ == 1 should be midpoint only!");
					}

					// extract current input and solution
					Vector<Double> extractedSolution;
					Vector<Double> extractedInput;
					ExtractSolutionAndInputForHystOperator(extractedSolution, extractedInput, operatorIndex, storageIndex, lpm, midpointOnly);

					steppingDirectionUsed = steppingDirection;

					// if steppingdirection is zero, take vector pointing in all directions
					if (steppingDirectionUsed.NormL2() == 0) {
						steppingDirectionUsed = Vector<Double>(dim_);
						steppingDirectionUsed.Init(1.0 / std::sqrt(dim_));
					}

					// alternate approach: use direction of last hyst ouput > may not contain all direction!
					//        if(steppingDirectionUsed.NormL2() == 0){
					//          // if no direction is specified, take direction of input
					//          steppingDirectionUsed = P_M_[storageIndex];
					//          if(steppingDirectionUsed.NormL2() == 0){
					//            // if it is still zero, take vector pointing in all directions
					//            steppingDirectionUsed = Vector<Double>(dim_);
					//            steppingDirectionUsed.Init(1.0/sqrt(dim_));
					//          }else{
					//            steppingDirectionUsed /= P_M_[storageIndex].NormL2();
					//          }
					//        }

					X_up = extractedInput;
					X_up.Add(scaling, steppingDirectionUsed);

					X_down = extractedInput;
					X_down.Add(-scaling, steppingDirectionUsed);

					// evaluate at corresponding operator index
					EvaluateHysteresisOperator(X_up, Y_up, operatorIndex, storageIndex, overwriteMemory, extractedSolution);

					EvaluateHysteresisOperator(X_down, Y_down, operatorIndex, storageIndex, overwriteMemory, extractedSolution);

					dX.Init();
					dX.Add(2 * scaling, steppingDirectionUsed);
					dY = Y_up;
					dY -= Y_down;

					CreateDeltaMatrix(dX, dY, deltaMat_estimated_[storageIndex], "Classic", storageIndex, false, false, false, steppingDirectionUsed);
				}

			} else {
				// only standard evaluation should end up here;
				// unlike function SetPreviousHystValues, extendend and full evaluation do not need the midpoint
				// Reason: we do not assemble at the midpoint (except for standard evaluation, where we only compute at
				//         midpoint and use this tensor for all integration points

				lp = Elem::shapes[el->type].midPointCoord;
				lpm.Set(lp, esm, 0.0);

				idxPoint = locPointIndices_[lpm.lp.number];
				if (XML_EvaluationDepth_ == 1) {
					// standard evaluation
					// storageIndex = operatorIndex = elementIndex
					storageIndex = idxElem;
					operatorIndex = idxElem;
				} else {
					EXCEPTION("EvaluationDepth_ == 2 and 3 should not need the midpoint");
				}

				// extract current input and solution
				Vector<Double> extractedSolution;
				Vector<Double> extractedInput;
				ExtractSolutionAndInputForHystOperator(extractedSolution, extractedInput, operatorIndex, storageIndex, lpm, midpointOnly);

				steppingDirectionUsed = steppingDirection;

				// if steppingdirection is zero, take vector pointing in all directions
				if (steppingDirectionUsed.NormL2() == 0) {
					steppingDirectionUsed = Vector<Double>(dim_);
					steppingDirectionUsed.Init(1.0 / std::sqrt(dim_));
				}

				// alternate approach: use direction of last hyst ouput > may not contain all direction!
				//        if(steppingDirectionUsed.NormL2() == 0){
				//          // if no direction is specified, take direction of input
				//          steppingDirectionUsed = P_M_[storageIndex];
				//          if(steppingDirectionUsed.NormL2() == 0){
				//            // if it is still zero, take vector pointing in all directions
				//            steppingDirectionUsed = Vector<Double>(dim_);
				//            steppingDirectionUsed.Init(1.0/sqrt(dim_));
				//          }else{
				//            steppingDirectionUsed /= P_M_[storageIndex].NormL2();
				//          }
				//        }

				X_up = extractedInput;
				X_up.Add(scaling, steppingDirectionUsed);

				X_down = extractedInput;
				X_down.Add(-scaling, steppingDirectionUsed);

				// evaluate at corresponding operator index
				EvaluateHysteresisOperator(X_up, Y_up, operatorIndex, storageIndex, overwriteMemory, extractedSolution);

				EvaluateHysteresisOperator(X_down, Y_down, operatorIndex, storageIndex, overwriteMemory, extractedSolution);

				dX.Init();
				dX.Add(2 * scaling, steppingDirectionUsed);
				dY = Y_up;
				dY -= Y_down;

				CreateDeltaMatrix(dX, dY, deltaMat_estimated_[storageIndex], "Classic", storageIndex, false, false, false, steppingDirectionUsed);

			}
		}

		/*
		 * 2.0 restore old state of evaluationPurpose
		 */
		RUN_evaluationPurpose_ = evalPurposeBackup;

	}

	void CoefFunctionHyst::SetPreviousHystVals(bool setLastTS, bool forceMemoryLock) {

		//if(XML_textOutputLevel_ == 2){
		std::cout << "++ SetPreviousHystVals ++" << std::endl;
		//std::cout << "lastTS? " << setLastTS << std::endl;
		// }

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
			// we only overwrite memorey for the midpoint (where the hystoperator is located)
			RUN_evaluationPurpose_ = 3;
		} else {
			// overwriteMemory = false
			RUN_evaluationPurpose_ = 2;
		}

		bool overwriteMemoryIntPoints = false;
		if (XML_EvaluationDepth_ == 3) {
			overwriteMemoryIntPoints = true;
		}
		bool overwriteMemoryMidPoint = OverwriteHystMemory();
		bool midpointOnly = EvaluateAtMidpointOnly();

		if (forceMemoryLock) {
			overwriteMemoryIntPoints = false;
			overwriteMemoryMidPoint = false;
		}

		//std::cout << "EvaluationDepth: " << XML_EvaluationDepth_ << std::endl;
		//std::cout << "RUN_evaluationPurpose_: " << RUN_evaluationPurpose_ << std::endl;
		//std::cout << "midpointOnly: " << midpointOnly << std::endl;
		//std::cout << "overwriteMemoryMidPoint: " << overwriteMemoryMidPoint << std::endl;

		/*
		 * 1. iterate over all elements
		 */
		EntityIterator it = SDList_->GetIterator();
		LocPoint lp;
		LocPointMapped lpm;
		const Elem * el;
		shared_ptr<ElemShapeMap> esm;
		Vector<Double> zeroVec = Vector<Double>(dim_);
		zeroVec.Init();

		for (it.Begin(); !it.IsEnd(); it++) {

			/*
			 * 1.0 get element and shape map
			 */
			el = it.GetElem();
			esm = it.GetGrid()->GetElemShapeMap(el, true);

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

					// extract solution and input for hystoperator
					ExtractSolutionAndInputForHystOperator(solution, input, operatorIdx, storageIdx, lpm, false);

					// evaluate at corresponding operator index
					EvaluateHysteresisOperator(input, output, operatorIdx, storageIdx, overwriteMemoryIntPoints, solution);

					// store at storage index
					if (setLastTS) {
						E_B_lastTS_[storageIdx] = solution;
						P_M_lastTS_[storageIdx] = output;
						E_H_lastTS_[storageIdx] = input;
						deltaMat_lastTS_[storageIdx] = deltaMat_[storageIdx];
						dY_sol[storageIdx] = zeroVec;
						X_low[storageIdx] = zeroVec;
						X_up[storageIdx] = zeroVec;
					} else {
						E_B_lastIt_[storageIdx] = solution;
						P_M_lastIt_[storageIdx] = output;
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

			// extract solution and input for hystoperator
			ExtractSolutionAndInputForHystOperator(solution, input, operatorIdx, storageIdx, lpm, true);

			// evaluate at corresponding operator index
			EvaluateHysteresisOperator(input, output, operatorIdx, storageIdx, overwriteMemoryMidPoint, solution);

			// store at storage index
			if (setLastTS) {
				E_B_lastTS_[storageIdx] = solution;
				P_M_lastTS_[storageIdx] = output;
				E_H_lastTS_[storageIdx] = input;
				deltaMat_lastTS_[storageIdx] = deltaMat_[storageIdx];
				dY_sol[storageIdx] = zeroVec;
				X_low[storageIdx] = zeroVec;
				X_up[storageIdx] = zeroVec;
			} else {
				E_B_lastIt_[storageIdx] = solution;
				P_M_lastIt_[storageIdx] = output;
				E_H_lastIt_[storageIdx] = input;
				deltaMat_lastIt_[storageIdx] = deltaMat_[storageIdx];
			}
		}

		if (setLastTS) {
			//EstimateCurrentSlope(Vector<Double>(dim_), 1e-4);
			//EstimateCurrentSlope(Vector<Double>(dim_), 1e-8);
			EstimateCurrentSlope(Vector<Double>(dim_), 1e-10);
		}


		/*
		 * 2.0 restore old state of evaluationPurpose
		 */
		RUN_evaluationPurpose_ = evalPurposeBackup;

	}

	void CoefFunctionHyst::GetScalar(Double& outputScalar, const LocPointMapped& lpm) {
		/*
		 * is this function needed at all?
		 * it shouldn't as we only need the Vector for RHS and Output
		 * and the Tensor for deltaMatrix computation
		 */
		EXCEPTION("GetScalar not implemented for coefFncHyst");
	}

	void CoefFunctionHyst::GetVector(Vector<Double>& outputVector, const LocPointMapped& lpm) {

		bool newImpl = true;
		if (newImpl) {
			bool invert = false;
			UInt timeLevel = 0;
			outputVector = GetOutputOfHysteresisOperator(lpm, timeLevel, invert);
		} else {




			if (XML_textOutputLevel_ == 2) {
				std::cout << "++ GetVector ++" << std::endl;
				//std::cout << "evaluationPurpose: " << RUN_evaluationPurpose_ << std::endl;
				//std::cout << "vectorToReturn: " << RUN_vectorToReturn_ << std::endl;
			}

			if (RUN_evaluationPurpose_ == 4) {
				/*
				 * for output computation, we should actually return the state of the hyst operator
				 */
				RUN_vectorToReturn_ = 1;
			}

			if (RUN_vectorToReturn_ == 2) {
				// return zero
				outputVector.Resize(dependCoef1_->GetVecSize());
				outputVector.Init();
				return;
			}

			/*
			 * 1. depending on evaluation depth and purpose, we
			 *    a) evaluate the hysteresis operator only at the midpoint
			 *          or
			 *       at the actual lpm
			 *    b) overwrite the hysteresis memory
			 *          or
			 *       work on temporal storage
			 */
			bool midpointOnly = EvaluateAtMidpointOnly();
			bool overwriteMemory = OverwriteHystMemory();

			/*
			 * in the current approach, we only write the hysteresis memory during
			 * SetPreviousHystValues
			 */
			assert(overwriteMemory == false);

			/*
			 *  2. extract solution and input to hysteresis operator
			 */
			Vector<Double> solution;
			Vector<Double> input;
			UInt operatorIdx;
			UInt storageIdx;

			ExtractSolutionAndInputForHystOperator(solution, input, operatorIdx, storageIdx, lpm, midpointOnly);

//			std::cout << "operatorIdx/storageIdx: " << operatorIdx << "/" << storageIdx << std::endl;
//			std::cout << "RUN_vectorToReturn_" << RUN_vectorToReturn_ << std::endl;

			if (RUN_vectorToReturn_ == 4) {
				// return value of lastTS
				// we have to call the extract function to get the storage index
				outputVector = P_M_lastTS_[storageIdx];
				return;
			}

			/*
			 *  3. evaluate hysteresis operator
			 */
			EvaluateHysteresisOperator(input, outputVector, operatorIdx, storageIdx, overwriteMemory, solution);
			//  std::cout << "evaluated values: " << std::endl;
			//  std::cout << "input: " << input.ToString() << std::endl;
			//  std::cout << "outputVector: " << outputVector.ToString() << std::endl;
			//  std::cout << "operatorIdx: " << operatorIdx << std::endl;
			//  std::cout << "storageIdx: " << storageIdx << std::endl;

			if (RUN_vectorToReturn_ == 3) {
				// subtract last state
				//std::cout << "subtract old state" << std::endl;
				outputVector -= P_M_lastTS_[storageIdx];
			}
			//std::cout << "done " << std::endl;
		}
	}

	void CoefFunctionHyst::GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm) {

		if (XML_textOutputLevel_ == 2) {
			std::cout << "-- GetTensor --" << std::endl;
		}

		if (RUN_evaluationPurpose_ == 4) {
			/*
			 * for output computation, we return only eps0/nu0
			 * > reason: D/H is computed during the output via eps0*E+P / nu*B - M
			 */
			RUN_tensorToReturn_ = 2;
		}

		if (RUN_tensorToReturn_ == 1) {
			if (XML_textOutputLevel_ == 2) {
				std::cout << "Return initial tensor" << std::endl;
			}
			outputTensor = MAT_initialTensor_;
			return;
		} else if (RUN_tensorToReturn_ == 2) {
			if (XML_textOutputLevel_ == 2) {
				std::cout << "Return free field tensor" << std::endl;
			}
			outputTensor = MAT_freeFieldTensor_;
			return;
		} else if (RUN_tensorToReturn_ == 4) {
			// special case; use ySat/xSat (for magnetics xSat/ySat)
			Double frac = 0;
			outputTensor = Matrix<Double>(dim_, dim_);
			if (PDEName_ == "Electromagnetics") {
				frac = MAT_xSat_ / MAT_ySat_;
			} else {
				frac = MAT_ySat_ / MAT_xSat_;
			}
			for (UInt i = 0; i < dim_; i++) {
				outputTensor[i][i] = frac;
			}
			return;
		}
		/*
		 * deltaMatrix case
		 */

		/*
		 * determine first if we need the deltaMatrix at the midpoint or
		 * at the actual integration point given by lpm
		 */
		bool midpointOnly = EvaluateAtMidpointOnly();

		/*
		 * check if memory shall get overwritten or if need to work on temporal
		 * storage
		 * In the current version, we ONLY work on the actual storage during
		 * setPreviousHystValues
		 */
		bool overwriteMemory = OverwriteHystMemory();

		assert(overwriteMemory == false);

		/*
		 * get current input for hyst operator
		 */
		Vector<Double> currentSolution;
		Vector<Double> currentInput;
		UInt storageIdx;
		UInt operatorIdx;

		ExtractSolutionAndInputForHystOperator(currentSolution, currentInput, operatorIdx, storageIdx, lpm, midpointOnly);

		//  std::cout << "Returned values: " << std::endl;
		//  std::cout << "currentInput: " << currentInput.ToString() << std::endl;
		//  std::cout << "operatorIdx: " << operatorIdx << std::endl;
		//  std::cout << "storageIdx: " << storageIdx << std::endl;

		if (RUN_tensorToReturn_ == 5) {
			outputTensor = deltaMat_estimated_[storageIdx];
			return;
		}


		/*
		 * TODO: rewrite comment
		 * different implementations selectable via flag deltaMatEvalVersion = xy
		 *
		 *  x > determines what values to use for deltaX and deltaY
		 *  y > determines what to with deltaX and deltaY
		 *
		 *  deltaMatEvalVersion = x0 (10,20,30,...)
		 *
		 *    deltaMat = addTensor +/- abs(deltaY)/abs(deltaX)
		 *
		 *    ->  see "DirectDivisionAbsNew" in CreateDeltaMatrix
		 *
		 *  deltaMatEvalVersion = x1 (11,21,31,...)
		 *
		 *    deltaMat =
		 *          addTensor +/- abs(deltaY)/abs(deltaX) (if deltaY != 0)
		 *
		 *          deltaMat_lastIt_[idx]/10;             (if deltaY == 0 && deltaMat_lastIt_[idx] > 10*rev_mat_fac_)
		 *
		 *    ->  see "DirectDivisionAbsSatStepDownNew" in CreateDeltaMatrix
		 *
		 *  deltaMatEvalVersion = 1y (10,11)
		 *
		 *    deltaX = currentInput - E_H_lastIt_ or
		 *           = currentInput - E_H_lastTS_ (depending on flag lastTS_)
		 *
		 *    deltaY = currentOutput - P_M_lastIt_ or
		 *           = currentOutput - P_M_lastTS_ (depending on flag lastTS_)
		 *
		 *  deltaMatEvalVersion = 2y (20,21)
		 *
		 *    a) to be done before evaluating the hysteresis operator
		 *    deltaX = currentInput - E_H_lastIt_ or
		 *           = currentInput - E_H_lastTS_ (depending on flag lastTS_)
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
		 *  deltaMatEvalVersion = 3< (30,31)
		 *
		 *    a) evaluate hystersis operator at currentInput and at a shiftedInput
		 *      with
		 *        shiftedInput = currentInput + minimalShift * (1,1,1)
		 *    b)
		 *      deltaX = minimalShift * (1,1,1)
		 *      deltaY = shiftedOutput - currentOutput
		 *
		 *    Advantage: no issues width division by 0
		 *               differential deltaMatrix computed
		 *    Disadvantage: hystOperator needs to be evaluated twice
		 *
		 */

		// prepare input vectors for function CreatteDeltaMatrix
		Vector<Double> deltaX = Vector<Double>(dim_);
		Vector<Double> deltaY = Vector<Double>(dim_);
		Vector<Double> currentOutput;

		Vector<Double> diff = currentInput;
		Vector<Double> toDiffX;
		if (RUN_useLastTS_) {
			toDiffX = E_H_lastTS_[storageIdx];
		} else {
			toDiffX = E_H_lastIt_[storageIdx];
		}
		diff -= toDiffX;

		bool outofSat = false;
		bool intoSat = false;
		bool satToSat = false;

		if (diff.NormL2() < -MAT_ampResolution_) {
			std::cout << "diff.NormL2 < MAT_ampResolution" << std::endl;
			/*
			 * the overall distance to the last value is close to zero
			 * > set both deltaX and deltaY to zero and let
			 *    function CreateDeltaMatrix do the rest
			 */
			deltaX.Init();
			deltaY.Init();

		} else {
			/*
			 * for deltaEvalVersion between 20 and 30, we check if single components
			 * of the difference vector are close to 0 (i.e. smaller than MAT_ampResolution_)
			 * in that case, the computation of the deltaMatrix could have troubles (division by values close to 0 is not a good idea)
			 * > try to avoid these troubles by shifting the input of these components to a distance of MAT_ampResolution before
			 *    evaluating the output
			 */
			if ((XML_deltaEvalVersion_ >= 20) && (XML_deltaEvalVersion_ < 30)) {
				for (UInt i = 0; i < dim_; i++) {
					if (abs(diff[i]) < MAT_ampResolution_) {
						if (diff[i] < 0) {
							/*
							 * new value is smaller than old value
							 */
							currentInput[i] = toDiffX[i] - MAT_ampResolution_;
						} else {
							/*
							 * new value is greater or equal old value
							 */
							currentInput[i] = toDiffX[i] + MAT_ampResolution_;
						}
					}
				}
			}

			//    std::cout << "operatorIdx: " << operatorIdx << std::endl;
			//    std::cout << "storageIdx: " << storageIdx << std::endl;
			// use current input (or shifted version) to evaluate hyst operator
			UInt evalResult = EvaluateHysteresisOperator(currentInput, currentOutput, operatorIdx, storageIdx, overwriteMemory, currentSolution);

			// evalResult can be either 0 or 1
			// 0: is only returned if currentInput != E_H_[storageIdx]
			//      > this value has not been used as input yet
			//      > deltaMatrix should be evaluated
			// 1: currentInput == E_H_[storageIdx]
			//      > this can only happen for standard evaluation where for each integration point
			//        the same hysteresis operator and the same storage is used
			//      > here we can basically reuse the current value of deltaMatrix as this should
			//        contain the computed value for this input

			// > does NOT work as expected; reason: before getTensor is called, getVector is used and thus sets E_H_ to the current input
			//    and thus prohibits an evaluation of deltaMat
			if (evalResult == 100) {
        // NOT USED > therefore = 100
				std::cout << "Reuse current tensor" << std::endl;
				outputTensor = deltaMat_[storageIdx];
				return;
			} else {

				// IMPORTANT: for magnetics we have to divide the deltaMat by B (not H!), i.e. use actual solution and not
				//              the retrieved input
				deltaX = currentSolution;
				if (RUN_useLastTS_) {
					toDiffX = E_B_lastTS_[storageIdx];
				} else {
					toDiffX = E_B_lastIt_[storageIdx];
				}
				deltaX -= toDiffX;

				//   std::cout << "E_B_lastIt_[storageIdx]: " << E_B_lastIt_[storageIdx].ToString() << std::endl;
				//   std::cout << "E_H_lastIt_[storageIdx]: " << E_H_lastIt_[storageIdx].ToString()<< std::endl;
				//   std::cout << "P_M_lastIt_[storageIdx]: " << P_M_lastIt_[storageIdx].ToString()<< std::endl;

				//   std::cout << "E_B_[storageIdx]: " << E_B_[storageIdx].ToString() << std::endl;
				//   std::cout << "E_H_[storageIdx]: " << E_H_[storageIdx].ToString()<< std::endl;
				//   std::cout << "P_M_[storageIdx]: " << P_M_[storageIdx].ToString()<< std::endl;

				//   std::cout << "currentSolution: " << currentSolution.ToString() << std::endl;
				//   std::cout << "currentInput: " << currentInput.ToString() << std::endl;
				//   std::cout << "currentOutput: " << currentOutput.ToString() << std::endl;


				Vector<Double> toDiffY;
				if (RUN_useLastTS_) {
					toDiffY = P_M_lastTS_[storageIdx];
				} else {
					toDiffY = P_M_lastIt_[storageIdx];
				}

				deltaY = currentOutput;
				deltaY -= toDiffY;

				//      std::cout << "Current Polarization: " << currentOutput.ToString() << std::endl;
				//      std::cout << "Previous Polarization: " << toDiffY.ToString() << std::endl;
				//      std::cout << "Current Polarization Norm: " << currentOutput.NormL2() << std::endl;
				//      std::cout << "Previous Polarization Norm: " << toDiffY.NormL2() << std::endl;
				//      std::cout << "deltaY: " << deltaY.ToString() << std::endl;
				//      std::cout << "deltaY.NormL2(): " << deltaY.NormL2() << std::endl;
				//      std::cout << "abs(deltaY.NormL2() - 2*MAT_ySat_): " << abs(deltaY.NormL2() - 2*MAT_ySat_) << std::endl;

				if (abs(currentOutput.NormL2() - MAT_ySat_) <= tol_) {
					// hyst operator is now in saturation
					if (abs(toDiffY.NormL2() - MAT_ySat_) > tol_) {
						// and was not in saturation the last time
						intoSat = true;
					} else {
						// and is in saturation
						// now we have to check if it stayed in the same saturation (i.e. that we did not
						// switch from positive to negative saturation!
						// allow much larger tolerance here
						if (abs(deltaY.NormL2() - 2 * MAT_ySat_) < 1e-8) {
							satToSat = true;
						}
					}
				} else {
					// hyst operator is not in saturation
					if (abs(toDiffY.NormL2() - MAT_ySat_) <= tol_) {
						// but was in saturation the last time
						outofSat = true;
					}
				}
			}
		}

		/*
		 * new we have deltaX and deltaY ready for the call to createDeltaMatrix
		 * we only have to determine what to do with deltaX and deltaY
		 */
		std::string evalMethodName;
		UInt lastDigit = XML_deltaEvalVersion_ % 10;

		if (lastDigit == 0) {
			evalMethodName = "Classic";
		} else if (lastDigit == 1) {
			evalMethodName = "DirectDivisionAbsSatStepDownNew";
		} else if (lastDigit == 2) {
			evalMethodName = "DirectDivisionAbsCut";
		} else if (lastDigit == 3) {
			evalMethodName = "GeometricMean";
		} else if (lastDigit == 4) {
			evalMethodName = "VeryNewApproach";
		} else {
			EXCEPTION("Evaluation approach not implemented yet");
		}
		//std::cout << "evalMethodName: " << evalMethodName << std::endl;

		CreateDeltaMatrix(deltaX, deltaY, outputTensor, evalMethodName, storageIdx, intoSat, outofSat, satToSat, currentInput);

		deltaMat_[storageIdx] = outputTensor;

	}

	void CoefFunctionHyst::CreateDeltaMatrix(Vector<Double>& dX, Vector<Double>& dY, Matrix<Double>& outputTensor, std::string evalMethod, UInt storageIdx, bool intoSat, bool outofSat, bool satToSat, Vector<Double>& X_current) {

		if (XML_textOutputLevel_ == 2) {
			std::cout << "### CreateDeltaMatrix ###" << std::endl;
			std::cout << "dX: " << dX.ToString() << std::endl;
			std::cout << "dY: " << dY.ToString() << std::endl;
			std::cout << "intoSat: " << intoSat << std::endl;
			std::cout << "outofSat: " << outofSat << std::endl;
			std::cout << "satToSat: " << satToSat << std::endl;
		}

		//  std::cout << "RUN_deltaMatComputedDuringCurrentIt_["<<idx<<"] = " << RUN_deltaMatComputedDuringCurrentIt_[idx] << std::endl;
		//  if(RUN_deltaMatComputedDuringCurrentIt_[idx] == false){
		//    // compute deltamat only once per index per iteration
		//    // RUN_cuttingAlreadyCounted_ has to be reset after each iteration
		//    RUN_deltaMatComputedDuringCurrentIt_[idx] = true;
		//  } else {
		//    // we only need to comput the deltaMat once for each idx per iteration
		//    // for full evaluation, each integration point correspo
		//    return;
		//  }

		Double sign;
		outputTensor = Matrix<Double>(dim_, dim_);

		if (PDEName_ == "Electromagnetics") {
			/*
			 * B = mu0*(H + M)
			 * -> H = nu0*B - M
			 *
			 * -> deltaMat = nu' - deltaM/deltaB
			 *                   ^
			 *             this sign is meant
			 */
			sign = -1.0;
		} else {
			/*
			 * D = eps0*E + P
			 *
			 * -> deltaMat = eps' + deltaP/deltaE
			 *                    ^
			 *             this sign is meant
			 */
			sign = 1.0;
		}

		if (evalMethod == "Classic") {
			// check if dX Vector is very close to zero > reuse old state
			if (dX.NormL2() < MAT_ampResolution_) {
				outputTensor = deltaMat_lastIt_[storageIdx];
			} else {
				Double addValue = 0;
				// go over single components
				for (UInt i = 0; i < dim_; i++) {

					if (RUN_tensorToAdd_ == 1) {
						addValue = MAT_initialTensor_[i][i];
					} else if (RUN_tensorToAdd_ == 2) {
						addValue = rev_mat_fac_;
					} else {
						EXCEPTION("Either addFreeFieldTensorToDeltaMat_ or addInitialTensorToDeltaMat_ should be true");
					}

					// check if single component is smaller than the tolerance
					// if so: set dX[i] to MAT_ampResolution_
					if (abs(dX[i]) == 0) {
						outputTensor[i][i] = deltaMat_lastIt_[storageIdx][i][i]; // dangerous
					} else {
						outputTensor[i][i] = addValue + sign * abs(dY[i] / dX[i]);
					}
					//        if(abs(dX[i]) < MAT_ampResolution_ ){
					//          outputTensor[i][i] = addValue + sign*abs(dY[i]/MAT_ampResolution_);
					//        } else {
					//          outputTensor[i][i] = addValue + sign*abs(dY[i]/dX[i]);
					//        }
				}
			}
		}
		else if (evalMethod == "VeryNewApproach") {
			Double addValue = 0.0;

			if (dY_sol[storageIdx].NormL2() > tol_) {
				// use bisection
				std::cout << "Use bisection" << std::endl;
				std::cout << "StorageIDX: " << storageIdx << std::endl;

				Vector<Double> dY_current = dY;
				dY_current += dX*rev_mat_fac_; //dY = deltaP whereas dY_sol = deltaD

				// compare dY_sol to dY
				if (dY_current.NormL2() > dY_sol[storageIdx].NormL2()) {
					// actual distance between starting point and current point is larger than target value
					// > stepped to far
					// > X_up = X_current
					std::cout << "dY_current.NormL2 > dY_sol.NormL2" << std::endl;
					X_up[storageIdx] = X_current;
				} else {
					// we stepped not far enough
					// > X_low = X_current
					std::cout << "dY_current.NormL2 <= dY_sol.NormL2" << std::endl;
					X_low[storageIdx] = X_current;
				}
				Vector<Double> diff = dY_current;
				diff -= dY_sol[storageIdx];
				std::cout << "dY_current: " << dY_current.ToString() << std::endl;
				std::cout << "dY_sol: " << dY_sol[storageIdx].ToString() << std::endl;
				std::cout << "dY - dY_sol: " << diff.ToString() << std::endl;

				std::cout << "dY_current.NormL2: " << dY_current.NormL2() << std::endl;
				std::cout << "dY_sol.NormL2: " << dY_sol[storageIdx].NormL2() << std::endl;

				Vector<Double> X_target = Vector<Double>(dim_);
				X_target.Add(0.5, X_low[storageIdx], 0.5, X_up[storageIdx]);

				std::cout << "X_low = " << X_low[storageIdx].ToString() << std::endl;
				std::cout << "X_up = " << X_up[storageIdx].ToString() << std::endl;
				std::cout << "X_target = " << X_target.ToString() << std::endl;
				std::cout << "X_current = " << X_current.ToString() << std::endl;

				for (UInt i = 0; i < dim_; i++) {
					outputTensor[i][i] = abs(dY_sol[storageIdx][i] / (X_target[i] - E_H_lastTS_[storageIdx][i]));
				}
			} else {

				// check if dX Vector is very close to zero > reuse old state
				if (dX.NormL2() < MAT_ampResolution_) {
					outputTensor = deltaMat_lastIt_[storageIdx];
				} else {
					// go over single components
					for (UInt i = 0; i < dim_; i++) {

						if (RUN_tensorToAdd_ == 1) {
							addValue = MAT_initialTensor_[i][i];
						} else if (RUN_tensorToAdd_ == 2) {
							addValue = rev_mat_fac_;
						} else {
							EXCEPTION("Either addFreeFieldTensorToDeltaMat_ or addInitialTensorToDeltaMat_ should be true");
						}

						std::cout << "AddValue: " << addValue << std::endl;

						// check if single component is smaller than the tolerance
						// if so: set dX[i] to MAT_ampResolution_
						if (abs(dX[i]) < MAT_ampResolution_) {
							outputTensor[i][i] = addValue + sign * abs(dY[i] / MAT_ampResolution_);
						} else {
							outputTensor[i][i] = addValue + sign * abs(dY[i] / dX[i]);
						}
					}
				}

				/*
				 * new check special cases
				 * 1 > go into saturation
				 *   > deltaMat may have to drop down significantly from e-6 or e-8 to e-12
				 * 2 > leave saturation
				 *   > deltaMat may have to increase significantly from e-12 to e-6 or e-8
				 * 3 > jump from saturation to opposite saturation
				 *   > extrem case of 2
				 */
				std::cout << "RUN_residualDecreased_ " << RUN_residualDecreased_ << std::endl;

				//      if((intoSat == true)&&(RUN_residualDecreased_ == false)){
				//        // std::cout << "TMP deltaMat: " << outputTensor.ToString() << std::endl;
				//        std::cout << "Going into staturation: Reduce slope by applying geometric mean with rev_mat_fac_" << std::endl;
				//        // you start with a quite large slope but the actual required slope may
				//        // be much smaller
				//        // > build geometric mean of current state and rev_mat_fac
				//        for(UInt i = 0; i < dim_; i++){
				//          outputTensor[i][i] = sqrt(outputTensor[i][i] * rev_mat_fac_);
				//        }
				//      }

				//      if(satToSat == true){
				//        // hard to solve case; try to solve by bisection; from now on (regardless of flag intoSat)
				//        // we will solve by bisection
				//
				//        // get actual dY_sol which should be achieved for this integration point/element
				//        //
				//        //    deltaX_new = (deltaMat_old)^-1*dY_sol
				//        //    > dY_sol = deltaMat_old*deltaX_new
				//        //
				//
				//        dY_sol[storageIdx] = deltaMat_lastIt_[storageIdx]*dX;
				//
				//        // setup points at ends of interval
				//        X_low[storageIdx] = E_H_lastTS_[storageIdx];
				//        X_up[storageIdx] = -E_H_lastTS_[storageIdx];; // dX = X_current - X_lastIt or X_lastTS
				//        // lets hope that X_sol acutally lies on the line between X_low and X_up
				//
				//        std::cout << "Actual dY to be achieved: " << dY_sol[storageIdx].ToString() << std::endl;
				//
				//        Vector<Double> X_target = Vector<Double>(dim_);
				//        X_target.Add(0.5,X_low[storageIdx],0.5,X_up[storageIdx]);
				//
				//        std::cout << "X_low = " << X_low[storageIdx].ToString() << std::endl;
				//        std::cout << "X_up = " << X_up[storageIdx].ToString() << std::endl;
				//        std::cout << "X_target = " << X_target.ToString() << std::endl;
				//
				//        for(UInt i = 0; i < dim_; i++){
				//          outputTensor[i][i] = abs(dY_sol[storageIdx][i]/( X_target[i] - E_H_lastTS_[storageIdx][i] ));
				//        }
				//
				//      }
			}
			//
			//    // important addition:
			//    // check if change to last outputTensor is larger than some tolerance
			//    // in that case use a geometric mean to smooth out abrupt change
			//    Matrix<Double> tmp = outputTensor;
			//    tmp -= deltaMat_lastIt_[storageIdx];
			//    std::cout << "Difference between current deltaMat ( " << std::endl;
			//    std::cout << outputTensor.ToString() << std::endl;
			//    std::cout << ") and matrix from previous iteration (" << std::endl;
			//    std::cout << deltaMat_lastIt_[storageIdx] << std::endl;
			//    std::cout << ") is " << tmp.NormL2() << std::endl;
			//    Matrix<Double> tmp2 = outputTensor;
			//    tmp2 += deltaMat_lastIt_[storageIdx];
			//
			//    std::cout << "Relative difference: " << std::endl;
			//    Double relDifference = tmp.NormL2()/tmp2.NormL2();
			//    std::cout << relDifference << std::endl;
			//
			//    if(relDifference > 0.10){
			//      std::cout << "Difference between current deltaMat ( " << std::endl;
			//      std::cout << outputTensor.ToString() << std::endl;
			//      std::cout << ") and matrix from previous iteration (" << std::endl;
			//      std::cout << deltaMat_lastIt_[storageIdx] << std::endl;
			//      std::cout << ") is larger then 1e5: create geometric mean = " << std::endl;
			//      for(UInt i = 0; i < dim_; i++){
			//        outputTensor[i][i] = sqrt(outputTensor[i][i] * deltaMat_lastIt_[storageIdx][i][i] );
			//      }
			//      std::cout << outputTensor.ToString() << std::endl;
			//    }
			//
			//  }else if(evalMethod == "DirectDivisionAbsSatStepDownNew"){
			//
			//    if(XML_textOutputLevel_ == 2){
			//      std::cout << "DirectDivisionAbsSatStepDownNew" << std::endl;
			//    }
			//
			//    Double addValue;
			//    for(UInt i = 0; i < dim_; i++){
			//      if(RUN_tensorToAdd_ == 1){
			//        addValue = MAT_initialTensor_[i][i];
			//
			//      } else if(RUN_tensorToAdd_ == 2){
			//        addValue = rev_mat_fac_;
			//
			//      } else {
			//        EXCEPTION("Either addFreeFieldTensorToDeltaMat_ or addInitialTensorToDeltaMat_ should be true");
			//      }
			//
			//      if( (dY[i] == 0) && (deltaMat_lastIt_[idx][i][i] > 10*rev_mat_fac_)) {
			//        /*
			//         * when we reach into saturation, we may step into it with a slope >> rev_mat_fac
			//         * but stepping back (in case that solution was overestimated) is done with rev_mat_fac
			//         * > solution update has to be extremely large and negative which might ruin the whole
			//         * solution process;
			//         * idea: slowly turn down slope
			//         */
			//
			//        std::cout << "Use step down " << std::endl;
			//
			//        deltaMat_[idx][i][i] = deltaMat_lastIt_[idx][i][i]/10;
			//      } else {
			//        if(dX[i] != 0){
			//          deltaMat_[idx][i][i] = addValue + sign*abs(dY[i]/dX[i]);
			//        } else {
			//          //stick with last value
			//          deltaMat_[idx][i][i] = deltaMat_lastIt_[idx][i][i];
			//        }
			//      }
			//    }
			//  } else if (evalMethod == "DirectDivisionAbsNew"){
			//
			//    if(XML_textOutputLevel_ == 2){
			//      std::cout << "DirectDivisionAbsNew" << std::endl;
			//    }
			//
			//    Double addValue;
			//    for(UInt i = 0; i < dim_; i++){
			//      if(RUN_tensorToAdd_ == 1){
			//        addValue = MAT_initialTensor_[i][i];
			//
			//      } else if(RUN_tensorToAdd_ == 2){
			//        addValue = rev_mat_fac_;
			//
			//      } else {
			//        EXCEPTION("Either addFreeFieldTensorToDeltaMat_ or addInitialTensorToDeltaMat_ should be true");
			//      }
			//
			//      if(dX[i] != 0){
			//        deltaMat_[idx][i][i] = addValue + sign*abs(dY[i]/dX[i]);
			//      } else {
			//        //stick with last value
			//        deltaMat_[idx][i][i] = deltaMat_lastIt_[idx][i][i];
			//      }
			//    }
			//  } else if (evalMethod == "DirectDivisionAbsCut"){
			//
			//    if(XML_textOutputLevel_ == 2){
			//      std::cout << "DirectDivisionAbsCut" << std::endl;
			//    }
			//
			//    Double addValue;
			//    for(UInt i = 0; i < dim_; i++){
			//      if(RUN_tensorToAdd_ == 1){
			//        addValue = MAT_initialTensor_[i][i];
			//
			//      } else if(RUN_tensorToAdd_ == 2){
			//        addValue = rev_mat_fac_;
			//
			//      } else {
			//        EXCEPTION("Either addFreeFieldTensorToDeltaMat_ or addInitialTensorToDeltaMat_ should be true");
			//      }
			//
			//
			//
			//
			//      if(dX[i] == 0){
			//        std::cout << "Use last value as dX["<<i<<"] is 0" << std::endl;
			//        //stick with last value
			//        std::cout << "idx: " << idx << std::endl;
			//        Double downStepping = std::pow(0.75,Double(cuttingApplied_[idx]));
			//
			//        std::cout << "Downstepping: " << downStepping << std::endl;
			//        Double tmpValue = deltaMat_lastIt_[idx][i][i]*downStepping;
			//        if(tmpValue < rev_mat_fac_){
			//          tmpValue = rev_mat_fac_;
			//        }
			//        deltaMat_[idx][i][i] = tmpValue;
			//
			//      } else {
			//        if( abs(dX[i]) >= 2*MAT_xSat_){
			//          std::cout << "abs(dX["<<i<<"]) is " << abs(dX[i]) <<">"<< 2*MAT_xSat_<< std::endl;
			//          deltaMat_[idx][i][i] = addValue + sign*abs(dY[i]/dX[i]/2);
			//        } else {
			//          std::cout << "use std division" << std::endl;
			//          deltaMat_[idx][i][i] = addValue + sign*abs(dY[i]/dX[i]);
			//        }
			//      }
			//    }
			//  } else if(evalMethod == "GeometricMean"){
			//    std::cout << "GeometricMean" << std::endl;
			//    std::cout << "RUN_cuttingAlreadyCounted_["<<idx<<"] = " << RUN_cuttingAlreadyCounted_[idx] << std::endl;
			//    if(RUN_cuttingAlreadyCounted_[idx] == 0){
			//      // compute deltamat only once per index per iteration
			//      // RUN_cuttingAlreadyCounted_ has to be reset after each iteration
			//      RUN_cuttingAlreadyCounted_[idx] = 1;
			//
			//      if(RUN_tensorToAdd_ == 1){
			//        std::cout << "Compute geometric mean from eps0 (nu0) and pSat/eSat (hSat/mSat)" << std::endl;
			//        // normally this flag means to add epsInitial
			//        // in this case we use it for the computation of an initial stepping
			//        // calculated from the geometric mean of rev_mat_fac and ySat/xSat
			//        if(PDEName_ == "Electromagnetics"){
			//          for(UInt i = 0; i < dim_; i++){
			//            deltaMat_[idx][i][i] = sqrt(rev_mat_fac_ * MAT_xSat_/MAT_ySat_);
			//          }
			//        } else {
			//          for(UInt i = 0; i < dim_; i++){
			//            deltaMat_[idx][i][i] = sqrt(rev_mat_fac_ * MAT_ySat_/MAT_xSat_);
			//          }
			//        }
			//      } else {
			//        std::cout << "Compute geometric mean from previous and current deltaMats" << std::endl;
			//        // compute current deltaMat with standard division approach
			//        Matrix<Double> current = Matrix<Double>(dim_,dim_);
			//        current.Init();
			//
			//        for(UInt i = 0; i < dim_; i++){
			//          if( abs(dX[i]) == 0){
			//            // reuse value if dX = 0
			//            current[i][i] = deltaMat_[idx][i][i];
			//          } else {
			//            //std::cout << "use std division" << std::endl;
			//            current[i][i] = rev_mat_fac_ + sign*abs(dY[i]/dX[i]);
			//          }
			//          // make geometric mean
			//          // note: older values will loose more and more impact as multiple
			//          //       root are take of them
			//          deltaMat_[idx][i][i] = sqrt(deltaMat_[idx][i][i] * current[i][i]);
			//        }
			//      }
			//
			//    } else {
			//      // reuse current state > deltaMat_ not changed
			//      return;
			//    }
			//
			//  } else if(evalMethod == "DirectDivisionCheckForSat"){
			//    std::cout << "DirectDivisionCheckForSat" << std::endl;
			//
			///*
			// * New approach:
			// *
			// *  case I - no saturation reached
			// *    > used DirectDivisionAbs > works well
			// *  case II - previous and current state in saturation
			// *    > works well with DirectDivisionAbs, too
			// *  case III - going into saturation (previousState not in saturation, current state in saturation)
			// *    > DirectDivisionAbs does not work well here as we grow with slope >> eps0/nu0 so that
			// *      computed increment will be small even though we would need a larger one
			// *    > Idea: check for this case
			// *
			// *
			// */
			//
			//    if(RUN_tensorToAdd_ == 1){
			//      std::cout << "Compute geometric mean from eps0 (nu0) and pSat/eSat (hSat/mSat)" << std::endl;
			//      // normally this flag means to add epsInitial
			//      // in this case we use it for the computation of an initial stepping
			//      // calculated from the geometric mean of rev_mat_fac and ySat/xSat
			//      if(PDEName_ == "Electromagnetics"){
			//        for(UInt i = 0; i < dim_; i++){
			//          deltaMat_[idx][i][i] = sqrt(rev_mat_fac_ * MAT_xSat_/MAT_ySat_);
			//        }
			//      } else {
			//        for(UInt i = 0; i < dim_; i++){
			//          deltaMat_[idx][i][i] = sqrt(rev_mat_fac_ * MAT_ySat_/MAT_xSat_);
			//        }
			//      }
			//    } else {
			//      std::cout << "Compute geometric mean from previous and current deltaMats" << std::endl;
			//      // compute current deltaMat with standard division approach
			//      Matrix<Double> current = Matrix<Double>(dim_,dim_);
			//      current.Init();
			//
			//      for(UInt i = 0; i < dim_; i++){
			//        if( abs(dX[i]) == 0){
			//          // reuse value if dX = 0
			//          current[i][i] = deltaMat_[idx][i][i];
			//        } else {
			//          //std::cout << "use std division" << std::endl;
			//          current[i][i] = rev_mat_fac_ + sign*abs(dY[i]/dX[i]);
			//        }
			//        // make geometric mean
			//        // note: older values will loose more and more impact as multiple
			//        //       root are take of them
			//        deltaMat_[idx][i][i] = sqrt(deltaMat_[idx][i][i] * current[i][i]);
			//      }
			//    }

		} else {
			EXCEPTION("Method " << evalMethod << " not implemented yet");
		}

		if (XML_textOutputLevel_ == 2) {
			std::cout << "deltaMat_: " << outputTensor.ToString() << std::endl;
		}

	}

	void CoefFunctionHyst::SetRuntimeDependentFlag(std::string flagName, UInt intState) {

		/*
		 * integer flags
		 */
		if (flagName == "evaluationPurpose") {
			/*
			 * 1 -> assemble
			 * 2 -> store iteration values; hyst memory locked
			 * 3 -> store ts values; hyst memory overwritten
			 * 4 -> output
			 *
			 * -> details: see .hh file
			 */
			if (intState < 1) {
				RUN_evaluationPurpose_ = 1;
			} else if (intState > 4) {
				RUN_evaluationPurpose_ = 4;
			} else {
				RUN_evaluationPurpose_ = intState;
			}
		} else if (flagName == "estimateSlope") {
			// here we do not acutally set a flag but instead execute the estimateCurrentSlope function
			Double fac = std::pow(10, -1.0 * Double(intState));
			//std::cout << "fac: " << fac << std::endl;
			EstimateCurrentSlope(Vector<Double>(dim_), fac);

		} else if (flagName == "vectorToReturn") {
			/*
			 *    1 -> evaluate Hysteresis operator and returns its value in rhs load
			 *          integrator
			 *    2 -> return zero vector (needed for implementations where we do not
			 *          put P/M to the rhs or only for the first two iteratons)
			 *    3 -> evaluate Hysteresis operator but return not its value, but the
			 *          difference to the value from the last TS
			 */
			if (intState < 1) {
				RUN_vectorToReturn_ = 1;
			} else if (intState > 4) {
				RUN_vectorToReturn_ = 4;
			} else {
				RUN_vectorToReturn_ = intState;
			}
		} else if (flagName == "tensorToReturn") {
			/*
			 *  1 -> initialTensor
			 *  2 -> freeFieldTensor
			 *  3 -> deltaMatrix
			 *  4 -> Pmax/Emax // Mmax/Bmax
			 *  5 -> estimatedState
			 */
			if (intState < 1) {
				RUN_tensorToReturn_ = 1;
			} else if (intState > 5) {
				RUN_tensorToReturn_ = 5;
			} else {
				RUN_tensorToReturn_ = intState;
			}

			if (RUN_tensorToReturn_ == 3) {
				// return deltaMatrix
				deltaMatActive_ = true;
			} else {
				deltaMatActive_ = false;
			}

		} else if (flagName == "tensorToAdd") {
			/*
			 *  1 -> add InitialTensor
			 *  2 -> add FreeFieldTensor
			 */
			if (intState < 1) {
				RUN_tensorToAdd_ = 1;
			} else if (intState > 2) {
				RUN_tensorToAdd_ = 2;
			} else {
				RUN_tensorToAdd_ = intState;
			}
		}
			/*
			 * boolean flags
			 */
		else if (flagName == "useLastTS") {
			/*
			 * flagName == "useLastTS"
			 * 1: old values Y_old and X_old used in deltaComputation
			 *      (see flags useDeltaY and useDeltaX) will be the
			 *      ones of the last timestep rather than from the last
			 *      iteration. I.e.
			 *       Y_old = YLastTS_ , YLastTSVEC_
			 *       X_old = XLastTS_ , XLastTSVEC_
			 * 0:
			 *     old values of Y_old and X_old will be taken from
			 *     last iteration
			 *       Y_old = YpreviousIt_ , YpreviousItVEC_
			 *       X_old = XpreviousIt_ , XpreviousItVEC_
			 */
			RUN_useLastTS_ = bool(intState);

		} else if (flagName == "allowBMP") {
			/*
			 * flagName == "allowBMP"
			 * 1: computeXY_vec outputs the overlaid switching and rotation state of
			 *      the vector Preisach model as a BMP image, but only if printOut was
			 *      set to true, too (printOut is to be set via material file)
			 * 0: disables output of bmp, even if printOut was set to true
			 */
			RUN_allowBMP_ = bool(intState);

		} else if (flagName == "residualDecreased") {

			RUN_residualDecreased_ = bool(intState);

		} else if (flagName == "overwriteDirection") {
			/*
			 * flagName = "overwriteDirection"
			 * 1: rotation and switching operator in vector model as set as usual
			 * 0: rotation state remains locked; only switching state may change
			 */
			RUN_overwriteDirection_ = bool(intState);
		} else if (flagName == "cuttingAlreadyCounted") {
			for (UInt i = 0; i < numHystOperators_; i++) {
				RUN_cuttingAlreadyCounted_[i] = intState;
			}
		} else if (flagName == "deltaMatComputedDuringCurrentIt") {
			for (UInt i = 0; i < numHystOperators_; i++) {
				RUN_deltaMatComputedDuringCurrentIt_[i] = (intState);
			}
		} else if (flagName == "deltaMatComputedDuringCurrentTS") {
			for (UInt i = 0; i < numHystOperators_; i++) {
				RUN_deltaMatComputedDuringCurrentTS_[i] = (intState);
			}


		} else if (flagName == "deltaForm") {
			/*
			 * RUN_deltaForm_ = -1 > delta = current - last ts
			 *                       rhs = last ts
			 * RUN_deltaForm_ = 0  > no delta mat
			 *                       rhs = current
			 * RUN_deltaForm_ = x  > delta between currernt and last xth iteration value
			 *                       rhs = last xth iteration value
			 */
			if (intState == 0) {
				RUN_deltaForm_ = 0;
				timeLevel_Mat_ = 0;
				timeLevel_RHS_ = 0;
				timeLevel_BC_ = 0;
				timeLevel_Output_ = 0;
			} else if (int(intState) == -1) {
				RUN_deltaForm_ = -1;
				timeLevel_Mat_ = -1;
				timeLevel_RHS_ = -1;
				timeLevel_BC_ = 0;
				timeLevel_Output_ = 0;
			} else {
				RUN_deltaForm_ = intState;
				timeLevel_Mat_ = intState;
				timeLevel_RHS_ = intState;
				timeLevel_BC_ = 0;
				timeLevel_Output_ = 0;
			}
		} else if (flagName == "resetReeval") {
//			std::cout << "Reset reeval flag" << std::endl;
			for (UInt i = 0; i < numStorageEntries_; i++) {
				requiresReeval_[i] = true;
				rotatedCouplingTensor_requiresReeval_[i] = true;
			}
		} else {
			WARN("flagName = " << flagName << " not known");
		}
	}

	void CoefFunctionHyst::SetInputDependentFlags(UInt multiDigitInteger) {
		/*
		 * Input:
		 *  integer of multiple digits; will be divided into single digits which
		 *  then are used to set corresponding flags
		 *
		 *  multiDigitInteger = e d cb a
		 *  a = evaluation depth
		 *  cb = evaluation version of delta matrix
		 *  d = approximation of H in case of magnetics
		 *  e = enable performance measurement
		 *  f = depth of text output
		 *  ... possible more
		 *
		 */

		if (XMLParameterSet_ == false) {
			UInt a = multiDigitInteger % 10;
			/*
			 * 1 -> standard evaluation
			 * 2 -> extended evaluation
			 * 3 -> full evaluation
			 *
			 * -> details: see .hh file
			 */
			if (a < 1) {
				XML_EvaluationDepth_ = 1;
			} else if (a > 3) {
				XML_EvaluationDepth_ = 3;
			} else {
				XML_EvaluationDepth_ = a;
			}

			multiDigitInteger = multiDigitInteger / 10;
			UInt bc = multiDigitInteger % 100;
			/*
			 * 1x -> standard, delta between current input/ouput and last input/output
			 * 2x -> standard but with slightly shifted input if diff is close to 0 in single
			 *        but not all components
			 *
			 * -> details: see .hh file
			 */
			if (bc < 10) {
				XML_deltaEvalVersion_ = 10;
			} else if (bc >= 30) {
				XML_deltaEvalVersion_ = 20;
			} else {
				XML_deltaEvalVersion_ = bc;
			}

			multiDigitInteger = multiDigitInteger / 100;
			UInt d = multiDigitInteger % 10;

			/*
			 *  1 -> H_new = nu0*B_new - M_old
			 *  2 -> H_new = deltaMat_lastIt_ * (B_new - B_old) + H_old
			 *  3 -> as 1 but with iterative refinement
			 */
			if (PDEName_ == "Electromagnetics") {
				if (d < 1) {
					XML_HApproxVersion_ = 1;
				} else if (d > 3) {
					XML_HApproxVersion_ = 3;
				} else {
					XML_HApproxVersion_ = d;
				}
			} else {
				// not needed in electrostatics
				XML_HApproxVersion_ = 0;
			}

			multiDigitInteger = multiDigitInteger / 10;
			UInt e = multiDigitInteger % 10;
			/*
			 * flagName == "XML_performanceMeasurement_"
			 * 2: do performance measurement also in Preisach operator
			 * 1: measure time it takes to evaluate hysteresis operator
			 *    (and maybe other functions)
			 * 0: do not measure time
			 */
			if (e > 2) {
				XML_performanceMeasurement_ = 2;
			} else {
				XML_performanceMeasurement_ = e;
			}

			multiDigitInteger = multiDigitInteger / 10;
			UInt f = multiDigitInteger % 10;
			/*
			 * flagName == "XML_textOutputLevel_"
			 * 0: no text output
			 * 1: info output
			 * 2: debug output
			 */
			if (f > 2) {
				XML_textOutputLevel_ = 2;
			} else {
				XML_textOutputLevel_ = f;
			}

			//std::cout << "XML_EvaluationDepth_: " << XML_EvaluationDepth_ << std::endl;
//			std::cout << "init storage" << std::endl;
			InitStorage();

			if (XML_textOutputLevel_ > 0) {
				// has to be called after InitStorage as otherwise some parameters have not been read in
				std::cout << this->ToString() << std::endl;
			}

			/*
			 * also set output flags for hysteresis operator
			 * (after InitStorage was called, as this sets up the hyst operators)
			 */
			multiDigitInteger = multiDigitInteger / 10;
			UInt mappingFlag = multiDigitInteger % 10;
			//std::cout << "multiDigitInteger_end: " << multiDigitInteger << std::endl;
			hyst_->setFlags(XML_performanceMeasurement_, XML_textOutputLevel_, mappingFlag);

			XMLParameterSet_ = true;
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

  void CoefFunctionHyst::TestInversion(Double eps_mu){
    // Tests inversion of VECTORPreisach Operator and returns 
    //  estimation for alpha (linesearch stepping paramater) that can be used
    //  for later computation
    // > can be used in actual simulations to setup an appropriate linesearch factor alpha
    //   to accelerate the actual computations later on!
    
    UInt numTests = 200;
    Vector<Double>* xIn = new Vector<Double>[numTests]; 
    
    /*
     * Testcase: Put material into full saturation in y-direction; from that state on
     * apply a decreasing triangular signal in x-direction; for initial period hold
     * signla in y-direction; during first real period decrease to 1/2 of value;
     * then to 1/4 and finally to 0
     */
    xIn[0] = Vector<Double>(dim_);
    xIn[0].Init(0.0);
    xIn[0][1] = -MAT_xSat_;
    
    UInt numStepsFirstIncrease = (UInt) ( (Double)numTests/8.0 ) +1;
    Double incr = MAT_xSat_/((Double) numStepsFirstIncrease-1);
    for(UInt i = 1; i < numStepsFirstIncrease; i++){
      xIn[i] = Vector<Double>(dim_);
      xIn[i].Init(0.0);
      xIn[i][0] = xIn[i-1][0] + incr; 
      xIn[i][1] = -MAT_xSat_;
    }
    UInt cnt = numStepsFirstIncrease;
    
    UInt remainingSteps = numTests - numStepsFirstIncrease;
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
        
        if(j==0){
          xIn[cnt][1] = -MAT_xSat_/2.0;
        } else if(j==1){
          xIn[cnt][1] = -MAT_xSat_/4.0;
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
        
    Vector<Double>* xRetrieved = new Vector<Double>[numTests];
    Vector<Double>* yRetrieved = new Vector<Double>[numTests];
    Vector<Double>* yError = new Vector<Double>[numTests];
    Vector<Double>* xError = new Vector<Double>[numTests];
    Vector<Double>* yIn = new Vector<Double>[numTests]; 
    
    bool isVirgin = true;
    // create one temporary VectorPreisach operator with the same
    // parameter as the later used ones, except, we just need it to store
    // entries for one single element / point
    
    if (MAT_methodType_ == SCALAR) {
      // this inversion test (as well as the overall inversion) is just used
      // for the VECTOR model (as there is a (hopefully) functioning inverse scalar
      // preisach model implemented
      WARN("Inversion test only required for vector model. Will do some tests nevertheless");
      // just take some random values to allow testing
      MAT_rotRes_ = 1.5;
      MAT_vecPreisachImplementationVersion_ = 2;
      MAT_angResistance_ = 10; // in degree
      MAT_angClipping_ = 0.0;
    }

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
    
//    float r = 2*static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX) - 1.0;
//    
//    xIn[0] = Vector<Double>(dim_);
//    xIn[0].Init(0.0);
//    xIn[0][0] = -r*MAT_xSat_;
//    
//    r = 2*static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX) - 1.0;
//    xIn[0][1] = -r*MAT_xSat_;
    
    std::ofstream testOutput;
    testOutput.open("ResultsOfHystInversionTest");
    testOutput << "## Results of TestInversion" << std::endl; 
    testOutput << "# Number of tests: " << numTests << std::endl;
    testOutput << "# eps_mu: " << eps_mu << std::endl;
    testOutput << "# xSat: " << MAT_xSat_ << std::endl;
    testOutput << "# ySat: " << MAT_ySat_ << std::endl;
    testOutput << "# Starting value of inversion " << xIn[0][0] << " " << xIn[0][1] << std::endl;
    testOutput << "# " << std::endl;
    testOutput << "# Testnumber    x_sol[0]    x_sol[1]    x_retrieved[0]    x_retrieved[1]    y_target[0]    y_target[1]    y_retrieved[0]    y_retrieved[1]    y_error[0]    y_error[1]" << std::endl;
    
    bool overwriteMemory = true;
    bool overwriteDirection = true;
    yIn[0] = Vector<Double>(dim_);
    yIn[0].Init(0.0);
    
    // Get polarization P
    yIn[0] = hystTMP->computeValue_vec(xIn[0], 0, overwriteMemory, overwriteDirection);
//    std::cout << "xIn[0]" << xIn[0].ToString() << std::endl;
//    std::cout << "yIn[0]" << yIn[0].ToString() << std::endl;
    // D = eps*E + P
    yIn[0].Add(eps_mu,xIn[0]);
       
    xRetrieved[0] = Vector<Double>(dim_);
    xRetrieved[0].Init(0.0);
    
    overwriteMemory = false;
    xRetrieved[0] = hystTMP->computeInput_vec(yIn[0], 0, eps_mu, alphaLinesearch_, overwriteMemory, overwriteDirection);
    
    // evaluate Preisach OP to get the actual output; DO NOT OVERWRITE MEMORY
    yRetrieved[0] = Vector<Double>(dim_);
    yRetrieved[0].Init(0.0);
    
    xError[0] = Vector<Double>(dim_);
    xError[0].Init(0.0);
    
    xError[0] = xRetrieved[0];
    xError[0] -= xIn[0];
    
    yError[0] = Vector<Double>(dim_);
    yError[0].Init(0.0);
    
    yRetrieved[0] = hystTMP->computeValue_vec(xRetrieved[0], 0, overwriteMemory, overwriteDirection);
    yError[0] = yRetrieved[0];
    yError[0] -= yIn[0];
    
    testOutput << "0    "<< xIn[0][0] <<"    "<< xIn[0][1] <<"    "<< xRetrieved[0][0] <<"    "<< xRetrieved[0][1] << "    "<<yIn[0][0]<<"    "<<yIn[0][1]<<"    "<<yRetrieved[0][0]<<"    "<<yRetrieved[0][1]<<"    "<<yError[0][0]<<"    "<<yError[0][1]<< std::endl;
    
//    Double inc = 2*MAT_xSat_/( (Double) numTests);
    
    for(UInt i = 1; i < numTests; i++){
      std::cout << "##### TEST NR " << i << "#####" << std::endl;
//      xIn[i] = Vector<Double>(dim_);
//      xIn[i][0] = xIn[i-1][0] + inc;
//      xIn[i][1] = xIn[i-1][0] + inc/2;
      overwriteMemory = true;
      yIn[i] = hystTMP->computeValue_vec(xIn[i], 0, overwriteMemory, overwriteDirection);
      yIn[i].Add(eps_mu,xIn[i]);
      
      overwriteMemory = false;
      xRetrieved[i] = hystTMP->computeInput_vec(yIn[i], 0, eps_mu, alphaLinesearch_, overwriteMemory, overwriteDirection);
      
      yRetrieved[i] = Vector<Double>(dim_);
      yRetrieved[i].Init(0.0);

      yError[i] = Vector<Double>(dim_);
      yError[i].Init(0.0);
      
      yRetrieved[i] = hystTMP->computeValue_vec(xRetrieved[i], 0, overwriteMemory, overwriteDirection);
      yError[i] = yRetrieved[i];
      yError[i] -= yIn[i];
      
      xError[i] = Vector<Double>(dim_);
      xError[i].Init(0.0);
      
      xError[i] = xRetrieved[i];
      xError[i] -= xIn[i];
    
      testOutput << i<<"    "<<xIn[i][0]<<"    "<<xIn[i][1]<<"    "<<xRetrieved[i][0]<<"    "<<xRetrieved[i][1]<<"    "<<xError[i][0]<<"    "<<xError[i][1]<<"    "<<yIn[i][0]<<"    "<<yIn[i][1]<<"    "<<yRetrieved[i][0]<<"    "<<yRetrieved[i][1]<<"    "<<yError[i][0]<<"    "<<yError[i][1]<< std::endl;
    
      std::cout << "##### TARGET Y-VECTOR #####" << std::endl;
      std::cout << yIn[i].ToString() << std::endl;
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
  
  
	std::string CoefFunctionHyst::ToString() const {

		static UInt initialOutput = 0;
		std::ostringstream oss;

		if ((XML_textOutputLevel_ >= 1) && (initialOutput == 0)) {
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
				<< "  deltaEvalVersion " << XML_deltaEvalVersion_ << "\n"
				<< "  H_approxVersion " << XML_HApproxVersion_ << "\n"
				<< "  performanceMeasurement " << XML_performanceMeasurement_ << "\n"
				<< "  textOutputLevel " << XML_textOutputLevel_ << "\n"
				<< "\n";
		}

		if (initialOutput > 0) {
			oss << "+++ CoefFunctionHyst +++\n";
			if (XML_textOutputLevel_ == 2) {
				oss << "++ Current state of runtime parameter: \n"
					<< "  evaluationPurpose " << RUN_evaluationPurpose_ << "\n"
					<< "  vectorToReturn " << RUN_vectorToReturn_ << "\n"
					<< "  tensorToReturn " << RUN_tensorToReturn_ << "\n"
					<< "  tensorToAdd " << RUN_tensorToAdd_ << "\n"
					<< "  allowBMP " << RUN_allowBMP_ << "\n"
					<< "  useLastTS " << RUN_useLastTS_ << "\n";
			}

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
