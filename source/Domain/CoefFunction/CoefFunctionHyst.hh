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
                    PtrCoefFct dependency,
          SubTensorType tensorType,
          MaterialType matType,
          shared_ptr<FeSpace> ptFeSpace);

  //! Destructor
  virtual ~CoefFunctionHyst();

//  //! Initialize with data
//  void Init( BaseMaterial* const material, shared_ptr<ElemList> actSDList);

  //! Return size of vector in case coefficient function is a vector
  virtual UInt GetVecSize() const {
    return dependCoef_->GetVecSize();
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

private:

  void InitStorage();

  void CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& outputTensor, std::string evalMethod,
                         UInt storageIdx, bool intoSat, bool outofSat, bool satToSat, Vector<Double>& X_current);

  void ExtractSolutionAndInputForHystOperator(Vector<Double>& extractedSolution, Vector<Double>& extractedInput, UInt& operatorIndex,
                                                 UInt& storageIndex, const LocPointMapped& lpm, bool midpointOnly);

  UInt EvaluateHysteresisOperator( Vector<Double> inputToHystOperator, Vector<Double>& outputOfHystOperator,
                                   UInt operatorIndex, UInt storageIndex, bool overwriteMemory, Vector<Double>& currentSolution );

  bool EvaluateAtMidpointOnly();
  bool OverwriteHystMemory();

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
   * determines whether the vector or the scalar model shall be usedl
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
  PtrCoefFct dependCoef_;

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
