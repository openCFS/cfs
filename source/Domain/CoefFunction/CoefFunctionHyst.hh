#ifndef FILE_COEFFUNCTION_HYST_HH
#define FILE_COEFFUNCTION_HYST_HH

#include <boost/tr1/type_traits.hpp>
#include "General/Environment.hh"
#include "CoefFunction.hh"
#include "Materials/BaseMaterial.hh"
#include "FeBasis/FeFunctions.hh"
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
					MaterialType matType);

  //! Destructor
  virtual ~CoefFunctionHyst();
  
  //! Initialize with data
  void Init( BaseMaterial* const material, shared_ptr<ElemList> actSDList);
  
  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, const LocPointMapped& lpm );

  void GetVector(Vector<Double>& coefScalar, const LocPointMapped& lpm);

  //! Return real-valued tensor at integration point
  void GetTensor(Matrix<Double>& tensor, const LocPointMapped& lpm );

  //!
  void SetPreviousHystVals(bool setNextToLastTS_too = false);

  //! Create for the vector case deltaMat from dX and dY
  void CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& deltaMat, UInt idxElem);

  //! \see CoefFunction::ToString
  std::string ToString() const;

  //! Return size of vector in case coefficient function is a vector
  virtual UInt GetVecSize() const {
	  return dependCoef_->GetVecSize();
  }

  //! Return row and columns size of tensor if coefficient function is a tensor
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      numRows = matTensorStart_.GetNumRows();
      numCols = matTensorStart_.GetNumCols();
  }

  void setOverwrite(bool overwrite_new){
    overwrite_ = overwrite_new;
  }

  void setOverwriteDirection(bool overwrite_new){
    overwriteDirection_ = overwrite_new;
  }

  void setUseNextToLastTS(bool useNextToLastTS_new){
    useNextToLastTS_ = useNextToLastTS_new;
  }

  void setDeltaComputation(bool deltaComputation_new){
    deltaComputation_ = deltaComputation_new;
  }

protected:
  
  //! compute X and Y value
  void ComputeXY( const LocPointMapped& lpm, Double& X, Double& Y, bool overwrite);

  //! compute X and Y vector for vector model
  void ComputeXY_vec( const LocPointMapped& lpm, Vector<Double>& X, Vector<Double>& Y, bool overwrite);

  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! hysteresis object
  Hysteresis * hyst_;

  //! previous Xval in hysteresis
  //! -> this is the value of the last timestep (if any)
  //! no longer needed -> instead XpreviousEval_ to check if this
  //! value was used during the last evaluation (there can be several evaluations
  //! during one iteration!)
  //Vector<Double> Xprevious_;
  //Vector<Double> Yprevious_;
  Vector<Double> XpreviousEval_;
  Vector<Double> YpreviousEval_;

  //! previous iteration value of hystersis
  //! -> this is the value of the last iteration (if any)
  //! -> needed to check if input has changed during iteration steps
  //! -> if not -> avoid recomputation and return YpreviousIt_
  Vector<Double> XpreviousIt_;
  Vector<Double> YpreviousIt_;

  Vector<Double> XnextToLastTS_;
  Vector<Double> YnextToLastTS_;

  //! delta material tensor
  //! new: extend to array of matrices so that we can store the old state of each element
  Matrix<Double>* matDeltaTensor_;

  //! non-delta matrial tensor (basically eps0, nu0 as diagonal matrix)
  Matrix<Double> matTensor_;

  //! for vector version
  //! XcurrentItVec_, YcurrentItVec_ store values of the current iteration; if function is called several times
  //                             with the same (or very close input), we avoid recomputations
  Vector<Double>* XcurrentItVEC_;
  Vector<Double>* YcurrentItVEC_;
  //! XpreviousItVec_, YcurrentItVec_ store values of the previous iteration; these vectors are needed to compute
  //                             the deltaMaterial deltaP/deltaE
  Vector<Double>* XpreviousItVEC_;
  Vector<Double>* YpreviousItVEC_;

  //!
  // values from the FIRST iteration of the last time step
  // > to estimate the deltaMatrxi for the first iteration of a new timestep
  //    it is not (always) a good idea to use the difference between current input
  //    and the previous iteration direction
  //    > reason: during the last iterations of a timestep, the steppings can be in
  //              totally different direction as the outer field is; if the delta matrix
  //              gets computed with such a difference, the first step of the new timestep
  //              (with possible new idbc values) might perform a huge step in the wrong
  //              direction
  //  > try to approximate the delta matrix for the first timestep by calculating
  //    delta between the current input (i.e. the solution of the last timestep) and
  //    the value at the first iteration of the last timestep (= solution of the next to
  //    last timestep)
  Vector<Double>* XnextToLastTSVEC_;
  Vector<Double>* YnextToLastTSVEC_;

  //! dXpreviousItVec_, dYcurrentItVec_ store values of dX and dY of the previous iteration;
  //                             we can compare these vectors to the current difference dX and dY to avoid
  //                             recomputation
  Vector<Double>* dXpreviousItVEC_;
  Vector<Double>* dYpreviousItVEC_;

  //! XpreviousEvalVec_, YpreviousEvalVec_ store the last used input/output of the hysteresis operator
  //                             that can avoid recomputations if the same input is used several times for example
  Vector<Double>* XpreviousEvalVEC_;
  Vector<Double>* YpreviousEvalVEC_;

  //! globale element to local element numbering
  std::map<UInt, UInt> globalElem2Local_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;

  //! type of material
   MaterialType matType_;

  //! type of subtensor
  SubTensorType tensorType_;

  //! direction to be traken from vector in case of sclalar hystersis
  Directions dirP_;

  //! initial material tensor
  Matrix<Double> matTensorStart_;

  //! list of all elements
  shared_ptr<ElemList> SDList_;

  //!
  BaseMaterial* material_;

  //! dim for vector model
  UInt dim_;

  //! this one is to distinguish between scalar and vector preisach
  //! do not confuse this with dimType_!
  CoefDimType methodType_;

  Double xSat_;
  Double ySat_;
  Double rotRes_;
  UInt printOut_;
  UInt bmpResolution_;
  Double rev_mat_fac_;

  UInt evalVersion_;
  UInt numRows_;
  Double tol_;

  //! flag needed in context of deltaStepping and linesearch;
  //! as different steps are tested during linesearch, we do not want them to have
  //! a permanent effect on the memory of the hysteresis operator
  //! if overwrite = false -> hysteresis operator applies changes only to a temporary copy
  bool overwrite_;
  bool overwriteDirection_;

  /*! New important flag: deltaComputation_
   *
   * Initial approach -> please read to end to get adapted approach
   *
   * if set to true, GetTensor will return
   *      [deltaY/deltaX] + rev_mat_factor
   *    = [(YcurrentItVEC_ - YpreviousItVEC_)/(XcurrentItVEC_ - XpreviousItVEC_)] + rev_mat_factor
   *
   * if set to false, GetTensor will return
   *       [YcurrentItVEC_/XcurrentItVEC_] + rev_mat_factor
   *
   * Notes:
   *  a) the "[ ]" denote a not yet determined way of evaluating this expressions (remember we have a division of vectors)
   *  b) for magnetics, we have to return the inverse of that expression (as we have to compute \nu instead of \mu)
   *
   * Background/Idea behind the two different return values:
   *
   *  1. Assume the following iteration approach (for the example of electrostatics):
   *      D_n+1 = D_n + deltaD
   *      E_n+1 = E_n + deltaE
   *      P_n+1 = P_n + deltaP
   *      rev_mat_factor = eps0
   *      n = time step
   *
   *      div (D_n+1) = div (D_n + deltaD)
   *                  = div (eps0 * E_n + P_n + eps0 * deltaE + deltaP)
   *                  = div ( [eps0 + P_n/E_n]*E_n + [eps_0 + deltaP/deltaE]*deltaE ) = f_n+1
   *      -> div( [eps_0 + deltaP/deltaE]*deltaE ) = f_n+1 - div( [eps0 + P_n/E_n]*E_n )
   *
   *      Here we can already see, that we have different "material relations"
   *      [eps_0 + deltaP/deltaE] vs. [eps0 + P_n/E_n] to be used in the same expression.
   *
   *      Note further, that we compute E not D during the iterations (in fact we use potential but that gives us E).
   *      To compute the correct value of D_n+1, we have to use the "material relation" from the rhs, i.e.
   *        D_n+1 = [eps0 + P_n+1/E_n+1]*E_n+1 = eps0*E_n+1 + P_n+1
   *      Using the relation from the lhs, we would get
   *        D_n+1 = [eps_0 + deltaP/deltaE]*E_n+1 = eps0*E_n+1 + (P_n+1 - P_n) * E_n+1/(E_n+1 - E_n)
   *        -> this does not fit!
   *
   *  2. As we cannot compute div( [eps_0 + deltaP/deltaE]*deltaE ) directly, as we do not know deltaP, deltaE yet,
   *      we use the following iteration scheme:
   *
   *      D_n+1^k+1 = D_n+1^k + deltaD^k+1
   *      E_n+1^k+1 = E_n+1^k + deltaE^k+1
   *      P_n+1^k+1 = P_n+1^k + deltaP^k+1
   *      with D_n+1^0 = D_n etc.
   *      n = time step; k = iteration
   *
   *      -> div( [eps_0 + deltaP^k/deltaE^k]*deltaE^k+1 ) = f_n+1 - div( [eps0 + P_n+1^k/E_n+1^k]*E_n+1^k )
   *
   *      D_n+1^k+1 = [eps0 + P_n+1^k+1/E_n+1^k+1]*E_n+1^k+1
   *
   *  3. Sidenote:
   *      If we use a delta-stepping for non-linear material curves (e.g. B-H-curve in magnetics) we do NOT need
   *      such a flag. This comes from the different approach that reads
   *        D_n+1^k+1 = D_n+1^k + deltaD^k+1
   *              = eps(E_n+1^k) [ E_n+1^k + deltaE^k+1 ]
   *
   *      As one hopefully can see, both terms E_n+1^k and deltaE^k+1 use the same material relation eps(E_n+1^k).
   *
   *      We could use a similar approach for the hysteresis case, too, by setting:
   *      div (D_n+1) = div (D_n + deltaD)
   *                  = div (eps0 * E_n + P_n + eps0 * deltaE + deltaP)
   *                  = div ( [eps0 * (P_n + deltaP)/(E_n + deltaE)] * (E_n + deltaE) )
   *                  = div ( [eps0 * P_n+1/E_n+1] * E_n+1 )
   *
   *      By doing so, we would always calculate the stepping in a fixed-point way as we would
   *      compute the slope deltaP/deltaE using deltaP = P_n+1 - 0 and deltaE = E_n+1 - 0.
   *
   *      -> Keep this alternative in mind in case the other approach does not work.
   *
   *  4. BIG ISSUE:
   *      Remanence points!
   *        if E -> 0 and P != 0, deltaP/deltaE might still work (as deltaP hopefully tends to 0 when deltaE does)
   *        however, P/E is a problem
   *          a) we could simply return eps0 in that case, but that would lead to a wrong rhs of 0 instead of P_remanence
   *          b) add P_remanence to rhs (as it was tried during a rhs-update strategy) and set material tensor for
   *              the rhs-modification to eps0
   *
   *     Adapted scheme:
   *
   *        div( [eps_0 + deltaP^k/deltaE^k]*deltaE^k+1 ) = f_n+1 - div( eps0*E_n+1^k ) - div( P_n+1^k )
   *
   *        in that sense, deltaComputation_ has the following purpose now:
   *
   * if deltaComputation_ is set to true, GetTensor will return
   *      [deltaY/deltaX] + rev_mat_factor
   *    = [(YcurrentItVEC_ - YpreviousItVEC_)/(XcurrentItVEC_ - XpreviousItVEC_)] + rev_mat_factor
   *
   * if deltaComputation_ is set to false, GetTensor will return
   *       rev_mat_factor
   *
   */
  bool deltaComputation_;
  /*
   * use flag to switch between the delta computation using
   * a) the previousIterationValues > useNextLastTS_ = false
   * b) the nextLastTimestepValues -> useNextLastTS_ = true
   */
  bool useNextToLastTS_;

  /*
   * for magnetics, we need the inverse deltaMatrix, as we not need \mu but \nu
   * -> indicate that by a flag
   */
  bool compute_inverse_;
};


}
#endif
