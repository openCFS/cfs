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

  void GetVector(Vector<Double>& coefScalar, const LocPointMapped& lpm );

  //! Return real-valued tensor at integration point
  void GetTensor(Matrix<Double>& tensor, const LocPointMapped& lpm );

  //!
  void SetPreviousHystVals();

  //! Create for the vector case deltaMat from dX and dY
  void CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& deltaMat, UInt version, UInt idxElem);

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

protected:
  
  //! compute X and Y value
  void ComputeXY( const LocPointMapped& lpm, Double& X, Double& Y);

  //! compute X and Y vector for vector model
  void ComputeXY_vec( const LocPointMapped& lpm, Vector<Double>& X, Vector<Double>& Y, bool overwrite);

  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! hysteresis object
  Hysteresis * hyst_;

  //! previous Xval in hysteresis
  //! -> this is the value of the last timestep (if any)
  Vector<Double> Xprevious_;

  //! previous iteration value of hystersis
  //! -> this is the value of the last iteration (if any)
  //! -> needed to check if input has changed during iteration steps
  //! -> if not -> avoid recomputation and return YpreviousIt_
  Vector<Double> XpreviousIt_;

  //! previous Yal in hysteresis
  Vector<Double> Yprevious_;

  Vector<Double> YpreviousIt_;

  //! delta material tensor
  //! new: extend to array of matrices so that we can store the old state of each element
  Matrix<Double>* matDeltaTensor_;

  //! for vector version
  Vector<Double>* XpreviousVEC_;
  Vector<Double>* XpreviousItVEC_;
  Vector<Double>* YpreviousVEC_;
  Vector<Double>* YpreviousItVEC_;

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
  Double rev_mat_fac_;

  UInt evalVersion_;
  UInt numRows_;
  Double tol_;

};


}
#endif
