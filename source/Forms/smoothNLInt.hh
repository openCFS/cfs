#ifndef FILE_SMOOTHNLINT
#define FILE_SMOOTHNLINT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <Materials/baseMaterial.hh>

namespace CoupledField
{
  
  
  /// base class for calculation of linear elasticity
class SmoothNLInt : public BDBInt
{
public:

  //! Constructor
  SmoothNLInt(BaseMaterial* matData, 
              Matrix<Double> & elastWeights,
              Matrix<Double> & couplingNodes,
              Double AMax,
              Double AMin,
              Double characteristicLength,
              Double expo,
              SubTensorType type = FULL,
              bool coordUpdate = false );
  
  //! Destructor
  virtual ~SmoothNLInt();

  //! Compute element matrix associated to ADB form
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
  void CalcMinMaxStrain( EntityIterator& ent1, 
                         EntityIterator& ent2 );
  void ResetMinMaxStrain(){
    elastFactorMax_=0.0;
    elastFactorMin_=2.0;
  };

protected:    

  Matrix<Double> & elastWeights_;
  Matrix<Double> & couplingNodes_;

  Double AMax_, AMin_;
  Double characteristicLength_;
  Double exponent_;
  Double elastFactorMax_, elastFactorMin_, eMean_;

  Double ComputeElastFactor( Double & e11, Double & e22, Double & e12 );

  //! returns B - matrix for BDB
  void calcBMat(Matrix<Double> & bMat, 
                UInt ip, 
                Matrix<Double> & ptCoord);

  //! returns D - matrix for BDB
  void calcDMat( Matrix<Double> &dMat );


  //! calculate the data-matrix for 2D plain-strain
  void calcDMat(Matrix<Double> & dMat, 
                UInt ip, 
                Matrix<Double> & ptCoord);

  //! returns dimension of D matrix
  UInt getDimD() {
    return dimD_; 
  };

  //! Query material type for \f$D\f$ tensor
  MaterialType getDMaterialType() { return MECH_STIFFNESS_TENSOR; }

  
  //! returns nr. of degrees of freedom
  UInt getNrDofs() {
    return nrDofs_;
  };

  private:

  //dimension of Dmatrix
  UInt dimD_;

  //! number of degrees 
  UInt nrDofs_;

  std::fstream file;

};
  
} //end namespace

#endif // FILE_SMOOTHINT
