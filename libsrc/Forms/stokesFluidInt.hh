#ifndef FILE_STOKESFLUIDINT
#define FILE_STOKESFLUIDINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a 3D Stokes operator
class StokesFluidInt : public BaseForm
{
public:
  /// Constructor
  StokesFluidInt(Double density, Double dynamicViscosity);

  /// 
  virtual ~StokesFluidInt();

  /// Calculation of stiffmess matrix
  virtual void CalcElementMatrix(Matrix<Double> & ptCoord, 
                                 Matrix<Double> & elemMat) {
    Error( "StokesFluidInt::CalcElementMatrix() not correctly overwritten!",
             __FILE__, __LINE__);
  };


  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

  
  virtual void SetActElemSol(Matrix<Double>& disp){};

protected: 

  void ResortElementMatrix(Matrix<Double> & sortedElemMat, 
                           const Matrix<Double> & elemMat,
                           const UInt & nrOfNodes, 
                           const UInt & dofPerNode);
  
  /// multiplicative value for laplace integration 
  Double density_;
  Double dynamicViscosity_;
private:
};


}

#endif // FILE_LAPLACEINT
