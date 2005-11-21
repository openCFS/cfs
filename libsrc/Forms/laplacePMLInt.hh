#ifndef FILE_LAPLACEPMLINT
#define FILE_LAPLACEPMLINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
  /// for PML formulation

class LaplacePMLInt : public BaseForm
{
public:
  //! Constructor
  //  LaplacePMLInt(BaseFE * aptelem, Double damp, DoublelaplVal, Boolean axi=FALSE);

  //! Constructor
  LaplacePMLInt(Double damp, Double layerThick, Boolean axi=FALSE);

  /// 
  virtual ~LaplacePMLInt();

  //! Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! set min/max of x,y,z coordinates form where PML starts
  virtual void SetPosPML(Matrix<Double> & pos) {
    minMax_ = pos;
  }

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

  
  virtual void SetActElemSol(Matrix<Double>& disp){};


private:

  //! calculates position and values
  void ComputeFactorPML(Vector<Double>& factorsPML, Matrix<Double> & ptCoord);

  //! damping factor 
  Double damping_;

  //!layer thickness
  Double layerThickness_;

  //! coordinates for inner box at wich PML starts
  Matrix<Double> minMax_;
};


}

#endif // FILE_LAPLACEINT
