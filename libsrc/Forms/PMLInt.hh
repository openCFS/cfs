#ifndef FILE_PMLINT
#define FILE_PMLINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness/mass matrix for PML formulation

class PMLInt : public BaseForm
{
public:

  //! Constructor
  PMLInt(std::string type, Double factor, std::string dampingTypePML, Double damp, 
	 Double layerThick, Boolean axi=FALSE);

  /// 
  virtual ~PMLInt();

  //! Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! set min/max of x,y,z coordinates form where PML starts
  void SetPosPML(Matrix<Double> & pos);

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

  
  virtual void SetActElemSol(Matrix<Double>& disp){};


private:

  //! Calculation of stiffmess matrix
  void CalcElementMatrixStiff(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! Calculation of mass matrix
  void CalcElementMatrixMass(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! calculates position and values
  void ComputeFactorPML(Vector<Double>& factorsPML, Matrix<Double> & ptCoord);

  //! type of bilinear form
  std::string formsType_;

  //! multiplicative factor for forms
  Double formsFactor_;

  //! type of PML damping
  std::string dampingTypePML_;

  //! damping factor 
  Double dampingFactor_;

  //!layer thickness
  Double layerThickness_;

  //! coordinates for inner box at which PML starts
  Double minX_, maxX_, minY_, maxY_, minZ_, maxZ_;
};


}

#endif // FILE_LAPLACEINT
