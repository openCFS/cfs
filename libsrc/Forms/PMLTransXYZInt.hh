#ifndef FILE_PMLTRANSXYZINT
#define FILE_PMLTRANSXYZINT

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  of element matrix for transient PML formulation

class PMLTransXYZInt : public BaseForm
{
public:

  //! Constructor
  PMLTransXYZInt(std::string type, Double factor, std::string dampingTypePML, Double damp, 
		 const UInt nrDofsPerNode=1, Boolean axi=FALSE);

  /// 
  virtual ~PMLTransXYZInt();

  //! Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! set min/max of x,y,z coordinates form where PML starts and ends
  void SetPosPML(Matrix<Double> & inner, Matrix<Double> & outer);

  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

  
  virtual void SetActElemSol(Matrix<Double>& disp){};


private:

  //! Calculation of stiffmess matrix
  void CalcElementMatrixStiff(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! Calculation of mass matrix
  void CalcElementMatrixMass(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);

  //! calculates position and values
  void ComputeFactorPML(Vector<Double>& factorsPML, Matrix<Double> & ptCoord);

  //! calculates the damping factor
  Double ComputeDampingFactor(Vector<Double>& pos, Directions dir);

  //! degrees of freedom per node
  UInt nrDofsPerNode_;   

  //! type of bilinear form
  std::string formsType_;

  //! multiplicative factor for forms
  Double formsFactor_;

  //! type of PML damping
  std::string dampingTypePML_;

  //! damping factor 
  Double dampingFactor_;

  //!layer thickness
  Matrix<Double> layerThickness_;

  //! coordinates for inner box at which PML starts
  Double minX_, maxX_, minY_, maxY_, minZ_, maxZ_;
};


}

#endif // FILE_LAPLACEINT
