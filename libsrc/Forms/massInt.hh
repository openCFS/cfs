#ifndef FILE_MASSINT_1
#define FILE_MASSINT_1

#include "baseForm.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix
class MassInt : public BaseForm
{
public:
  /// Constructor
  MassInt(BaseFE * aptelemt, Double aDensity, Boolean axi=FALSE);

  /// Constructor
  MassInt(const Double aDensity, const Integer nrDofsPerNode=1, Boolean axi=FALSE);

  /// Constructor
  MassInt(const Double aDensity, const Integer nrDofsPerNode, Integer dofzero, Boolean axi);

  /// Destructor
  virtual ~MassInt();

  /// Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMa);

  ///
  virtual void SetFracDamping() 
  {isFracDamping_ = TRUE;};
  
  //
  virtual void SetFactor(Double factor) 
  {density_ = factor;};
   
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected:
  //! generates a multi-dof-matrix with similar entries for all dofs
  virtual void MassMultiDof(Matrix<Double>& massMultDof, const Matrix<Double>& massMatSingleDof,  
			    const Integer nrDofs);

  virtual void MassMultiDofZero(Matrix<Double>& massMultDofZero, 
				const Matrix<Double>& massMatSingleDof);
  

private:
  // multiplicative value for mass integrator
  Double density_;
  Integer nrDofsPerNode_;
  Integer dofzero_;
  
};





}

#endif // FILE_MASSINT
