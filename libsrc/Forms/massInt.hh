#ifndef FILE_MASSINT_1
#define FILE_MASSINT_1

#include "baseForm.hh"

namespace CoupledField
{

  //! Class for calculation  element mass matrix
class MassInt : public BaseForm
{
public:
  // Constructor
  MassInt(BaseFE * aptelemt, Double aDensity, Boolean axi=FALSE);

  // Constructor
  MassInt(const Double aDensity, const Integer nrDofsPerNode=1, Boolean axi=FALSE);

  // Constructor
  MassInt(const Double aDensity, const Integer nrDofsPerNode, Integer dofzero, Boolean axi);

  // Destructor
  virtual ~MassInt();

  //! Calculation of stiffmess matrix
  void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMa);

  //! for fractional damping model, is called in AcousticPDE::DefineIntegrators
  void SetFracDamping() 
  {isFracDamping_ = TRUE;};

  //! for fractional damping model, is called in AcousticPDE::DefineIntegrators
  void UnsetFracDamping() 
  {isFracDamping_ = FALSE;};
  
  //! set additional multiplicative factor of mass matrix
  void SetFactor(Double aFactor) 
  {factor_ = aFactor;};
   
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected:
  //! generates a multi-dof-matrix with similar entries for all dofs
  virtual void MassMultiDof(Matrix<Double>& massMultDof, const Matrix<Double>& massMatSingleDof,  
			    const Integer nrDofs);

  virtual void MassMultiDofZero(Matrix<Double>& massMultDofZero, 
				const Matrix<Double>& massMatSingleDof);
  

private:

  Double density_;          //!< multiplicative value for mass integrator
  Double factor_;           //!< yet another multiplicative value for mass integrator
  Integer nrDofsPerNode_;   //!< degrees of freedom per node
  Integer dofzero_;
  
};





}

#endif // FILE_MASSINT
