// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
*		\file		SurfaceNormalFluxDensityOperator.hh
*		\brief	Collection of Operators acting on the surfaces of volume elements
*           These operators assume always a local point defined on a surface
*           and a volume BaseFE. If Not, they will throw an exception
*
*		\date		Feb 7, 2015
*		\author	Kirill Shaposhnikov (kirill.shaposhnikov@tuwien.ac.at)
*/
//================================================================================================
#ifndef SOURCE_DIRECTORY__SOURCE_FORMS_OPERATORS_SURFACENORMALFLUXDENSITYOPERATOR_HH_
#define SOURCE_DIRECTORY__SOURCE_FORMS_OPERATORS_SURFACENORMALFLUXDENSITYOPERATOR_HH_

#include "BaseBOperator.hh"
#include "GradientOperator.hh"

namespace CoupledField
{

template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceNormalFluxDensityOperator : public BaseBOperator
{
public:

  // ------------------
  //  STATIC CONSTANTS
  // ------------------
  //@{
  //! \name Static constants

  //! Order of differentiation
  static const UInt ORDER_DIFF = 1;

  //! Number of components of the problem (scalar, vector)
  static const UInt DIM_DOF = D_DOF;

  //! Dimension of the underlying domain / space
  static const UInt DIM_SPACE = D;

  //! Dimension of the finite element
  static const UInt DIM_ELEM = D;

  //! Dimension of the related material
  static const UInt DIM_D_MAT = D_DOF;
  //@}

  SurfaceNormalFluxDensityOperator(const std::string& subType)
  {
   if (subType == "2.5d")
   {
     gradOp_ = new GradientOperator2p5D<FeH1, D, D_DOF, TYPE>();
   }
   else
   {
     gradOp_ = new GradientOperator<FeH1, D, D_DOF, TYPE>();
   }

   this->name_ = "surfNormFluxOp";
   gradOp_->SetOperator2SurfOperator();
  }

  SurfaceNormalFluxDensityOperator(const std::string& subType, PtrCoefFct baseOpCoef)
  {
    if (subType == "2.5d")
    {
      gradOp_ = new ScaledGradientOperator2p5D<FeH1, D, TYPE>();
    }
    else
    {
      gradOp_ = new ScaledGradientOperator<FeH1, D, TYPE>();
    }

    this->name_ = "surfNormFluxOp";
    gradOp_->SetOperator2SurfOperator();
    gradOp_->SetCoefFunction(baseOpCoef);
  }

  //! Copy constructor
  SurfaceNormalFluxDensityOperator(const SurfaceNormalFluxDensityOperator & other)
     : BaseBOperator(other){
    this->name_ = "surfNormFluxOp";
    this->gradOp_ = other.gradOp_->Clone();
    this->gradOp_->SetOperator2SurfOperator();
  }

  //! \copydoc BaseBOperator::Clone()
  virtual SurfaceNormalFluxDensityOperator * Clone(){
    return new SurfaceNormalFluxDensityOperator(*this);
  }

  virtual ~SurfaceNormalFluxDensityOperator() { return; }

  virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  virtual void CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  virtual void CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  // ===============
  //  QUERY METHODS
  // ===============
  //@{ \name Query Methods
  //! \copydoc BaseBOperator::GetDiffOrder
  virtual UInt GetDiffOrder() const
  {
    return ORDER_DIFF;
  }

  //! \copydoc BaseBOperator::GetDimDof()
  virtual UInt GetDimDof() const
  {
    return DIM_DOF;
  }

  //! \copydoc BaseBOperator::GetDimSpace()
  virtual UInt GetDimSpace() const
  {
    return DIM_SPACE;
  }

  //! \copydoc BaseBOperator::GetDimElem()
  virtual UInt GetDimElem() const
  {
    return DIM_ELEM;
  }

  //! \copydoc BaseBOperator::GetDimDMat()
  virtual UInt GetDimDMat() const
  {
    return gradOp_->GetDimDMat();
  }
  //@}

protected:

  //! computes the normal operator
  virtual void ComputeNormalOperator(Matrix<Double>& normalOp, const LocPointMapped& lp)
  {
    UInt dimD = gradOp_->GetDimDMat();
    normalOp.Resize(dimD, DIM_DOF);
    normalOp.Init();

    if (DIM_SPACE == 2)
    {
      //2D and 2.5D cases
      normalOp[0][0] = lp.normal[0];
      normalOp[1][0] = lp.normal[1];
    }
    else if (DIM_SPACE == 3)
    {
      // 3D case
      normalOp[0][0] = lp.normal[0];
      normalOp[1][0] = lp.normal[1];
      normalOp[2][0] = lp.normal[2];
    }
    else
    {
      EXCEPTION("Unknown formulation for 'DIM_SPACE' = " << DIM_SPACE);
    }
  }

  //! computes the normal operator
  virtual void ComputeNormalOperator(Matrix<Complex>& normalOp, const LocPointMapped& lp)
  {
    UInt dimD = gradOp_->GetDimDMat();
    normalOp.Resize(dimD, DIM_DOF);
    normalOp.Init();

    Vector<Complex> coeffs;
    PtrCoefFct gradOpCoef = gradOp_->GetCoefFunction();
    if (gradOpCoef == NULL)
    {
      coeffs.Resize(DIM_SPACE, Complex(1.0, 0.0));
    }
    else
    {
      gradOpCoef->GetVector(coeffs, lp);
    }

    if (DIM_SPACE == 2)
    {
      //2D and 2.5D cases
      normalOp[0][0] = lp.normal[0]*coeffs[0];
      normalOp[1][0] = lp.normal[1]*coeffs[1];
    }
    else if (DIM_SPACE == 3)
    {
      // 3D case
      normalOp[0][0] = lp.normal[0]*coeffs[0];
      normalOp[1][0] = lp.normal[1]*coeffs[1];
      normalOp[2][0] = lp.normal[2]*coeffs[2];
    }
    else
    {
      EXCEPTION("Unknown formulation for 'DIM_SPACE' = " << DIM_SPACE);
    }
  }

  //! gradient operator
  BaseBOperator* gradOp_;
};

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalFluxDensityOperator<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Double> tmpMat;
  this->CalcOpMatTransposed(tmpMat, lp, ptFe);
  bMat = Transpose(tmpMat);
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalFluxDensityOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(this->coef_);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize(numFncs*D_DOF, D_DOF);
  bMat.InitValue(0.0);

  //obtain external field
  Matrix<Double> matTensor;
  this->coef_->GetTensor(matTensor, lp);

  // Get B-Operator of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Double> gradMat;
  FE *fe = (static_cast<FE*>(ptFe));
  gradOp_->CalcOpMatTransposed(gradMat, lp, fe);

  // get normal opertaor
  Matrix<Double> normalOp;
  ComputeNormalOperator(normalOp, lp);

  bMat = gradMat*matTensor*normalOp;
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalFluxDensityOperator<FE, D, D_DOF, TYPE>::CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  Matrix<Complex> tmpMat;
  this->CalcOpMatTransposed(tmpMat, lp, ptFe);
  bMat = Transpose(tmpMat);
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalFluxDensityOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);
  assert(this->coef_);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize(numFncs*D_DOF, D_DOF);
  bMat.InitValue(0.0);

  //obtain external field
  Matrix<Complex> matTensor;
  this->coef_->GetTensor(matTensor, lp);

  // Get B-Operator of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Complex> gradMat;
  FE *fe = (static_cast<FE*>(ptFe));
  gradOp_->CalcOpMatTransposed(gradMat, lp, fe);

  // get normal opertaor
  Matrix<Complex> normalOp;
  ComputeNormalOperator(normalOp, lp);

  bMat = gradMat*matTensor*normalOp;
}

template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceNormalPiezoStrainOperator : public SurfaceNormalFluxDensityOperator<FE, D, D_DOF, TYPE>
{

public:

  // ------------------
  //  STATIC CONSTANTS
  // ------------------
  //@{
  //! \name Static constants

  //! Order of differentiation
  static const UInt ORDER_DIFF = 1;

  //! Number of components of the problem (scalar, vector)
  static const UInt DIM_DOF = D_DOF;

  //! Dimension of the underlying domain / space
  static const UInt DIM_SPACE = D;

  //! Dimension of the finite element
  static const UInt DIM_ELEM = D;

  //! Dimension of the related material
  static const UInt DIM_D_MAT = D_DOF;
  //@}

  SurfaceNormalPiezoStrainOperator(const std::string& subType)
    : SurfaceNormalFluxDensityOperator<FE, D, D_DOF, TYPE>(subType), dimNorm_(0)
  {
    this->name_ = "surfNormPiezoStrainOp";

    if (subType == "axi")
      dimNorm_ = 4;
    else if (subType == "planeStrain" || subType == "planeStress")
      dimNorm_ = 3;
    else if (subType == "3d" || subType == "2.5d")
      dimNorm_ = 6;
    else
      EXCEPTION("Subtype '" << subType << "' in SurfaceNormalStressOperator");
  }

  SurfaceNormalPiezoStrainOperator(const std::string& subType, PtrCoefFct baseOpCoef)
    : SurfaceNormalFluxDensityOperator<FE, D, D_DOF, TYPE>(subType, baseOpCoef), dimNorm_(0)
  {
    this->name_ = "surfNormPiezoStrainOp";

    if (subType == "axi")
      dimNorm_ = 4;
    else if (subType == "planeStrain" || subType == "planeStress")
      dimNorm_ = 3;
    else if (subType == "3d" || subType == "2.5d")
      dimNorm_ = 6;
    else
      EXCEPTION("Subtype '" << subType << "' in SurfaceNormalStressOperator");
  }

  //! Copy constructor
  SurfaceNormalPiezoStrainOperator(const SurfaceNormalPiezoStrainOperator & other)
     : SurfaceNormalFluxDensityOperator<FE, D, D_DOF, TYPE>(other) {
    this->name_ = other.name_;
  }

  //! \copydoc BaseBOperator::Clone()
  virtual SurfaceNormalPiezoStrainOperator * Clone(){
    return new SurfaceNormalPiezoStrainOperator(*this);
  }

  virtual ~SurfaceNormalPiezoStrainOperator() { return; }

  virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

  virtual void CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

protected:

  //! computes the normal operator
  virtual void ComputeNormalOperator(Matrix<Double>& normalOp, const LocPointMapped& lp)
  {
    //normalOp.Resize(dimNorm_, DIM_SPACE);
    normalOp.Resize(dimNorm_, this->gradOp_->GetDimDMat());
    normalOp.Init();

    if (dimNorm_ == 3)
    {
      //2D case
      normalOp[0][0] = lp.normal[0];
      normalOp[1][1] = lp.normal[1];
      normalOp[2][0] = lp.normal[1];
      normalOp[2][1] = lp.normal[0];
    }
    else if (dimNorm_ == 4)
    {
      //2D axisymmetric case
      normalOp[0][0] = lp.normal[0];
      normalOp[1][1] = lp.normal[1];
      normalOp[2][0] = lp.normal[1];
      normalOp[2][1] = lp.normal[0];
    }
    else if (dimNorm_ == 6)
    {
      if (DIM_SPACE == 3)
      {
        // 3D case
        normalOp[0][0] = lp.normal[0];
        normalOp[1][1] = lp.normal[1];
        normalOp[2][2] = lp.normal[2];
        normalOp[3][1] = lp.normal[2];
        normalOp[3][2] = lp.normal[1];
        normalOp[4][0] = lp.normal[2];
        normalOp[4][2] = lp.normal[0];
        normalOp[5][0] = lp.normal[1];
        normalOp[5][1] = lp.normal[0];
      }
      else if (DIM_SPACE == 2)
      {
        // 2.5D case
        normalOp[0][0] = lp.normal[0];
        normalOp[1][1] = lp.normal[1];
        normalOp[3][2] = lp.normal[1];
        normalOp[4][2] = lp.normal[0];
        normalOp[5][0] = lp.normal[1];
        normalOp[5][1] = lp.normal[0];
      }
      else
      {
        EXCEPTION("Unknown formulation: DIM_DOF = " << DIM_DOF << " and DIM_SPACE = " << DIM_SPACE);
      }
    }
  }

  //! computes the normal operator
  virtual void ComputeNormalOperator(Matrix<Complex>& normalOp, const LocPointMapped& lp)
  {
    //normalOp.Resize(dimNorm_, DIM_SPACE);
    normalOp.Resize(dimNorm_, this->gradOp_->GetDimDMat());
    normalOp.Init();

    Vector<Complex> coeffs;
    PtrCoefFct strainOpCoef = this->gradOp_->GetCoefFunction();
    if (strainOpCoef == NULL)
    {
      coeffs.Resize(3, Complex(1.0, 0.0));
    }
    else
    {
      strainOpCoef->GetVector(coeffs, lp);
    }

    if (dimNorm_ == 3)
    {
      //2D case
      normalOp[0][0] = lp.normal[0]*coeffs[0];
      normalOp[1][1] = lp.normal[1]*coeffs[1];
      normalOp[2][0] = lp.normal[1]*coeffs[1];
      normalOp[2][1] = lp.normal[0]*coeffs[0];
    }
    else if (dimNorm_ == 4)
    {
      //2D axisymmetric case
      normalOp[0][0] = lp.normal[0]*coeffs[0];
      normalOp[1][1] = lp.normal[1]*coeffs[1];
      normalOp[2][0] = lp.normal[1]*coeffs[1];
      normalOp[2][1] = lp.normal[0]*coeffs[0];
    }
    else if (dimNorm_ == 6)
    {
      if (DIM_SPACE == 3)
      {
        // 3D case
        normalOp[0][0] = lp.normal[0]*coeffs[0];
        normalOp[1][1] = lp.normal[1]*coeffs[1];
        normalOp[2][2] = lp.normal[2]*coeffs[2];
        normalOp[3][1] = lp.normal[2]*coeffs[2];
        normalOp[3][2] = lp.normal[1]*coeffs[1];
        normalOp[4][0] = lp.normal[2]*coeffs[2];
        normalOp[4][2] = lp.normal[0]*coeffs[0];
        normalOp[5][0] = lp.normal[1]*coeffs[1];
        normalOp[5][1] = lp.normal[0]*coeffs[0];
      }
      else if (DIM_SPACE == 2)
      {
        // 2.5D case
        normalOp[0][0] = lp.normal[0]*coeffs[0];
        normalOp[1][1] = lp.normal[1]*coeffs[1];
        normalOp[3][2] = lp.normal[1]*coeffs[1];
        normalOp[4][2] = lp.normal[0]*coeffs[0];
        normalOp[5][0] = lp.normal[1]*coeffs[1];
        normalOp[5][1] = lp.normal[0]*coeffs[0];
      }
      else
      {
        EXCEPTION("Unknown formulation: DIM_DOF = " << DIM_DOF << " and DIM_SPACE = " << DIM_SPACE);
      }
    }
  }

  //! Dimension for normal operator
  UInt dimNorm_;

};

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalPiezoStrainOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize(numFncs*D_DOF, this->gradOp_->GetDimDMat());
  bMat.InitValue(0.0);

  //obtain external field
  Matrix<Double> matTensor, matTensTrans;
  this->coef_->GetTensor(matTensor, lp);
  matTensor.Transpose(matTensTrans);

  // Get B-Operator of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Double> gradMat;
  FE *fe = (static_cast<FE*>(ptFe));
  this->gradOp_->CalcOpMatTransposed(gradMat, lp, fe);

  // get normal opertaor
  Matrix<Double> normalOp;
  ComputeNormalOperator(normalOp, lp);

  bMat = gradMat*matTensTrans*normalOp;
}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalPiezoStrainOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize(numFncs*D_DOF, this->gradOp_->GetDimDMat());
  bMat.InitValue(0.0);

  //obtain external field
  Matrix<Complex> matTensor, matTensTrans;
  this->coef_->GetTensor(matTensor, lp);
  matTensor.Transpose(matTensTrans);

  // Get B-Operator of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Complex> gradMat;
  FE *fe = (static_cast<FE*>(ptFe));
  this->gradOp_->CalcOpMatTransposed(gradMat, lp, fe);

  // get normal opertaor
  Matrix<Complex> normalOp;
  ComputeNormalOperator(normalOp, lp);

  bMat = gradMat*matTensTrans*normalOp;
}

}
#endif /* SOURCE_DIRECTORY__SOURCE_FORMS_OPERATORS_SURFACENORMALFLUXDENSITYOPERATOR_HH_ */
