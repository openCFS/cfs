// ================================================================================================
/*!
 *       \file     SurfaceNormalStressOperator.hh
 *       \brief    Collection of Operators acting on the surfaces of volume elements
 *                 These operators assume always a local point defined on a surface
 *                 and a volume BaseFE. If Not, they will throw an exception
 *
 *       \date     September, 2013
 *       \author   Manfred Kaltenbacher (manfred.kaltenbacher@tuwien.ac.at)
 */
//================================================================================================

#ifndef SURFACENORMALSTRESSOPERATORS_HH
#define SURFACENORMALSTRESSOPERATORS_HH

#include "BaseBOperator.hh"
#include "StrainOperator.hh"

namespace CoupledField{

template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
class SurfaceNormalStressOperator : public BaseBOperator{

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

   SurfaceNormalStressOperator(const std::string& subType, bool icModes){
     if( subType == "axi" ) {
       strainOp_ = new StrainOperatorAxi<FeH1,TYPE>(icModes);
     } else if( subType == "planeStrain" ) {
       strainOp_ = new StrainOperator2D<FeH1,TYPE>(icModes);
     } else if( subType == "planeStress" ) {
       strainOp_ = new StrainOperator2D<FeH1,TYPE>(icModes);
     } else if( subType == "3d") {
       strainOp_ = new StrainOperator3D<FeH1,TYPE>(icModes);
     } else if ( subType == "2.5d") {
       strainOp_ = new StrainOperator2p5D<FeH1,TYPE>(icModes);
     }else {
       EXCEPTION( "Subtype '" << subType << "' in SurfaceNormalStressOperator" );
     }

     this->name_ = "surfNormStressOp";
     strainOp_->SetOperator2SurfOperator();
    }

   //! Copy constructor
   SurfaceNormalStressOperator(const SurfaceNormalStressOperator & other)
      : BaseBOperator(other){
     this->name_ = other.name_;
     this->strainOp_ = other.strainOp_->Clone();
     strainOp_->SetOperator2SurfOperator();
   }

   //! \copydoc BaseBOperator::Clone()
   virtual SurfaceNormalStressOperator * Clone(){
     return new SurfaceNormalStressOperator(*this);
   }

    SurfaceNormalStressOperator(const std::string& subType, PtrCoefFct baseOpCoef, bool icModes){
     if( subType == "axi" ) {
       //TODO: there's no scaled strain operator for 'axi' subtype
       strainOp_ = new StrainOperatorAxi<FeH1,TYPE>(icModes);
     } else if( subType == "planeStrain" ) {
       strainOp_ = new ScaledStrainOperator2D<FeH1,TYPE>(icModes);
     } else if( subType == "planeStress" ) {
       strainOp_ = new ScaledStrainOperator2D<FeH1,TYPE>(icModes);
     } else if( subType == "3d") {
       strainOp_ = new ScaledStrainOperator3D<FeH1,TYPE>(icModes);
     } else if ( subType == "2.5d") {
       strainOp_ = new ScaledStrainOperator2p5D<FeH1,TYPE>(icModes);
     }else {
       EXCEPTION( "Subtype '" << subType << "' in SurfaceNormalStressOperator" );
     }

     this->name_ = "surfNormStressOp";
     strainOp_->SetOperator2SurfOperator();
     strainOp_->SetCoefFunction(baseOpCoef);
    }

    virtual ~SurfaceNormalStressOperator(){
      return;
    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe );

    virtual void CalcOpMat(Matrix<Complex> & bMat, const LocPointMapped& lp, BaseFE* ptFe);

    virtual void CalcOpMatTransposed(Matrix<Complex> & bMat, const LocPointMapped& lp, BaseFE* ptFe);

    // ===============
    //  QUERY METHODS
    // ===============
    //@{ \name Query Methods
    //! \copydoc BaseBOperator::GetDiffOrder
    virtual UInt GetDiffOrder() const {
      return ORDER_DIFF;
    }

    //! \copydoc BaseBOperator::GetDimDof()
    virtual UInt GetDimDof() const {
      return DIM_DOF;
    }

    //! \copydoc BaseBOperator::GetDimSpace()
    virtual UInt GetDimSpace() const {
      return DIM_SPACE;
    }

    //! \copydoc BaseBOperator::GetDimElem()
    virtual UInt GetDimElem() const {
      return DIM_ELEM;
    }

    //! \copydoc BaseBOperator::GetDimDMat()
    virtual UInt GetDimDMat() const {
      return strainOp_->GetDimDMat();
    }
    //@}

  protected:
    //! computes the normal operator
    virtual void ComputeNormalOperator( Matrix<Double>& normalOp,
                                       const LocPointMapped& lp ) {

      UInt dimD = strainOp_->GetDimDMat();
      // kirill: In the case of 2.5D we have 6x3 matrix. This means DIM_SPACE != DIM_DOF
      // In order take account of such a formulation, Resize is done with DIM_DOF instead of DIM_SPASE
      normalOp.Resize(dimD, DIM_DOF);
      normalOp.Init();

      if ( dimD == 3 ){
        //2D case
        normalOp[0][0] = lp.normal[0];
        normalOp[1][1] = lp.normal[1];
        normalOp[2][0] = lp.normal[1];
        normalOp[2][1] = lp.normal[0];
      }
      else if ( dimD == 4) {
        //2D axisymmetric case
        normalOp[0][0] = lp.normal[0];
        normalOp[1][1] = lp.normal[1];
        normalOp[2][0] = lp.normal[1];
        normalOp[2][1] = lp.normal[0];
      }
      else if ( dimD == 6 ) {
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
      UInt dimD = strainOp_->GetDimDMat();
      // kirill: In the case of 2.5D we have 6x3 matrix. This means DIM_SPACE != DIM_DOF
      // In order take account of such a formulation, Resize is done with DIM_DOF instead of DIM_SPASE
      normalOp.Resize(dimD, DIM_DOF);
      normalOp.Init();

      Vector<Complex> coeffs;
      PtrCoefFct strainOpCoef = strainOp_->GetCoefFunction();
      if (strainOpCoef == NULL)
      {
        coeffs.Resize(3, Complex(1.0, 0.0));
      }
      else
      {
        strainOpCoef->GetVector(coeffs, lp);
      }

      if ( dimD == 3 ){
        //2D case
        normalOp[0][0] = lp.normal[0]*coeffs[0];
        normalOp[1][1] = lp.normal[1]*coeffs[1];
        normalOp[2][0] = lp.normal[1]*coeffs[1];
        normalOp[2][1] = lp.normal[0]*coeffs[0];
      }
      else if ( dimD == 4) {
        //2D axisymmetric case
        normalOp[0][0] = lp.normal[0]*coeffs[0];
        normalOp[1][1] = lp.normal[1]*coeffs[1];
        normalOp[2][0] = lp.normal[1]*coeffs[1];
        normalOp[2][1] = lp.normal[0]*coeffs[0];
      }
      else if ( dimD == 6 ) {
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

    //! strain operator
    BaseBOperator* strainOp_;

};

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalStressOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                              const LocPointMapped& lp,
                                              BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  Matrix<Double> tmpMat;
  this->CalcOpMatTransposed(tmpMat,lp,ptFe);
  bMat = Transpose(tmpMat);

}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalStressOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                                       const LocPointMapped& lp,
                                                                       BaseFE* ptFe ){
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize( numFncs*D_DOF, D_DOF );
  bMat.InitValue(0.0);

  //obtain external field
  Matrix<Double> matTensor;
  this->coef_->GetTensor(matTensor,lp);

  // Get B-Operator of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Double> strainMat;
  FE *fe = (static_cast<FE*>(ptFe));
  strainOp_->CalcOpMatTransposed(strainMat, lp, fe);

  bMat = strainMat * matTensor;

  // get normal opertaor
  Matrix<Double> normalOp;
  ComputeNormalOperator( normalOp, lp );

  bMat *= normalOp;

}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalStressOperator<FE,D,D_DOF,TYPE>::CalcOpMat(Matrix<Complex> & bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  Matrix<Complex> tmpMat;
  this->CalcOpMatTransposed(tmpMat, lp, ptFe);
  bMat = Transpose(tmpMat);
}

template<class FE,  UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalStressOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(Matrix<Complex> & bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize(numFncs*D_DOF, D_DOF);
  bMat.InitValue(0.0);

  //obtain external field
  Matrix<Complex> matTensor;
  this->coef_->GetTensor(matTensor, lp);

  // Get B-Operator of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Complex> strainMat;
  FE *fe = (static_cast<FE*>(ptFe));
  strainOp_->CalcOpMatTransposed(strainMat, lp, fe);

  // get normal opertaor
  Matrix<Complex> normalOp;
  ComputeNormalOperator(normalOp, lp);

  bMat = strainMat*matTensor*normalOp;
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
class SurfaceNormalPiezoFluxOperator : public SurfaceNormalStressOperator<FE, D, D_DOF, TYPE>
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

  SurfaceNormalPiezoFluxOperator(const std::string& subType)
    : SurfaceNormalStressOperator<FE, D, D_DOF, TYPE>(subType, false)
  {
    this->name_ = "surfNormPiezoFluxOp";
  }

  SurfaceNormalPiezoFluxOperator(const std::string& subType, PtrCoefFct baseOpCoef)
    : SurfaceNormalStressOperator<FE, D, D_DOF, TYPE>(subType, baseOpCoef, false)
  {
    this->name_ = "surfNormPiezoFluxOp";
  }

  //! copy constructor
  SurfaceNormalPiezoFluxOperator(const SurfaceNormalPiezoFluxOperator & other)
     : SurfaceNormalStressOperator<FE, D, D_DOF, TYPE>(other){
    this->name_ = other.name_;
  }

  //! \copydoc baseboperator::clone()
  virtual SurfaceNormalPiezoFluxOperator * Clone(){
    return new SurfaceNormalPiezoFluxOperator(*this);
  }

  virtual ~SurfaceNormalPiezoFluxOperator() { }

  virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe );

  virtual void CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

protected:

  //! computes the normal operator
  virtual void ComputeNormalOperator(Matrix<Double>& normalOp, const LocPointMapped& lp)
  {
    normalOp.Resize(DIM_DOF, DIM_DOF);
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
    normalOp.Resize(DIM_DOF, DIM_DOF);
    normalOp.Init();

    Vector<Complex> coeffs;
    PtrCoefFct strainOpCoef = this->strainOp_->GetCoefFunction();
    if (strainOpCoef == NULL)
    {
      coeffs.Resize(3, Complex(1.0, 0.0));
    }
    else
    {
      strainOpCoef->GetVector(coeffs, lp);
    }

    if (DIM_SPACE == 2)
    {
      //2D and 2.5D cases
      normalOp[0][0] = lp.normal[0]*coeffs[0];
      normalOp[1][1] = lp.normal[1]*coeffs[1];
    }
    else if (DIM_SPACE == 3)
    {
      // 3D case
      normalOp[0][0] = lp.normal[0]*coeffs[0];
      normalOp[1][1] = lp.normal[1]*coeffs[1];
      normalOp[2][2] = lp.normal[2]*coeffs[2];
    }
    else
    {
      EXCEPTION("Unknown formulation for 'DIM_SPACE' = " << DIM_SPACE);
    }
  }

};

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalPiezoFluxOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize(numFncs*D_DOF, DIM_D_MAT);
  bMat.InitValue(0.0);

  //obtain external field
  Matrix<Double> matTensor;
  this->coef_->GetTensor(matTensor,lp);

  // Get B-Operator of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Double> strainMat;
  FE* fe = (static_cast<FE*>(ptFe));
  this->strainOp_->CalcOpMatTransposed(strainMat, lp, fe);

  bMat = strainMat*matTensor;

  // get normal opertaor
  Matrix<Double> normalOp;
  ComputeNormalOperator( normalOp, lp );

  bMat *= normalOp;
}

template<class FE, UInt D, UInt D_DOF, class TYPE>
void SurfaceNormalPiezoFluxOperator<FE, D, D_DOF, TYPE>::CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
{
  //check if lp is surface and ptFe is volume
  assert(lp.isSurface);
  assert(D == ptFe->shape_.dim);

  UInt numFncs = ptFe->GetNumFncs();
  // Set correct size of matrix B and initialise with zeros
  bMat.Resize(numFncs*D_DOF, DIM_D_MAT);
  bMat.InitValue(0.0);

  //obtain external field
  Matrix<Complex> matTensor;
  this->coef_->GetTensor(matTensor, lp);

  // Get B-Operator of local shape functions with respect to global
  // coords (format: nrNodes x spaceDim)
  Matrix<Complex> strainMat;
  FE* fe = (static_cast<FE*>(ptFe));
  this->strainOp_->CalcOpMatTransposed(strainMat, lp, fe);

  // get normal opertaor
  Matrix<Complex> normalOp;
  ComputeNormalOperator(normalOp, lp);

  bMat = strainMat*matTensor*normalOp;
}

}
#endif
