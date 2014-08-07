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

   SurfaceNormalStressOperator(std::string subType, bool icModes){
     if( subType == "axi" ) {
       strainOp_ = new StrainOperatorAxi<FeH1,TYPE>(icModes);
     } else if( subType == "planeStrain" ) {
       strainOp_ = new StrainOperator2D<FeH1,TYPE>(icModes);
     } else if( subType == "planeStress" ) {
       strainOp_ = new StrainOperator2D<FeH1,TYPE>(icModes);
     } else if( subType == "3d") {
       strainOp_ = new StrainOperator3D<FeH1,TYPE>(icModes);
     } else {
       EXCEPTION( "Subtype '" << subType << "' in SurfaceNormalStressOperator" );
     }
     strainOp_->SetOperator2SurfOperator();
   }

    virtual ~SurfaceNormalStressOperator(){
      return;
    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe );

    //avoid reimplementation of complex operator by making the base class function
    //available
    using BaseBOperator::CalcOpMat;

    using BaseBOperator::CalcOpMatTransposed;

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
      normalOp.Resize( dimD, DIM_SPACE);
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

}
#endif
