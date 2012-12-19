#ifndef GRADIENTOP_HH
#define GRADIENTOP_HH

#include "BaseBOperator.hh"

namespace CoupledField{
  //! Calculate the gradient of the shape functions
  //!    / N_1x N_2x ...\
  //! b =| N_1y N_2y ...|
  //!    \ N_1z N_2z .../
  //!  here N_1x denotes the x-derivative of the first
  //!  shape function at a given local point
  //! \tparam FE Type of Finite Element used
  //! \tparam D Dimension of the problem space
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, UInt D, class TYPE = Double>
  class GradientOperator : public BaseBOperator{
    public:

    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1; // m=1 (Skalar), m=2,3 Vektor

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D; // n=2,3

    // Gradient normal: n x m, DIM_SPACE x DIM_DOF
    // Gradient transposed: m x n, DIM_DOF x DIM_SPACE <- Normale Schreibweise für Gradient

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D; // Dimension von Referenzelement
    
    //! Dimension of the related material 
    static const UInt DIM_D_MAT = D; 
    //@}


    GradientOperator(){
      this->name_ = "GradientOperator";
    }

    virtual ~GradientOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, 
                                     BaseFE* ptFe );

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
      return DIM_D_MAT;
    }
    //@}
    
    protected:

  };

  template<class FE, UInt D, class TYPE>
  void GradientOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                                const LocPointMapped& lp, 
                                                BaseFE* ptFe ){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( DIM_SPACE, numFncs );

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    bMat= Transpose(xiDx);
  }

  template<class FE, UInt D, class TYPE>
  void GradientOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                          const LocPointMapped& lp, 
                                                          BaseFE* ptFe ){
    const UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize(numFncs , DIM_SPACE );

    // Get derivatives of local shape functions with respect to global
    // coords (format: spaceDim x nrNodes )
    FE *fe = (static_cast<FE*>(ptFe));
    fe->GetGlobDerivShFnc( bMat, lp, lp.shapeMap->GetElem() , 1 );


  }

  //! Calculate the gradient of the shape functions scaled by a
  //! vector valued coefficient function.
  //!    / (s_x N_1x) (s_x N_2x) ...\
  //! b =| (s_y N_1y) (s_y N_2y) ...|
  //!    \ (s_z N_1z) (s_z N_2z) .../
  //!  here N_1x denotes the x-derivative of the first
  //!  shape function at a given local point
  //! \tparam FE Type of Finite Element used
  //! \tparam D Dimension of the problem space
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
  template<class FE, UInt D, class TYPE = Double>
  class ScaledGradientOperator : public GradientOperator<FE,D,TYPE>{
    public:

    // ------------------
    //  STATIC CONSTANTS
    // ------------------
    //@{
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D;

    //! Dimension of the related material
    static const UInt DIM_D_MAT = D;
    //@}


    ScaledGradientOperator(){
      this->name_ = "ScaledGradientOperator";
    }

    virtual ~ScaledGradientOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp,
                                     BaseFE* ptFe );

    virtual void CalcOpMat(Matrix<Complex> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,
                                     const LocPointMapped& lp,
                                     BaseFE* ptFe );

    protected:
    
  };

  template<class FE, UInt D, class TYPE>
  void ScaledGradientOperator<FE,D,TYPE>::CalcOpMat(Matrix<Complex> & bMat,
                                                const LocPointMapped& lp,
                                                BaseFE* ptFe ){
    assert(this->coef_ != NULL);
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( DIM_SPACE, numFncs * DIM_DOF );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx, xiDxTmp, rotMat;
    FE *fe = (static_cast<FE*>(ptFe));
    
    if (this->coef_->GetCoordinateSystem() ) {
      fe->GetGlobDerivShFnc( xiDxTmp, lp, lp.shapeMap->GetElem() , 1 );

      // If coordinate system is set at the coefficient function, rotate B-matrix
      Vector<Double> globPoint;
      lp.shapeMap->Local2Global(globPoint,lp.lp.coord);
      this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
      xiDx = xiDxTmp * rotMat;
    } else {
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    }


    Vector<Complex> coefs;
    this->coef_->GetVector(coefs,lp);

    if(DIM_DOF == 1){
      for(UInt i=0;i<numFncs;++i){
        for(UInt d = 0;d<DIM_SPACE;++d){
          bMat[d][i] = coefs[d] * xiDx[i][d];
        }
      }
    }else{

      for(UInt i = 0; i< numFncs ; ++i){
        for(UInt sDim = 0; sDim < DIM_SPACE; ++sDim){
          bMat[sDim][i*DIM_DOF + sDim] = xiDx[i][sDim] * coefs[sDim];
        }
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void ScaledGradientOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Complex> & bMat,
                                                          const LocPointMapped& lp,
                                                          BaseFE* ptFe ){
    assert(this->coef_ != NULL);
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( numFncs * DIM_DOF, DIM_SPACE  );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx, xiDxTmp, rotMat;
    FE *fe = (static_cast<FE*>(ptFe));

    if (this->coef_->GetCoordinateSystem() ) {
      fe->GetGlobDerivShFnc( xiDxTmp, lp, lp.shapeMap->GetElem() , 1 );

      // If coordinate system is set at the coefficient function, rotate B-matrix
      Vector<Double> globPoint;
      lp.shapeMap->Local2Global(globPoint,lp.lp.coord);
      this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
      xiDx = xiDxTmp * rotMat;
    } else {
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    }

    Vector<Complex> coefs;
    this->coef_->GetVector(coefs,lp);

    if(DIM_DOF == 1){
      for(UInt i=0;i<numFncs;++i){
        for(UInt d = 0;d<DIM_SPACE;++d){
          bMat[i][d] = coefs[d] * xiDx[i][d];
        }
      }
    }else{
      assert(DIM_SPACE == DIM_DOF);

      for(UInt i = 0; i< numFncs ; ++i){
        for(UInt sDim = 0; sDim < DIM_SPACE; ++sDim){
          bMat[i*DIM_DOF + sDim][sDim] = xiDx[i][sDim] * coefs[sDim];
        }
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void ScaledGradientOperator<FE,D,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                                const LocPointMapped& lp,
                                                BaseFE* ptFe ){
    assert(this->coef_ != NULL);

    GradientOperator<FE,D,TYPE>::CalcOpMat(bMat,lp,ptFe);
    Vector<Double> coefs;
    this->coef_->GetVector(coefs,lp);
    for(UInt i=0;i<bMat.GetNumCols();++i){
      for(UInt d = 0;d<bMat.GetNumRows();++d){
        bMat[d][i] *= coefs[d];
      }
    }
  }

  template<class FE, UInt D, class TYPE>
  void ScaledGradientOperator<FE,D,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                          const LocPointMapped& lp,
                                                          BaseFE* ptFe ){
    assert(this->coef_ != NULL);

    GradientOperator<FE,D,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
    Vector<Double> coefs;
    this->coef_->GetVector(coefs,lp);
    for(UInt i=0;i<bMat.GetNumRows();++i){
      for(UInt d = 0;d<bMat.GetNumCols();++d){
        bMat[i][d] *= coefs[d];
      }
    }
  }


}
#endif
