// =====================================================================================
// 
//       Filename:  DaeOperator.hh
// 
//    Description:  This class implements operators which are necessary for 
//                  the construction of differential-algebraic equations (DAE)
// 
//        Version:  1.0
//        Created:  03/09/2022
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Samuel Kvasnicka, samuel.kvasnicka@tugraz.at
//        Company:  TU Graz
// 
// =====================================================================================

#ifndef DAEOPERATOR_HH
#define DAEOPERATOR_HH
#include "BaseBOperator.hh"
#include "FeBasis/HCurl/HCurlElems.hh"

namespace CoupledField{
  
  //! This class implements the identity operator, which gets projected in 
  //! normal direction. The element matrix is computed as 
  //!
  //! b = ( N_1*n_x + N_1*n_y + N_1*n_z  N_2*n_x + N_2*n_y + N_2*n_z  ..  ) = scalarproduct(\vec{N}, \vec{n})
  //!
  //! and is of size (1 x Number of functions (DOF)).
  template<class FE, UInt D , UInt D_DOF = 1, class TYPE = Double>
  class DaeIdentityOperatorNormal : public BaseBOperator{

  public:
    
    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D_DOF;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 1; 
    //@}


    DaeIdentityOperatorNormal(){
      return;
    }

    //! Copy constructor
    DaeIdentityOperatorNormal(const DaeIdentityOperatorNormal & other)
       : BaseBOperator(other){
    }

    //! \copydoc BaseBOperator::Clone()
    virtual DaeIdentityOperatorNormal * Clone(){
      return new DaeIdentityOperatorNormal(*this);
    }

    virtual ~DaeIdentityOperatorNormal(){
      return;
    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe ){
      EXCEPTION("Operator only implemented for HCurl")
    }

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe ){
      EXCEPTION("Operator only implemented for HCurl")
    }

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
      return DIM_D_MAT;
    }
    //@}
    
  protected:

};

  // ===============================================
  //  DaeIdentityOperatorNormal FOR H-CURL ELEMENTS
  // ===============================================
  //! Partial specialized class for HCurl elements

  //! This class implements the identity operator, which gets projected in 
  //! normal direction. The element matrix is computed as 
  //!
  //! b = ( N_1*n_x + N_1*n_y + N_1*n_z  N_2*n_x + N_2*n_y + N_2*n_z  ..  ) = scalarproduct(\vec{N}, \vec{n})
  //!
  //! and is of size (1 x Number of functions (DOF)).
  template<UInt D , UInt D_DOF, class TYPE>
  class DaeIdentityOperatorNormal<FeHCurl, D, D_DOF, TYPE> : public BaseBOperator{
  public:
    
    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = D_DOF;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 1; 
    //@}


    DaeIdentityOperatorNormal(){
      return;
    }

    //! Copy constructor
    DaeIdentityOperatorNormal(const DaeIdentityOperatorNormal & other)
       : BaseBOperator(other){
    }

    //! \copydoc BaseBOperator::Clone()
    virtual DaeIdentityOperatorNormal * Clone(){
      return new DaeIdentityOperatorNormal(*this);
    }

    virtual ~DaeIdentityOperatorNormal(){
      return;
    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe );

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe );

    //avoid reimplementation of complex operator by making the bas class function
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
      return DIM_D_MAT;
    }
    //@}
    
  protected:

};
  
  template<UInt D, UInt D_DOF, class TYPE>
  void DaeIdentityOperatorNormal<FeHCurl,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){

    assert(D == ptFe->shape_.dim); // ptFe should be a volume
    
    const UInt numFncs = ptFe->GetNumFncs();

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( 1, numFncs);
    bMat.Init();

    Matrix<Double> shFncMat;
    FeHCurl *fe = (static_cast<FeHCurl*>(ptFe));
    // Get HCurl shape function evaluated in lp.lp
    // Remark: last parameter "1" (component selection) does not have 
    //         any influence (currently not yet implemented).
    
    fe->GetShFnc( shFncMat, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem());

    for(UInt sh = 0; sh < numFncs; sh++){
      for(UInt d = 0; d < DIM_SPACE; d++){
    	  bMat[0][sh] += shFncMat[d][sh] * lp.normal[d];
      }
    }

    std::cout << std::endl;
    std::cout << "shFncMat: " << std::endl;
    for (int i = 0; i < DIM_SPACE; i++)
    {
      for (int j = 0; j < numFncs; j++)
      {
        std::cout << shFncMat[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "normal: " << std::endl;
    for (int i = 0; i < DIM_SPACE; i++)
    {
      std::cout << lp.normal[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "bMat: " << std::endl;
    for (int i = 0; i < 1; i++)
    {
      for (int j = 0; j < numFncs; j++)
      {
        std::cout << bMat[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  template<UInt D, UInt D_DOF, class TYPE>
  void DaeIdentityOperatorNormal<FeHCurl,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    // ensure, that the surface information (i.e. normal direction)
    // is set at the mapped local point 
    //assert(lp.isSurface);
    assert(D == ptFe->shape_.dim); // ptFe should be a volume
    
    Matrix<Double> tmpMat;
    this->CalcOpMat(tmpMat,lp,ptFe);
    bMat = Transpose(tmpMat);
  }

} // end of namespace
#endif