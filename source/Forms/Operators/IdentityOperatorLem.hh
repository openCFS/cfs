// =====================================================================================
// 
//       Filename:  identityOperatorLem.hh
// 
//    Description:  This class implements the identity operator just
//                  consisting of the predefined values for network elements
//                      / 1  -1 \
//                  b = \-1   1 /
//                      
//
//        Version:  1.0
//        Created:  27/03/2025 09:10:00 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Dominik Mayrhofer (dmayrhofer), dominik.mayrhofer@tugraz.at
//        Company:  TU Graz
// 
// =====================================================================================

#ifndef IDENTITYOPLEM_HH
#define IDENTITYOPLEM_HH
#include "BaseBOperator.hh"


namespace CoupledField{
  
  template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
  class IdentityOperatorLem : public BaseBOperator{

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


    IdentityOperatorLem(){
      return;
    }

    //! Copy constructor
    IdentityOperatorLem(const IdentityOperatorLem & other)
       : BaseBOperator(other){
    }

    //! \copydoc BaseBOperator::Clone()
    virtual IdentityOperatorLem * Clone(){
      return new IdentityOperatorLem(*this);
    }

    virtual ~IdentityOperatorLem(){
      return;
    }
.
    //! \param bMat (out) the output matrix
    //! \param lp (in) the local point mapped
    //! \param ptFe (in) pointer to the Fe function space
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
  
  template<class FE,  UInt D, UInt D_DOF, class TYPE>
  void IdentityOperatorLem<FE,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( 2, 2 );
    bMat.Init();

    // multiplying
    //               / 1  -1\
    // 1/sqrt(2)    |        |
    //               \-1   1/
    // with its transposed gives us the desired
    //  1  -1
    // -1   1
    //matrix
    bMat[0,0] = 1/sqrt(2);
    bMat[0,1] = -1/sqrt(2);
    bMat[1,0] = -1/sqrt(2);
    bMat[1,1] = 1/sqrt(2);
   }

  template<class FE,  UInt D, UInt D_DOF, class TYPE>
    void IdentityOperatorLem<FE,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( 2, 2 );
    bMat.Init();

    // multiplying
    //               / 1  -1\
    // 1/sqrt(2)    |        |
    //               \-1   1/
    // with its transposed gives us the desired
    //  1  -1
    // -1   1
    //matrix
    bMat[0,0] = 1/sqrt(2);
    bMat[0,1] = -1/sqrt(2);
    bMat[1,0] = -1/sqrt(2);
    bMat[1,1] = 1/sqrt(2);

  }
 
} // end of namespace
#endif
