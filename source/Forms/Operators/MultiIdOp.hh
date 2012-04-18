// =====================================================================================
// 
//       Filename:  MultiIdOp.hh
// 
//    Description:  This class implements the identity operator just
//                  consisting of the elements shape functions evaulated at the 
//                  given local coordinate
//                      / N_1  N_2 .. \
//                  b = | N_1  N_2
//                      \ N_1  N_2 .. /
//                        ...  ... .. 
//                  for the basic shape functions N_1 N_2
//                  for scalar unknowns the operator is just a vector
// 
//        Version:  1.0
//        Created:  10/04/2011 09:10:00 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef MULTIIDENTITYOP_HH
#define MULTIIDENTITYOP_HH
#include "BaseBOperator.hh"


#define USE_BLAS_VERSION

namespace CoupledField{
  
  template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
  class MultiIdOp : public BaseBOperator{

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


    MultiIdOp(){
      return;
    }

    virtual ~MultiIdOp(){
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
  
  template<class FE,  UInt D, UInt D_DOF, class TYPE>
  void MultiIdOp<FE,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){
     UInt numFncs = ptFe->GetNumFncs();
     
     // Set correct size of matrix B and initialize with zeros
     bMat.Resize( DIM_SPACE, numFncs * DIM_DOF );
     bMat.Init();

     Vector<Double> s;
     FE *fe = (static_cast<FE*>(ptFe));
     for(UInt d = 0; d < DIM_SPACE ; d ++){
       fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
       for(UInt sh = 0; sh < numFncs; sh ++){
         bMat[d][sh] = s[sh];
       }
     }

   }

  template<class FE,  UInt D, UInt D_DOF, class TYPE>
    void MultiIdOp<FE,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_SPACE, numFncs * DIM_DOF );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_SPACE ; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[sh][d] = s[sh];
      }
    }

  }

} // end of namespace
#endif
