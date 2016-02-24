// =====================================================================================
// 
//       Filename:  identityOperator.hh
// 
//    Description:  This class implements the identity operator with a constant projection
//                  (just usefull for linear elements)
//                  Used in the stabilzation of flow computations applying linear basis 
//                  functions both for velocity and pressure
//                  consisting of the elements shape functions evaulated at the 
//                  given local coordinate
//                      / N_1x - const       0        0             N_2_x - const. ... \
//                  b = |  0            N_1y - const  0             0    .. |
//                      \  0                 0        N_1z-const.   0    .. /
//
//                  for the basic shape functions N_1x = N_1y = N_1z
//                  for scalar unknowns the operator is just a vector
// 
//         Author:  Manfred Kaltenbacher (TUW), manfred.kaltenbacher@tuwien.ac.at
// 
// =====================================================================================

#ifndef IDENTITYOPPROJECTED_HH
#define IDENTITYOPPROJECTED_HH
#include "BaseBOperator.hh"


#define USE_BLAS_VERSION

namespace CoupledField{
  
  template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
  class IdentityOperatorProjected : public BaseBOperator{

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


    IdentityOperatorProjected(){
      return;
    }

    //! Copy constructor
    IdentityOperatorProjected(const IdentityOperatorProjected & other)
       : BaseBOperator(other){
    }

    //! \copydoc BaseBOperator::Clone()
    virtual IdentityOperatorProjected * Clone(){
      return new IdentityOperatorProjected(*this);
    }

    virtual ~IdentityOperatorProjected(){
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
  void IdentityOperatorProjected<FE,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){
     UInt numFncs = ptFe->GetNumFncs();

     Double meanVal = 1.0 / (Double)numFncs;

     // Set correct size of matrix B and initialize with zeros
     bMat.Resize( DIM_DOF, numFncs * DIM_DOF );
     bMat.Init();

     Vector<Double> s;
     FE *fe = (static_cast<FE*>(ptFe));
     for(UInt d = 0; d < DIM_DOF ; d ++){
       fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
       for(UInt sh = 0; sh < numFncs; sh ++){
         bMat[d][sh*DIM_DOF + d] = s[sh] - meanVal;
       }
     }

   }

  template<class FE,  UInt D, UInt D_DOF, class TYPE>
    void IdentityOperatorProjected<FE,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    UInt numFncs = ptFe->GetNumFncs();
    Double meanVal = 1.0 / (Double)numFncs;

    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( numFncs * DIM_DOF , DIM_DOF );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < DIM_DOF ; d ++){
      fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < numFncs; sh ++){
        bMat[sh*DIM_DOF + d][d] = s[sh] - meanVal;
      }
    }

  }
  
   
 
} // end of namespace
#endif
