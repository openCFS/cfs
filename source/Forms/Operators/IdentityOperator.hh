// =====================================================================================
// 
//       Filename:  identityOperator.hh
// 
//    Description:  This class implements the identity operator just
//                  consisting of the elements shape functions evaulated at the 
//                  given local coordinate
//                      / N_1x  0    0   N_2_x .. \
//                  b = |  0   N_1y  0    0    .. |
//                      \  0    0   N_1z  0    .. /
//
//                  for the basic shape functions N_1x = N_1y = N_1z
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

#ifndef IDENTITYOP_HH
#define IDENTITYOP_HH
#include "BaseBOperator.hh"
#include "FeBasis/HCurl/HCurlElems.hh"


namespace CoupledField{
  
  template<class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double>
  class IdentityOperator : public BaseBOperator{

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


    IdentityOperator(){
      return;
    }

    //! Copy constructor
    IdentityOperator(const IdentityOperator & other)
       : BaseBOperator(other){
    }

    //! \copydoc BaseBOperator::Clone()
    virtual IdentityOperator * Clone(){
      return new IdentityOperator(*this);
    }

    virtual ~IdentityOperator(){
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
  void IdentityOperator<FE,D,D_DOF,TYPE>::
  CalcOpMat(Matrix<Double> & bMat,
            const LocPointMapped& lp, BaseFE* ptFe){
     UInt numFncs = ptFe->GetNumFncs();
     
     // Set correct size of matrix B and initialize with zeros
     bMat.Resize( DIM_DOF, numFncs * DIM_DOF );
     bMat.Init();

     Vector<Double> s;
     FE *fe = (static_cast<FE*>(ptFe));
     for(UInt d = 0; d < DIM_DOF ; d ++){
       fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem() , d );
       for(UInt sh = 0; sh < numFncs; sh ++){
         bMat[d][sh*DIM_DOF + d] = s[sh];
       }
     }
   }

  template<class FE,  UInt D, UInt D_DOF, class TYPE>
    void IdentityOperator<FE,D,D_DOF,TYPE>::
  CalcOpMatTransposed(Matrix<Double> & bMat,
                      const LocPointMapped& lp, BaseFE* ptFe){
    UInt numFncs = ptFe->GetNumFncs();
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
        bMat[sh*DIM_DOF + d][d] = s[sh];
      }
    }

  }
  
   
   // =======================================================================
  
  //! Specialization for edge curl problems
  template<UInt D, class TYPE >
  class IdentityOperator<FeHCurl, D,1, TYPE> : public BaseBOperator{
    public:
    // =============================================
    //  STATIC CONSTANTS 
    // =============================================
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 1; 
    //@}
    
      IdentityOperator(){
        return;
      }

      //! Copy constructor
      IdentityOperator(const IdentityOperator & other)
         : BaseBOperator(other){
      }

      //! \copydoc BaseBOperator::Clone()
      virtual IdentityOperator * Clone(){
        return new IdentityOperator(*this);
      }

      ~IdentityOperator(){
        return;
      }

      virtual void CalcOpMat(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
        FeHCurl* fe = static_cast<FeHCurl*>(ptFe);
        fe->GetShFnc( bMat, lp, lp.shapeMap->GetElem(), 0);
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                       const LocPointMapped& lp, BaseFE* ptFe ){
        FeHCurl* fe = static_cast<FeHCurl*>(ptFe);
        Matrix<Double> xiDx;
        fe->GetShFnc( xiDx, lp, lp.shapeMap->GetElem(), 0);
        xiDx.Transpose(bMat);
      }
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

  //! Identity operator, which gets scaled by the average edge length
  
  //! This integrator gets gets used e.g. for the mass integrator within the
  //! magneticEdge PDE.
  template<class FE, UInt D, class TYPE>
  class ScaledByEdgeIdentityOperator : public IdentityOperator<FE,D,1,TYPE>{
    public:
      ScaledByEdgeIdentityOperator(){
        return;
      }

      //! Copy constructor
      ScaledByEdgeIdentityOperator(const ScaledByEdgeIdentityOperator & other)
         : IdentityOperator<FE,D,1,TYPE>(other){
      }

      //! \copydoc BaseBOperator::Clone()
      virtual ScaledByEdgeIdentityOperator * Clone(){
        return new ScaledByEdgeIdentityOperator(*this);
      }

      virtual ~ScaledByEdgeIdentityOperator(){
        return;
      }

      virtual void CalcOpMat(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
        IdentityOperator<FE,D,1,TYPE>::CalcOpMat(bMat,lp,ptFe);
        //scale By Edge
        Double minE,maxE;
        lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
        bMat /= maxE;
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                       const LocPointMapped& lp, BaseFE* ptFe ){
        IdentityOperator<FE,D,1,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
        //scale By Edge
        Double minE,maxE;
        lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
        bMat /= maxE;
      }

  };

  //! Identity operator, which takes account piola trasformation

  //! This integrator gets gets used e.g. for the mass integrator within the
  //! magneticEdge PDE.
  template< class FE, UInt D = 1, UInt D_DOF = 1, class TYPE = Double >
  class IdentityOperatorPiola : public IdentityOperator<FE,D,D_DOF,TYPE>{
    public:
      IdentityOperatorPiola(){
        return;
      }

      //! Copy constructor
      IdentityOperatorPiola(const IdentityOperatorPiola & other)
         : IdentityOperator<FE,D,D_DOF,TYPE>(other){
      }

      //! \copydoc BaseBOperator::Clone()
      virtual IdentityOperatorPiola * Clone(){
        return new IdentityOperatorPiola(*this);
      }

      virtual ~IdentityOperatorPiola(){
        return;
      }

      virtual void CalcOpMat(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
        Matrix<Double> bMatInitial;
        IdentityOperator<FE,D,D_DOF,TYPE>::CalcOpMat(bMatInitial,lp,ptFe);
        bMat.Resize(bMatInitial.GetNumRows(),bMatInitial.GetNumCols());
        bMat.Init();
        //now apply piola transform
        //in case of NC_SURF_ELEMs we evaluate the piola matrix on the volume element

        if(lp.isSurface){
#ifdef NDEBUG
          Double jacDetInv = (1.0/lp.lpmVol->jacDet);
          lp.lpmVol->jac.Mult_Blas(bMatInitial,bMat,false,false,jacDetInv,0.0);
#else
          bMat = lp.lpmVol->jac * bMatInitial;
          bMat *= (1.0/lp.lpmVol->jacDet);
#endif
        }else{
#ifdef NDEBUG
          Double jacDetInv = (1.0/lp.jacDet);
          lp.jac.Mult_Blas(bMatInitial,bMat,false,false,jacDetInv,0.0);
#else
          bMat = lp.jac * bMatInitial;
          bMat *= (1.0/lp.jacDet);
#endif
        }

      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                       const LocPointMapped& lp, BaseFE* ptFe ){
        Matrix<Double> bMatInitial;
        IdentityOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed(bMatInitial,lp,ptFe);
        //now apply piola transform
        bMat.Resize(bMatInitial.GetNumRows(),bMatInitial.GetNumCols());
        bMat.Init();
        if(lp.isSurface){
#ifdef NDEBUG
          Double jacDetInv = (1.0/lp.lpmVol->jacDet);
          bMatInitial.Mult_Blas(lp.lpmVol->jac,bMat,false,true,jacDetInv,0.0);
#else
          Matrix<Double> jacTmp;
          lp.lpmVol->jac.Transpose(jacTmp);
          bMat = bMatInitial * jacTmp;
          bMat *= (1.0/lp.lpmVol->jacDet);
#endif
        }else{
#ifdef NDEBUG
          Double jacDetInv = (1.0/lp.jacDet);
          bMatInitial.Mult_Blas(lp.jac,bMat,false,true,jacDetInv,0.0);
#else
          Matrix<Double> jacTmp;
          lp.jac.Transpose(jacTmp);
          bMat = bMatInitial * jacTmp;
          bMat *= (1.0/lp.jacDet);
#endif
        }
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available
      using IdentityOperator<FE,D,D_DOF,TYPE>::CalcOpMat;

      using IdentityOperator<FE,D,D_DOF,TYPE>::CalcOpMatTransposed;

  };

 
} // end of namespace
#endif
