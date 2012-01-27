// =====================================================================================
// 
//       Filename:  curlOp.hh
// 
//    Description:  This operator implements the curl of the shape function
//                  Its operator-matrix has the structure
//                      / d N_1z/dy - d N_1y/dz ...\
//                  B = | d N_1x/dz - d N_1z/dx ...|
//                      \ d N_1y/dx - d N_1x/dy .../
//                      
//        Version:  1.0
//        Created:  10/28/2011 05:35:08 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#include "BaseBOperator.hh"
#include "FeBasis/H1/H1Elems.hh"

namespace CoupledField{
  template<class FE, UInt D, class TYPE>
  class CurlOperator : public BaseBOperator<FE,TYPE>{
    public:
    
    CurlOperator();

      ~CurlOperator();

      virtual void CalcOpMat(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
        EXCEPTION("Curl Operator not implemented for general functions yet")
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                       const LocPointMapped& lp, BaseFE* ptFe ){
        EXCEPTION("Curl Operator not implemented for general functions yet")
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available, template caused
      using BaseBOperator<FE,TYPE>::CalcOpMat;

      using BaseBOperator<FE,TYPE>::CalcOpMatTransposed;

    protected:

  };

  
  // ===================================
  //  CURL-OPERATOR FOR H-CURL ELEMENTS
  // ===================================
  //!Specialized class for HCurl elements
  template<UInt D, class TYPE >
  class CurlOperator<FeHCurl,D,TYPE> : public BaseBOperator<FeHCurl,TYPE>{
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

    CurlOperator(){
      this->name_ = "CurlOperator";

    }

    ~CurlOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe ){
      FeHCurl *fe = (static_cast<FeHCurl*>(ptFe));
      fe->GetCurlShFnc(bMat, lp, lp.shapeMap->GetElem(), 1);
    }

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe ){
      FeHCurl *fe = (static_cast<FeHCurl*>(ptFe));
      Matrix<Double> xiDx;
      fe->GetCurlShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
      xiDx.Transpose(bMat);
    }

    //avoid reimplementation of complex operator by making the bas class function
    //available, template caused
    using BaseBOperator<FeHCurl,TYPE>::CalcOpMat;

    using BaseBOperator<FeHCurl,TYPE>::CalcOpMatTransposed;

    protected:

  };

  // ===========================================================
  //  CURL-OPERATOR FOR H-CURL ELEMENTS (Scaled by Edge length)
  // ===========================================================
  
  //!Specialized class for curl elements
  template<UInt D, class TYPE>
  class ScaledByEdgeCurlOperator : public CurlOperator<FeHCurl,D,TYPE>{
  public:

    ScaledByEdgeCurlOperator(){

    }

    ~ScaledByEdgeCurlOperator(){

    }

    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp, BaseFE* ptFe ){
      CurlOperator<FeHCurl,D,TYPE>::CalcOpMat(bMat,lp,ptFe);
      Double minE,maxE;
      lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
      bMat /= maxE;
    }

    virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                     const LocPointMapped& lp, BaseFE* ptFe ){
      CurlOperator<FeHCurl,D,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
      Double minE,maxE;
      lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
      bMat /= maxE;
    }

    //avoid reimplementation of complex operator by making the bas class function
    //available, template caused
    using BaseBOperator<FeHCurl,TYPE>::CalcOpMat;

    using BaseBOperator<FeHCurl,TYPE>::CalcOpMatTransposed;

  protected:

  };
  
  // =====================================
  //  3D - CURL-OPERATOR FOR H-1 ELEMENTS
  // =====================================
  
  template<class TYPE>
  class CurlOperator<FeH1, 3, TYPE> : public BaseBOperator<FeH1,TYPE>{
  public:

    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 3;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = 3;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 3;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 3; 
    //@}
    
    CurlOperator(){
      this->name_ = "CurlOperator3D";

    }

    ~CurlOperator(){

    }
    
    void CalcOpMat(Matrix<Double> & bMat,
                   const LocPointMapped& lp, BaseFE* ptFe ){

      UInt numFncs = ptFe->GetNumFncs();
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( 3, numFncs * 3 );
      bMat.Init();

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

      // see Kaltenbacher, 2nd edition, p. 124
      UInt i = 0;
      UInt pos = 0;
      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[0][pos+1] = -xiDx[i][2];
        bMat[0][pos+2] =  xiDx[i][1];
      }

      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[1][pos+0] =  xiDx[i][2];
        bMat[1][pos+2] = -xiDx[i][0];
      }

      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[2][pos+0] = -xiDx[i][1];
        bMat[2][pos+1] =  xiDx[i][0];
      }
    }

    void CalcOpMatTransposed(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
      UInt numFncs = ptFe->GetNumFncs();
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( numFncs * 3, 3  );
      bMat.Init();

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

      // see Kaltenbacher, 2nd edition, p. 124
      UInt i = 0;
      UInt pos = 0;
      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[pos+0][1] =  xiDx[i][2];
        bMat[pos+0][2] = -xiDx[i][1];
      }

      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[pos+1][0] = -xiDx[i][2];
        bMat[pos+1][2] =  xiDx[i][0];
      }

      for( i = 0, pos = 0; i < numFncs; ++i, pos+=3 ) {
        bMat[pos+2][1] = -xiDx[i][0];
        bMat[pos+2][0] =  xiDx[i][1];
      }
    }

    //avoid reimplementation of complex operator by making the bas class function
    //available, template caused
    using BaseBOperator<FeH1,TYPE>::CalcOpMat;

    using BaseBOperator<FeH1,TYPE>::CalcOpMatTransposed;

  protected:

  };

  // =====================================
  //  2D - CURL-OPERATOR FOR H-1 ELEMENTS
  // =====================================
  
  //! Curl operator for 2D fields
  
  //! Note: We assume, that the field itself is a vector with only 
  //! a z-component. 
  template<class TYPE>
  class CurlOperator<FeH1, 2, TYPE> : public BaseBOperator<FeH1,TYPE>{
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
    static const UInt DIM_SPACE = 2;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 2;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 2; 
    //@}
    
    CurlOperator(){
      this->name_ = "CurlOperator2D";

    }

    ~CurlOperator(){

    }

    void CalcOpMat(Matrix<Double> & bMat,
                   const LocPointMapped& lp, BaseFE* ptFe ){
      UInt numFncs = ptFe->GetNumFncs();
      
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( 2, numFncs );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[0][iFunc] = -xiDx[iFunc][1];
      }
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[1][iFunc] =  xiDx[iFunc][0];
      }
    }

    void CalcOpMatTransposed(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
      
      UInt numFncs = ptFe->GetNumFncs();

      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( numFncs, 2 );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[iFunc][0] =  xiDx[iFunc][1];
        bMat[iFunc][1] = -xiDx[iFunc][0];
      }
    }

    //avoid reimplementation of complex operator by making the base class function
    //available, template caused
    using BaseBOperator<FeH1,TYPE>::CalcOpMat;

    using BaseBOperator<FeH1,TYPE>::CalcOpMatTransposed;

  protected:
  };
  
  
  // ====================================================
  //  2D - AXISYMMETRIC - CURL-OPERATOR FOR H-1 ELEMENTS
  // ====================================================
   template<class TYPE>
   class CurlOperatorAxi : public BaseBOperator<FeH1,TYPE>{
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
    static const UInt DIM_SPACE = 2;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 2;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 2; 
    //@}
    
    CurlOperatorAxi(){
      this->name_ = "CurlOperatorAxi";

    }

    ~CurlOperatorAxi(){

    }

    void CalcOpMat(Matrix<Double> & bMat,
                   const LocPointMapped& lp, BaseFE* ptFe ){
      UInt numFncs = ptFe->GetNumFncs();
      
      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( 2, numFncs );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );

      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[0][iFunc] = -xiDx[iFunc][1];
      }
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[1][iFunc] =  xiDx[iFunc][0];
      }
    }

    void CalcOpMatTransposed(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
      
      UInt numFncs = ptFe->GetNumFncs();

      // Set correct size of matrix B and initialize with zeros
      bMat.Resize( numFncs, 2 );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      FeH1 *fe = (static_cast<FeH1*>(ptFe));
      fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
      UInt iFunc = 0;
      for( iFunc = 0; iFunc < numFncs; ++iFunc ) {
        bMat[iFunc][0] = -xiDx[iFunc][1];
        bMat[iFunc][1] =  xiDx[iFunc][0];
      }
    }

    //avoid reimplementation of complex operator by making the base class function
    //available, template caused
    using BaseBOperator<FeH1,TYPE>::CalcOpMat;

    using BaseBOperator<FeH1,TYPE>::CalcOpMatTransposed;

  protected:  
   };
}
