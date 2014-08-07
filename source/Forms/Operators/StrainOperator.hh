#ifndef FILE_STRAINOP_HH
#define FILE_STRAINOP_HH

#include "BaseBOperator.hh"

namespace CoupledField{
  // ============================
  //  2D STRAIN OPERATOR (PLANE) 
  // ============================
  //! Strain-like differential operator in 2D
  template<class FE, class TYPE = Double >
  class StrainOperator2D : public BaseBOperator{

  public:

    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 2;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = 2;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 2;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 3; 
    //@}

    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    StrainOperator2D( bool useICModes = false )
    :  useICModes_(useICModes) {
      this->name_ = "StrainOperator2D";
    }

    //! Destructor
    virtual ~StrainOperator2D(){

    }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    //! Calculate transposed operator matrix
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
    
    //! Flag, if incompatible modes are used
    const bool useICModes_;

  };

  template<class FE, class TYPE>
  void StrainOperator2D<FE,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                            const LocPointMapped& lp, 
                                            BaseFE* ptFe ){
    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if( useICModes_ ) {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFncICModes(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFncICModes( xiDx, lp, lp.shapeMap->GetElem() , 1 );

    } else {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFnc(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    }
    
    const UInt numFncs = xiDx.GetNumRows();
    
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_D_MAT, numFncs * DIM_SPACE );
    bMat.Init();

    
    UInt iFunc = 0;
    UInt pos = 0;
    for( ; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[0][pos+0] = xiDx[iFunc][0];
    }

    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[1][pos+1] = xiDx[iFunc][1];
    }
    
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[2][pos+0] = xiDx[iFunc][1];
      bMat[2][pos+1] = xiDx[iFunc][0];
    }
  }

  template<class FE, class TYPE>
  void StrainOperator2D<FE,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                      const LocPointMapped& lp, 
                                                      BaseFE* ptFe ){
    // Get derivatives of local shape functions with respect to global
    // coords (format: spaceDim x nrNodes )
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    // query const variable, should be pretty much optimized away
    if( useICModes_ ) {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFncICModes(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFncICModes( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    } else {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFnc(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
    }

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize(numFncs * DIM_SPACE , DIM_D_MAT );
    bMat.Init();

    UInt iFunc = 0;
    UInt pos = 0;
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[pos+0][0] = xiDx[iFunc][0];
      bMat[pos+0][2] = xiDx[iFunc][1];
    }
    
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[pos+1][1] = xiDx[iFunc][1];
      bMat[pos+1][2] = xiDx[iFunc][0];
    }
  }

  // =====================================
  //  2D-AXIAL-SYMMETFRIC STRAIN OPERATOR  
  // =====================================
  //! Strain-like differential operator for axial symmetry 
  template<class FE, class TYPE = Double >
  class StrainOperatorAxi : public BaseBOperator{

  public:

    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 2;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = 2;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 2;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 4; 
    //@}

    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    StrainOperatorAxi( bool useICModes = false )
    :  useICModes_(useICModes) {
      this->name_ = "StrainOperatorAxi";
    }

    //! Destructor
    ~StrainOperatorAxi(){

    }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double> & bMat,
                           const LocPointMapped& lp,
                           BaseFE* ptFe );

    //! Calculate transposed operator matrix
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

    //! Flag, if incompatible modes are used
    const bool useICModes_;

  };

  template<class FE, class TYPE>
  void StrainOperatorAxi<FE,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                            const LocPointMapped& lpm, 
                                            BaseFE* ptFe ){

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if( useICModes_ ) {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFncICModes(  xiDx, *lpm.lpmVol, lpm.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFncICModes( xiDx, lpm, lpm.shapeMap->GetElem() , 1 );
    } else {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFnc(  xiDx, *lpm.lpmVol, lpm.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFnc( xiDx, lpm, lpm.shapeMap->GetElem() , 1 );
    }
    const UInt numFncs = xiDx.GetNumRows();
    
    // Set correct size of matrix B and initialize with zeros
    bMat.Resize( DIM_D_MAT, numFncs * DIM_SPACE );
    bMat.Init();
        
    // Calculate phi-phi component
    Vector<Double> shape;
    if( useICModes_ ) {
//      if ( isSurfOpt_ )
//        fe->GetShFncICModes( shape, *lpm.lpmVol, lpm.lpmVol->shapeMap->GetElem() );
//      else
        fe->GetShFncICModes( shape, lpm.lp, lpm.shapeMap->GetElem() );
    } else {
//      if ( isSurfOpt_ )
//        fe->GetShFnc( shape, *lpm.lpmVol, lpm.lpmVol->shapeMap->GetElem() );
//      else
        fe->GetShFnc( shape, lpm.lp, lpm.shapeMap->GetElem() );
    }
    Vector<Double> globPoint;
    lpm.shapeMap->Local2Global(globPoint, lpm.lp);
    const Double oneOverR = 1.0 /  globPoint[0];
    
    UInt iFunc = 0;
    UInt pos = 0;
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) { 
      bMat[0][pos+0] = xiDx[iFunc][0];
    }
    
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[1][pos+1] = xiDx[iFunc][1];
    }
    
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[2][pos+0] = xiDx[iFunc][1];
      bMat[2][pos+1] = xiDx[iFunc][0];
    }
      
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      // phi-phi component: N_i / r
      // The formula is taken from:
      // Zienkiewicz, The FEM - Vol.1 - The Basis, 5th ed., p. 114
      bMat[3][pos+0] = shape[iFunc] * oneOverR;
    }
  }

  template<class FE, class TYPE>
  void StrainOperatorAxi<FE,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                      const LocPointMapped& lpm, 
                                                      BaseFE* ptFe ){

    // Get derivatives of local shape functions with respect to global
    // coords (format: spaceDim x nrNodes )
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if( useICModes_ ) {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFncICModes(  xiDx, *lpm.lpmVol, lpm.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFncICModes( xiDx, lpm, lpm.shapeMap->GetElem() , 1 );
    } else {
      if ( isSurfOpt_ )
        fe->GetGlobDerivShFnc(  xiDx, *lpm.lpmVol, lpm.lpmVol->shapeMap->GetElem() , 1 );
      else
        fe->GetGlobDerivShFnc( xiDx, lpm, lpm.shapeMap->GetElem() , 1 );
    }
    
    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize(numFncs * DIM_SPACE , DIM_D_MAT );
    bMat.Init();

    // Calculate phi-phi component
    Vector<Double> shape;
    if( useICModes_ ) {
      if ( isSurfOpt_ )
        fe->GetShFncICModes( shape, lpm.lpmVol->lp, lpm.lpmVol->shapeMap->GetElem() );
      else
        fe->GetShFncICModes( shape, lpm.lp, lpm.shapeMap->GetElem() );
    } else {
        fe->GetShFnc( shape, lpm.lp, lpm.shapeMap->GetElem() );
    }
    Vector<Double> globPoint;
    lpm.shapeMap->Local2Global(globPoint, lpm.lp);
    const Double oneOverR = 1.0 /  globPoint[0];

    UInt iFunc = 0;
    UInt pos = 0;
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[pos+0][0] = xiDx[iFunc][0];
      bMat[pos+0][2] = xiDx[iFunc][1];
      
      // phi-phi component: N_i / r
      // The formula is taken from:
      // Zienkiewicz, The FEM - Vol.1 - The Basis, 5th ed., p. 114
      bMat[pos+0][3] = shape[iFunc] * oneOverR;
    }
    
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
      bMat[pos+1][1] = xiDx[iFunc][1];
      bMat[pos+1][2] = xiDx[iFunc][0];
    }
  }
  
  // ============================
   //  3D STRAIN OPERATOR 
   // ============================
   //! Strain-like differential operator in 3D
   template<class FE, class TYPE = Double >
   class StrainOperator3D : public BaseBOperator  {

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
     static const UInt DIM_D_MAT = 6; 
     //@}

     //! Constructor
     //! \param useIcModes Use incompatible modes shape functions
     StrainOperator3D( bool useICModes = false)
     :  useICModes_(useICModes) {
       this->name_ = "StrainOperator3D";
     }

     //! Destructor
     ~StrainOperator3D(){

     }

     //! Calculate operator matrix
     virtual void CalcOpMat(Matrix<Double> & bMat,
                            const LocPointMapped& lp,
                            BaseFE* ptFe );

     //! Calculate transposed operator matrix
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

     //! Flag, if incompatible modes are used
     const bool useICModes_;

   };

   template<class FE, class TYPE>
   void StrainOperator3D<FE,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                             const LocPointMapped& lp, 
                                             BaseFE* ptFe ){
     // Get derivatives of local shape functions with respect to global
     // coords (format: nrNodes x spaceDim)
     Matrix<Double> xiDx;
     FE *fe = (static_cast<FE*>(ptFe));
     if( useICModes_ ) {
       if ( isSurfOpt_ )
         fe->GetGlobDerivShFncICModes(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
       else
         fe->GetGlobDerivShFncICModes( xiDx, lp, lp.shapeMap->GetElem() , 1 );
     } else {
       if ( isSurfOpt_ )
         fe->GetGlobDerivShFnc(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
       else
         fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
     }
     
     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialize with zeros
     bMat.Resize( DIM_D_MAT, numFncs * DIM_SPACE );
     bMat.Init();
         
     UInt iFunc = 0;
     UInt pos = 0;
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[0][pos+0] = xiDx[iFunc][0];
     }

     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) { 
       bMat[1][pos+1] = xiDx[iFunc][1];
     }
     
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) { 
       bMat[2][pos+2] = xiDx[iFunc][2];
     }
     
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) { 
       bMat[3][pos+1] = xiDx[iFunc][2];
       bMat[3][pos+2] = xiDx[iFunc][1];
     }
     
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) { 
       bMat[4][pos+0] = xiDx[iFunc][2];
       bMat[4][pos+2] = xiDx[iFunc][0];
     }
     
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) { 
       bMat[5][pos+0] = xiDx[iFunc][1];
       bMat[5][pos+1] = xiDx[iFunc][0];
     }
   }

   template<class FE, class TYPE>
   void StrainOperator3D<FE,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                       const LocPointMapped& lp, 
                                                       BaseFE* ptFe ){
     // Get derivatives of local shape functions with respect to global
     // coords (format: spaceDim x nrNodes )
     Matrix<Double> xiDx;
     FE *fe = (static_cast<FE*>(ptFe));
     if( useICModes_ ) {
       if ( isSurfOpt_ )
         fe->GetGlobDerivShFncICModes(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
       else
         fe->GetGlobDerivShFncICModes( xiDx, lp, lp.shapeMap->GetElem() , 1 );
     } else {
       if ( isSurfOpt_ )
         fe->GetGlobDerivShFnc(  xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem() , 1 );
       else
         fe->GetGlobDerivShFnc( xiDx, lp, lp.shapeMap->GetElem() , 1 );
     }
     
     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialise with zeros
     bMat.Resize(numFncs * DIM_SPACE , DIM_D_MAT );
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[pos+0][0] = xiDx[iFunc][0];
       bMat[pos+0][4] = xiDx[iFunc][2];
       bMat[pos+0][5] = xiDx[iFunc][1];
     }

     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) { 
       bMat[pos+1][1] = xiDx[iFunc][1];
       bMat[pos+1][3] = xiDx[iFunc][2];
       bMat[pos+1][5] = xiDx[iFunc][0];
     }

     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) { 
       bMat[pos+2][2] = xiDx[iFunc][2];
       bMat[pos+2][3] = xiDx[iFunc][1];
       bMat[pos+2][4] = xiDx[iFunc][0];
     }
   }
}
#endif
