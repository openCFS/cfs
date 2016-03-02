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

    //! Copy constructor
    StrainOperator2D(const StrainOperator2D & other)
       :  BaseBOperator(other),
          useICModes_(other.useICModes_) {
    }

    //! \copydoc BaseBOperator::Clone()
    virtual StrainOperator2D * Clone(){
      return new StrainOperator2D(*this);
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
    bool useICModes_;

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

    //! Copy constructor
    StrainOperatorAxi(const StrainOperatorAxi & other)
       :  BaseBOperator(other), useICModes_(other.useICModes_) {
    }

    //! \copydoc BaseBOperator::Clone()
    virtual StrainOperatorAxi * Clone(){
      return new StrainOperatorAxi(*this);
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

    //! Copy constructor
    StrainOperator3D(const StrainOperator3D & other)
       : BaseBOperator(other),
         useICModes_(other.useICModes_){
      this->name_ = other.name_;
    }

    //! \copydoc BaseBOperator::Clone()
    virtual StrainOperator3D * Clone(){
      return new StrainOperator3D(*this);
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

  // ============================
  //  2.5D STRAIN OPERATOR
  // ============================
  //! Strain-like 3D differential operator in 2D
  template<class FE, class TYPE = Double>
  class StrainOperator2p5D : public BaseBOperator  {

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
    static const UInt DIM_SPACE = 2;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 2;

    //! Dimension of the related material
    static const UInt DIM_D_MAT = 6;
    //@}

    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    StrainOperator2p5D(bool useICModes = false)
    :  useICModes_(useICModes) {
      this->name_ = "StrainOperator2.5D";
    }

    //! Copy constructor
    StrainOperator2p5D(const StrainOperator2p5D & other)
       : BaseBOperator(other),
         useICModes_(other.useICModes_){
      this->name_ = other.name_;
    }

    //! \copydoc BaseBOperator::Clone()
    virtual StrainOperator2p5D * Clone(){
      return new StrainOperator2p5D(*this);
    }

    //! Destructor
    ~StrainOperator2p5D(){

    }

    //! Calculate operator matrix
    virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

    //! Calculate transposed operator matrix
    virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

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
  void StrainOperator2p5D<FE, TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if(useICModes_)
    {
      if (isSurfOpt_)
        fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
    }
    else
    {
      if (isSurfOpt_)
        fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
    }

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialize with zeros
    // In 2.5D case DIM_DOF = 3 is to be used instead of conventional DIM_SPACE = 2
    bMat.Resize(DIM_D_MAT, numFncs*DIM_DOF);
    bMat.Init();

    UInt iFunc = 0;
    UInt pos = 0;
    for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
    {
      bMat[0][pos+0] = xiDx[iFunc][0];
    }

    for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
    {
      bMat[1][pos+1] = xiDx[iFunc][1];
    }

    for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
    {
      bMat[3][pos+2] = xiDx[iFunc][1];
    }

    for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
    {
      bMat[4][pos+2] = xiDx[iFunc][0];
    }

    for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
    {
      bMat[5][pos+0] = xiDx[iFunc][1];
      bMat[5][pos+1] = xiDx[iFunc][0];
    }
  }

  template<class FE, class TYPE>
  void StrainOperator2p5D<FE, TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
  {
    // Get derivatives of local shape functions with respect to global
    // coords (format: spaceDim x nrNodes )
    Matrix<Double> xiDx;
    FE *fe = (static_cast<FE*>(ptFe));
    if (useICModes_)
    {
      if (isSurfOpt_)
        fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
    }
    else
    {
      if (isSurfOpt_)
        fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
      else
        fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
    }

    const UInt numFncs = xiDx.GetNumRows();
    // Set correct size of matrix B and initialise with zeros
    // In 2.5D case DIM_DOF = 3 is to be used instead of conventional DIM_SPACE = 2
    bMat.Resize(numFncs*DIM_DOF , DIM_D_MAT);
    bMat.Init();

    UInt iFunc = 0;
    UInt pos = 0;
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF ) {
      bMat[pos+0][0] = xiDx[iFunc][0];
      bMat[pos+0][5] = xiDx[iFunc][1];
    }

    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF ) {
      bMat[pos+1][1] = xiDx[iFunc][1];
      bMat[pos+1][5] = xiDx[iFunc][0];
    }

    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF ) {
      bMat[pos+2][3] = xiDx[iFunc][1];
      bMat[pos+2][4] = xiDx[iFunc][0];
    }
  }


   /** This is the StrainOperator2D Variant for Bloch mode eigenfrequency analysis with additional imaginary part for B */
   template<class FE, class TYPE = Complex >
   class StrainOperatorBloch2D : public StrainOperator2D<FE, Double>
   {
   public:

     /** The wave vector needs to be set!
      * @see SetWaveVector() */
     StrainOperatorBloch2D(bool useICModes = false)
     {
       assert(useICModes == false);

       this->name_ = "StrainOperatorBloch2D";
       this->wave_vector_ = NULL; // to be set
     }

     //! Copy constructor
     StrainOperatorBloch2D(const StrainOperatorBloch2D & other)
        : StrainOperator2D<FE, Double>(other){
       this->name_ = other.name_;
       //deep copy
       this->wave_vector_ = new Vector<double>(other.wave_vector_->GetSize());
       //copy values
       for(UInt i=0;i<this->wave_vector_->GetSize();i++){
         (*this->wave_vector_)[i] = (*other.wave_vector_)[i];
       }
     }

     //! \copydoc BaseBOperator::Clone()
     virtual StrainOperatorBloch2D * Clone(){
       return new StrainOperatorBloch2D(*this);
     }

     //! Destructor
     virtual ~StrainOperatorBloch2D(){ }

     /** reference to the always up to date wave vector. Comes from EigenFrequencyDrive::GetCurrentWaveVector() */
     void SetWaveVector(Vector<double>& current_wave_vector) {
       wave_vector_ = &current_wave_vector;
     };


     //! Calculate operator matrix
     void CalcOpMat(Matrix<Complex> & bMat, const LocPointMapped& lp, BaseFE* ptFe );

     //! Calculate transposed operator matrix
     void CalcOpMatTransposed(Matrix<Complex> & bMat, const LocPointMapped& lp, BaseFE* ptFe );

     static const UInt ORDER_DIFF = 1;
     static const UInt DIM_DOF = 2;
     static const UInt DIM_SPACE = 2;
     static const UInt DIM_ELEM = 2;
     static const UInt DIM_D_MAT = 3;
     UInt GetDiffOrder() const { return ORDER_DIFF; }
     UInt GetDimDof() const { return DIM_DOF; }
     UInt GetDimSpace() const { return DIM_SPACE; }
     UInt GetDimElem() const { return DIM_ELEM; }
     UInt GetDimDMat() const { return DIM_D_MAT; }


   private:
     Vector<double>* wave_vector_;
   };

   template<class FE, class TYPE>
   void StrainOperatorBloch2D<FE,TYPE>::CalcOpMat(Matrix<Complex> & bMat, const LocPointMapped& lpm, BaseFE* ptFe )
   {
     unsigned int numFncs = ptFe->GetNumFncs();
     // the implementation follows Hussein; Reduced Bloch mode expansion for periodic media band structure calculations
     assert(wave_vector_ != NULL && wave_vector_->GetSize() >= 2);
     Vector<double>& wv = *wave_vector_;

     // the real part
     Matrix<double> real;
     StrainOperator2D<FE,double>::CalcOpMat(real, lpm, ptFe);
     assert(real.GetNumRows() == DIM_D_MAT);
     assert(real.GetNumCols() == numFncs * DIM_SPACE);

     // std::cout << "SOB2d:COM k=" << wv.ToString() << std::endl;

     // assemble with the imaginary part
     bMat.Resize(real.GetNumRows(), real.GetNumCols());
     bMat.Init(); // reset to 0 as we only fill partially
     Vector<double> shape; // N_a
     static_cast<FE*>(ptFe)->GetShFnc(shape, lpm.lp, lpm.shapeMap->GetElem());

     unsigned int iFunc = 0;
     unsigned int pos = 0;
     for(iFunc = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[0][pos+0] = Complex(real[0][pos+0], shape[iFunc] * wv[0]);  // N_a,x + j k_x N_a
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[1][pos+1] = Complex(real[1][pos+1], shape[iFunc] * wv[1]);  // N_a,y + j k_y N_a
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[2][pos+0] = Complex(real[2][pos+0], shape[iFunc] * wv[1]);  // N_a,y + j k_y N_a
       bMat[2][pos+1] = Complex(real[2][pos+1], shape[iFunc] * wv[0]);  // N_a,x + j k_x N_a
     }
     // std::cout << "SOB2D::COM e=" << lpm.ptEl->elemNum << " wv=" << wv.ToString() << " bMat -> " << bMat.ToString() << std::endl;
   }

   template<class FE, class TYPE>
   void StrainOperatorBloch2D<FE,TYPE>::CalcOpMatTransposed(Matrix<Complex> & bMat, const LocPointMapped& lpm, BaseFE* ptFe )
   {
     // is typically not called but BDBInt::CalcElementMatrix calls CalcOpMat and BLAS transposes it.

     Matrix<Complex> tmp;
     CalcOpMat(tmp, lpm, ptFe);
     tmp.Transpose(bMat);
     // tmp.TransposeConjugate(bMat);
     // std::cout << "SOB2D::COM^T e=" << lpm.ptEl->elemNum << " bMat -> " << bMat.ToString() << std::endl;
   }


   template<class FE, class TYPE = Complex >
   class StrainOperatorBloch3D : public StrainOperator3D<FE, Double>
   {
   public:

     /** The wave vector needs to be set!
      * @see SetWaveVector() */
     StrainOperatorBloch3D(bool useICModes = false)
     {
       assert(useICModes == false);

       this->name_ = "StrainOperatorBloch3D";
       this->wave_vector_ = NULL; // to be set
     }

     //! Copy constructor
     StrainOperatorBloch3D(const StrainOperatorBloch3D & other)
        : StrainOperator3D<FE, Double>(other){
       this->wave_vector_ = new Vector<double>(other.wave_vector_->GetSize());
       //copy values
       for(UInt i=0;i<this->wave_vector_->GetSize();i++){
         (*this->wave_vector_)[i] = (*other.wave_vector_)[i];
       }
     }

     //! \copydoc BaseBOperator::Clone()
     StrainOperatorBloch3D * Clone(){
       return new StrainOperatorBloch3D(*this);
     }

     //! Destructor
     ~StrainOperatorBloch3D(){ }

     /** reference to the always up to date wave vector. Comes from EigenFrequencyDrive::GetCurrentWaveVector() */
     void SetWaveVector(Vector<double>& current_wave_vector) {
       wave_vector_ = &current_wave_vector;
     };


     //! Calculate operator matrix
     void CalcOpMat(Matrix<Complex> & bMat, const LocPointMapped& lp, BaseFE* ptFe );

     //! Calculate transposed operator matrix
     void CalcOpMatTransposed(Matrix<Complex> & bMat, const LocPointMapped& lp, BaseFE* ptFe );

     static const UInt ORDER_DIFF = 1;
     static const UInt DIM_DOF = 3;
     static const UInt DIM_SPACE = 3;
     static const UInt DIM_ELEM = 3;
     static const UInt DIM_D_MAT = 6;
     UInt GetDiffOrder() const { return ORDER_DIFF; }
     UInt GetDimDof() const { return DIM_DOF; }
     UInt GetDimSpace() const { return DIM_SPACE; }
     UInt GetDimElem() const { return DIM_ELEM; }
     UInt GetDimDMat() const { return DIM_D_MAT; }


   private:
     Vector<double>* wave_vector_;
   };

   template<class FE, class TYPE>
   void StrainOperatorBloch3D<FE,TYPE>::CalcOpMat(Matrix<Complex> & bMat, const LocPointMapped& lpm, BaseFE* ptFe )
   {
     unsigned int numFncs = ptFe->GetNumFncs();
     // the implementation follows Hussein; Reduced Bloch mode expansion for periodic media band structure calculations -> (A5)
     assert(wave_vector_ != NULL && wave_vector_->GetSize() == 3);
     Vector<double>& wv = *wave_vector_;

     // the real part
     Matrix<double> real;
     StrainOperator3D<FE,double>::CalcOpMat(real, lpm, ptFe);
     assert(real.GetNumRows() == DIM_D_MAT);
     assert(real.GetNumCols() == numFncs * DIM_SPACE);

     // std::cout << "SOB3D:COM k=" << wv.ToString() << std::endl;

     // assemble with the imaginary part
     bMat.Resize(real.GetNumRows(), real.GetNumCols());
     bMat.Init(); // reset to 0 as we only fill partially
     Vector<double> shape; // N_a
     static_cast<FE*>(ptFe)->GetShFnc(shape, lpm.lp, lpm.shapeMap->GetElem());

     unsigned int iFunc = 0;
     unsigned int pos = 0;
     for(iFunc = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[0][pos+0] = Complex(real[0][pos+0], shape[iFunc] * wv[0]);  // N_a,x + j k_x N_a
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[1][pos+1] = Complex(real[1][pos+1], shape[iFunc] * wv[1]);  // N_a,y + j k_y N_a
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[2][pos+2] = Complex(real[2][pos+2], shape[iFunc] * wv[2]);  // N_a,z + j k_z N_a
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[3][pos+1] = Complex(real[3][pos+1], shape[iFunc] * wv[2]);  // N_a,z + j k_z N_a
       bMat[3][pos+2] = Complex(real[3][pos+2], shape[iFunc] * wv[1]);  // N_a,y + j k_y N_a
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[4][pos+0] = Complex(real[4][pos+0], shape[iFunc] * wv[2]);  // N_a,z + j k_z N_a
       bMat[4][pos+2] = Complex(real[4][pos+2], shape[iFunc] * wv[0]);  // N_a,x + j k_x N_a
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE ) {
       bMat[5][pos+0] = Complex(real[5][pos+0], shape[iFunc] * wv[1]);  // N_a,y + j k_y N_a
       bMat[5][pos+1] = Complex(real[5][pos+1], shape[iFunc] * wv[0]);  // N_a,x + j k_x N_a
     }
     //std::cout << bMat.ToString() << std::endl;
   }

   template<class FE, class TYPE>
   void StrainOperatorBloch3D<FE,TYPE>::CalcOpMatTransposed(Matrix<Complex> & bMat, const LocPointMapped& lpm, BaseFE* ptFe )
   {
     // is typically not called but BDBInt::CalcElementMatrix calls CalcOpMat and BLAS transposes it.

     Matrix<Complex> tmp;
     CalcOpMat(tmp, lpm, ptFe);
     tmp.Transpose(bMat);
     // Conjugate(bMat);
   }

   // ============================
   //  2D SCALED STRAIN OPERATOR (PLANE)
   // ============================
   //! Scaled strain-like differential operator in 2D
   template<class FE, class TYPE = Double >
   class ScaledStrainOperator2D : public StrainOperator2D<FE, TYPE>
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
     ScaledStrainOperator2D(bool useICModes = false) : StrainOperator2D<FE, TYPE>(useICModes)
     {
       this->name_ = "ScaledStrainOperator2D";
     }

     //! Copy constructor
     ScaledStrainOperator2D(const ScaledStrainOperator2D & other)
        : StrainOperator2D<FE, TYPE>(other){
       this->name_ = other.name_;
     }

     //! \copydoc BaseBOperator::Clone()
     virtual ScaledStrainOperator2D * Clone(){
       return new ScaledStrainOperator2D(*this);
     }

     //! Destructor
     virtual ~ScaledStrainOperator2D() { }

     //! Calculate operator matrix
     virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate transposed operator matrix
     virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate operator matrix
     virtual void CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate transposed operator matrix
     virtual void CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     // ===============
     //  QUERY METHODS
     // ===============
     //@{ \name Query Methods
     //! \copydoc BaseBOperator::GetDiffOrder
     virtual UInt GetDiffOrder() const
     {
       return ORDER_DIFF;
     }

     //! \copydoc BaseBOperator::GetDimDof()
     virtual UInt GetDimDof() const
     {
       return DIM_DOF;
     }

     //! \copydoc BaseBOperator::GetDimSpace()
     virtual UInt GetDimSpace() const
     {
       return DIM_SPACE;
     }

     //! \copydoc BaseBOperator::GetDimElem()
     virtual UInt GetDimElem() const
     {
       return DIM_ELEM;
     }

     //! \copydoc BaseBOperator::GetDimDMat()
     virtual UInt GetDimDMat() const
     {
       return DIM_D_MAT;
     }
     //@}

   };

   template<class FE, class TYPE>
   void ScaledStrainOperator2D<FE,TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Double> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: nrNodes x spaceDim)
     Matrix<Double> xiDx;
     FE *fe = (static_cast<FE*>(ptFe));
     if (this->useICModes_)
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
     else
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialize with zeros
     bMat.Resize(DIM_D_MAT, numFncs*DIM_SPACE);
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0 ; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[0][pos+0] = xiDx[iFunc][0]*coefs[0];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[1][pos+1] = xiDx[iFunc][1]*coefs[1];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[2][pos+0] = xiDx[iFunc][1]*coefs[1];
       bMat[2][pos+1] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator2D<FE,TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Double> coefs;
     this->coef_->GetVector(coefs, lp);
     // Get derivatives of local shape functions with respect to global
     // coords (format: spaceDim x nrNodes )
     Matrix<Double> xiDx;
     FE *fe = (static_cast<FE*>(ptFe));
     // query const variable, should be pretty much optimized away
     if (this->useICModes_)
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
     else
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialise with zeros
     bMat.Resize(numFncs*DIM_SPACE , DIM_D_MAT);
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+0][0] = xiDx[iFunc][0]*coefs[0];
       bMat[pos+0][2] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+1][1] = xiDx[iFunc][1]*coefs[1];
       bMat[pos+1][2] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator2D<FE,TYPE>::CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Complex> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: nrNodes x spaceDim)
     Matrix<Double> xiDx, xiDxTmp, rotMat;
     FE *fe = (static_cast<FE*>(ptFe));

     if (this->coef_->GetCoordinateSystem())
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);

       // If coordinate system is set at the coefficient function, rotate B-matrix
       Vector<Double> globPoint;
       lp.shapeMap->Local2Global(globPoint, lp.lp.coord);
       this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
       xiDx = xiDxTmp*rotMat;
     }
     else
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
     }

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialize with zeros
     bMat.Resize(DIM_D_MAT, numFncs*DIM_SPACE);
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0 ; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[0][pos+0] = xiDx[iFunc][0]*coefs[0];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[1][pos+1] = xiDx[iFunc][1]*coefs[1];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[2][pos+0] = xiDx[iFunc][1]*coefs[1];
       bMat[2][pos+1] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator2D<FE,TYPE>::CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Complex> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: spaceDim x nrNodes )
     Matrix<Double> xiDx, xiDxTmp, rotMat;
     FE *fe = (static_cast<FE*>(ptFe));
     // query const variable, should be pretty much optimized away
     if (this->coef_->GetCoordinateSystem())
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);

       // If coordinate system is set at the coefficient function, rotate B-matrix
       Vector<Double> globPoint;
       lp.shapeMap->Local2Global(globPoint, lp.lp.coord);
       this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
       xiDx = xiDxTmp*rotMat;
     }
     else
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
     }

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialise with zeros
     bMat.Resize(numFncs*DIM_SPACE , DIM_D_MAT);
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+0][0] = xiDx[iFunc][0]*coefs[0];
       bMat[pos+0][2] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+1][1] = xiDx[iFunc][1]*coefs[1];
       bMat[pos+1][2] = xiDx[iFunc][0]*coefs[0];
     }
   }

   // ============================
   //  3D SCALED STRAIN OPERATOR
   // ============================
   //! Scaled strain-like differential operator in 3D
   template<class FE, class TYPE = Double >
   class ScaledStrainOperator3D : public StrainOperator3D<FE, TYPE>
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
     ScaledStrainOperator3D(bool useICModes = false) : StrainOperator3D<FE, TYPE>(useICModes)
     {
       this->name_ = "ScaledStrainOperator3D";
     }

     //! Copy constructor
     ScaledStrainOperator3D(const ScaledStrainOperator3D & other)
        : StrainOperator3D<FE, TYPE>(other){
       this->name_ = other.name_;
     }

     //! \copydoc BaseBOperator::Clone()
     virtual ScaledStrainOperator3D * Clone(){
       return new ScaledStrainOperator3D(*this);
     }

     //! Destructor
     virtual ~ScaledStrainOperator3D() { }

     //! Calculate operator matrix
     virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate transposed operator matrix
     virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate operator matrix
     virtual void CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate transposed operator matrix
     virtual void CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     // ===============
     //  QUERY METHODS
     // ===============
     //@{ \name Query Methods
     //! \copydoc BaseBOperator::GetDiffOrder
     virtual UInt GetDiffOrder() const
     {
       return ORDER_DIFF;
     }

     //! \copydoc BaseBOperator::GetDimDof()
     virtual UInt GetDimDof() const
     {
       return DIM_DOF;
     }

     //! \copydoc BaseBOperator::GetDimSpace()
     virtual UInt GetDimSpace() const
     {
       return DIM_SPACE;
     }

     //! \copydoc BaseBOperator::GetDimElem()
     virtual UInt GetDimElem() const
     {
       return DIM_ELEM;
     }

     //! \copydoc BaseBOperator::GetDimDMat()
     virtual UInt GetDimDMat() const
     {
       return DIM_D_MAT;
     }
     //@}

   };

   template<class FE, class TYPE>
   void ScaledStrainOperator3D<FE,TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Double> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: nrNodes x spaceDim)
     Matrix<Double> xiDx;
     FE *fe = (static_cast<FE*>(ptFe));
     if (this->isSurfOpt_)
       fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
     else
       fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialize with zeros
     bMat.Resize(DIM_D_MAT, numFncs*DIM_SPACE);
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[0][pos+0] = xiDx[iFunc][0]*coefs[0];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[1][pos+1] = xiDx[iFunc][1]*coefs[1];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[2][pos+2] = xiDx[iFunc][2]*coefs[2];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[3][pos+1] = xiDx[iFunc][2]*coefs[2];
       bMat[3][pos+2] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[4][pos+0] = xiDx[iFunc][2]*coefs[2];
       bMat[4][pos+2] = xiDx[iFunc][0]*coefs[0];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[5][pos+0] = xiDx[iFunc][1]*coefs[1];
       bMat[5][pos+1] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator3D<FE,TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Double> coefs;
     this->coef_->GetVector(coefs, lp);
     // Get derivatives of local shape functions with respect to global
     // coords (format: spaceDim x nrNodes )
     Matrix<Double> xiDx;
     FE *fe = (static_cast<FE*>(ptFe));
     // query const variable, should be pretty much optimized away
     if (this->useICModes_)
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
     else
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialise with zeros
     bMat.Resize(numFncs*DIM_SPACE, DIM_D_MAT);
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+0][0] = xiDx[iFunc][0]*coefs[0];
       bMat[pos+0][4] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+0][5] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+1][1] = xiDx[iFunc][1]*coefs[1];
       bMat[pos+1][3] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+1][5] = xiDx[iFunc][0]*coefs[0];
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+2][2] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+2][3] = xiDx[iFunc][1]*coefs[1];
       bMat[pos+2][4] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator3D<FE,TYPE>::CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Complex> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: nrNodes x spaceDim)
     Matrix<Double> xiDx, xiDxTmp, rotMat;
     FE *fe = (static_cast<FE*>(ptFe));

     if (this->coef_->GetCoordinateSystem())
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);

       // If coordinate system is set at the coefficient function, rotate B-matrix
       Vector<Double> globPoint;
       lp.shapeMap->Local2Global(globPoint, lp.lp.coord);
       this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
       xiDx = xiDxTmp*rotMat;
     }
     else
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
     }

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialize with zeros
     bMat.Resize(DIM_D_MAT, numFncs*DIM_SPACE);
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[0][pos+0] = xiDx[iFunc][0]*coefs[0];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[1][pos+1] = xiDx[iFunc][1]*coefs[1];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
       bMat[2][pos+2] = xiDx[iFunc][2]*coefs[2];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[3][pos+1] = xiDx[iFunc][2]*coefs[2];
       bMat[3][pos+2] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[4][pos+0] = xiDx[iFunc][2]*coefs[2];
       bMat[4][pos+2] = xiDx[iFunc][0]*coefs[0];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[5][pos+0] = xiDx[iFunc][1]*coefs[1];
       bMat[5][pos+1] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator3D<FE,TYPE>::CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Complex> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: spaceDim x nrNodes )
     Matrix<Double> xiDx, xiDxTmp, rotMat;
     FE *fe = (static_cast<FE*>(ptFe));
     // query const variable, should be pretty much optimized away
     if (this->coef_->GetCoordinateSystem())
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);

       // If coordinate system is set at the coefficient function, rotate B-matrix
       Vector<Double> globPoint;
       lp.shapeMap->Local2Global(globPoint, lp.lp.coord);
       this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
       xiDx = xiDxTmp*rotMat;
     }
     else
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
     }

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialise with zeros
     bMat.Resize(numFncs*DIM_SPACE , DIM_D_MAT );
     bMat.Init();

     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+0][0] = xiDx[iFunc][0]*coefs[0];
       bMat[pos+0][4] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+0][5] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+1][1] = xiDx[iFunc][1]*coefs[1];
       bMat[pos+1][3] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+1][5] = xiDx[iFunc][0]*coefs[0];
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     {
       bMat[pos+2][2] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+2][3] = xiDx[iFunc][1]*coefs[1];
       bMat[pos+2][4] = xiDx[iFunc][0]*coefs[0];
     }
   }

   // ============================
   //  2.5D SCALED STRAIN OPERATOR
   // ============================
   //! Scaled strain-like 3D differential operator in 2D
   template<class FE, class TYPE = Double >
   class ScaledStrainOperator2p5D : public StrainOperator2p5D<FE, TYPE>
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
     static const UInt DIM_DOF = 3;

     //! Dimension of the underlying domain / space
     static const UInt DIM_SPACE = 2;

     //! Dimension of the finite element
     static const UInt DIM_ELEM = 2;

     //! Dimension of the related material
     static const UInt DIM_D_MAT = 6;
     //@}

     //! Constructor
     //! \param useIcModes Use incompatible modes shape functions
     ScaledStrainOperator2p5D(bool useICModes = false) : StrainOperator2p5D<FE, TYPE>(useICModes)
     {
       this->name_ = "ScaledStrainOperator2p5D";
     }

     //! Copy constructor
     ScaledStrainOperator2p5D(const ScaledStrainOperator2p5D & other)
        : StrainOperator2p5D<FE, TYPE>(other){
       this->name_ = other.name_;
     }

     //! \copydoc BaseBOperator::Clone()
     virtual ScaledStrainOperator2p5D * Clone(){
       return new ScaledStrainOperator2p5D(*this);
     }


     //! Destructor
     virtual ~ScaledStrainOperator2p5D() { }

     //! Calculate operator matrix
     virtual void CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate transposed operator matrix
     virtual void CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate operator matrix
     virtual void CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     //! Calculate transposed operator matrix
     virtual void CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe);

     // ===============
     //  QUERY METHODS
     // ===============
     //@{ \name Query Methods
     //! \copydoc BaseBOperator::GetDiffOrder
     virtual UInt GetDiffOrder() const
     {
       return ORDER_DIFF;
     }

     //! \copydoc BaseBOperator::GetDimDof()
     virtual UInt GetDimDof() const
     {
       return DIM_DOF;
     }

     //! \copydoc BaseBOperator::GetDimSpace()
     virtual UInt GetDimSpace() const
     {
       return DIM_SPACE;
     }

     //! \copydoc BaseBOperator::GetDimElem()
     virtual UInt GetDimElem() const
     {
       return DIM_ELEM;
     }

     //! \copydoc BaseBOperator::GetDimDMat()
     virtual UInt GetDimDMat() const
     {
       return DIM_D_MAT;
     }
     //@}

   };

   template<class FE, class TYPE>
   void ScaledStrainOperator2p5D<FE,TYPE>::CalcOpMat(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Double> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: nrNodes x spaceDim)
     Matrix<Double> xiDx;
     FE *fe = (static_cast<FE*>(ptFe));
     if (this->useICModes_)
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
     else
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialize with zeros
     // In 2.5D case DIM_DOF = 3 is to be used instead of conventional DIM_SPACE = 2
     bMat.Resize(DIM_D_MAT, numFncs*DIM_DOF);
     bMat.Init();

     // commented entries are zeros in 2.5D case
     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
       bMat[0][pos+0] = xiDx[iFunc][0]*coefs[0];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
       bMat[1][pos+1] = xiDx[iFunc][1]*coefs[1];


     //for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     //  bMat[2][pos+2] = xiDx[iFunc][2]*coefs[2];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
     //  bMat[3][pos+1] = xiDx[iFunc][2]*coefs[2];
       bMat[3][pos+2] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
     //  bMat[4][pos+0] = xiDx[iFunc][2]*coefs[2];
       bMat[4][pos+2] = xiDx[iFunc][0]*coefs[0];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
       bMat[5][pos+0] = xiDx[iFunc][1]*coefs[1];
       bMat[5][pos+1] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator2p5D<FE,TYPE>::CalcOpMatTransposed(Matrix<Double>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Double> coefs;
     this->coef_->GetVector(coefs, lp);
     // Get derivatives of local shape functions with respect to global
     // coords (format: spaceDim x nrNodes )
     Matrix<Double> xiDx;
     FE *fe = (static_cast<FE*>(ptFe));
     // query const variable, should be pretty much optimized away
     if (this->useICModes_)
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
     else
       if (this->isSurfOpt_)
         fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
       else
         fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialise with zeros
     // In 2.5D case DIM_DOF = 3 is to be used instead of conventional DIM_SPACE = 2
     bMat.Resize(numFncs*DIM_DOF, DIM_D_MAT);
     bMat.Init();

     // commented entries are zeros in 2.5D case
     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
       bMat[pos+0][0] = xiDx[iFunc][0]*coefs[0];
     //  bMat[pos+0][4] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+0][5] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
       bMat[pos+1][1] = xiDx[iFunc][1]*coefs[1];
     //  bMat[pos+1][3] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+1][5] = xiDx[iFunc][0]*coefs[0];
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
     //  bMat[pos+2][2] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+2][3] = xiDx[iFunc][1]*coefs[1];
       bMat[pos+2][4] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator2p5D<FE,TYPE>::CalcOpMat(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Complex> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: nrNodes x spaceDim)
     Matrix<Double> xiDx, xiDxTmp, rotMat;
     FE *fe = (static_cast<FE*>(ptFe));

     if (this->coef_->GetCoordinateSystem())
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);

       // If coordinate system is set at the coefficient function, rotate B-matrix
       Vector<Double> globPoint;
       lp.shapeMap->Local2Global(globPoint, lp.lp.coord);
       this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
       xiDx = xiDxTmp*rotMat;
     }
     else
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
     }

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialize with zeros
     // In 2.5D case DIM_DOF = 3 is to be used instead of conventional DIM_SPACE = 2
     bMat.Resize(DIM_D_MAT, numFncs*DIM_DOF);
     bMat.Init();

     // commented entries are zeros in 2.5D case
     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
       bMat[0][pos+0] = xiDx[iFunc][0]*coefs[0];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
       bMat[1][pos+1] = xiDx[iFunc][1]*coefs[1];

     //for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_SPACE)
     //  bMat[2][pos+2] = xiDx[iFunc][2]*coefs[2];

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
     //  bMat[3][pos+1] = xiDx[iFunc][2]*coefs[2];
       bMat[3][pos+2] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
     //  bMat[4][pos+0] = xiDx[iFunc][2]*coefs[2];
       bMat[4][pos+2] = xiDx[iFunc][0]*coefs[0];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
       bMat[5][pos+0] = xiDx[iFunc][1]*coefs[1];
       bMat[5][pos+1] = xiDx[iFunc][0]*coefs[0];
     }
   }

   template<class FE, class TYPE>
   void ScaledStrainOperator2p5D<FE,TYPE>::CalcOpMatTransposed(Matrix<Complex>& bMat, const LocPointMapped& lp, BaseFE* ptFe)
   {
     assert(this->coef_ != NULL);
     Vector<Complex> coefs;
     this->coef_->GetVector(coefs, lp);

     // Get derivatives of local shape functions with respect to global
     // coords (format: spaceDim x nrNodes )
     Matrix<Double> xiDx, xiDxTmp, rotMat;
     FE *fe = (static_cast<FE*>(ptFe));
     // query const variable, should be pretty much optimized away
     if (this->coef_->GetCoordinateSystem())
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDxTmp, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDxTmp, lp, lp.shapeMap->GetElem(), 1);

       // If coordinate system is set at the coefficient function, rotate B-matrix
       Vector<Double> globPoint;
       lp.shapeMap->Local2Global(globPoint, lp.lp.coord);
       this->coef_->GetCoordinateSystem()->GetGlobRotationMatrix(rotMat, globPoint);
       xiDx = xiDxTmp*rotMat;
     }
     else
     {
       if (this->useICModes_)
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFncICModes(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFncICModes(xiDx, lp, lp.shapeMap->GetElem(), 1);
       else
         if (this->isSurfOpt_)
           fe->GetGlobDerivShFnc(xiDx, *lp.lpmVol, lp.lpmVol->shapeMap->GetElem(), 1);
         else
           fe->GetGlobDerivShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
     }

     const UInt numFncs = xiDx.GetNumRows();
     // Set correct size of matrix B and initialise with zeros
     // In 2.5D case DIM_DOF = 3 is to be used instead of conventional DIM_SPACE = 2
     bMat.Resize(numFncs*DIM_DOF, DIM_D_MAT);
     bMat.Init();

     // commented entries are zeros in 2.5D case
     UInt iFunc = 0;
     UInt pos = 0;
     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
       bMat[pos+0][0] = xiDx[iFunc][0]*coefs[0];
     //  bMat[pos+0][4] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+0][5] = xiDx[iFunc][1]*coefs[1];
     }

     for (iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
       bMat[pos+1][1] = xiDx[iFunc][1]*coefs[1];
     //  bMat[pos+1][3] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+1][5] = xiDx[iFunc][0]*coefs[0];
     }

     for(iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos+=DIM_DOF)
     {
     //  bMat[pos+2][2] = xiDx[iFunc][2]*coefs[2];
       bMat[pos+2][3] = xiDx[iFunc][1]*coefs[1];
       bMat[pos+2][4] = xiDx[iFunc][0]*coefs[0];
     }
   }

}
#endif
