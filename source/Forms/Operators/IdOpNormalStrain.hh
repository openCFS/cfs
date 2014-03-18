#ifndef FILE_ID_OP_NORMAL_STRAIN_HH
#define FILE_ID_OP_NORMAL_STRAIN_HH

#include "BaseBOperator.hh"

namespace CoupledField{

//! Strain-like normal mapping operator

//! Perform the normal projection of tensorial quantity in Voigt notation:
//! 
//!    / N_1*n_x    0       0    ... \
//!    |    0    N_1*n_y    0    ... |
//! b =|    0       0    N_1*n_z ... |
//!    |    0    N_1*n_z N_1*n_y ... |
//!    | N_1*n_z    0    N_1*n_x ... |     
//!    \ N_1*n_y N_1*n_x    0    ... /
//!
//!  Here n denotes the normal vector (n_x, n_y, n_z). 
//!  shape function at a given local point
//! \tparam FE Type of Finite Element used
//! \tparam DIM Dimension of the problem space
//! \tparam TYPE Data type (DOUBLE, COMPLEX)
template<class FE, UInt DIM, class TYPE = Double >
  class IdNormalStrainOperator : public BaseBOperator{

    IdNormalStrainOperator();
    ~IdNormalStrainOperator();
  };

  // ===========================
  //  Specialization in 2D  
  // ===========================
  //! Strain-like normal mapping operator for 2D
  
  //! Perform the normal projection of tensorial quantity in Voigt notation
  //! for the 2D-case:
  //! 
  //!     / N_1*n_x    0     ... \
  //! b = |    0    N_1*n_y  ... |
  //!     \ N_1*n_y N_1*n_x  ... /
  //!
  //!  Here n denotes the normal vector (n_x, n_y). Note that the operator
  //!  is only defined on surface elements 
  //! \tparam FE Type of Finite Element used
  //! \tparam D Dimension of the problem space
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)

  template<class FE, class TYPE>
  class IdNormalStrainOperator<FE, 2, TYPE> : public BaseBOperator{

  public:

    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 0;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 2;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = 2;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = 1; // only defined on surfaces

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = 3; 
    //@}

    //! Constructor
    //! \param useIcModes Use incompatible modes shape functions
    IdNormalStrainOperator( ) {
      this->name_ = "IdNormalStrainOperator";
    }

    //! Destructor
    virtual ~IdNormalStrainOperator(){

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
    
  };

  template<class FE, class TYPE>
  void IdNormalStrainOperator<FE, 2, TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                                    const LocPointMapped& lp, 
                                                    BaseFE* ptFe ) {

    // Set correct size of matrix B and initialize with zeros
    const UInt numFncs = ptFe->GetNumFncs();
    bMat.Resize( DIM_D_MAT, numFncs * DIM_SPACE );
    bMat.Init();

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    
    UInt iFunc = 0;
    UInt pos = 0;
    fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 0 );
    for( ; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) {
      bMat[0][pos+0] = s[iFunc] * lp.normal[0];
      bMat[2][pos+0] = s[iFunc] * lp.normal[1];
    }

    fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 1 );
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) {
      bMat[1][pos+1] = s[iFunc] * lp.normal[1];
      bMat[2][pos+1] = s[iFunc] * lp.normal[0];
    }
  }

  template<class FE, class TYPE>
  void IdNormalStrainOperator<FE, 2, TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                              const LocPointMapped& lp, 
                                                              BaseFE* ptFe ) {

    // Set correct size of matrix B and initialise with zeros
    const UInt numFncs = ptFe->GetNumFncs();
    bMat.Resize(numFncs * DIM_SPACE , DIM_D_MAT );
    bMat.Init();

    Vector<Double> s;
    FE *fe = (static_cast<FE*>(ptFe));
    
    UInt iFunc = 0;
    UInt pos = 0;
    fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 0 );
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) {
      bMat[pos+0][0] = s[iFunc] * lp.normal[0];
      bMat[pos+0][2] = s[iFunc] * lp.normal[1];
    }
    
    fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 1 );
    for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) {
      bMat[pos+1][1] = s[iFunc] * lp.normal[1];
      bMat[pos+1][2] = s[iFunc] * lp.normal[0];
    }
  }

  
  // ============================
  //  3D NORMAL STRAIN OPERATOR 
  // ============================
  //! Strain-like normal mapping operator for 3D
  
  //! Perform the normal projection of tensorial quantity in Voigt notation
  //! for the 2D-case:
  //! 
  //!    / N_1*n_x    0       0    ... \
  //!    |    0    N_1*n_y    0    ... |
  //! b =|    0       0    N_1*n_z ... |
  //!    |    0    N_1*n_z N_1*n_y ... |
  //!    | N_1*n_z    0    N_1*n_x ... |     
  //!    \ N_1*n_y N_1*n_x    0    ... /
  //!
  //!  Here n denotes the normal vector (n_x, n_y, n_z). 
  //!  shape function at a given local point
  //! \tparam FE Type of Finite Element used
  //! \tparam D Dimension of the problem space
  //! \tparam TYPE Data type (DOUBLE, COMPLEX)
   template<class FE, class TYPE>
   class IdNormalStrainOperator<FE, 3, TYPE> : public BaseBOperator  {

   public:

     // ------------------
     //  STATIC CONSTANTS 
     // ------------------
     //@{ 
     //! \name Static constants

     //! Order of differentiation
     static const UInt ORDER_DIFF = 0;

     //! Number of components of the problem (scalar, vector)
     static const UInt DIM_DOF = 3;

     //! Dimension of the underlying domain / space
     static const UInt DIM_SPACE = 3;

     //! Dimension of the finite element
     static const UInt DIM_ELEM = 2; // only defined on surfaces

     //! Dimension of the related material 
     static const UInt DIM_D_MAT = 6; 
     //@}

     //! Constructor
     //! \param useIcModes Use incompatible modes shape functions
     IdNormalStrainOperator( ) {
       this->name_ = "IdNormalStrainOperator";
     }

     //! Destructor
     ~IdNormalStrainOperator(){

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

   };

   template<class FE, class TYPE>
   void IdNormalStrainOperator<FE,3,TYPE>::CalcOpMat(Matrix<Double> & bMat,
                                                   const LocPointMapped& lp, 
                                                   BaseFE* ptFe ) {
     
     // Set correct size of matrix B and initialize with zeros
     const UInt numFncs = ptFe->GetNumFncs();
     bMat.Resize( DIM_D_MAT, numFncs * DIM_SPACE );
     bMat.Init();

     Vector<Double> s;
     FE *fe = (static_cast<FE*>(ptFe));
     
     UInt iFunc = 0;
     UInt pos = 0;
     fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 0 );
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) {
       bMat[0][pos+0] = s[iFunc] * lp.normal[0];
       bMat[4][pos+0] = s[iFunc] * lp.normal[2];
       bMat[5][pos+0] = s[iFunc] * lp.normal[1];
     }

     fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 1 );
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) { 
       bMat[1][pos+1] = s[iFunc] * lp.normal[1];
       bMat[3][pos+1] = s[iFunc] * lp.normal[2];
       bMat[5][pos+1] = s[iFunc] * lp.normal[0];
     }
     
     fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 2 );
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) { 
       bMat[2][pos+2] = s[iFunc] * lp.normal[2];
       bMat[3][pos+2] = s[iFunc] * lp.normal[1];
       bMat[4][pos+2] = s[iFunc] * lp.normal[0];
     }
     
   }

   template<class FE, class TYPE>
   void IdNormalStrainOperator<FE,3,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,
                                                             const LocPointMapped& lp, 
                                                             BaseFE* ptFe ) {
     // Set correct size of matrix B and initialise with zeros
     const UInt numFncs = ptFe->GetNumFncs();
     bMat.Resize(numFncs * DIM_SPACE , DIM_D_MAT );
     bMat.Init();
     
     Vector<Double> s;
     FE *fe = (static_cast<FE*>(ptFe));
     
     UInt iFunc = 0;
     UInt pos = 0;
     fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 0 );
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) {
       bMat[pos+0][0] = s[iFunc] * lp.normal[0];
       bMat[pos+0][4] = s[iFunc] * lp.normal[2];
       bMat[pos+0][5] = s[iFunc] * lp.normal[1];
     }

     fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 1 );
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) {
       bMat[pos+1][1] = s[iFunc] * lp.normal[1];
       bMat[pos+1][3] = s[iFunc] * lp.normal[2];
       bMat[pos+1][5] = s[iFunc] * lp.normal[0];
     }
     
     fe->GetShFnc( s, lp.lp, lp.shapeMap->GetElem(), 2 );
     for( iFunc = 0, pos = 0; iFunc < numFncs; iFunc++, pos += DIM_SPACE ) { 
       bMat[pos+2][2] = s[iFunc] * lp.normal[2];
       bMat[pos+2][3] = s[iFunc] * lp.normal[1];
       bMat[pos+2][4] = s[iFunc] * lp.normal[0];
     }
   }
}
#endif
