#ifndef FILE_ABINT
#define FILE_ABINT

#include "BBInt.hh"
#include "FeBasis/HCurl/HCurlElemsHi.hh"


namespace CoupledField {
  //! General class for calculation of AB-Forms
  //! \tparam COEF_DAATA_TYPE Data type of the material tensor  
  //! \tparam B_DATA_TYPE Data type of the differential operator
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
    class ABInt : public BBInt<COEF_DATA_TYPE, B_DATA_TYPE> {
  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

    //! Constructor with pointer to BaseElem
    ABInt( BaseBOperator * aOp, BaseBOperator * bOp,
           PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
           bool coordUpdate = false );

    //! Copy Constructor
    ABInt(const ABInt& right):BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(right){
      //here we would also need to create a new operator
      this->aOperator_ = right.aOperator_->Clone();
      this->solDependent_ = right.solDependent_;
      this->aMat_ = right.aMat_;
    }

    //! \copydoc BiLinearForm::Clone
    virtual ABInt* Clone(){
      return new ABInt( *this );
    }

    //! Destructor
    virtual ~ABInt(){
      delete aOperator_;
    }
    
    //! Compute element matrix associated to AB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );

    //! Set Coefficient Function of A operator
    virtual void SetBCoefFunctionOpA(PtrCoefFct coef){
      this->aOperator_->SetCoefFunction(coef);
    }

    //! Set Coefficient Function of B operator
    virtual void SetBCoefFunctionOpB(PtrCoefFct coef){
      this->bOperator_->SetCoefFunction(coef);
    }

    //! \copydoc BiLinearForm::IsSolDependent
    virtual void SetSolDependent() {
      solDependent_ = true;
    }

    //! \copydoc BiLinearForm::IsSolDependent
    virtual bool IsSolDependent() {
      return solDependent_;
    }

  protected:

    //! First differential operator
    BaseBOperator* aOperator_;
    
    //! Store intermediate A-matrix
    Matrix<MAT_DATA_TYPE> aMat_;
    
    bool solDependent_;
  };


  //! general class for calculation of AB-Forms
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
    class ABIntLem : public ABInt<COEF_DATA_TYPE,B_DATA_TYPE>{
  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;
    
    //! Constructor with pointer to CoefFunction for surface itself
    ABIntLem( BaseBOperator * aOp, BaseBOperator * bOp,
                  PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                  bool coordUpdate = false);

    //! Constructor with CoefFunctions for a number of volume regions
    ABIntLem( BaseBOperator * aOp, BaseBOperator * bOp,
                  const std::map< RegionIdType, PtrCoefFct >& regionCoefs,
                  MAT_DATA_TYPE factor,
                  bool coordUpdate = false);

    //! Copy constructor
    ABIntLem(const ABIntLem& right)
      : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(right){
      //here we would also need to create a new operator
      this->regionCoefs_ = right.regionCoefs_;
    }

    //! \copydoc BiLinearForm::Clone
    virtual ABIntLem* Clone(){
      return new ABIntLem( *this );
    }

    //! Destructor
    virtual ~ABIntLem() {}

    //! Compute element matrix associated to AB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );
  protected:
    //! Map containing all coefficient functions for volume regions for operator A
    std::map< RegionIdType, PtrCoefFct > regionCoefs_;
  };


  //! general class for calculation of AB-Forms
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
    class SurfaceABInt : public ABInt<COEF_DATA_TYPE,B_DATA_TYPE>{
  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;
    
    //! Constructor with pointer to CoefFunction for surface itself
    SurfaceABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                  PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                  const std::set<RegionIdType>& volRegions, bool coordUpdate = false);

    //! Constructor with CoefFunctions for a number of volume regions
    SurfaceABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                  const std::map< RegionIdType, PtrCoefFct >& regionCoefs,
                  MAT_DATA_TYPE factor,
                  const std::set<RegionIdType>& volRegions,
                  bool coordUpdate = false);

    //! Copy constructor
    SurfaceABInt(const SurfaceABInt& right)
      : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(right){
      //here we would also need to create a new operator
      this->volRegions_ = right.volRegions_;
      this->regionCoefs_ = right.regionCoefs_;
    }

    //! \copydoc BiLinearForm::Clone
    virtual SurfaceABInt* Clone(){
      return new SurfaceABInt( *this );
    }

    //! Destructor
    virtual ~SurfaceABInt() {}

    //! Compute element matrix associated to AB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );
  protected:
    //! Set containing all volume regions for surface integrators
    std::set<RegionIdType> volRegions_;

    //! Map containing all coefficient functions for volume regions for operator A
    std::map< RegionIdType, PtrCoefFct > regionCoefs_;
  };


  //! class for calculation of bilinear forms for Mortar non-conforming coupling
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
    class SurfaceMortarABInt : public ABInt<COEF_DATA_TYPE,B_DATA_TYPE> {

  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

    //! Constructor with pointer to CoefFunction for surface itself
    SurfaceMortarABInt(BaseBOperator * aOp,
                       BaseBOperator * bOp,
                       PtrCoefFct scalCoef, 
                       MAT_DATA_TYPE factor,
                       RegionIdType masterVolRegion,
                       RegionIdType slaveVolRegion,
                       bool coplanar,
                       bool coordUpdate = false,
                       BiLinearForm::CouplingDirection cplDirection = BiLinearForm::PRIM_SEC);

    //! Copy Constructor
    SurfaceMortarABInt(const SurfaceMortarABInt& right): 
      ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(right) {
      this->ptFeSpaceLM_ = right.ptFeSpaceLM_;
      this->ptFeSpaceField_ = right.ptFeSpaceField_;
      this->primaryVolRegion_ = right.primaryVolRegion_;
      this->secondaryVolRegion_ = right.secondaryVolRegion_;
      this->volRegions_ = right.volRegions_;
      this->isCoplanar_ = right.isCoplanar_;
      this->cplDirection_ = right.cplDirection_;
    }

    //! \copydoc BiLinearForm::Clone
    virtual SurfaceMortarABInt* Clone(){
      return new SurfaceMortarABInt( *this );
    }

    //! Destructor
    virtual ~SurfaceMortarABInt() {};

    //! Compute element matrix associated to AB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );
    
    //! Set finite element space in cases of mixed spaces
    void SetFeSpace( shared_ptr<FeSpace> feSpace1,
                     shared_ptr<FeSpace> feSpace2);

    //! Set coupling direction
    void SetCoupling(BiLinearForm::CouplingDirection cplDir) { cplDirection_ = cplDir; };

  protected:
    
    //! pointer to the FeSpace of the Lagrange multiplier
    shared_ptr<FeSpace> ptFeSpaceLM_ = nullptr;
    
    //! pointer to the FeSpace of the field variable to be coupled
    shared_ptr<FeSpace> ptFeSpaceField_ = nullptr;
    
    // RegionId of volume region on primary side of the ncInterface
    RegionIdType primaryVolRegion_;

    // RegionId of volume region on secondary side of the ncInterface
    RegionIdType secondaryVolRegion_;

    //! Set containing all volume regions for surface integrators
    std::set<RegionIdType> volRegions_;
    
    //! Is the interface coplanar?
    bool isCoplanar_;

    //! Coupling direction indicates the entities the operators are defined on
    BiLinearForm::CouplingDirection cplDirection_;
  };


  //! class for calculation of Nitsche-type Mortar non-conforming coupling
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
    class SurfaceNitscheABInt : public ABInt<COEF_DATA_TYPE,B_DATA_TYPE> {
  public:
    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

    //! Constructor with pointer to CoefFunction for surface itself
    SurfaceNitscheABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                         PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                         BiLinearForm::CouplingDirection cplDir,
                         bool coordUpdate = false, bool isSym = false, bool isPenalty=false);

    //! Copy constructor
    SurfaceNitscheABInt(const SurfaceNitscheABInt& right):
      ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(right) {
      this->cplDirection_ = right.cplDirection_;
      this->isPenalty_ = right.isPenalty_;
      this->volRegions_ = right.volRegions_;
    }

    //! \copydoc BiLinearForm::Clone
    virtual SurfaceNitscheABInt* Clone(){
      return new SurfaceNitscheABInt( *this );
    }

    //! Destructor
    virtual ~SurfaceNitscheABInt() {}

    //! Compute element matrix associated to BDB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );

  protected:
    //! Set containing all volume regions for surface integrators
    BiLinearForm::CouplingDirection cplDirection_;

    //! flag if a penalty term is assembled
    bool isPenalty_;

    //! Set containing all volume regions for surface integrators
    std::set<RegionIdType> volRegions_;

  private:
    //! Sets the passed volume elems, fe functions, and bools according to which part of the Nitsche form is constructed
    //! @param ptSurfPrimary (in) pointer to the primary surface element
    //! @param ptSurfSecondary (in) pointer to the secondary surface element
    //! @param ptVolElem1 (out) pointer to the first volume element
    //! @param ptVolElem2 (out) pointer to the second volume element
    //! @param ptFe1 (out) pointer to the first fe function
    //! @param ptFe1 (out) pointer to the second fe function
    //! @param usePrimary1 (out) bool indicating if the test function sits on the primary (true) or secondary (false) side
    //! @param usePrimary2 (out) bool indicating if the unknown sits on the primary (true) or secondary (false) side
    void SetCouplingDirection(const SurfElem* ptSurfPrimary, const SurfElem* ptSurfSecondary, Elem*& ptVolElem1, Elem*& ptVolElem2, 
                              BaseFE*& ptFe1, BaseFE*& ptFe2, bool& usePrimary1, bool& usePrimary2);

    //! Compute penalty factor for Nitsche coupling
    //! @param ptSurfPrimary (in) pointer to the primary surface element
    //! @param ptSurfSecondary (in) pointer to the secondary surface element
    //! @param ent1 (in) corresponding entity iterator
    //! @param penaltyFactor (out) penalty factor
    void ComputePenalty(const SurfElem* ptSurfPrimary, const SurfElem* ptSurfSecondary, EntityIterator& ent1, MAT_DATA_TYPE& penaltyFactor);
  };


  //! class for calculation of bilinear forms for non-conforming coupling
  //! in case of mechanical-Acoustic coupling or LinFlow-Acoustic coupling
  //! (could be extended to use on other inter-PDE couplings without Lagrange multiplier)
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
    class SurfaceMortarABIntMA : public SurfaceABInt<COEF_DATA_TYPE,B_DATA_TYPE> {
  public:
    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

    //! Constructor with pointer to CoefFunction for surface itself
    //! This one gets a single CoefFunction (scalCoef)
    SurfaceMortarABIntMA(BaseBOperator * aOp, BaseBOperator * bOp,
                         PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                         bool coplanar, bool coordUpdate = false);

    //! Constructor with CoefFunctions for a number of volume regions
    //! This one gets a map of CoefFunctions for individual regions (regionCoefs)
    SurfaceMortarABIntMA(BaseBOperator * aOp, BaseBOperator * bOp,
                         const std::map< RegionIdType, PtrCoefFct >& regionCoefs,
                         MAT_DATA_TYPE factor, bool coplanar, bool coordUpdate = false);

    //! Copy Constructor
    SurfaceMortarABIntMA(const SurfaceMortarABIntMA& right): 
      SurfaceABInt<COEF_DATA_TYPE,B_DATA_TYPE>(right) {
      this->ptPrimaryOp_ = right.ptPrimaryOp_;
      this->ptSecondaryOp_  = right.ptSecondaryOp_;
      this->ptFeSpaceSecondary_ = right.ptFeSpaceSecondary_;
      this->ptFeSpacePrimary_ = right.ptFeSpacePrimary_;
      this->isCoplanar_ = right.isCoplanar_;
      doTranspose_ = right.doTranspose_;
    }

    //! \copydoc BiLinearForm::Clone
    virtual SurfaceMortarABIntMA* Clone() {
      return new SurfaceMortarABIntMA( *this );
    }

    //! Destructor
    virtual ~SurfaceMortarABIntMA(){};

    //! Compute element matrix associated to AB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );

    //! Set finite element space in cases of mixed spaces
    //! Here we use SetFeSpace to assure that the correct coupling direction is set
    //! Hence, the function sets the ptFeSpacePrimary_ and ptFeSpaceSecondary_ pointers
    //! and the ptPrimaryOp_ and ptSecondaryOp_ pointers depending on which PDE is present 
    //! in the aOperator_ and bOperator_ pointers
    //! SetFeSpace is usually called via the SetFeFunctions() function in the BiLinearForm class
    //! when the BiLinearFormContext is defined
    void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2);

  protected:

    //! Pointer to the differential operator for the master side
    BaseBOperator* ptPrimaryOp_ = nullptr;

    //! Pointer to the differential operator for the slave side
    BaseBOperator* ptSecondaryOp_ = nullptr;

    //! pointer to the FeSpace of the secondary side
    shared_ptr<FeSpace> ptFeSpaceSecondary_ = nullptr;

    //! pointer to the FeSpace of the primary side
    shared_ptr<FeSpace> ptFeSpacePrimary_ = nullptr;

    //! Is the interface coplanar?
    bool isCoplanar_;

    //! element matrix needs to be transposed
    bool doTranspose_ = false;
  };
}
#endif
