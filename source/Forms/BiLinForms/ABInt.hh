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
    ABInt(const ABInt& right)
     : BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(right){
      //here we would also need to create a new operator
      this->aOperator_ = right.aOperator_->Clone();
      this->solDependent_ = right.solDependent_;
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

    //! Set Coefficient Function of B operator
    virtual void SetBCoefFunctionOpA(PtrCoefFct coef){
      this->aOperator_->SetCoefFunction(coef);
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

  //! class for calculation of bilinear forms for non-conforming coupling
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
  class SurfaceMortarABInt : public ABInt<COEF_DATA_TYPE,B_DATA_TYPE> {

  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

    //! Constructor with pointer to CoefFunction for surface itself
    SurfaceMortarABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                        PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                        RegionIdType masterVolRegion,
                        RegionIdType slaveVolRegion,
                        bool coplanar,
                        bool coordUpdate = false,
                        BiLinearForm::CouplingDirection cplDirection = BiLinearForm::MASTER_SLAVE);

    //! Copy Constructor
    SurfaceMortarABInt(const SurfaceMortarABInt& right)
     : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(right){
      this->ptFeSpaceLM_ = right.ptFeSpaceLM_;
      this->ptFeSpaceField_ = right.ptFeSpaceField_;
      this->masterVolRegion_ = right.masterVolRegion_;
      this->slaveVolRegion_ = right.slaveVolRegion_;
      this->volRegions_ = right.volRegions_;
      this->isCoplanar_ = right.isCoplanar_;
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
    shared_ptr<FeSpace> ptFeSpaceLM_;
    
    //! pointer to the FeSpace of the field variable to be coupled
    shared_ptr<FeSpace> ptFeSpaceField_;
    
    // RegionId of volume region on master side of the ncInterface
    RegionIdType masterVolRegion_;

    // RegionId of volume region on slave side of the ncInterface
    RegionIdType slaveVolRegion_;

    //! Set containing all volume regions for surface integrators
    std::set<RegionIdType> volRegions_;
    
    //! Is the interface coplanar?
    bool isCoplanar_;

    //! Coupling direction indicates the entities the operators are defined on
    BiLinearForm::CouplingDirection cplDirection_;
  };


  //! general class for calculation of AB-Forms
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
  class SurfaceNitscheABInt : public ABInt<COEF_DATA_TYPE,B_DATA_TYPE>{
    
  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

    //! Constructor with pointer to CoefFunction for surface itself
    SurfaceNitscheABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                         PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                         BiLinearForm::CouplingDirection cplDir,
                         bool coordUpdate = false, bool isSym = false, bool isPenalty=false);

    //! Copy constructor
    SurfaceNitscheABInt(const SurfaceNitscheABInt& right)
     : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(right){
      this->myDirection_ = right.myDirection_;
      this->isPenalty_ = right.isPenalty_;
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
    
    //! Returns the two volume elements necessary to calculate the operators
    void GetVolFromSurfElem(bool & uMaster1, bool & uMaster2);

    //! Set containing all volume regions for surface integrators
    BiLinearForm::CouplingDirection myDirection_;

    bool isPenalty_;

    //TODO: for future purpose it would be helpful to be able
    //to add two and more volume regions to one master or slave side
  };
}

#endif
