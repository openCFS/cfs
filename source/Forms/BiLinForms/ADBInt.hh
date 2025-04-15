#ifndef FILE_ADBINT_HH
#define FILE_ADBINT_HH


#include "BDBInt.hh"
#include "FeBasis/BaseFE.hh"


namespace CoupledField {


  //! General class for calculating coupling integrators of form ADB
  //! \tparam COEF_DATA_TYPE Data type of the material tensor  
  //! \tparam B_DATA_TYPE Data type of the differential operator
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
  class ADBInt : public BDBInt<COEF_DATA_TYPE, B_DATA_TYPE> {
  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

    //! Constructor 
    ADBInt( BaseBOperator * aOp, BaseBOperator * bOp,
            PtrCoefFct dData, MAT_DATA_TYPE factor,
            bool coordUpdate = false );

    //! Copy COnstructor
    ADBInt(const ADBInt& right)
     : BDBInt<COEF_DATA_TYPE, B_DATA_TYPE>(right){
      //here we would also need to create a new operator
      this->aOperator_ = right.aOperator_->Clone();
    }

    //! \copydoc BiLinearForm::Clone
    virtual ADBInt* Clone(){
      return new ADBInt( *this );
    }


    //! Destructor
    virtual ~ADBInt();

    //! Compute element matrix associated to ADB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );

    //! Compute element matrix associated to BDB form for a specific lpm
    void CalcElementMatrixLpm( Matrix<MAT_DATA_TYPE>& elemMat,
                               BaseFE* ptFe,
                               const LocPointMapped& lp, 
                               bool overrideIsSurfOpt ){
      EXCEPTION("CalcElementMatrixLpm is not implemented for ADB type integrators!");
    };

    //@{
    void ApplyATransMat( Vector<Double>&ret, 
                         const Vector<Double>& sol,
                         const LocPointMapped& lpm );

    void ApplyATransMat( Vector<Complex>&ret, 
                         const Vector<Complex>& sol,
                         const LocPointMapped& lpm );
    //@}

    //@{
    void ApplydATransMat( Vector<Double>&ret, 
                          const Vector<Double>& sol,
                          const LocPointMapped& lpm );
    void ApplydATransMat( Vector<Complex>&ret, 
                          const Vector<Complex>& sol,
                          const LocPointMapped& lpm );
    //@}

    //! Calculate integration kernel, i.e. A*d*B without integration
    void CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
                     const LocPointMapped& lpm );

    //! Set Coefficient Function of A operator
    virtual void SetBCoefFunctionOpA(PtrCoefFct coef){
      this->aOperator_->SetCoefFunction(coef);
    }

  protected:

    //! First differential operator
    BaseBOperator* aOperator_;

    //! Auxiliary matrices
    Matrix<MAT_DATA_TYPE> aMat_, dAMat_;
  };

}
#endif
