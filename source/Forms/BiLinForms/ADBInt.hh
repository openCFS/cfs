#ifndef FILE_ADBINT_HH
#define FILE_ADBINT_HH

#include <boost/tr1/type_traits.hpp>

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

    //! Destructor
    virtual ~ADBInt();

    //! Compute element matrix associated to ADB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );

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

    //! Set Coefficient Function of B operator
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
