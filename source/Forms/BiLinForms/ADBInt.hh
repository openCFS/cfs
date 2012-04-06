

#ifndef FILE_ADBINT_HH
#define FILE_ADBINT_HH

#include <boost/tr1/type_traits.hpp>

#include "BDBInt.hh"

#include "FeBasis/BaseFE.hh"


namespace CoupledField {


  //! General class for calculating coupling integrators of form ADB
  template< class A_OP,
            class B_OP,
            class MAT_DATA_TYPE=Double,
            class COEF_DATA_TYPE=Double >
  class ADBInt : public BDBInt<B_OP, MAT_DATA_TYPE, COEF_DATA_TYPE> {
    public:

      //! Constructor 
      ADBInt(PtrCoefFct dData, MAT_DATA_TYPE factor,
             bool coordUpdate = false );

      //! Destructor
      virtual ~ADBInt();

      //! Compute element matrix associated to ADB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                              EntityIterator& ent1,
                              EntityIterator& ent2 );

      //! Apply A-Trans-Operator on vector
      void ApplyATransMat( Vector<MAT_DATA_TYPE>&ret, 
                           const Vector<MAT_DATA_TYPE>& sol,
                           const LocPointMapped& lpm );

      //! Apply dATrans-Operator on vector
      void ApplydATransMat( Vector<MAT_DATA_TYPE>&ret, 
                            const Vector<MAT_DATA_TYPE>& sol,
                            const LocPointMapped& lpm );
      
      //! Calculate integration kernel, i.e. A*d*B without integration
      void CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
                       const LocPointMapped& lpm );
      
      //! Set Coefficient Function of B operator
      virtual void SetBCoefFunctionOpA(PtrCoefFct coef){
        this->aOperator_.SetCoefFunction(coef);
      }

    protected:
      
      //! First differential operator
      A_OP aOperator_;
  };

}

//Include template definition file
#include "ADBInt.cc"

#endif
