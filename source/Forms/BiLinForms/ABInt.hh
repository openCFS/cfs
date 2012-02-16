#ifndef FILE_ABINT
#define FILE_ABINT

#include "BBInt.hh"
#include "FeBasis/HCurl/HCurlElemsHi.hh"


namespace CoupledField {

  //! general class for calculation of AB-Forms
  template<class A_OP, class B_OP, class MAT_DATA_TYPE=Double>
  class ABInt : public BBInt<B_OP, MAT_DATA_TYPE> {
    public:

      //! Constructor with pointer to BaseElem
      ABInt( shared_ptr<CoefFunction> scalCoef, MAT_DATA_TYPE factor,
             bool coordUpdate = false );

      //! Destructor
      ~ABInt(){

      }
      //! Compute element matrix associated to AB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );

      //! Set Coefficient Function of B operator
      virtual void SetBCoefFunctionOpA(shared_ptr<CoefFunction> coef){
        this->aOperator_.SetCoefFunction(coef);
      }

    protected:
      
      //! First differential operator
      A_OP aOperator_;
  };
}

//Include template definition file
#include "ABInt.cc"
#endif
