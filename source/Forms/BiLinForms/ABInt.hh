#ifndef FILE_ABINT
#define FILE_ABINT

#include "BBInt.hh"
#include "FeBasis/HCurl/HCurlElemsHi.hh"


namespace CoupledField {

  //! general class for calculation of AB-Forms
  template<class A_OP, class B_OP, class MAT_DATA_TYPE=Double>
  class ABInt : public BBInt<B_OP, MAT_DATA_TYPE,MAT_DATA_TYPE> {
    public:

      //! Constructor with pointer to BaseElem
      ABInt( PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
             bool coordUpdate = false );

      //! Destructor
      virtual ~ABInt(){

      }
      //! Compute element matrix associated to AB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );

      //! Set Coefficient Function of B operator
      virtual void SetBCoefFunctionOpA(PtrCoefFct coef){
        this->aOperator_.SetCoefFunction(coef);
      }

    protected:
      
      //! First differential operator
      A_OP aOperator_;
  };

  //! general class for calculation of AB-Forms
  template<class A_OP, class B_OP, class MAT_DATA_TYPE=Double>
  class SurfaceABInt : public ABInt<A_OP,B_OP, MAT_DATA_TYPE>{
  public:
      //! Constructor with pointer to BaseElem
      SurfaceABInt(PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                   const std::set<RegionIdType>& volRegions, bool coordUpdate = false);

      //! Constructor with CoefFunctions for a number of volume regions
      SurfaceABInt(const std::map< RegionIdType, PtrCoefFct >& regionCoefs,
                   MAT_DATA_TYPE factor,
                   const std::set<RegionIdType>& volRegions,
                   bool coordUpdate = false);

      //! Destructor
      ~SurfaceABInt() {}

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );
    protected:
      //! Set containing all volume regions for surface integrators
      std::set<RegionIdType> volRegions_;

      //! Map containing all coefficient functions for volume regions for operator A
      std::map< RegionIdType, PtrCoefFct > regionCoefs_;    
  };
}

//Include template definition file
#include "ABInt.cc"
#endif
