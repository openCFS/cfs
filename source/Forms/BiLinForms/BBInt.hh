// =====================================================================================
// 
//       Filename:  BBInt.hh
// 
//    Description:  This class implements a general symmetric integrator without any 
//                  material factors. The other opprtunity, to pass a id-Matrix to the
//                  general BDB-Integrator is not preferable due to overhead in
//                  computing an identity matrix of correct size in the coefficient
//                  function and an unnecessay matrix-matrix calculation
// 
//        Version:  1.0
//        Created:  10/29/2011 02:39:39 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef FILE_BBINT
#define FILE_BBINT

#include "BDBInt.hh"
#include "FeBasis/BaseFE.hh"
#include <boost/tr1/type_traits.hpp>
#include "FeBasis/HCurl/HCurlElemsHi.hh"


namespace CoupledField {

  //! general class for calculation of bb forms
  template<class B_OP, 
           class MAT_DATA_TYPE=Double,
           class COEF_DATA_TYPE=Double>
  class BBInt : public BaseBDBInt {
    public:

      //! Constructor with pointer to BaseElem
      BBInt( shared_ptr<CoefFunction> scalCoef, MAT_DATA_TYPE factor,
             bool coordUpdate = false);

      //! Destructor
      ~BBInt(){

      }
      
      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );

      //! Multiply element matrix with vector
      void ApplyElemMat( Vector<MAT_DATA_TYPE>&ret, 
                         const Vector<Double>& sol,
                         EntityIterator& ent1,
                         EntityIterator& ent2 );


      //! Calculate integration kernel, i.e. B*d*B without integration
      void CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
                       const LocPointMapped& lpm );

      //! Apply B-Operator on vector 
      void ApplyBMat( Vector<MAT_DATA_TYPE>&ret, 
                      const Vector<MAT_DATA_TYPE>& sol,
                      const LocPointMapped& lpm );

      //! Apply dB-Operator on vector
      void ApplydBMat( Vector<MAT_DATA_TYPE>&ret, 
                       const Vector<MAT_DATA_TYPE>& sol,
                       const LocPointMapped& lpm );

      bool IsComplex(){
        return std::tr1::is_same<MAT_DATA_TYPE,Complex>::value;
      }
      
      //! \copydoc BiLinearForm::IsSolDependent
      bool IsSolDependent() {
        return coefScalar_->GetDependency() == CoefFunction::SOLUTION;
      }
            
      void SetFeSpace( shared_ptr<FeSpace> feSpace ) {
        this->ptFeSpace1_ = feSpace;
        this->intScheme_ = ptFeSpace1_->GetIntScheme();
      }

      virtual void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2) {
        this->ptFeSpace1_ = feSpace1;
        this->ptFeSpace2_ = feSpace2;
        this->intScheme_ = ptFeSpace1_->GetIntScheme();
      }
      //! Set Coefficient Function of B operator
      virtual void SetBCoefFunctionOpB(shared_ptr<CoefFunction> coef){
        this->bOperator_.SetCoefFunction(coef);
      }

    protected:
      
      //! Differential operator
      B_OP bOperator_;

      //! A constant factor for multiplication with the element matrix
      MAT_DATA_TYPE factor_;

      
      //! Pointer to coefficient function for scalar values
      shared_ptr<CoefFunction > coefScalar_;
  };

  //! general class for calculation of bb forms
  template<class B_OP, class MAT_DATA_TYPE=Double>
  class BBIntMassEdge : public BBInt<B_OP, MAT_DATA_TYPE> {
    public:

      //! Constructor with pointer to BaseElem
      BBIntMassEdge(shared_ptr<CoefFunction> scalCoef, MAT_DATA_TYPE factor,
                    bool coordUpdate = false):
        BBInt<B_OP,MAT_DATA_TYPE>(scalCoef, factor, coordUpdate ){
        this->name_ = "BBIntMassEdge";
      }

      //! Destructor
      ~BBIntMassEdge(){

      }

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );
  };

  //! general class for calculation of bb forms on surfaces
  template<class B_OP, class MAT_DATA_TYPE=Double,
      class COEF_DATA_TYPE=Double>
  class SurfaceBBInt : public BBInt<B_OP, MAT_DATA_TYPE,COEF_DATA_TYPE> {
    public:

      //! Constructor with pointer to BaseElem
    SurfaceBBInt(shared_ptr<CoefFunction> scalCoef, MAT_DATA_TYPE factor,
                 const std::set<RegionIdType>& volRegions, bool coordUpdate = false):
        BBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>(scalCoef, factor, coordUpdate ){
        this->name_ = "SurfaceBBInt";
        volRegions_ = volRegions;
        this->isSymmetric_ = false;
      }

      //! Destructor
      ~SurfaceBBInt(){

      }

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );
    protected:
      //! set containing all volume regions for surface integrators
      std::set<RegionIdType> volRegions_;
  };
}

//Include template definition file
#include "BBInt.cc"
#endif
