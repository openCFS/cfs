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
#include "MatVec/promote.hh"
#include "FeBasis/HCurl/HCurlElemsHi.hh"


namespace CoupledField {

  //! general class for calculation of bb forms
  template< class COEF_DATA_TYPE=Double,
            class B_DATA_TYPE=Double >
  class BBInt : public BaseBDBInt {
    public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;
    
      //! Constructor with pointer to BaseElem
      BBInt( BaseBOperator * bOp,
             PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
             bool coordUpdate = false);

      //! Copy Constructor
      BBInt(const BBInt& right)
        : BaseBDBInt(right){
        //here we would also need to create a new operator
        this->bOperator_ = right.bOperator_->Clone();
        this->factor_ = right.factor_;
        this->coefScalar_ = right.coefScalar_;
        this->bMat_  = right.bMat_;
      }

      virtual BBInt* Clone(){
        return new BBInt( *this );
      }


      //! Destructor
      ~BBInt(){
        delete this->bOperator_;

      }
      //! \copydoc BaseBDBInt::GetBOp
      virtual BaseBOperator* GetBOp() {
        return bOperator_;
      }
      
      //! \copydoc BaseBDBInt::GetCoef
      virtual PtrCoefFct GetCoef() {
        return coefScalar_;
      }

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                              EntityIterator& ent1,
                              EntityIterator& ent2 );

      //! Compute element matrix associated to BDB form for a specific lpm
      virtual void CalcElementMatrixLpm( Matrix<MAT_DATA_TYPE>& elemMat,
                                 BaseFE* ptFe,
                                 const LocPointMapped& lp, 
                                 bool overrideIsSurfOpt );

      //@{
      void ApplyElemMat( Vector<Double>&ret, 
                         const Vector<Double>& sol,
                         EntityIterator& ent1,
                         EntityIterator& ent2 );
      
      void ApplyElemMat( Vector<Complex>&ret, 
                         const Vector<Complex>& sol,
                         EntityIterator& ent1,
                         EntityIterator& ent2 );

      //@}


      //! Calculate integration kernel, i.e. B*d*B without integration
      void CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
                       const LocPointMapped& lpm );

      //@{
      void ApplyBMat( Vector<Double>&ret, 
                      const Vector<Double>& sol,
                      const LocPointMapped& lpm );
      
      void ApplyBMat( Vector<Complex>&ret, 
                      const Vector<Complex>& sol,
                      const LocPointMapped& lpm );
      //@}

      //@{
      void ApplydBMat( Vector<Double>&ret, 
                       const Vector<Double>& sol,
                       const LocPointMapped& lpm );
      
      void ApplydBMat( Vector<Complex>&ret, 
                       const Vector<Complex>& sol,
                       const LocPointMapped& lpm );
      //@}
      
      //@{
      void ApplyATransMat( Vector<Double>&ret, 
                           const Vector<Double>& sol,
                           const LocPointMapped& lpm ) {
        EXCEPTION("Not implemented");
      }
      void ApplyATransMat( Vector<Complex>&ret, 
                           const Vector<Complex>& sol,
                           const LocPointMapped& lpm ) {
        EXCEPTION("Not implemented");
      }
      //@}

      //@{
      void ApplydATransMat( Vector<Double>&ret, 
                            const Vector<Double>& sol,
                            const LocPointMapped& lpm ) {
        EXCEPTION("Not implemented");
      }
      void ApplydATransMat( Vector<Complex>&ret, 
                            const Vector<Complex>& sol,
                            const LocPointMapped& lpm ) {
        EXCEPTION("Not implemented");
      }
      //@}

      bool IsComplex() const {
        return std::is_same<MAT_DATA_TYPE,Complex>::value;
      }
      
      //! \copydoc BiLinearForm::IsSolDependent
      virtual bool IsSolDependent() {
        return (coefScalar_->GetDependency() == CoefFunction::SOLUTION) || (coefScalar_->GetDependency() == CoefFunction::SPACE);
      }
            
      //! Set Finite Element Space
      void SetFeSpace( shared_ptr<FeSpace> feSpace ) {
        this->ptFeSpace1_ = feSpace;
        this->intScheme_ = ptFeSpace1_->GetIntScheme();
      }

      //! Set finite element space in cases of mixed spaces
      virtual void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2) {
        this->ptFeSpace1_ = feSpace1;
        this->ptFeSpace2_ = feSpace2;
        this->intScheme_ = ptFeSpace1_->GetIntScheme();
      }
      //! Set Coefficient Function of B operator
      virtual void SetBCoefFunctionOpB(PtrCoefFct coef){
        this->bOperator_->SetCoefFunction(coef);
      }

    protected:
      
      //! Differential operator
      BaseBOperator * bOperator_;

      //! Store intermediate operator matrix for B
      Matrix<MAT_DATA_TYPE> bMat_;
      
      //! A constant factor for multiplication with the element matrix
      MAT_DATA_TYPE factor_;

      //! Pointer to coefficient function for scalar values
      shared_ptr<CoefFunction > coefScalar_;
  };
  
  // ========================================================================

  //! Specialized class for calculating the mass bilienarform for edge elements
  template<class COEF_DATA_TYPE=Double,
      class B_DATA_TYPE=Double >
  class BBIntMassEdge : public BBInt<COEF_DATA_TYPE, B_DATA_TYPE> {
    public:

    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;
    
      //! Constructor with pointer to BaseElem
      BBIntMassEdge(BaseBOperator * bOp,
                    PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                    bool coordUpdate = false):
        BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(bOp, scalCoef, factor, coordUpdate ){
        this->name_ = "BBIntMassEdge";
      }

      BBIntMassEdge(const BBIntMassEdge & right)
      : BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(right){

      }

      virtual BBIntMassEdge* Clone(){
        return new BBIntMassEdge( *this );
      }

      //! Destructor
      virtual ~BBIntMassEdge(){

      }

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );
  };

  
  // ========================================================================
  //! general class for calculation of bb forms on surfaces
  template<class COEF_DATA_TYPE=Double,
        class B_DATA_TYPE=Double >
  class SurfaceBBInt : public BBInt<COEF_DATA_TYPE, B_DATA_TYPE> {
    public:
    
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;

      //! Constructor with pointer to BaseElem
    SurfaceBBInt(BaseBOperator * bOp,
                 PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                 const std::set<RegionIdType>& volRegions, bool coordUpdate = false):
        BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(bOp, scalCoef, factor, coordUpdate ){
        this->name_ = "SurfaceBBInt";
        volRegions_ = volRegions;
        this->isSymmetric_ = true;
      }

     SurfaceBBInt(const SurfaceBBInt & right)
     : BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(right){
        this->volRegions_ = right.volRegions_;
      }

    virtual SurfaceBBInt* Clone(){
      return new SurfaceBBInt( *this );
    }

      //! Destructor
      virtual ~SurfaceBBInt(){

      }

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );

    protected:
      //! set containing all volume regions for surface integrators
      std::set<RegionIdType> volRegions_;
      
      Matrix<MAT_DATA_TYPE> bMatT_;
  };
}

#endif
