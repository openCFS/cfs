// =====================================================================================
// 
//       Filename:  bdbInt.hh
// 
//    Description:  New implementation of the BDB integrator class
//                  Takes as a template parameter the operator it should evaluate
//                  new implementation to avoid old structures
// 
//        Version:  1.0
//        Created:  10/04/2011 09:28:38 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef FILE_BDBINT_NEW
#define FILE_BDBINT_NEW

#include "BiLinearForm.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "MatVec/promote.hh"
#include "FeBasis/BaseFE.hh"
#include "Domain/Domain.hh"

namespace CoupledField {


  //! Base class for BDB integrator

  //! This is the base class for all BDB-integrators. It mainly serves as
  //! a common hook to store BDB-integrators, independent of the type of
  //! finite element functions and the coefficient data, by which the derived
  //! classes are templatized.
  //! This interface mainly provides abstract methods to apply the B-operator
  //! and the dB-operator to vectors, which is used e.g. in postprocessing to
  //! evaluate flux and energy quantities.
class BaseBDBInt : public BiLinearForm {
public:

  //! Constructor
  BaseBDBInt( bool coordUpate = false) :
    BiLinearForm(coordUpate) {

  }

  //! Copy constructor
  BaseBDBInt(const BaseBDBInt& right) :
    BiLinearForm(right){
  }

  //! \copydoc BiLinearForm::Clone
  virtual BaseBDBInt* Clone()=0;

  //! Destructor
  virtual ~BaseBDBInt() {
  }

  //! Obtain differential operator
  virtual BaseBOperator* GetBOp() = 0;
  
  //! Obtain coefficient function
  virtual PtrCoefFct GetCoef() = 0;
  
  //@{
  //! Apply / multiply element matrix to vector

  //! This method multiplies the element matrix with a vector. This can be used e.g.
  //! to calculate the energy of a field on the element level. 
  virtual void ApplyElemMat( Vector<Double>&ret, const Vector<Double>& sol,
                             EntityIterator& ent1,
                             EntityIterator& ent2 ) = 0;
  virtual void ApplyElemMat( Vector<Complex>&ret, const Vector<Complex>& sol,
                             EntityIterator& ent1,
                             EntityIterator& ent2 ) = 0;
  //@}

  //@{
  //! Apply B-Operator on vector
  virtual void ApplyBMat( Vector<Double>&ret, const Vector<Double>& sol,
                          const LocPointMapped& lpm ) = 0;
  virtual void ApplyBMat( Vector<Complex>&ret, const Vector<Complex>& sol,
                          const LocPointMapped& lpm )  = 0;
  //@}


  //@{
  //! Apply dB-Operator on vector
  virtual void ApplydBMat( Vector<Double>&ret, const Vector<Double>& sol,
                           const LocPointMapped& lpm ) = 0;

  virtual void ApplydBMat( Vector<Complex>&ret, const Vector<Complex>& sol,
                           const LocPointMapped& lpm ) = 0;
  //@}
  
  //! Apply transposed A-Operator on vector
  
  //! This method is only properly overwritten for the derived ADB-Integrator
  //! and allows to apply the transposed A-Matrix (i.e. first differential
  //! operator) to a vector.
   virtual void ApplyATransMat( Vector<Double>&ret, const Vector<Double>& sol,
                                const LocPointMapped& lpm ) = 0;
   
   virtual void ApplyATransMat( Vector<Complex>&ret, const Vector<Complex>& sol,
                           const LocPointMapped& lpm ) = 0;
   //@}


   //@{
   //! Apply dA-Operator on vector
   virtual void ApplydATransMat( Vector<Double>&ret, const Vector<Double>& sol,
                                 const LocPointMapped& lpm ) = 0;

   virtual void ApplydATransMat( Vector<Complex>&ret, const Vector<Complex>& sol,
                                 const LocPointMapped& lpm ) = 0;
   //@}


  //@{
  //! Calculate integration kernel, i.e. B*d*B without integration
  virtual void CalcKernel( Matrix<Double>& kernel,
                           const LocPointMapped& lpm ) {
    EXCEPTION("BaseBDBInt: CalcKernel not implemented");
  }

  virtual void CalcKernel( Matrix<Complex>& kernel,
                           const LocPointMapped& lpm ) {
    EXCEPTION("BaseBDBInt: CalcKernel not implemented");
  }
  //@}

  
};


  //! General class for calculating BdB bilinearforms
  
  //! This class calculates the general BdB-form, where B
  //! denotes an arbitrary differential operator and d
  //! denotes the material tensor
  //! \tparam COEF_DATA_TYPE Data type of the material tensor  
  //! \tparam B_DATA_TYPE Data type of the differential operator
  template<class COEF_DATA_TYPE=Double, class B_DATA_TYPE=Double>
  class BDBInt : public BaseBDBInt {
  public:

    //! Define data type for matrix entries, derived by type trait
    typedef PROMOTE(B_DATA_TYPE, COEF_DATA_TYPE) MAT_DATA_TYPE;
    
    //! Constructor with pointer to BaseElem
    BDBInt( BaseBOperator * bOp,
            PtrCoefFct dData, MAT_DATA_TYPE factor,
            bool coordUpdate = false );

    //! Copy constructor
    BDBInt(const BDBInt& right)
      : BaseBDBInt(right)
    {
      //here we would also need to create a new operator
      this->bOperator_ = right.bOperator_->Clone();
      this->factor_ = right.factor_;
      this->dData_ = right.dData_;
    }

    //! \copydoc BiLinearForm::Clone
    virtual BDBInt* Clone(){
      return new BDBInt( *this );
    }

    //! Destructor
    virtual ~BDBInt();

    //! \copydoc BaseBDBInt::GetBOp
    virtual BaseBOperator* GetBOp() {
      return bOperator_;
    }

    //! \copydoc BaseBDBInt::GetCoef
    virtual PtrCoefFct GetCoef() {
      return dData_;
    }

    //! Compute element matrix associated to BDB form
    void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );

    //! Compute element matrix associated to BDB form for a specific lpm
    virtual void CalcElementMatrixLpm( Matrix<MAT_DATA_TYPE>& elemMat,
                            BaseFE* ptFe,
                            const LocPointMapped& lp );

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
    virtual void ApplydBMat( Vector<Double>&ret, const Vector<Double>& sol,
                             const LocPointMapped& lpm );
    virtual void ApplydBMat( Vector<Complex>&ret, const Vector<Complex>& sol,
                             const LocPointMapped& lpm );
    //@}


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


    bool IsComplex() const {
      return std::is_same<MAT_DATA_TYPE,Complex>::value;
    }

    //! \copydoc BiLinearForm::IsSolDependent
    virtual bool IsSolDependent() {
      return isSolDependent_;
    }

    void SetFeSpace( shared_ptr<FeSpace> feSpace ) {
      this->ptFeSpace1_ = feSpace;
      this->ptFeSpace2_ = feSpace;
      this->intScheme_ = ptFeSpace1_->GetIntScheme();
    }

    virtual void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2) {
      this->ptFeSpace1_ = feSpace1;
      this->ptFeSpace2_ = feSpace2;
      this->intScheme_ = ptFeSpace1_->GetIntScheme();
    }


    //! Set Coefficient Function of B operator
    virtual void SetBCoefFunctionOpA(PtrCoefFct coef){
      this->bOperator_->SetCoefFunction(coef);
    }
    void __Instantiate();


  protected:

    //! Differential operator
    BaseBOperator* bOperator_;

    //! set a constant factor for multiplication with the element matrix
    MAT_DATA_TYPE factor_;

    //! Pointer to coefficient function computing the d-matrix of the BDB Integrator
    shared_ptr<CoefFunction> dData_;

    //! Store intermediate operator matrix for B
    Matrix<MAT_DATA_TYPE> bMat_;

    //! Store intermediate material matrix c
    Matrix<MAT_DATA_TYPE> dMat_;

    //! Store intermediate matrix
    Matrix<MAT_DATA_TYPE> dbMat_;
      
  };

}

#endif
