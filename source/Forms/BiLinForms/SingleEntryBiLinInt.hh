#ifndef CFS_SINGLEENTRY_BILIN_INT_HH
#define CFS_SINGLEENTRY_BILIN_INT_HH

#include "BiLinearForm.hh" 
#include "MatVec/Vector.hh"

namespace CoupledField{

  //! (Biinearform, which sets just one single entry 
  class SingleEntryBiLinInt : public BiLinearForm {

  public:

    //! Constructor
    //! \param val Coefficient function (vector valued)
    SingleEntryBiLinInt( UInt numDofs, PtrCoefFct& val );

    //! Constructor for scalar function and dof value
    SingleEntryBiLinInt( UInt numDofs, const std::string& val, UInt dof,
                         MathParser* mp );
    
    //! Constructor for scalar function and dof value
    SingleEntryBiLinInt( UInt numDofs, const std::string& real, 
                         const std::string& imag, UInt dof,
                         MathParser* mp );

    //! Copy constructor
    SingleEntryBiLinInt(const SingleEntryBiLinInt& right)
     : BiLinearForm(right) {
      //here we would also need to create a new operator
      //! Number of unknowns
      this->numDofs_ = right.numDofs_;
      this->val_ = right.val_;
    }

    //! \copydoc BiLinearForm::Clone
    virtual SingleEntryBiLinInt* Clone(){
      return new SingleEntryBiLinInt( *this );
    }
    
    //! Destructor
    virtual ~SingleEntryBiLinInt();

    //! Calculation of element 'matrix' (real case )
    void CalcElementMatrix( Matrix<Double>& stiffMat,
                            EntityIterator& ent1, EntityIterator& ent2);

    //! Calculation of element 'matrix' (complex case )
    void CalcElementMatrix( Matrix<Complex>& stiffMat,
                            EntityIterator& ent1, EntityIterator& ent2);

    //! \copydoc LinearForm::IsSolDependent
    bool IsSolDependent() {
      return val_->GetDependency() == CoefFunction::SOLUTION;
    }
    
    //! \see LinearForm::IsComplex
    bool IsComplex() const {
      return val_->IsComplex();
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
    
  protected:

    //! Number of unknowns 
    UInt numDofs_;
    
    //! Coefficient function
    PtrCoefFct val_;
  };
}

#endif
