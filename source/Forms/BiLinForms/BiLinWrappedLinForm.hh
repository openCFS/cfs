#ifndef CFS_BILIN_WRAPPED_LINFORM_HH
#define CFS_BILIN_WRAPPED_LINFORM_HH

#include "BiLinearForm.hh"

namespace CoupledField {


// forward class declarations
class LinearForm;

//! Encapsulates a LinearForm for calculating row/column element "matrix"

//! This class basically wraps a LinearForm, which normally only calculates 
//! a row-vector. By encapsulating this form it can deliver either a 
//! (numFcts x 1 ) column matrix or a (1 x numFctrs) matrix (if transposed is set).
class BiLinWrappedLinForm : public BiLinearForm {
  
public:
  
  //! Constructor
  BiLinWrappedLinForm(LinearForm* linForm,
                      bool assembleTranposed = false);

  //! Destructor
  virtual ~BiLinWrappedLinForm();

  //! \copydoc BiLinearForm::CalcElementMatrix
  virtual void CalcElementMatrix( Matrix<Double>& elemMat,
                                  EntityIterator& ent1,
                                  EntityIterator& ent2);

  //! \copydoc BiLinearForm::CalcElementMatrix
  virtual void CalcElementMatrix( Matrix<Complex>& elemMat,
                                  EntityIterator& ent1,
                                  EntityIterator& ent2);

  virtual bool IsSolDependent();

  virtual bool IsComplex();

  virtual void SetFeSpace( shared_ptr<FeSpace> feSpace );

  virtual void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2 );

protected:
  
  //! Pointer to LinearForm
  LinearForm* linForm_;
  
  //! Flag if transposed element vector is to be assembled
  bool assembleTransposed_;
};

}

#endif
