#include "MaterialTensor.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Enum.hh"
#include "MatVec/Matrix.hh"

DEFINE_LOG(mt, "materialTensor")

namespace CoupledField
{

template <class TYPE>
MaterialTensor<TYPE>::MaterialTensor(MaterialTensorNotation notation, const Matrix<TYPE>& m):
  notation_(notation)
{
  matrix_ = new Matrix<TYPE>(m);
}

template <class TYPE>
MaterialTensor<TYPE>::MaterialTensor(MaterialTensorNotation notation, Matrix<TYPE>* m, bool copy):
  notation_(notation)
{
  if (!copy) {
    matrix_ = m;
    delmem_ = false;
  } else
    matrix_ = new Matrix<TYPE>(*m);
}

template <class TYPE>
MaterialTensor<TYPE>::MaterialTensor(MaterialTensorNotation notation):
  notation_(notation)
{
  matrix_ = new Matrix<TYPE>();
}

template <class TYPE>
MaterialTensor<TYPE>::~MaterialTensor()
{
  if (delmem_) {
    delete matrix_;
    matrix_ = NULL;
  }
}

template <class TYPE>
void MaterialTensor<TYPE>::Replace(MaterialTensorNotation notation, const Matrix<TYPE>& m)
{
  assert(notation_ == notation);
  *matrix_ = m; // copy constructor fun
}

// in-place conversion of tensor to Voigt notation
// throws an exception if validate=true and notation before conversion is already Voigt or None
template <class TYPE>
void MaterialTensor<TYPE>::ToVoigt(bool validate)
{
  if (validate && notation_ != HILL_MANDEL)
   throw Exception("Converting material tensor to Voigt notation is not allowed! Current notation: " + tensorNotation.ToString(notation_));

  if (notation_ != VOIGT)
    HillMandelToVoigt();

  notation_ = VOIGT;
}

// in-place conversion of tensor to Voigt notation
// throws an exception if validate=true and notation before conversion is already Hill-Mandel or None
template<class TYPE>
void MaterialTensor<TYPE>::ToHillMandel(bool validate)
{
  if (validate && notation_ != VOIGT)
    throw Exception("Converting material tensor to Hill-Mandel notation is not allowed! Current notation: " + tensorNotation.ToString(notation_));

  if (notation_ != HILL_MANDEL)
    VoigtToHillMandel();
  notation_ = HILL_MANDEL;
}

template <class TYPE>
MaterialTensor<TYPE> MaterialTensor<TYPE>::AsVoigt(bool validate) const
{
  MaterialTensor<TYPE> ret(HILL_MANDEL, matrix_, true);
  // ToVoigt() does the validation check
  ret.ToVoigt(validate);
  return ret;
}

template <class TYPE>
MaterialTensor<TYPE> MaterialTensor<TYPE>::AsHillMandel(bool validate) const
{
  MaterialTensor<TYPE> ret(VOIGT, matrix_, true);
  // ToHillMandel() does the validation check
  ret.ToHillMandel(validate);
  return ret;
}

/** Private conversion functions **/

template<class TYPE>
void MaterialTensor<TYPE>::VoigtToHillMandel()
{
  assert(matrix_->IsQuadratic());
  assert(notation_ == VOIGT);
  // matrix is quadratic -> rows == cols
  unsigned int rows = matrix_->GetNumRows();
  assert(rows == 3 || rows == 6);

  double sqrt2 = sqrt(2);
  if (rows == 3) {
    for(unsigned int i = 0; i < rows-1; i++)
    {
      (*matrix_)[i][rows-1] *= sqrt2;
      (*matrix_)[rows-1][i] *= sqrt2;
    }
    (*matrix_)[rows-1][ rows-1] *= 2.0;
  } else {
    for(unsigned int i = 0; i < rows; i++) {
      for (unsigned int j = i; j < rows; j++) {
        if (i > 2 || j > 2) {
          (*matrix_)[i][j] *= sqrt2;
          (*matrix_)[j][i] *= sqrt2;
        } else if (i == j && i > 2)
          (*matrix_)[i][j] *= 2.0;
      }
    }
  }
}

/* Convert from Hill-Mandel to Voigt notation */
template<class TYPE>
void MaterialTensor<TYPE>::HillMandelToVoigt()
{
  // based on mativ_rot.py
  assert(matrix_->IsQuadratic());
  assert(notation_ == HILL_MANDEL);
  // matrix is quadratic -> rows == cols
  unsigned int rows = matrix_->GetNumRows();

  assert(rows == 3 || rows == 6);

  double inv_sqrt2 = 1/sqrt(2);
  if (rows == 3) {
    for(unsigned int i = 0; i < rows-1; i++) {
      (*matrix_)[i][rows-1] *= inv_sqrt2;
      (*matrix_)[rows-1][i] *= inv_sqrt2;
    }
    (*matrix_)[rows-1][rows-1] *= 0.5;
  } else {
    for(unsigned int i = 0; i < rows; i++) {
      for (unsigned int j = i; j < rows; j++) {
        if (i > 2 || j > 2) {
          (*matrix_)[i][j] *= inv_sqrt2;
          (*matrix_)[j][i] *= inv_sqrt2;
        } else if (i == j && i > 2)
          (*matrix_)[i][j] *= 0.5;
      }
    }
  }
}

template class MaterialTensor<double>;
template class MaterialTensor<Complex>;

} // end of namespace
