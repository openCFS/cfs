#ifndef MATERIALTENSOR_HH_
#define MATERIALTENSOR_HH_

#include "General/Enum.hh"
#include "MatVec/Matrix.hh"

namespace CoupledField {

template<class TYPE>
class MaterialTensor
{
  public:
    // constructors
    MaterialTensor(MaterialTensorNotation notation, const Matrix<TYPE>& m);
    MaterialTensor(MaterialTensorNotation notation, Matrix<TYPE>* m, bool copy);
    MaterialTensor(MaterialTensorNotation notation);

    //! Destructor
    virtual ~MaterialTensor();

    /** in-place conversion to Voigt / H-M notation **/
    // validate: throw exception if, before conversion, tensor is already in Voigt notation
    void ToVoigt(bool validate = true);
    // validate: check if conversion is allowed - see ToVoigt()
    void ToHillMandel(bool validate = true);

    /** returns new objects in Voigt / H-M notation **/
    // validate: see ToVoigt()
    MaterialTensor<TYPE> AsVoigt(bool validate = true) const;
    // validate: see ToVoigt()
    MaterialTensor<TYPE> AsHillMandel(bool validate = true) const;

    // service functions for casting
//    inline MaterialTensor<double> AsDouble() {
//      return static_cast<MaterialTensor<double> >(*this);
//    }

    /** returns the matrix in the current notation.
     * Parameter notation is only used for checks in debug. matrix is NOT
     * returned in notation "notation"! For this use AsVoigt() or AsHillMandel().
     */
    Matrix<TYPE>& GetMatrix(MaterialTensorNotation notation) {
      assert(notation == notation_);
      return *matrix_;
    };

    MaterialTensorNotation GetNotation() const { return notation_; }

  private:
    MaterialTensorNotation notation_ = NO_NOTATION;
    Matrix<TYPE>* matrix_ = NULL;

    bool delmem_ = true;

    /** In-place conversion from Voigt to Hill Mandel */
    void VoigtToHillMandel();

    /** In-place conversion from Hill-Mandel to Voigt Notation */
    void HillMandelToVoigt();
};

}// namespace

#endif /*DESIGN_SPACE_HH_*/
