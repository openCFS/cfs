#ifndef FILE_LINELASTINT
#define FILE_LINELASTINT

#include <Elements/basefe.hh>
#include <Forms/baseForm.hh>

namespace CoupledField
{

  //! class for calculation of mechanical plain strain condition
  class mechPlainStrainInt : public BaseForm
  {
  public:
    //! Constructor
    mechPlainStrainInt(BaseFE * aptelem);
    
    //! Deconstructor
    virtual ~mechPlainStrainInt();
    
    //! Virtual function, implemented in derived classes
    virtual void CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & StiffMat);
    
    //! Prints the bilinear form
    virtual void Print(std::ostream * out, const Matrix<Double> Result) const;
    
  protected:    
  };
  
} //end namespace

#endif // FILE_BASERFORM
