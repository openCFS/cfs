#ifndef FILE_LINEARFORM
#define FILE_LINEARFORM

namespace CoupledField
{

/// class for calculation right hand side
template <class Dim>
class LinearForm : public BaseForm<Dim> 
{
public:
  ///
  LinearForm(BaseElem * aptelem);

  ///
  virtual ~LinearForm();

  /// Calculation of vector of right hand side 
  void CalcElemMatrix(Dim * ptCoord, Vector<Double> & Result);

protected:

private:

};

template class  LinearForm <Point2D>;
template class  LinearForm <Point3D>;

}

#endif // FILE_LINEARFORM
