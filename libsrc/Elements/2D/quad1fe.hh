#ifndef FILE_QUAD1FE_2003
#define FILE_QUAD1FE_2003

#include "Elements/basefe.hh"
#include "Elements/2D/rectanglefe.hh"

namespace CoupledField
{
  //! Quadrilateral finite element with four nodes (linear interpolation function)
  
class Quad1FE : public RectangleFE
{
public:
  
  //! Constructor with type of integration rule
  Quad1FE();
  
  //! Destructor
  virtual ~Quad1FE();
  

protected:

  //! Initialize Quadrilateral element
  virtual void Init();

  //! Set local corner coordinates
  virtual void SetCornerCoords();

  //! calculates the shape functions at an arbitrary local point
  virtual void CalcShapeFnc(std::vector<Double> & LShape, 
			    const std::vector<Double> & LCoord);
  
  //! calculates the local derivatives of shape functions at an arbitrary local point
  virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				      const std::vector<Double> & LCoord);

private:
};

} // end of namespace

#endif // FILE_QUAD1FE_2003
