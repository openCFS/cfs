#ifndef FILE_QUAD1FE_2003
#define FILE_QUAD1FE_2003

#include <Elements/basefe.hh>
#include <Elements/2D/rectanglefe.hh>

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
  /*!
    \param Shape (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void CalcShapeFnc(std::vector<Double> & LShape, 
			    const std::vector<Double> & LCoord);
  
  //! calculates the local derivatives of shape functions at an arbitrary local point
  /*!
    \param LDeriv (output) Matrix with local derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
                                  N_{2,d\xi} & N_{2,d\eta} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				      const std::vector<Double> & LCoord);
  
   //! Calculates a measure for the geometric distortion of an element
  /*!
    \param cornerCoords (input) Corner coordinates of the element
    \param size (input) Absolute size of element in all dimensions
    \param displacement (input) Displacement of the corner points (same ordering as CornerCoords!!)
  */
  virtual Double CalcMeanStrain(Matrix<Double> &cornerCoords, Matrix<Double> &displacements);

private:
};

} // end of namespace

#endif // FILE_QUAD1FE_2003
