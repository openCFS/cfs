#ifndef FILE_LINE1FE_2003
#define FILE_LINE1FE_2003

#include <Elements/basefe.hh>
#include <Elements/1D/gelinefe.hh>

namespace CoupledField
{
  //! 1D line element with two nodes (linear interpolation function)
  
class Line1FE : public GeLineFE
{
public:
  
  //! Constructor with type of integration rule
  Line1FE();
  
  //! Destructor
  virtual ~Line1FE();
  

protected:

  //! Initialize line element
  virtual void Init();

  //! Set local corner coordinates
  virtual void SetCornerCoords();

  //! calculates the shape functions at an arbitrary local point
  /*!
    \param Shape (output) Vector of shape fnc values \f$ (N_{1},N_{2})^T \f$
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void CalcShapeFnc(std::vector<Double> & LShape, 
			    const std::vector<Double> & LCoord);
  
  //! calculates the local derivatives of shape functions at an arbitrary local point
  /*!
    \param LDeriv (output) Matrix with local derivatives of all shape functions
    \f[ \left( \begin{array}{cc} N_{1,d\xi}  \\
                                  N_{2,d\xi} \end{array}\right) \f]
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				      const std::vector<Double> & LCoord);


private:
};

} // end of namespace

#endif // FILE_LINE1FE_2003
