#ifndef FILE_RECTANGLEFE_2003
#define FILE_RECTANGLEFE_2003

#include <Elements/basefe.hh>

namespace CoupledField
{
//! Class with general procedures for quadrilateral finite elements
/*! This class is derived from BaseFE. It stores general procedures for each type of finite element on quadralateral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights */
  
class RectangleFE : public BaseFE
{
public:

  //! Constructor with type of integration rule
  RectangleFE(); 
  
  //! Deconstructor
  virtual ~RectangleFE();

  //! return FE-Type for CLA++
  virtual Integer feType() { return QUAD;}

protected:


   //! Set integration points
  virtual void SetIntPoints();
  
  //! Set value of shape fnc at integration points
  virtual void SetShapeFncAtIp();
  
  //! Set value of shape fnc derivatives at integration points
  virtual void SetShapeFncDerivAtIp();

  //! Calculation of Jacobian determinant at arbitrary local point
  virtual Double CalcJacobianDet(const Vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords);

  //! Calculation of Jacobian determinant at integration point ip
  virtual Double CalcJacobianDetAtIp(const Integer ip, 
				     const Matrix<Double> & CornerCoors);

  //! calculates the Jacobian Matrix at an arbitrary local point
  virtual void CalcJacobian(Matrix<Double> & J, 
			    const Vector<Double> & LCoord, 
			    const Matrix<Double> & CornerCoords);
  
   //! Calculates the Jacobian Matrix at integration point ip
  virtual void CalcJacobianAtIp(Matrix<Double> & J, 
				const Integer ip, 
				const Matrix<Double> & CornerCoords);

  //! calculates the Inverse Jacobian Matrix at an arbitrary local point
  virtual void CalcInvJacobian(Matrix<Double> & JInv,
			       const Vector<Double> & LCoord,
			       const Matrix<Double> & CornerCoords);
  
  //! Calculates the Inverse Jacobian Matrix at integration point ip
  virtual void CalcInvJacobianAtIp(Matrix<Double> & JInv,
				   const Integer ip,
				   const Matrix<Double> & CornerCoords);

};

} // end of namespace

#endif // FILE_RECTANGLEFE_2003
