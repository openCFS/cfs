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
  /*! 
    \param LCoord (input) Local Coordinates of evaluation point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
  */
  virtual Double CalcJacobianDet(const std::vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords);

  //! Calculation of Jacobian determinant at integration point ip
  /*! 
    \param ip (input) Integration point
    \param CornerCoords (input)
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
  */
  virtual Double CalcJacobianDetAtIp(const Integer ip, 
				     const Matrix<Double> & CornerCoors);

  //! calculates the Jacobian Matrix at an arbitrary local point
   /*!
    \param J (output) Jacobian Matrix
    \f[ J = \left( \begin{array}{ccc} x_{\xi} & x_{\eta} & x_{\zeta} \\
                                    y_{\xi} & y_{\eta} & y_{\zeta}\\
                                    z_{\xi} & z_{\eta} & z_{\zeta} \end{array}\right)\f] 
    \param LCoord (input) Local coorindates of evaluation point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f] 
  */
  virtual void CalcJacobian(Matrix<Double> & J, 
			    const std::vector<Double> & LCoord, 
			    const Matrix<Double> & CornerCoords);
  
   //! Calculates the Jacobian Matrix at integration point ip
   /*!
    \param J (output) Jacobian Matrix
    \f[ J = \left( \begin{array}{ccc} x_{\xi} & x_{\eta} & x_{\zeta} \\
                                    y_{\xi} & y_{\eta} & y_{\zeta}\\
                                    z_{\xi} & z_{\eta} & z_{\zeta} \end{array}\right)\f]  
    \param ip (input) Integration point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f] 
   */
  virtual void CalcJacobianAtIp(Matrix<Double> & J, 
				const Integer ip, 
				const Matrix<Double> & CornerCoords);

  //! calculates the Inverse Jacobian Matrix at an arbitrary local point
  /*!
    \param JInv (output) Inverse Jacobian Matrix 
    \f[ J^{-1} = \left( \begin{array}{ccc} \xi_{x} & \xi_{y} & \xi_{z} \\
                                    \eta_{x} & \eta_{y} & \eta_{z}\\
                                    \zeta_{x} & \zeta_{y} & \zeta_{z} \end{array}\right)\f]
    \param LCoord (input) Local coorindates of evaluation point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f] 
  */
  virtual void CalcInvJacobian(Matrix<Double> & JInv,
			       const std::vector<Double> & LCoord,
			       const Matrix<Double> & CornerCoords);
  
  //! Calculates the Inverse Jacobian Matrix at integration point ip
  /*!
    \param JInv (output) Inverse Jacobian Matrix 
    \f[ J^{-1} = \left( \begin{array}{ccc} \xi_{x} & \xi_{y} & \xi_{z} \\
                                    \eta_{x} & \eta_{y} & \eta_{z}\\
                                    \zeta_{x} & \zeta_{y} & \zeta_{z} \end{array}\right)\f]
    \param ip (input) Integration point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f] 
  */
  virtual void CalcInvJacobianAtIp(Matrix<Double> & JInv,
				   const Integer ip,
				   const Matrix<Double> & CornerCoords);

};

} // end of namespace

#endif // FILE_RECTANGLEFE_2003
