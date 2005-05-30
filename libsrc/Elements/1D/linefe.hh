#ifndef FILE_GEOMETRY_LINEFE_2003
#define FILE_GEOMETRY_LINEFE_2003

#include <Elements/basefe.hh>
#include <Elements/jacobian.hh>

namespace CoupledField
{

  //! Class with general procedures for finite elements: line
  /*! This class is derived from BaseElem. It stores general procedures for
    each type of finite element on lines, such as calculation of Jacobian
    matrix and information about integration points and integration weights.
  */

  class LineFE : public BaseFE {

  public:

    //! Constructor with type of integration rule
    LineFE();

    //! Deconstructor
    virtual ~LineFE();

    //! return FE-Type
    virtual FEType feType() {
      return LINE;
    };


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
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \end{array} \right) \f]
    */
    virtual Double CalcJacobianDet(const Vector<Double> & LCoord,
                                   const Matrix<Double> & CornerCoords);

    //! Calculation of Jacobian determinant at integration point ip
    /*! 
      \param ip (input) Integration point
      \param CornerCoords (input)
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \end{array} \right) \f]
    */
    virtual Double CalcJacobianDetAtIp(const UInt ip, 
                                       const Matrix<Double> & CornerCoors);

    //! calculates the Jacobian Matrix at an arbitrary local point
    /*!
      \param J (output) Jacobian Matrix
      \f[ J = \left( \begin{array}{c} x_{\xi} \end{array}\right)\f] 
      \param LCoord (input) Local coorindates of evaluation point
      \param CornerCoords (input) Coordinates of element corners
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots  \end{array} \right) \f] 
    */
    virtual void CalcJacobian(Matrix<Double> & J, 
                              const Vector<Double> & LCoord, 
                              const Matrix<Double> & CornerCoords);
  
    //! Calculates the Jacobian Matrix at integration point ip
    /*!
      \param J (output) Jacobian Matrix
      \f[ J = \left( \begin{array}{c} x_{\xi}\end{array}\right)\f]  
      \param ip (input) Integration point
      \param CornerCoords (input) Coordinates of element corners
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \end{array} \right) \f] 
    */
    virtual void CalcJacobianAtIp(Matrix<Double> & J, 
                                  const UInt ip, 
                                  const Matrix<Double> & CornerCoords);

    //! calculates the Inverse Jacobian Matrix at an arbitrary local point
    /*!
      \param JInv (output) Inverse Jacobian Matrix 
      \f[ J^{-1} = \left( \begin{array}{cc} \xi_{x} \end{array}\right)\f]
      \param LCoord (input) Local coorindates of evaluation point
      \param CornerCoords (input) Coordinates of element corners
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \end{array} \right) \f] 
    */
    virtual void CalcInvJacobian(Matrix<Double> & JInv,
                                 const Vector<Double> & LCoord,
                                 const Matrix<Double> & CornerCoords);
  
    //! Calculates the Inverse Jacobian Matrix at integration point ip
    /*!
      \param JInv (output) Inverse Jacobian Matrix 
      \f[ J^{-1} = \left( \begin{array}{c} \xi_{x} \end{array}\right)\f]
      \param ip (input) Integration point
      \param CornerCoords (input) Coordinates of element corners
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots  \end{array} \right) \f] 
    */
    virtual void CalcInvJacobianAtIp(Matrix<Double> & JInv,
                                     const UInt ip,
                                     const Matrix<Double> & CornerCoords);
  };

}

#endif // FILE_LINE
