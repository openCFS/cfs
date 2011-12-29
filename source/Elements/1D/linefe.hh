// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GEOMETRY_LINEFE_2003
#define FILE_GEOMETRY_LINEFE_2003

#include "Elements/basefe.hh"

#include "Domain/elem.hh"
#include "General/defs.hh"
#include "MatVec/vector.hh"

namespace CoupledField {
template <class TYPE> class Matrix;
template <class TYPE> class StdVector;
}  // namespace CoupledField

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
    virtual Elem::FEType feType() {
      return Elem::UNDEF;
    };


  protected:

   /** the childs fill here the integration points map via AddIntegrationPoints() */    
    virtual void FillIntegrationPoints();
  
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
                                   const Matrix<Double> & CornerCoords,
                                   const Elem* elem );

    //! Calculation of Jacobian determinant at integration point ip
    /*! 
      \param ip (input) Integration point
      \param CornerCoords (input)
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \end{array} \right) \f]
    */
    virtual Double CalcJacobianDetAtIp(const UInt ip, 
                                       const Matrix<Double> & CornerCoors,
                                       const Elem* elem );

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
                              const Matrix<Double> & CornerCoords,
                              const Elem* elem );
  
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
                                  const Matrix<Double> & CornerCoords,
                                  const Elem* elem );

    //! Returns wether a given local coordinate is
    //! within this element
    //! \param localCoords (input) local coordinates
    //! \param coordsInside (output) which local coordinates are inside
    virtual void CoordsInsideElem(const Matrix<Double> & localCoords,
                                  const Double tolerance,
                                  StdVector<bool> & coordsInside) const;

  };

}

#endif // FILE_LINE
