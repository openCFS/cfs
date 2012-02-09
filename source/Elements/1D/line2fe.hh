// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LINE2FE_2004
#define FILE_LINE2FE_2004

#include "Elements/1D/linefe.hh"

#include "Domain/ansatzFct.hh"
#include "Domain/elem.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"

namespace CoupledField {
template <class TYPE> class Matrix;
}  // namespace CoupledField

namespace CoupledField
{
  //! 1D line element with three nodes (quadratic interpolation function)
  
  class Line2FE : public LineFE
  {
  public:
  
    //! Constructor with type of integration rule
    Line2FE();
  
    //! Destructor
    virtual ~Line2FE();
  
    //! return FE-Type
    virtual Elem::FEType feType() {
      return Elem::LINE3;
    };

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
    virtual void CalcShapeFnc( Vector<Double> & LShape, 
                               const Vector<Double> & LCoord,
                               const Elem* elem , UInt dof,
                               AnsatzFct::FctEntityType );
    
    //! calculates the local derivatives of shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{cc} N_{1,d\xi}  \\
      N_{2,d\xi} \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv, 
                                         const Vector<Double> & LCoord,
                                         const Elem* elem , UInt dof,
                                         AnsatzFct::FctEntityType );

    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 5; // 2+3 are same -> avoid msg
    }

    /** Sets the default reduced integration */ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 3; 
    }

    //! Get the local coordinates for given global ones
    //! \param localCoords (output) local coordinates
    //! \param globalCoords (input) global coordinates
    //! \param coordMat (input) global corner coordinates of element
    //!                         (spaceDim \f$\times\f$ nrNodes)
    virtual void Global2LocalCoords(Matrix<Double> & localCoords,
                                    const Matrix<Double> & globalCoords,
                                    const Matrix<Double> & coordMat);

  private:
  };

} // end of namespace

#endif // FILE_LINE2FE_2004
