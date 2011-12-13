// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_TRIANGLE2FE_2004
#define FILE_TRIANGLE2FE_2004

#include "Elements/2D/trianglefe.hh"

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
  //!  Triangle finite element with six nodes (quadratic interpolation function)
  
  class Triangle2FE : public TriangleFE
  {
  public:

    //! Constructor with type of integration rule
    Triangle2FE();
  
    //! Destructor
    virtual ~Triangle2FE();
  
    //! return FE-Type
    virtual Elem::FEType feType() {
      return Elem::TRIA6;
    };

  protected:

    //! Initialize Trianglerilateral element
    virtual void Init();

    //! Set local corner coordinates
    virtual void SetCornerCoords();

    //! calculates the shape functions at an arbitrary local point
    /*!
      \param Shape (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param LCoord (input) Local coordinates of evalutation point 
    */
    void CalcShapeFnc( Vector<Double> & LShape, 
                       const Vector<Double> & LCoord,
                       const Elem* elem , UInt dof,
                       AnsatzFct::FctEntityType );
    
    //! calculates the local derivatives of shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
      N_{2,d\xi} & N_{2,d\eta} & \cdots \\
      \cdots     & \cdots      & \cdots \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    void CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv, 
                                 const Vector<Double> & LCoord,
                                 const Elem* elem , UInt dof,
                                 AnsatzFct::FctEntityType);
    
    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 5;
    }

    /** Sets the default reduced integration*/ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 3;
    }


  private:
  };

} // end of namespace

#endif // FILE_TRIANGLE2FE_2004
