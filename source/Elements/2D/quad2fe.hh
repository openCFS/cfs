// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_QUAD2FE_2004
#define FILE_QUAD2FE_2004

#include <Elements/basefe.hh>
#include <Elements/2D/rectanglefe.hh>

namespace CoupledField
{
  //! Quadrilateral finite element with eight nodes (quadratic interpolation function)
  
  class Quad2FE : public RectangleFE
  {
  public:
  
    //! Constructor with type of integration rule
    Quad2FE();
  
    //! Destructor
    virtual ~Quad2FE();
  
    //! return FE-Type
    virtual FEType feType() {
      return ET_QUAD8;
    };

  protected:

    //! Initialize Quadrilateral element
    virtual void Init();

    //! Set local nodal coordinates
    virtual void SetCornerCoords();

    //! calculates the shape functions at an arbitrary local point
    /*!
      \param Shape (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcShapeFnc(Vector<Double> & LShape, 
                              const Vector<Double> & LCoord,
                              const Elem*, UInt dof,
                              AnsatzFct::FctEntityType );
  
    //! calculates the local derivatives of shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
      N_{2,d\xi} & N_{2,d\eta} & \cdots \\
      \cdots     & \cdots      & \cdots \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                        const Vector<Double> & LCoord,
                                        const Elem*, UInt dof,
                                        AnsatzFct::FctEntityType );
  
    //! NOT YET IMPLEMENTED FOR QUADRATIC ELEMENTS!!!Calculates a measure for the geometric distortion of an element
    /*!
      \param cornerCoords (input) Corner coordinates of the element
      \param size (input) Absolute size of element in all dimensions
      \param displacement (input) Displacement of the corner points (same ordering as CornerCoords!!)
    */
    virtual Double CalcMeanStrain(Matrix<Double> &cornerCoords, Matrix<Double> &displacements);

    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 5; // actually 2 but 2+3 is same and we avoid the msg
    }

    /** Sets the default reduced integration  */ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 3; 
    }

  private:
  };

} // end of namespace

#endif // FILE_QUAD2FE_2004
