#ifndef FILE_TRIANGLE1FE_2003
#define FILE_TRIANGLE1FE_2003

#include <Elements/basefe.hh>
#include <Elements/2D/trianglefe.hh>

namespace CoupledField
{
  //!  Triangle finite element with three nodes (linear interpolation function)
  
  class Triangle1FE : public TriangleFE
  {
  public:

    /** Constructor with optional integraption paramters
     * used only to construct higher order wedge elements
     * @param method leave blank as normal user! Then the default or XML stuff is loaded
     * @param oder see method */
    Triangle1FE(IntegrationMethod method = UNDEFINED, int order=0);
  
    //! Destructor
    virtual ~Triangle1FE();
  

  protected:

    //! Initialize Trianglerilateral element
    virtual void Init(IntegrationMethod method, int order);

    //! Set local corner coordinates
    virtual void SetCornerCoords();

    //! Set local edge indices
    void SetEdgeIndices();

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
                                 AnsatzFct::FctEntityType );
    
    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 2;
    }

    /** Sets the default reduced integration  */ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 1;
    }


  private:
  };

} // end of namespace

#endif // FILE_TRIANGLE1FE_2003
