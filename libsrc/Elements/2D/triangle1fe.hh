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

    //! Constructor with type of integration rule
    Triangle1FE();
  
    //! Destructor
    virtual ~Triangle1FE();
  

  protected:

    //! Initialize Trianglerilateral element
    virtual void Init();

    //! Set local corner coordinates
    virtual void SetCornerCoords();

    //! Set local edge indices
    void SetEdgeIndices();

    //! calculates the shape functions at an arbitrary local point
    /*!
      \param Shape (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcShapeFnc(Vector<Double> & LShape, 
                              const Vector<Double> & LCoord);
  
    //! calculates the local derivatives of shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
      N_{2,d\xi} & N_{2,d\eta} & \cdots \\
      \cdots     & \cdots      & \cdots \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                        const Vector<Double> & LCoord);

  private:
  };

} // end of namespace

#endif // FILE_TRIANGLE1FE_2003
