#ifndef FILE_LINE2FE_2004
#define FILE_LINE2FE_2004

#include <Elements/basefe.hh>
#include <Elements/1D/linefe.hh>

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
    virtual void CalcShapeFnc(Vector<Double> & LShape, 
                              const Vector<Double> & LCoord);
  
    //! calculates the local derivatives of shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{cc} N_{1,d\xi}  \\
      N_{2,d\xi} \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                        const Vector<Double> & LCoord);


  private:
  };

} // end of namespace

#endif // FILE_LINE2FE_2004
