#ifndef FILE_HEXA2FE_2004
#define FILE_HEXA2FE_2004

#include "hexaFE.hh"

namespace CoupledField
{

  //! Class with  description of hexahedral element

  class Hexa2FE : public HexaFE
  {
  public:
    //! Constructor with type of integration rule
    Hexa2FE();

    //! Deconstructor
    virtual ~Hexa2FE();

  protected:
    //! Define variables of this class
    virtual void Init();

    //! Set local corner coordinates
    virtual void SetCornerCoords();


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

    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 5; // 4+5 is same -> avoid warning
    }

    /** Sets the default reduced integration  */ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 5; 
    }


  private:
   
  };

}
#endif //
