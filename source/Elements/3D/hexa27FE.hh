// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HEXA27FE_2008
#define FILE_HEXA27FE_2008

#include "Domain/ansatzFct.hh"
#include "Domain/elem.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "hexaFE.hh"

namespace CoupledField {
template <class TYPE> class Matrix;
}  // namespace CoupledField

namespace CoupledField
{

  // WARNING!!!
  // THIS CLASS IS JUST COPY OF THE 20 NODE HEXAHEDRON AT THE MOMENT

  //! Class with  description of hexahedral element

  class Hexa27FE : public HexaFE
  {
  public:
    //! Constructor with type of integration rule
    Hexa27FE();

    //! Deconstructor
    virtual ~Hexa27FE();

    //! return FE-Type
    virtual Elem::FEType feType() {
      return Elem::HEXA27;
    };

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
    virtual void CalcShapeFnc( Vector<Double> & LShape, 
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
    virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                        const Vector<Double> & LCoord,
                                        const Elem* elem , UInt dof,
                                        AnsatzFct::FctEntityType );

    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 5; // 4+5 is same -> avoid warning
    }

    /** Sets the default reduced integration  */ 
    void SetDefaultReducedIntegration()
    {
      IntegMethod = SPECIAL;
      IntegOrder  = 1; 
      SetIntPoints(IntegMethod,IntegOrder);
    }


  private:
   
  };

}
#endif //
