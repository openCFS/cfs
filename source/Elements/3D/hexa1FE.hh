// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HEXA1FE_2003
#define FILE_HEXA1FE_2003

#include <bitset>

#include "hexaFE.hh"

namespace CoupledField
{

  //! Class with  description of hexahedral element

  class Hexa1FE : public HexaFE
  {
  public:
    //! Constructor with type of integration rule
    Hexa1FE();

    //! Deconstructor
    virtual ~Hexa1FE();

    //! return FE-Type
    virtual FEType feType() {
      return ET_HEXA8;
    };

  protected:
    //! Define variables of this class
    virtual void Init();

    //! Set local corner coordinates
    virtual void SetCornerCoords();

    //! Set local edge indices
    void SetEdgeIndices();

    //! Set local face indices
    void SetFaceIndices();

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

    //! Calculates the local derivatives of incompatible mode shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
      N_{2,d\xi} & N_{2,d\eta} & \cdots \\
      \cdots     & \cdots      & \cdots \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalICModesDerivShapeFnc(Matrix<Double> & LDeriv, 
					       const Vector<Double> & LCoord,
					       const Elem* elem, UInt dof,
					       AnsatzFct::FctEntityType = AnsatzFct::ALL );

    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = CLASSICAL;
        IntegOrder  = 1; // 4+5 is same -> avoid warning
        //        IntegMethod = ECONOMICAL;
        //        IntegOrder  = 5; // 4+5 is same -> avoid warning
    }

    /** Sets the default reduced integration */ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 1; 
    }


    void SetAnsatzFct( shared_ptr<AnsatzFct>& actFct,
                       bool setIntPoints = true );

    void GetNumFncs(Vector<UInt>& numFcns, 
                    const shared_ptr<AnsatzFct>& fcnType, 
                    AnsatzFct::FctEntityType fctEntityType, 
                    UInt dof = 1);


    UInt GetNumFncs( const shared_ptr<AnsatzFct>& fncType );

  private:
   
  };

}
#endif //
