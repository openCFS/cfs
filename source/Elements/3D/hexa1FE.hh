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
    virtual Elem::FEType feType() {
      return Elem::HEXA8;
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
        IntegOrder  = 3; 
        //        IntegMethod = ECONOMICAL;
        //        IntegOrder  = 5; // 4+5 is same -> avoid warning
    }

    /** Sets the default reduced integration */
    void SetDefaultReducedIntegration()
    {
      IntegMethod = CLASSICAL;
      IntegOrder  = 1; // 4+5 is same -> avoid warning
        
    }


    void SetAnsatzFct( shared_ptr<AnsatzFct>& actFct,
                       bool setIntPoints = true );

    void GetNumFncs(StdVector<UInt>& numFcns,
                    const shared_ptr<AnsatzFct>& fcnType,
                    AnsatzFct::FctEntityType fctEntityType,
                    UInt dof = 1);


    UInt GetNumFncs( const shared_ptr<AnsatzFct>& fncType );

    // ============================= methods for edge elements =======================
    // ===============================================================================

    //! calculates the vectorial shape functions at an arbitrary local point
    /*!
      \param shape (output) Matrix of shape function values 
      \f[ \left( \begin{array}{c} E_1 \\ E_2 \\ \cdots \end{array} \right) = 
      \left( \begin{array}{ccc} N_{1\xi} & N_{1\eta} & N_{1\zeta} \\
      N_{2\xi} & N_{2\eta} & N_{2\zeta} \\
      \cdots     & \cdots      & \cdots 
      \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcEdgeShapeFnc(Matrix<Double> & shape, 
                                  const Vector<Double> & LCoord, 
                                  const Matrix<Double> & CornerCoords,
                                  const Elem* elem);
  
    //! calculates the local derivatives of the edge shape functions at an arbitrary local point
    /*!
      \param deriv (output) Vector of matrices with local derivatives of all shape functions.
      Every matrix stores a complete set of global derivations.
      \f[ deriv[edge1] = \left( \begin{array}{ccc} N_{1\xi,d\xi} & N_{1\xi,d\eta} & N_{1\xi,d\zeta}\\
      N_{1\eta,d\xi} & N_{1\eta,d\eta} & N_{1\eta,d\zeta}\\
      N_{1\zeta,d\xi} & N_{1\zeta,d\eta} & N_{1\zeta,d\zeta}
      \end{array}\right) \f]
      \param lCoord (input) Local coordinates of evalutation point 
    */
    virtual void GetEdgeGlobalDerivShapeFnc(StdVector<Matrix<Double> > & deriv, 
                                            const Vector<Double> & lCoord,
                                            const Matrix<Double> & CornerCoords,
                                            const Elem* elem);

    void Global2LocalCoords(Matrix<Double> & localCoords,
                            const Matrix<Double> & globalCoords,
                            const Matrix<Double> & coordMat );
  private:

    virtual void CalcSpectralShFct( Vector<Double> & Shape,
                                    const Vector<Double> & LCoord,
                                    const Elem* elem, UInt dof,
                                    AnsatzFct::FctEntityType type );

    virtual void CalcSpectralDerivFct( Matrix<Double> & LDeriv,
                                       const Vector<Double> & LCoord,
                                       const Elem* elem, UInt dof,
                                       AnsatzFct::FctEntityType type);

    //! Set value of shape fnc at integration points
    virtual void SetShapeFnc2AtIp();


    //! Set value of shape fnc derivatives at integration points
    virtual void SetShapeFnc2DerivAtIp();

    //! calculates the shape functions used for edge elements at an arbitrary local point
    /*!
    \param Shape (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcShapeFnc2( Vector<Double> & LShape, 
                                const Vector<Double> & LCoord);


    //! calculates the local derivatives of shape functions for edge elements at an arbitrary local point
    /*!
    \param LDeriv (output) Matrix with local derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
    N_{2,d\xi} & N_{2,d\eta} & \cdots \\
    \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalDerivShapeFnc2( Matrix<Double> & LDeriv, 
                                          const Vector<Double> & LCoord,
                                          const Elem* elem , UInt dof);

    void GetGlobDerivShFnc2(Matrix<Double> & Deriv, 
                            const Vector<Double> & LCoord,
                            const Matrix<Double> & CornerCoords,
                            const Elem * elem, 
                            UInt dof );

    Vector<Double> * ShFnc2AtIp_; //!< shape functions for edge elements at IPs
    Matrix<Double> * ShFnc2DerivAtIp_; //!< deivative of  shape functions for edge elements at IPs
  };

}
#endif //
