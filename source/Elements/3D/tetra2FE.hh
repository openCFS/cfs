// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_TETRA2FE
#define FILE_TETRA2FE

#include "Domain/ansatzFct.hh"
#include "Domain/elem.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "tetraFE.hh"

namespace CoupledField {
template <class TYPE> class Matrix;
}  // namespace CoupledField

namespace CoupledField
{

//! Class with  description of tetrahedral element

class Tetra2FE : public TetraFE
{
public:
   //! Constructor with type of integration rule
   Tetra2FE();

   //! Deconstructor
   virtual ~Tetra2FE();

    //! return FE-Type
    virtual Elem::FEType feType() {
      return Elem::TET10;
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
  virtual void CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv, 
                                       const Vector<Double> & LCoord,
                                       const Elem* elem , UInt dof,
                                        AnsatzFct::FctEntityType );
  
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
				const Matrix<Double> & CornerCoords);
  
  

  
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
  virtual void GetEdgeGlobalDerivShapeFnc(StdVector<Matrix<Double>* > & deriv, 
					  const Vector<Double> & lCoord,
					  const Matrix<Double> & CornerCoords);

    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 5;
    }

    /** Sets the default reduced integration  */ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 3;
    }


protected:
  /// defines the connected nodes with every edge 
  void SetEdgeVertices();   
};

}
#endif //
