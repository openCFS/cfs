#ifndef FILE_PYRA2FE
#define FILE_PYRA2FE

#include "pyraFE.hh"

namespace CoupledField
{

//! Class with  description of pyramidal element

class Pyra2FE : public PyraFE
{
public:
   //! Constructor with type of integration rule
   Pyra2FE();

   //! Deconstructor
   virtual ~Pyra2FE();

    //! return FE-Type
    virtual FEType feType() {
      return ET_PYRA13;
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



  
  



protected:

};

}
#endif //
