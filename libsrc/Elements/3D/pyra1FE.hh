#ifndef FILE_PYRA1FE
#define FILE_PYRA1FE

#include "pyraFE.hh"

namespace CoupledField
{

//! Class with  description of pyramidal element

class Pyra1FE : public PyraFE
{
public:
   //! Constructor with type of integration rule
   Pyra1FE();

   //! Deconstructor
   virtual ~Pyra1FE();

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
  virtual void CalcShapeFnc(std::vector<Double> & LShape, 
			    const std::vector<Double> & LCoord);


  
  //! calculates the local derivatives of shape functions at an arbitrary local point
  /*!
    \param LDeriv (output) Matrix with local derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
                                  N_{2,d\xi} & N_{2,d\eta} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				      const std::vector<Double> & LCoord);




  // ============================= methods for edge elements =======================
  // ===============================================================================



  
  



protected:

};

}
#endif //
