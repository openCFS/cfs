#ifndef FILE_TETRA1FE
#define FILE_TETRA1FE

#include "tetraFE.hh"

namespace CoupledField
{

//! Class with  description of tetrahedral element

class Tetra1FE : public TetraFE
{
public:
   //! Constructor with type of integration rule
   Tetra1FE();

   //! Deconstructor
   virtual ~Tetra1FE();

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
   
};

}
#endif //
