#ifndef FILE_HEXA1FE_2003
#define FILE_HEXA1FE_2003

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


  //! Calculates corresponding volume point of neighbouring surface

  //! For a given surface element and a neighbouring volume element this
  //! mehtod calculates the local volume-coordinates out of the given
  //! local surface-coordinates, which have one less dimension.
  //! This can be used to get the corrsponding volume coordinates of 
  //! the integration points of a surface. Therefore it calculates 
  //! on which side of the volume element the surface elemente lies
  //! and creates the according volume point.
  /*!
    \param surfConnect (input) Node numbers of surface element
    \param volConnect (input) Node numbers of colume element
    \param surfIntPoint (input) Surface integration point, which gets mapped
    onto the volume element
    \param volIntPoint (output) Corresponding volume integration point
  */
  void GetLocalIntPoints4Surface(const StdVector<Integer> & surfConnect,
				 const StdVector<Integer> & volConnect,
				 const Vector<Double> & surfIntPoint,
				 Vector<Double> & volIntPoint);


private:
   
};

}
#endif //
