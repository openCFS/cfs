#ifndef FILE_EDGETETRA1FE
#define FILE_EDGETETRA1FE

#include "tetra1FE.hh"

namespace CoupledField
{

//! Class with  description of tetrahedral element

class EdgeTetra1FE : public Tetra1FE
{
public:
   //! Constructor with type of integration rule
   EdgeTetra1FE();

   //! Deconstructor
   virtual ~EdgeTetra1FE();

   //! Define variables of this class
   virtual void Init();




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
				const std::vector<Double> & LCoord);
  
  

  
  //! calculates the local derivatives of shape functions at an arbitrary local point
  /*!
    \param LDeriv (output) Matrix with local derivatives of all shape functions.
                  Every matrix stores a complete set of global derivations
    \f[ LDeriv_1 = \left( \begin{array}{ccc} N_{1\xi,d\xi} & N_{1\xi,d\eta} & N_{1\xi,d\zeta}\\
                                             N_{1\eta,d\xi} & N_{1\eta,d\eta} & N_{1\eta,d\zeta}\\
					     N_{1\zeta,d\xi} & N_{1\zeta,d\eta} & N_{1\zeta,d\zeta}
					     \end{array}\right) \f]
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void CalcEdgeGlobalDerivShapeFnc(std::vector<Matrix<Double> > & LDeriv, 
					  const std::vector<Double> & LCoord);

protected:
  /// defines the connected nodes with every edge 
  void SetEdgeVertices();
  

  /// Matrix with connected nodes to the proper edge
  /*!
    \f[ \left( \begin{array}{c} E_1 \\ E_2 \\ \cdots \end{array} \right) = 
    \left( \begin{array}{cc} node_1^{E1} & node_2^{E1} \\
                             node_1^{E2} & node_2^{E2}  \\
                              \cdots & \cdots 
			      \end{array} \right) \f]
  */
  Matrix<Integer> edgeVertices_;
  
  

  /// number of edges of the element
  Integer NumEdges_;
};

}
#endif //
