#ifndef FILE_FE_NODAL_HH
#define FILE_FE_NODAL_HH

#include "BaseFE.hh"

namespace CoupledField {

//! Base class for all nodal / non-hierarchic Finite Elements

//! This class encapsulates all functionality needed for nodal based
//! finite elements where the shape functions are evaluate at discrete 
//! supporting points, i.e. Lagrangian finite elements. In this case
//! the shape functions and their local derivatives can be pre-computed.
class FeNodal {
public:
  
  //! Constructor
  FeNodal();

  //! Destructor
  virtual ~FeNodal();
  
  // ----------------------------------------------------------------------
  //  Polynomial Evaluation
  // ----------------------------------------------------------------------

  //! Evaluate Lagrange polynomial 
  void EvaluateLagrangePolynomial( Vector<Double> & shape, Double coord,
                                   UInt order);

  //! Evaluate derivative of Lagrange polynomials 
  void EvaluateDerivLagrangePolynomial( Vector<Double> & deriv, Double coord,
                                        UInt order);
     
  //! Pre-compute all shape functions and their derivative / curl / divergence

  //! This method pre-computes the shape functions and their derivative /curl/ 
  //! divergence for a list of local points. This list typically contains all 
  //! the points at which an integration point is defined by the IntScheme.
  virtual void SetFunctionsAtIp( const StdVector<LocPoint>& iPoints ) = 0;

  //! Pre-compute shape functions and their derivatives at specific points
  
  //! Overloaded method to be able to set only specific points
  virtual void SetFunctionsAtIp( const std::map<Integer, LocPoint >& 
                                 iPoints) = 0;


  //! Get Coordinates of local degrees of freedom
  virtual void GetLocalDOFCoordinates(Matrix<Double> & coordMat){
    Exception("GetLocalDOFCoordinates: Not implemented here...");
  }

  //! @copydoc BaseFE::ComputeMonomialCoefficients
  //! Overloaded method for lagrange Elements
  virtual void ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C){
    Exception("ComputeMonomialCoefficients: Not implemented here...");
  }

protected:
  
  //! Calculate the location of unknowns for a line up to given order
  //! Righ now the calculation is done here but has to move into the 
  //! Integration Scheme class!
  //! This method gets reimplemented in the Spectral classes to ask the 
  //! Integration scheme to provide the points
  void CalcAllSupportingPoints(UInt maxOrder);

  //! Calculate the location of unknowns for a given order
  void CalcSupportingPoints(UInt order);

  //! Calculates the Gauss Lobatto Points for a given order
  //! Is to be moved into Integration Scheme Class
  StdVector<Double> CalcGaussLobattoPoints(UInt order );
  
  //! Stores the Locations of the Element DOFs for a line for every order
  std::map<UInt,StdVector<Double> > supportingPoints_;

}; // class FeNodal

} // end of namespace

#endif // header guard

