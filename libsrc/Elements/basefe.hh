#ifndef FILE_BASEFE_2003
#define FILE_BASEFE_2003

#include <vector>
#include <Matrix/matrix.hh>


namespace CoupledField
{
  //! Base class for description of elements
  /*! From this class all special finite elements as triangles, quadrilaterals,
    tetraheader, etc. are derived. Some Functions are pure virtual.
  */

class BaseFE
{
public:

  //! constructor (does nothing)
  BaseFE();

  //! decstructor 
  virtual ~BaseFE();
 
  //! Get value of all shape fnc at arbitrary local point
  virtual void GetShFnc(std::vector<double> & S, 
			const std::vector<Double> & LCoord);

  
  //! Get value of all shape fnc at integration point ip
  virtual void GetShFncAtIp(std::vector<double> & S, 
			    const Integer ip);

  //! Get global derivatives of all shape fnc at arbitrary local point
  virtual  void GetGlobDerivShFnc(Matrix<Double> & Deriv, 
				  const std::vector<Double> & LCoord,
				  const Matrix<Double> & CornerCoords);

  //! Get global derivatives of all shape fnc at integration point ip
  virtual  void GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				      const Integer ip,
				      const Matrix<Double> & CornerCoords);

  //! Calculation of Jacobian determinant at arbitrary local point
  virtual Double CalcJacobianDet(const std::vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords) = 0;

  //! Calculation of Jacobian determinant at integration point ip
  virtual Double CalcJacobianDetAtIp(const Integer ip, 
				     const Matrix<Double> & CornerCoords) = 0;

  //! Return number of nodes   
  ShortInt GetDim() {return Dim_;}
 
  //! Return number of integration pointes
  ShortInt GetNumNodes() {return NumNodes_;}

  //! Return number of edges
  ShortInt GetNumEdges() {return NumEdges_;}

  //! Return number of faces
  ShortInt GetNumFaces() {return NumFaces_;}

  //! Returns number of integration points
  ShortInt GetNumIntPoints() {return NumIntPoints_;}
  
  //! Return FE-Type for LAS++
  virtual Integer feType()=0;

protected:
  
  //! Calculates the shape functions at an arbitrary local point
  virtual void CalcShapeFnc(std::vector<Double> & LDeriv, 
			    const std::vector<Double> & LCoord) = 0;
  
  //! Calculates the local derivatives of shape functions at an arbitrary local point
  virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				      const std::vector<Double> & LCoord) = 0;

  //! Calculates the Jacobian Matrix at an arbitrary local point
  virtual void CalcJacobian(Matrix<Double> & J, 
				   const std::vector<Double> & LCoord, 
				   const Matrix<Double> & CornerCoords) = 0;

  //! Calculates the Jacobian Matrix at integration point ip
  virtual void CalcJacobianAtIp(Matrix<Double> & J, 
				const Integer ip, 
				const Matrix<Double> & CornerCoords) = 0;

  //! Calculates the Inverse Jacobian Matrix at an arbitrary local point
  virtual void CalcInvJacobian(Matrix<Double> & JInv,
			       const std::vector<Double> & LCoord,
			       const Matrix<Double> & CornerCoords) = 0;
  
  //! Calculates the Inverse Jacobian Matrix at integration point ip
  virtual void CalcInvJacobianAtIp(Matrix<Double> & JInv,
				   const Integer ip,
				   const Matrix<Double> & CornerCoords) = 0;


  ShortInt Dim_;                 //!< space dimension
  ShortInt NumNodes_;            //!< number of nodes
  ShortInt NumEdges_;            //!< number of edges 
  ShortInt NumFaces_;            //!< number of faces
  ShortInt NumIntPoints_;        //!< number of integration points
  ShortInt DegreeInteg_;         //!< numerical integration order
  Matrix<Double> LCornerCoords_;  //!< Matrix of local corner coordinates (x: Ecken, y: Dimension)
  std::vector<Double> * ShFncAtIp_;     //!< Matrix of shape function values at IPs (x:local Dim, y:Number)
  Matrix<Double> * ShFncDerivAtIp_; //!< Array of local derivatives in each integration point
  std::vector<Double> * IntPoints_;     //!< integration points
  std::vector<Double> IntWeights_;  //!< integration weights
  
  enum IntegrationType IntegType;
  
  //! Converts the string used for the integration type to an integer
  enum IntegrationType String2EnumIntegrationType(const Char * inttype);

private:
 
};
 
}

#endif // FILE_BASEFE_2003
