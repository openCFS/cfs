#ifndef FILE_BASEFE_2003
#define FILE_BASEFE_2003

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
  /*!
    \param S (output) Vector with ShapeFncValues  \f$ (N_{1},\cdots,N_{NumNodes})^T \f$
    \param LCoord (input) Vector with local evalulation Coordinates
  */
  virtual void GetShFnc(Vector<double> & S, 
			const Vector<Double> & LCoord);

  
  //! Get value of all shape fnc at integration point ip
  /*!
    \param S (output) Vector with ShapeFncValues  \f$ (N_{1},..,N_{NumNodes})^T \f$
    \param ip (input) Number of integration point for evaluattion
  */
  virtual void GetShFncAtIp(Vector<double> & S, 
			    const Integer ip);

  //! Get global derivatives of all shape fnc at arbitrary local point
  virtual  void GetGlobDerivShFnc(Matrix<Double> & Deriv, 
				  const Vector<Double> & LCoord,
				  const Matrix<Double> & CornerCoords);

  //! Get global derivatives of all shape fnc at integration point ip
  virtual  void GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				      const Integer ip,
				      const Matrix<Double> & CornerCoords);

  //! Calculation of Jacobian determinant at arbitrary local point
  virtual Double CalcJacobianDet(const Vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords) = 0;

  //! Calculation of Jacobian determinant at integration point ip
  virtual Double CalcJacobianDetAtIp(const Integer ip, 
				     const Matrix<Double> & CornerCoords) = 0;

  //! Return number of nodes   
  ShortInt GetDim() {return Dim;}
 
  //! Return number of integration pointes
  ShortInt GetNumNodes() {return NumNodes;}

  //! Return number of edges
  ShortInt GetNumEdges() {return NumEdges;}

  //! Return number of faces
  ShortInt GetNumFaces() {return NumFaces;}

  //! Returns number of integration points
  ShortInt GetNumIntPoints() {return NumIntPoints;}
  
  //! Return FE-Type for LAS++
  virtual Integer feType()=0;

protected:
  
  //! Calculates the shape functions at an arbitrary local point
  virtual void CalcShapeFnc(Vector<Double> & LDeriv, 
			    const Vector<Double> & LCoord) = 0;
  
  //! Calculates the local derivatives of shape functions at an arbitrary local point
  virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				      const Vector<Double> & LCoord) = 0;

  //! Calculates the Jacobian Matrix at an arbitrary local point
  virtual void CalcJacobian(Matrix<Double> & J, 
				   const Vector<Double> & LCoord, 
				   const Matrix<Double> & CornerCoords) = 0;

  //! Calculates the Jacobian Matrix at integration point ip
  virtual void CalcJacobianAtIp(Matrix<Double> & J, 
				const Integer ip, 
				const Matrix<Double> & CornerCoords) = 0;

  //! Calculates the Inverse Jacobian Matrix at an arbitrary local point
  virtual void CalcInvJacobian(Matrix<Double> & JInv,
			       const Vector<Double> & LCoord,
			       const Matrix<Double> & CornerCoords) = 0;
  
  //! Calculates the Inverse Jacobian Matrix at integration point ip
  virtual void CalcInvJacobianAtIp(Matrix<Double> & JInv,
				   const Integer ip,
				   const Matrix<Double> & CornerCoords) = 0;


  ShortInt Dim;                 //!< space dimension
  ShortInt NumNodes;            //!< number of nodes
  ShortInt NumEdges;            //!< number of edges 
  ShortInt NumFaces;            //!< number of faces
  ShortInt NumIntPoints;        //!< number of integration points
  ShortInt DegreeInteg;         //!< numerical integration order
  Matrix<Double> LCornerCoords;  //!< Matrix of local corner coordinates (x: Ecken, y: Dimension)
  Vector<Double> * ShFncAtIp;     //!< Matrix of shape function values at IPs (x:local Dim, y:Number)
  Matrix<Double> * ShFncDerivAtIp; //!< Array of local derivatives in each integration point
  Vector<Double> * IntPoints;     //!< integration points
  Vector<Double> IntWeights;  //!< integration weights
  
  enum IntegrationType IntegType;
  
  //! Converts the string used for the integration type to an integer
  enum IntegrationType String2EnumIntegrationType(const Char * inttype);

private:
 
};
 
}

#endif // FILE_BASEFE_2003
