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
  /*! 
    \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param LCoord (input) Local coordinates of evalutation point 
  */
 virtual void GetShFnc(std::vector<double> & S, 
			const std::vector<Double> & LCoord);

  
  //! Get value of all shape fnc at integration point ip
  /*! 
    \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param ip (input) Integration point
  */
  virtual void GetShFncAtIp(std::vector<double> & S, 
			    const Integer ip);

  //! Get global derivatives of all shape fnc at arbitrary local point
  /*! 
    \param S (output) Matrix with global derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,dx} & N_{1,dy} & \cdots \\
                                  N_{2,dx} & N_{2,dy} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param LCoord (input) Local Coordinates of evalutaion point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]       
  */
  virtual  void GetGlobDerivShFnc(Matrix<Double> & Deriv, 
				  const std::vector<Double> & LCoord,
				  const Matrix<Double> & CornerCoords);

  //! Get global derivatives of all shape fnc at integration point ip
  /*! 
    \param S (output) Matrix with global derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
                                  N_{2,d\xi} & N_{2,d\eta} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param ip(input) Integration point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
     \param jacDet jacobian determinant
  */
  virtual  void GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				      const Integer ip,
				      const Matrix<Double> & CornerCoords,
				      Double & jacDet);


  //! Get global derivatives of all shape fnc at integration point ip
  /*! 
    \param S (output) Matrix with global derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
                                  N_{2,d\xi} & N_{2,d\eta} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param ip(input) Integration point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
  */
  virtual  void GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				      const Integer ip,
				      const Matrix<Double> & CornerCoords);


  //! Calculation of Jacobian determinant at arbitrary local point
  /*! 
    \param LCoord (input) Local Coordinates of evaluation point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
  */
  virtual Double CalcJacobianDet(const std::vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords);

  //! Calculation of Jacobian determinant at integration point ip
  /*! 
    \param ip (input) Integration point
    \param CornerCoords (input)
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
  */
  virtual Double CalcJacobianDetAtIp(const Integer ip, 
				     const Matrix<Double> & CornerCoords);

  //! Return number of nodes   
  ShortInt GetDim() const {return Dim_;}
 
  //! Return number of integration pointes
  ShortInt GetNumNodes() const {return NumNodes_;}

  //! Return number of edges
  ShortInt GetNumEdges() const {return NumEdges_;}

  //! Return number of faces
  ShortInt GetNumFaces() const {return NumFaces_;}

  //! Returns number of integration points
  ShortInt GetNumIntPoints() const {return NumIntPoints_;}
  
  //! Return FE-Type for LAS++
  virtual Integer feType()=0;

  /// Return weightings of integration points
  std::vector<Double> GetIntWeights() const {return IntWeights_;};
  
  
protected:
  
  //! Calculates the shape functions at an arbitrary local point
  /*!
    \param Shape (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void CalcShapeFnc(std::vector<Double> & Shape, 
			    const std::vector<Double> & LCoord) = 0;
  
  //! Calculates the local derivatives of shape functions at an arbitrary local point
  /*!
    \param LDeriv (output) Matrix with local derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
                                  N_{2,d\xi} & N_{2,d\eta} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				      const std::vector<Double> & LCoord) = 0;

  //! Calculates the Jacobian Matrix at an arbitrary local point
  /*!
    \param J (output) Jacobian Matrix
    \f[ J = \left( \begin{array}{ccc} x_{\xi} & x_{\eta} & x_{\zeta} \\
                                    y_{\xi} & y_{\eta} & y_{\zeta}\\
                                    z_{\xi} & z_{\eta} & z_{\zeta} \end{array}\right)\f] 
    \param LCoord (input) Local coorindates of evaluation point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f] 
  */
  virtual void CalcJacobian(Matrix<Double> & J, 
				   const std::vector<Double> & LCoord, 
				   const Matrix<Double> & CornerCoords);

  //! Calculates the Jacobian Matrix at integration point ip
  /*!
    \param J (output) Jacobian Matrix
    \f[ J = \left( \begin{array}{ccc} x_{\xi} & x_{\eta} & x_{\zeta} \\
                                    y_{\xi} & y_{\eta} & y_{\zeta}\\
                                    z_{\xi} & z_{\eta} & z_{\zeta} \end{array}\right)\f]  
    \param ip (input) Integration point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f] 
  */
  virtual void CalcJacobianAtIp(Matrix<Double> & J, 
				const Integer ip, 
				const Matrix<Double> & CornerCoords);

  //! Calculates the Inverse Jacobian Matrix at an arbitrary local point
  /*!
    \param JInv (output) Inverse Jacobian Matrix 
    \f[ J^{-1} = \left( \begin{array}{ccc} \xi_{x} & \xi_{y} & \xi_{z} \\
                                    \eta_{x} & \eta_{y} & \eta_{z}\\
                                    \zeta_{x} & \zeta_{y} & \zeta_{z} \end{array}\right)\f]
    \param LCoord (input) Local coorindates of evaluation point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f] 
  */
  virtual void CalcInvJacobian(Matrix<Double> & JInv,
			       const std::vector<Double> & LCoord,
			       const Matrix<Double> & CornerCoords);
  
  //! Calculates the Inverse Jacobian Matrix at integration point ip
  /*!
    \param JInv (output) Inverse Jacobian Matrix 
    \f[ J^{-1} = \left( \begin{array}{ccc} \xi_{x} & \xi_{y} & \xi_{z} \\
                                    \eta_{x} & \eta_{y} & \eta_{z}\\
                                    \zeta_{x} & \zeta_{y} & \zeta_{z} \end{array}\right)\f]
    \param ip (input) Integration point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f] 
  */
  virtual void CalcInvJacobianAtIp(Matrix<Double> & JInv,
				   const Integer ip,
				   const Matrix<Double> & CornerCoords);


 //! Set value of shape fnc at integration points
  virtual void SetShapeFncAtIp();

  
  //! Set value of shape fnc derivatives at integration points
  virtual void SetShapeFncDerivAtIp();

  //! Set integration points
  virtual void SetIntPoints()=0;


  ShortInt Dim_;                    //!< space dimension
  ShortInt NumNodes_;               //!< number of nodes
  ShortInt NumEdges_;               //!< number of edges 
  ShortInt NumFaces_;               //!< number of faces
  ShortInt NumIntPoints_;           //!< number of integration points
  ShortInt DegreeInteg_;            //!< numerical integration order
  Matrix<Double> LCornerCoords_;    //!< Matrix of local corner coordinates (x:number, y:Dim)
  std::vector<Double> * ShFncAtIp_; //!< Array of vectors of function values at IPs (x:local Dim, y:Number)
  Matrix<Double> * ShFncDerivAtIp_; //!< Array of local derivatives in each integration point
  std::vector<Double> * IntPoints_; //!< integration points
  std::vector<Double> IntWeights_;  //!< integration weights
  
  enum IntegrationType IntegType;
  
  //! Converts the string used for the integration type to an integer
  enum IntegrationType String2EnumIntegrationType(const Char * inttype);

private:
 
};
 
}

#endif // FILE_BASEFE_2003
