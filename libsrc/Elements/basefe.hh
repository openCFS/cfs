#ifndef FILE_BASEFE_2003
#define FILE_BASEFE_2003

#include <vector>
#include <Matrix/matrix.hh>


#ifdef USE_OLAS
#include <olas.hh>
#else
#include <multigrid.hh>
#endif

namespace CoupledField
{
  //! Base class for description of elements
  /*! From this class all special finite elements as triangles, quadrilaterals,
    tetraheader, etc. are derived. Some Functions are pure virtual.
  */

class BaseFE
{
public:
  // Due to (Re-) initialization of elements (e.g. in reduced integration)
  friend class Grid;
  
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


  //! Get local coordinates of element corners 
  /*! 
    \param lCornerCoords (output) local coordinates of element corners
  */
  virtual void GetLocalCornerCoords(Matrix<Double>& lCornerCoords)
  { lCornerCoords = LCornerCoords_;};
  
  

  
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
    \param deriv (output) Matrix with global derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,dx} & N_{1,dx} & \cdots \\
                                  N_{2,dx} & N_{2,dy} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param ip(input) Integration point
    \param cornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
     \param jacDet jacobian determinant
  */
  virtual  void GetGlobDerivShFncAtIp(Matrix<Double> & deriv, 
				      const Integer ip,
				      const Matrix<Double> & cornerCoords,
				      Double & jacDet);



  //! Get global derivatives of all shape fnc at integration point ip
  /*! 
    \param S (output) Matrix with global derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,dx} & N_{1,dy} & \cdots \\
                                  N_{2,dx} & N_{2,dy} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param ip(input) Integration point
    \param CornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
  */
  virtual  void GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				      const Integer ip,
				      const Matrix<Double> & CornerCoords);

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
  //! Calculates a measure for the geometric distortion of an element
  /*!
    \param cornerCoords (input) Corner coordinates of the element
    \param size (input) Absolute size of element in all dimensions
    \param displacement (input) Displacement of the corner points (same ordering as CornerCoords!!)
  */
  virtual Double CalcMeanStrain(Matrix<Double> &cornerCoords,
				Array<Double> &displacements)
  {
    Error("Not implemented",__FILE__,__LINE__);
    return 0;
  }


  //! Return space dimension
  ShortInt GetDim() const {return Dim_;}
 
  //! Get integration points
  std::vector<Double> * GetIntPoints() {return IntPoints_;}

  //! Return number of nodes
  ShortInt GetNumNodes() const {return NumNodes_;}

  //! Return number of edges
  ShortInt GetNumEdges() const {return NumEdges_;}

  //! Return number of faces
  ShortInt GetNumFaces() const {return NumFaces_;}

  //! Returns number of integration points
  ShortInt GetNumIntPoints() const {return NumIntPoints_;}
  
  //! Return FE-Type for LAS++
#ifdef USE_OLAS
  virtual FEType feType()=0;
#else
  virtual Integer feType()=0;
#endif

  /// Return weightings of integration points
  std::vector<Double> GetIntWeights() const {return IntWeights_;};
  
  // return number of childs in refinement
  virtual Integer getNumChilds() const { return numChilds_;}

  // ============================= methods for edge elements =======================
  // =================================================================================

 //! Calculates the local derivatives of the edge shape functions at an arbitrary local point
  /*!
    \param deriv (output) Vector of matrices with local derivatives of all shape functions.
                  Every matrix stores a complete set of global derivations
    \f[ LDeriv_1 = \left( \begin{array}{ccc} N_{1\xi,d\xi} & N_{1\xi,d\eta} & N_{1\xi,d\zeta}\\
                                             N_{1\eta,d\xi} & N_{1\eta,d\eta} & N_{1\eta,d\zeta}\\
					     N_{1\zeta,d\xi} & N_{1\zeta,d\eta} & N_{1\zeta,d\zeta}
					     \end{array}\right) \f]
    \param LCoord (input) Local coordinates of evalutation point 
  */
  virtual void GetEdgeGlobalDerivShapeFnc(std::vector<Matrix<Double>* > & deriv, 
					  const std::vector<Double> & LCoord,
					  const Matrix<Double> & CornerCoords)
  { Error("GetEdgeGlobDerivShFnc called for non edge element! ",__FILE__,__LINE__);};



  //! Get global derivatives of all edge shape functions at integration point ip
  /*! 
    \param deriv (output) Matrix with global derivatives of all shape functions
    \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
                                  N_{2,d\xi} & N_{2,d\eta} & \cdots \\
                                  \cdots     & \cdots      & \cdots \end{array}\right) \f]
    \param ip(input) Integration point
    \param cornerCoords (input) Coordinates of element corners
    \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
                                  \cdots & \cdots & \cdots \end{array} \right) \f]
  */
  virtual void GetEdgeGlobDerivShFncAtIp(std::vector< Matrix<Double> *> & deriv, 
					 const Integer ip,
					 const Matrix<Double> & cornerCoords);
  



  //! Calculates the edge shape functions at an arbitrary local point
  /*!
    \param shape (output) Matrix of shape shape values (size: nrEdges x nrSpaceDim) 
  \f[ \left( \begin{array}{c} E_1 \\ E_2 \\ \cdots \end{array} \right) = 
   \left( \begin{array}{ccc} N_{1\xi} & N_{1\eta} & N_{1\zeta} \\
                             N_{2\xi} & N_{2\eta} & N_{2\zeta} \\
                             \cdots     & \cdots      & \cdots 
	  \end{array}\right) \f]
    \param lCoord (input) Local coordinates of evalutation point 
  */
 virtual void CalcEdgeShapeFnc(Matrix<Double> & shape, 
				const std::vector<Double> & LCoord, 
				const Matrix<Double> & CornerCoords)
  { Error("CalcEdgeShapeFnc called for non edge element! ",__FILE__,__LINE__);};
  
  

  //! Calculates the edge shape functions at an integration point
  /*!
    \param shape (output) Matrix of shape shape values (size: nrEdges x nrSpaceDim) 
    \f[ \left( \begin{array}{c} E_1 \\ E_2 \\ \cdots \end{array} \right) = 
    \left( \begin{array}{ccc} N_{1\xi} & N_{1\eta} & N_{1\zeta} \\
                              N_{2\xi} & N_{2\eta} & N_{2\zeta} \\
                              \cdots     & \cdots      & \cdots 
   	  \end{array}\right) \f]
    \param ip (input) Number of desired integration point
    \param cornerCoords (input) Coordinates of edge points 
  */
  virtual void CalcEdgeShapeFncAtIp(Matrix<Double> & shape, 
				    const Integer ip,
				    const Matrix<Double> & cornerCoords);
  

  //! Get nodes belonging to one edge
  /*! 
    \param edges (output) Matrix of assignment of edge index to according nodes 
    (dimension: nrEdges x 2 (every edge has two nodes))
  */
 virtual void GetEdgeVertices(Matrix<Integer> & edges){edges = edgeVertices_;};


  //! Get global edge numbers
  /*! 
    \param edges (output) Vector of global edge numbers
    \param pDENodes (input) Global index of nodes belonging to one element
    \param algsys (input) Pointer to the algebraic system
  */
  virtual void GetGlobalEdgeIndices(std::vector<Integer>& edges, Integer * pDENodes, BaseSystem * algsys);


  //! Get global coordinates based on local element coordinates
  /*! 
    \param globCoord (output) Vector of global coordinates
    \param ip (input) Integeration point at which global coord has to be calculated
    \param cornerCoords (input) Matrix of corner coordinates
  */
  virtual void GetGlobalEdgeIndicesAtIP( std::vector<Double> & globCoord,
					 Integer ip,
					 const Matrix<Double> & cornerCoords);
  

  //! Set integration type (GaussOrder1, GaussOrder2, ...)
  virtual void SetIntegrationType(const IntegrationType aIntType)
  {IntegType = aIntType;};
  
    
  //! Get integration type (GaussOrder1, GaussOrder2, ...)
  virtual IntegrationType GetIntType() {return IntegType;};



protected:
   //! Define variables of this class
  virtual void Init()
  { Error("Init not implemented for this element type! ",__FILE__,__LINE__);};

  /// defines the connected nodes with every edge 
  void SetEdgeVertices()
  { Error("SetEdgeVertices called for non edge element! ",__FILE__,__LINE__);};
  
  
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

 //! Set value of shape fnc at integration points
  virtual void SetShapeFncAtIp();

  
  //! Set value of shape fnc derivatives at integration points
  virtual void SetShapeFncDerivAtIp();


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
  Integer             numChilds_;   //!< number of children for element in refinement

  enum IntegrationType IntegType;
  
  //! Converts the string used for the integration type to an integer
  enum IntegrationType String2EnumIntegrationType(const Char * inttype);


  /// Matrix with connected nodes to the proper edge
  /*!
    \f[ \left( \begin{array}{c} E_1 \\ E_2 \\ \cdots \end{array} \right) = 
    \left( \begin{array}{cc} node_1^{E1} & node_2^{E1} \\
                             node_1^{E2} & node_2^{E2}  \\
                              \cdots & \cdots 
			      \end{array} \right) \f]
  */
  Matrix<Integer> edgeVertices_; 
};
 
}

#endif // FILE_BASEFE_2003
