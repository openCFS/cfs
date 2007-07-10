// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEFE_2003
#define FILE_BASEFE_2003

#include "Matrix/matrix.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Utils/StdVector.hh"
#include "Domain/ansatzFct.hh"
#include <map>

namespace CoupledField
{
  // Forward class declaration
  class BaseSystem;
  class Elem;

  //! Base class for description of elements

  //! From this class all special finite elements as triangles, quadrilaterals,
  //! tetraheader, etc. are derived. Some methods are purely virtual.
  class BaseFE {

  public:
    // Due to (Re-) initialization of elements (e.g. in reduced integration)
    friend class Grid;
#ifdef INTEGLIB
    friend class ElemIntegr;
#endif
  
    //! constructor (does nothing)
    BaseFE();

    //! decstructor 
    virtual ~BaseFE();

    //! Get value of all shape fnc at arbitrary local point
    /*! 
      \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void GetShFnc(Vector<Double> & S, 
                          const Vector<Double> & LCoord, 
                          const Elem* elem,
                          UInt dof = 1 );


    //! Get local coordinates of element corners 
    /*! 
      \param lCornerCoords (output) local coordinates of element corners
    */
    virtual void GetLocalCornerCoords(Matrix<Double>& lCornerCoords)
    { lCornerCoords = LCornerCoords_;};
  
  
    //! Get the global coordinate for a given local one
    //! \param globCoord (output) global coordinate
    //! \param locCoord (input) local coordinate
    //! \param coordMat (input) global corner coordinates of element
    //!                         (spaceDim \f$\times\f$ nrNodes)
    virtual void Local2GlobalCoord(Vector<Double> & globCoord,
                                   const Vector<Double> & locCoord,
                                   const Matrix<Double> & coordMat,
                                   const Elem* elem );

  
    //! Get the local coordinates for given global ones
    //! \param localCoords (output) local coordinates
    //! \param globalCoords (input) global coordinates
    //! \param coordMat (input) global corner coordinates of element
    //!                         (spaceDim \f$\times\f$ nrNodes)
    virtual void Global2LocalCoords(Matrix<Double> & localCoords,
                                    const Matrix<Double> & globalCoords,
                                    const Matrix<Double> & coordMat);

    //! Returns wether a given local coordinate is
    //! within this element
    //! \param localCoords (input) local coordinates
    //! \param coordsInside (output) which local coordinates are inside
    virtual void CoordsInsideElem(const Matrix<Double> & localCoords,
                                  const Double tolerance,
                                  StdVector<bool> & coordsInside) const
    { Error( "Not implemented", __FILE__, __LINE__); }


    //! Get value of all shape fnc at integration point ip
    /*! 
      \param S (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param ip (input) Integration point
    */
    virtual void GetShFncAtIp(Vector<Double> & S, 
                              const UInt ip, 
                              const Elem * elem, 
                              UInt dof = 1);
    

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
                                    const Vector<Double> & LCoord,
                                    const Matrix<Double> & CornerCoords,
                                    const Elem * elem, 
                                    UInt dof = 1);


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
                                        const UInt ip,
                                        const Matrix<Double> & cornerCoords,
                                        Double & jacDet,
                                        const Elem * elem, 
                                        UInt dof = 1);



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
                                        const UInt ip,
                                        const Matrix<Double> & CornerCoords,
                                        const Elem * elem, 
                                        UInt dof = 1 );

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
                              const Vector<Double> & LCoord, 
                              const Matrix<Double> & CornerCoords,
                              const Elem* elem );
    
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
                                  const UInt ip, 
                                  const Matrix<Double> & CornerCoords,
                                  const Elem* elem );

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
                                 const Vector<Double> & LCoord,
                                 const Matrix<Double> & CornerCoords,
                                 const Elem* elem);
  
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
                                     const UInt ip,
                                     const Matrix<Double> & CornerCoords,
                                     const Elem* elem);


    //! Calculation of Jacobian determinant at arbitrary local point
    /*! 
      \param LCoord (input) Local Coordinates of evaluation point
      \param CornerCoords (input) Coordinates of element corners
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
      \cdots & \cdots & \cdots \end{array} \right) \f]
    */
    virtual Double CalcJacobianDet(const Vector<Double> & LCoord,
                                   const Matrix<Double> & CornerCoords,
                                   const Elem* elem);

    //! Calculation of Jacobian determinant at integration point ip
    /*! 
      \param ip (input) Integration point
      \param CornerCoords (input)
      \f[ \left( \begin{array}{ccc} x_{1} & x_{2} & \cdots \\ y_{1} & y_{2} & \cdots \\
      \cdots & \cdots & \cdots \end{array} \right) \f]
    */
    virtual Double CalcJacobianDetAtIp(const UInt ip, 
                                       const Matrix<Double> & CornerCoords,
                                       const Elem* elem);

    //! Calculation the volume (area, length) of an element
    /*! 
      \param CornerCoords (input)
      \param isaxi (true: axisymmetric formulation)
    */
    virtual Double CalcVolume(const Matrix<Double> & CornerCoords, const bool isaxi);
    
    /** Calculates the barycenter of the element. For 3D and 2D (z=0 then) 
     * This is the average of every coordinate. 
     * @param coords might be 2D and 3D. See Grid::GetElemNodesCoord()
     * @param barycenter the output. In 2D z=0. If you want 2D then overload */
    static void CalcBarycenter(const Matrix<Double>& coords, Point& barycenter);


    //! Calculates a measure for the geometric distortion of an element
    /*!
      \param cornerCoords (input) Corner coordinates of the element
      \param size (input) Absolute size of element in all dimensions
      \param displacement (input) Displacement of the corner points (same ordering as CornerCoords!!)
    */
    virtual Double CalcMeanStrain(Matrix<Double> &cornerCoords,
                                  Matrix<Double> &displacements)
    {
      Error("Not implemented",__FILE__,__LINE__);
      return 0;
    }


    //! Calculates corresponding volume point of neighbouring surfaces
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
    virtual void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                           const StdVector<UInt> & volConnect,
                                           const Vector<Double> & surfIntPoint,
                                           Vector<Double> & volIntPoint)
    {
      Error( "BaseFE::GetLocalIntPoints4Surface: Not implemented here", 
             __FILE__, __LINE__);}

    virtual void CoordTrans( const Matrix<Double> &ptCoord, 
			     Matrix<Double> &TransMat, 
                             Matrix<Double> &ShellCoord );

    virtual void CoordTrans2D( const Matrix<Double> &ptCoord, 
			     Matrix<Double> &TransMat, 
                             Matrix<Double> &ShellCoord );

    //! Get total number of dofs for particular dof
    virtual UInt GetNumFncs( const shared_ptr<AnsatzFct>& fncType );

    //! Return space dimension
    UInt GetDim() const {return Dim_;}
 
    //! Get integration points
    Vector<Double> * GetIntPoints() {return IntPoints_;}

    //! Return number of nodes
    UInt GetNumNodes() const {return NumNodes_;}

    //! Return number of edges
    UInt GetNumEdges() const {return NumEdges_;}

    //! Return number of faces
    UInt GetNumFaces() const {return NumFaces_;}

    //! Return number of corners
    UInt GetNumCorners() const {return NumCorners_;}

    //! Returns number of integration points
    UInt GetNumIntPoints() const {return NumIntPoints_;}
  
    //! 
    void GetCoordMidPoint(Vector<Double> & coord) 
    {coord = MidPoint_;};

    //! compute length of edge with maximal/minimal size
    virtual void GetMaxMinEdgeLength( Matrix<Double> &ptCoord, Double &Lmax, Double &Lmin ) 
    { 
      Error("GetMaxMinEdgeLength not implemented", __FILE__, __LINE__);
    }

    //! Return FE-Type
    virtual FEType feType() = 0;

    /// Return weightings of integration points
    Vector<Double> GetIntWeights() const {return IntWeights_;};
  
    // return number of childs in refinement
    virtual UInt getNumChilds() const { return numChilds_;}


    // ========================================================================
    // ======================== methods for edge elements =====================
    // ========================================================================

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
    virtual void GetEdgeGlobalDerivShapeFnc(StdVector<Matrix<Double>* > & deriv, 
                                            const Vector<Double> & LCoord,
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
    virtual void GetEdgeGlobDerivShFncAtIp(StdVector< Matrix<Double> *> & deriv, 
                                           const UInt ip,
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
                                  const Vector<Double> & LCoord, 
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
                                      const UInt ip,
                                      const Matrix<Double> & cornerCoords);
  

    //! Get nodes belonging to one edge
    /*! 
      \param edges (output) Matrix of assignment of edge index to according nodes 
      (dimension: nrEdges x 2 (every edge has two nodes))
    */
    virtual void GetEdgeVertices(Matrix<UInt> & edges){edges = edgeVertices_;};


    //! Get the indices in the connect array belonging to given edge
    virtual void GetEdgeIndices( StdVector<UInt>& indices, UInt edgeNr ) {
      indices = edgeIndices_[edgeNr]; }

    //! Get the indices in the connect array belonging to given face
    virtual void GetFaceIndices( StdVector<UInt>& indices, UInt faceNr ) {
      indices = faceIndices_[faceNr]; }

    //! Get global edge numbers
    /*! 
      \param edges (output) Vector of global edge numbers
      \param pDENodes (input) Global index of nodes belonging to one element
      \param algsys (input) Pointer to the algebraic system
    */
    virtual void GetGlobalEdgeIndices(StdVector<UInt>& edges, UInt * pDENodes, BaseSystem * algsys);


    //! Get global coordinates based on local element coordinates
    /*! 
      \param globCoord (output) Vector of global coordinates
      \param ip (input) UIntation point at which global coord has to be calculated
      \param cornerCoords (input) Matrix of corner coordinates
    */
    virtual void GetGlobalEdgeIndicesAtIP( Vector<Double> & globCoord,
                                           UInt ip,
                                           const Matrix<Double> & cornerCoords);
  

    /** sets the element to default reduced integration mode
     * See the comment of SetStandardIntegration */
    void SetReducedIntegration()
    {
       ENTER_FCN( "BaseFE:SetReducedIntegration():" );
       SetDefaultReducedIntegration();
       //   SetIntPoints();
    }   

    /** sets the element to standard integration 
     * It would be smarter to save the old value in SetReducedIntegration()
     * and do a switch back here! */
    void SetStandardIntegration()
    {
       ENTER_FCN( "BaseFE:SetStandardIntegration():" );
       SetDefaultIntegration();
       SetIntPoints(IntegMethod, IntegOrder);
    }   
    
    // set incompatible modes cabaility to yes
    void UseICModes() {
      ICModes_ = true;
    } 

    //! compute now just with incompatible modes
    void SetCalcICModes() {
      CalcICModes_ = true;
    }

    //! compute now with standard basis functions
    void ResetCalcICModes() {
      CalcICModes_ = false;
    }

    //! a public helper method that dumps the current content of the map
    void DumpIntegrationPointsMap();


    // =======================================================================
    // L E G E N D R E    P A R T
    // =======================================================================

    //! Get number of functions for this element and this dof
    void virtual GetNumFncs(Vector<UInt>& numFcns, 
                    const shared_ptr<AnsatzFct>& fcnType, 
                    AnsatzFct::FctEntityType fctEntityType, 
                    UInt dof = 1);

    //! Evaluate polynom and its derivative using Honer's algorithm
    void EvalPolynom( Double& value, Double& deriv,
                      const UInt order, const Double* coeff, 
                      const Double xVal );

    //! Coefficients of 1D Legendre coefficients up to order 9
    static Double lCoeff_[9][10];

    //! Set current ansatz-fct object
    virtual void SetAnsatzFct( shared_ptr<AnsatzFct>& actFct,
                               bool setIntPoints = true ) {
    }

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
    virtual void CalcShapeFnc(Vector<Double> & Shape, 
                              const Vector<Double> & LCoord,
                              const Elem* elem , UInt dof,
                              AnsatzFct::FctEntityType = AnsatzFct::ALL ) = 0;
  
    //! Calculates the local derivatives of shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
      N_{2,d\xi} & N_{2,d\eta} & \cdots \\
      \cdots     & \cdots      & \cdots \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                        const Vector<Double> & LCoord,
                                        const Elem* elem, UInt dof,
                                        AnsatzFct::FctEntityType = AnsatzFct::ALL ) = 0;

    //! Calculates the local derivatives of incompatible mode shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
      N_{2,d\xi} & N_{2,d\eta} & \cdots \\
      \cdots     & \cdots      & \cdots \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalICModesDerivShapeFnc(Matrix<Double> & LDeriv, 
					       const Vector<Double> & LCoord,
					       const Elem* elem, UInt dof,
					       AnsatzFct::FctEntityType = AnsatzFct::ALL ) {
      Error("CalcLocalICModesDerivShapeFnc not implemented",__FILE__,__LINE__);
    };

    //! Set value of shape fnc at integration points
    virtual void SetShapeFncAtIp();

  
    //! Set value of shape fnc derivatives at integration points
    virtual void SetShapeFncDerivAtIp();

    //! Set local corner coordinates
    virtual void SetCornerCoords() = 0;

    /** returns the shape name as used for integration type in the XML file based on
     * the abstract feType() */
    const char* GetShapeName(); 

    /** the childs fill here the integration points map via AddIntegrationPoints() */    
    virtual void FillIntegrationPoints() = 0;

    /** Every Element has to set the preferred default integration mode */
    virtual void SetDefaultIntegration() = 0;

    /** Every Element has also a default reduced integration set -> see the tex documentation for details */
    virtual void SetDefaultReducedIntegration() = 0;

    /** Expicitly set and load the integration type and order */
    void SetIntPoints(IntegrationMethod method, int order);

    /** The child classes add here their integration point data to the elements map. Called in the constructors.
     * Too avoid errors, any key (type+order) must be unique, otherwise exit!
     * @param type e.g. ECONOMICAL
     * @param order the polynomial order
     * @param numberOfRows the first dimension of data
     * @param data actually [][2] for 1D (coord + weight) and [][3] for 2D and  [][4] for 3D */
    void AddIntegrationPoints(IntegrationMethod type, int order, int numberOfRows, Double* data);
     

    //! Helper function for printing a coordinate matrix in a string
    std::string CoordMatrix2String(const Matrix<Double> & coordMat);

    /** This is the common Init/Constructor for all element classes! 
     *  It calls the (virtual) inititialization methods.
     *  To be called AFTER all individual attributes (Dim_, ...) are set! 
     *  The parametes are only used in cartesian integration points mode, leave blank as normal user.
     * BaseFE::CommonInit() sets BaseFE::IntegMethod and BaseFE::IntegOrder via the 
     * XML file (see above) or from the virtual SetDefaultIntegration(). These two 
     * attributes should always represent the current setting (NumIntPoints_, 
     * IntPoints_ and IntWeights_) - even when a GetIntegrationPoints() choosed 
     * alternative values.
     * @param method only for LineFE in cartesian usage - leave to default in any other case
     * @param order same as for parameter method applies. */
    void CommonInit(IntegrationMethod method, int order);  

    void CommonInit()
    {
        CommonInit(UNDEFINED, 0);
    }  
    
    /** Creates the integration points by cartesian product of 1D for rectangle and cube.
     * The creation is quite expensive, but the results are cached in THIS element!
     * For rectangle, hexaeder (cube) and wedge practically arbitrary high integration orders
     * can be constructed via the cartesian product rule. For the hexaeder the ECONOMICAL
     * ones are much more efficient! For anisotropic elements, RectangeFE::SetCartesianIntegration()
     * and HexaFE::SetCartesianIntegration() creates the elements and stores them in
     * the map or simple loads it from the map. The internal method type is CARTESIAN and 
     * the order is zyx where z is the order of the coordinate x3 with z in [1;9], ...
     * This can also be set for testing in the XML file, always set for three coordinates (order > 111).
     * 
     * When debugging or altering code be carefull, the calling sequence of
     * SetCartesianIntegration() is not really straight forward. :(
     * 
     * @param order3 if smaller 1 this dimension is omitted (rectangle)
     * @param create_only false or leave to default for any user, true only for internal purpose (recursive call) */
    void SetCartesianInteg(int order1, int order2, int order3, bool create_only);



    UInt Dim_;                    //!< space dimension
    UInt NumNodes_;               //!< number of nodes
    UInt NumEdges_;               //!< number of edges 
    UInt NumFaces_;               //!< number of faces
    UInt NumCorners_;             //!< number of corners
    UInt NumIntPoints_;           //!< number of integration points
    Vector<Double> MidPoint_;         //!< coordinate of midpoint (for 1. order integration)
    Matrix<Double> LCornerCoords_;    //!< Matrix of local corner coordinates (x:number, y:Dim)
    Vector<Double> * ShFncAtIp_;      //!< Array of vectors of function values at IPs (x:local Dim, y:Number)
    Matrix<Double> * ShFncDerivAtIp_; //!< Array of local derivatives in each integration point
    Matrix<Double> * ShFncICModesDerivAtIp_; //!< Array of local derivatives of incomp. modes in each integration point
    Vector<Double> * IntPoints_;      //!< integration points
    Vector<Double> IntWeights_;       //!< integration weights
    UInt numChilds_;               //!< number of children for element in refinement

    bool ICModes_; //yes, if we use incompatible modes
    bool CalcICModes_; //yes, then we do the computations just with incompatible modes

    //! AddIntegrationPoints stores here all availabe integrations
    std::map<const std::string, StdVector<Double*>*> IntegrationPointsMap_; 

    //! Vectors with node indices of each edge
    StdVector<UInt> * edgeIndices_;

    //! Vectors with node indices of each face
    StdVector<UInt> * faceIndices_;

    enum IntegrationMethod IntegMethod; //! set in XML and SetStandard/ReducedIntegration

    UInt IntegOrder;                    //! accompanies IntegMethod 
  
    /// Matrix with connected nodes to the proper edge
    /*!
      \f[ \left( \begin{array}{c} E_1 \\ E_2 \\ \cdots \end{array} \right) = 
      \left( \begin{array}{cc} node_1^{E1} & node_2^{E1} \\
      node_1^{E2} & node_2^{E2}  \\
      \cdots & \cdots 
      \end{array} \right) \f]
    */
    Matrix<UInt> edgeVertices_; 
    

    //! Current ansatz function ( which is used for current ansatz functions)
    shared_ptr<AnsatzFct> actFct_;

    //! Actual number of ansatz functions
    Integer actNumFcns_;
    
    /** Helper method that generates the key for the map */
    void MakeKey(IntegrationMethod type, int order, std::string &out);
    
    /** generate a key for a cartesian integration
     * @param line3 optional, use < 1 for 2D */
    void MakeKey(int order1, int order2, int order3, std::string &out);

    /** Helper function for SetIntPoints. Makes a comfortable search and searches for default when fallback is true
     * @param fallback caller do true, false for internal recursive search
     * @return when fallback is true not null but exit with message */
    StdVector<Double*>* GetIntegrationPoints(IntegrationMethod type, int order, bool search_upwards, bool search_downwards, bool fallback);
    
    /** Creates the integration points via cartesian product for 2D and 3D elements. The result is stored in the map.
      *  @param element1 IntegMethod and IntegOrder shall be set
      *  @param element2 IntegMethod and IntegOrder shall be set
      *  @param element3 the only optional element */
    void CreateCartesianIntegration(BaseFE* line1, BaseFE* line2, BaseFE* line3);
    
   
    /** Private Helper Method that dumps the currently set integration points */
    void DumpIntegrationPoints();    

    /** Sets the integration points for the data filled by the virtual FillIntegrationPoints.
     * This mehtod reads the data from the map with the help of the dimension */
    void SetIntPoints();

    
    /** encodes the orders of the cartesian product 
     * @param order3 is ignored if < 1 or Dim_ != 3 
     * @return zyx where z = x3 in [1..9], ... */
    int EncodeCartesianOrder(int order1, int order2, int order3);

    /** private order decoder
     * @param order3 0 is written here if Dim_ != 3 */
    void DecodeCartesianOrder(int encoded_order, int* order1, int* order2, int* order3);


    // === Some helper variables to improve performance ===
    Matrix<Double> LDeriv, JInv, J;
  };
 
}

#endif // FILE_BASEFE_2003
