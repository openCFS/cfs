// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_QUAD1FE_2003
#define FILE_QUAD1FE_2003

#include <Elements/basefe.hh>
#include <Elements/2D/rectanglefe.hh>

namespace CoupledField
{
  //! Quadrilateral finite element with four nodes (linear interpolation function)
  
  class Quad1FE : public RectangleFE
  {
  public:
  
    //! Constructor with type of integration rule
    Quad1FE();
  
    //! Destructor
    virtual ~Quad1FE();
  
    //! return FE-Type
    virtual FEType feType() {
      return ET_QUAD4;
    };

  protected:

    //! Initialize Quadrilateral element
    virtual void Init();

    //! Set local corner coordinates
    virtual void SetCornerCoords();

    //! Set local edge indices
    void SetEdgeIndices();

    //! Set local face indices
    void SetFaceIndices();

    //! calculates the shape functions at an arbitrary local point
    /*!
      \param Shape (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcShapeFnc(Vector<Double> & LShape, 
                              const Vector<Double> & LCoord,
                              const Elem* elem,
                              UInt dof,
                              AnsatzFct::FctEntityType );
  
   //! Calculates the shape functions of incompatible modes at an arbitrary local point
    /*!
      \param Shape (output) Vector of shape fnc values \f$ (N_{1},\cdots\,N_{NumNodes})^T \f$
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcShapeFncICModes(Vector<Double> & Shape, 
                                     const Vector<Double> & LCoord,
                                     const Elem* elem , UInt dof,
                                     AnsatzFct::FctEntityType = AnsatzFct::ALL );

    //! calculates the local derivatives of shape functions at an arbitrary local point
    /*!
      \param LDeriv (output) Matrix with local derivatives of all shape functions
      \f[ \left( \begin{array}{ccc} N_{1,d\xi} & N_{1,d\eta} & \cdots \\
      N_{2,d\xi} & N_{2,d\eta} & \cdots \\
      \cdots     & \cdots      & \cdots \end{array}\right) \f]
      \param LCoord (input) Local coordinates of evalutation point 
    */
    virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                        const Vector<Double> & LCoord,
                                        const Elem* elem,
                                        UInt dof,
                                        AnsatzFct::FctEntityType);
  
    virtual void CalcLocal2ndDerivShapeFnc(Matrix<Double> & L2ndDeriv, 
                                           const Vector<Double> & LCoord,
                                           const Elem* elem,
                                           UInt dof,
                                           AnsatzFct::FctEntityType);

    
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
					       AnsatzFct::FctEntityType = AnsatzFct::ALL );

    //! Calculates a measure for the geometric distortion of an element
    /*!
      \param cornerCoords (input) Corner coordinates of the element
      \param size (input) Absolute size of element in all dimensions
      \param displacement (input) Displacement of the corner points (same ordering as CornerCoords!!)
    */
    virtual Double CalcMeanStrain(Matrix<Double> &cornerCoords, 
                                  Matrix<Double> &displacements);

    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 3; // actually 2 but 2+3 is same and we avoid the msg
    }

    /** Sets the default reduced integration */ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = ECONOMICAL;
        IntegOrder  = 1; 
    }

    void SetAnsatzFct( shared_ptr<AnsatzFct>& actFct,
                       bool setIntPoints = true);

    void GetNumFncs(Vector<UInt>& numFcns, 
                    const shared_ptr<AnsatzFct>& fcnType, 
                    AnsatzFct::FctEntityType fctEntityType, 
                    UInt dof = 1);


    UInt GetNumFncs( const shared_ptr<AnsatzFct>& fncType );

  private:
    virtual void CalcSpectralShFct( Vector<Double> & Shape, 
                                    const Vector<Double> & LCoord,
                                    const Elem* elem, UInt dof,
                                    AnsatzFct::FctEntityType type );

    virtual void CalcSpectralDerivFct( Matrix<Double> & LDeriv, 
                                       const Vector<Double> & LCoord,
                                       const Elem* elem, UInt dof,
                                       AnsatzFct::FctEntityType type);

    //! 1D Lagrange functions at IPs for spectral mode
    Vector<Double> * sShFcnAtIp_;

    //! 1D derivatives of Lagrange fncs at ips for spectral Mode
    Vector<Double> * sDerivAtIp_;

    //! 1D Legendre functions at IPs
    Vector<Double> * lShFcnAtIp_;
    
    //! 1D derivatives of legendre functions at IP
    Matrix<Double> * lDerivAtIp_;
  };

} // end of namespace

#endif // FILE_QUAD1FE_2003
