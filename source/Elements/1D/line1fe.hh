// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LINE1FE_2003
#define FILE_LINE1FE_2003

#include <Elements/basefe.hh>
#include <Elements/1D/linefe.hh>
#include <Domain/elem.hh>
#include "MatVec/vector.hh"

namespace CoupledField
{
  //! 1D line element with two nodes (linear interpolation function)

  class Line1FE : public LineFE
  {
  public:

    /** Constructor with type of integration rule for cartesian product rule.
    * Leave blank in any other case
    * @param method leave to default - usage only internally for product rule integratio points
    * @param order  same as for parameter method applies*/
    Line1FE(IntegrationMethod method = UNDEFINED, int order=0);

    //! Destructor
    virtual ~Line1FE();

    //! return FE-Type
    virtual Elem::FEType feType() {
      return Elem::LINE2;
    };

  protected:

    //! Initialize line element
    virtual void Init(IntegrationMethod method, int order);

    //! Set local corner coordinates
    virtual void SetCornerCoords();

    //! Set local edge indices
    void SetEdgeIndices();

    //! calculates the shape functions at an arbitrary local point
    /*!
    \param Shape (output) Vector of shape fnc values \f$ (N_{1},N_{2})^T \f$
    \param LCoord (input) Local coordinates of evalutation point
    */
    virtual void CalcShapeFnc(Vector<Double> & LShape,
                              const Vector<Double> & LCoord,
                              const Elem* elem , UInt dof,
                              AnsatzFct::FctEntityType );

    //! calculates the local derivatives of shape functions at an arbitrary local point
    /*!
    \param LDeriv (output) Matrix with local derivatives of all shape functions
    \f[ \left( \begin{array}{cc} N_{1,d\xi}  \\
    N_{2,d\xi} \end{array}\right) \f]
    \param LCoord (input) Local coordinates of evalutation point
    */
    virtual void CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
                                        const Vector<Double> & LCoord,
                                        const Elem* elem , UInt dof,
                                        AnsatzFct::FctEntityType );

    /** Sets the default numerical integration - can be overwritten in XML with integRules */
    void SetDefaultIntegration()
    {
      IntegMethod = ECONOMICAL;
      IntegOrder  = 3; // NOT confirmed :( - Fabian
    }

    /** Sets the default reduced numerical integration */
    void SetDefaultReducedIntegration()
    {
      IntegMethod = ECONOMICAL;
      IntegOrder  = 1; // NOT confirmed :( - Fabian
    }

    //! Get the local coordinates for given global ones
    //! \param localCoords (output) local coordinates
    //! \param globalCoords (input) global coordinates
    //! \param coordMat (input) global corner coordinates of element
    //!                         (spaceDim \f$\times\f$ nrNodes)
    virtual void Global2LocalCoords(Matrix<Double> & localCoords,
                                    const Matrix<Double> & globalCoords,
                                    const Matrix<Double> & coordMat);

    virtual void GetNumFncs(StdVector<UInt>& numFcns,
                            const shared_ptr<AnsatzFct>& fcnType,
                            AnsatzFct::FctEntityType fctEntityType,
                            UInt dof = 1);


    virtual UInt GetNumFncs( const shared_ptr<AnsatzFct>& fncType );

    virtual void SetAnsatzFct( shared_ptr<AnsatzFct>& actFct, bool setIntPoints ) ;

  private:
    void CalcSpectralShFct( Vector<Double> & Shape,
                            const Vector<Double> & LCoord,
                            const Elem* elem, UInt dof,
                            AnsatzFct::FctEntityType type );

    void CalcSpectralDerivFct( Matrix<Double> & LDeriv,
                               const Vector<Double> & LCoord,
                               const Elem* elem, UInt dof,
                               AnsatzFct::FctEntityType type);

    //! 1D Lagrange functions at IPs for spectral mode
    Vector<Double> * sShFcnAtIp_;

    //! 1D derivatives of Lagrange fncs at ips for spectral Mode
    Vector<Double> * sDerivAtIp_;

  };

} // end of namespace

#endif // FILE_LINE1FE_2003
