#ifndef FILE_BASEFORM_2001
#define FILE_BASEFORM_2001

#include <Elements/baseelem.hh>

namespace CoupledField
{

typedef Double (*RHS)(Double,Double,Double);

//! Base class for calculation of elements matrixes
template<Integer dim>
class BaseForm 
{
public:
  //! Constructor
  BaseForm(BaseElem * aptelem);

  //! Deconstructor
  virtual ~BaseForm();

  //! Virtual function, implemented in derived classes
  virtual void CalcElemMatrix(Point<dim> * ptCoord, Matrix<Double> & StiffMat);

  //! Virtual function, implemented in derived classes
  virtual void CalcElemVector(Point<dim> * ptCoord, Vector<Double> & StiffMat);

  //! Virtual function, implemented in derived classes
  virtual void CalcElemVector4InterpolatedFnc(Point<dim> * ptCoord, const Integer aComponent, Vector<Double> & aValueAtNodePoints, Vector<Double> & Result);

  //! Calculate element matrix for multiplication of derivatives respect to x of shape functions
  void CalcShapeFncDxShapeFncDx(Point<dim> * ptCoord, Matrix<Double> & Result);

  //! Calculate element matrix for multiplication of derivatives respect to y of shape functions
  void CalcShapeFncDyShapeFncDy(Point<dim> * ptCoord, Matrix<Double> & Result);

  //! Calculate element matrix for multiplication of shape functions
  void CalcShapeFncShapeFnc(Point<dim> * ptCoord, Matrix<Double> & Result);

  //! Calculation of vector of RHS Dipole having a FlowSrc term (Lighthill's Tensor)
  /*! 
    Calculation of term in the following form:
    \f$\int_{\Gamma}\frac{\partial P}{\partial \bar n} N_i d\Gamma,\f$
    where \f$P\f$ is fluid pressure fluctuation.
    \param ptCoord pointer to array with coordinates
    \param connecth connection array for the element
    \param Result (output) returned element vector
    \param gradN_x_P \f$\frac{\partial P}{\partial x}\f$ at the center of element
  */
  virtual void CalcElemVector4FlowSrcDip(Point<dim> * ptCoord,const Vector<Integer> & connecth, Vector<Double> & Result, const std::vector<Double> gradN_x_P)
  { Error(" This method is not implemented for this class",__FILE__,__LINE__);}

  //!Calculation of vector of RHS Quadrupole having a FlowSrc term (Lighthill's Tensor)
  /*!
    Calculation of the following term:
    \f $ \int_{\Omega} \nabla N_k \nabla \cdot T_{ij} d\Omega$,
    where \f $T_{ij}$ is Lighthill's Tensor, \f $N_k$ is a shape function 

    \param ptCoord pointer to array with coordinates
    \param connecth connection array for the element
    \param Result (output) returned element vector
    \param FlowData (input) matrix with $(P,U_x,U_y)$ data at the nodes of the element
  */
  virtual void CalcElemVector4FlowSrcQuad(Point<dim> * ptCoord,const Vector<Integer> & connecth,const Matrix<Double> & FlowData, Vector<Double> & Result)
  { Error(" This method is not implemented for this class",__FILE__,__LINE__);}

  //!
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected:

  //! 
  BaseElem  * ptelem;

  //! Calculation of function in integration points
  Double FuncAtIP(const ShortInt iIP, RHS f);  

};

#ifdef __GNUC__
template class BaseForm<3>;
template class BaseForm<2>;
#endif

} //end namespace

#endif // FILE_BASERFORM
