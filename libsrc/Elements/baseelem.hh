#ifndef FILE_BASEELEM_2001
#define FILE_BASEELEM_2001

namespace CoupledField
{

  //! Class for calculation Jacobian 
  /*! In this class we store Jacobian of transformation element in standart, det of Jacobian and calculate inverse of Jacobian matrix */

template<class Dim>
class Jacobian
{
  public:
    Double          detJ;          //!< determinant of Jacobian
    Matrix<Double>  J;            //!< Jacobian
    Matrix<Double>  Jinv;         //!< inverse Jacobian

    //! Constructor 
    Jacobian();
    
    //! Return colomn with derivatives respect to x in inverse Jacobian
    void GetJinvX(Vector<Double> & JinvX);

    //! Return colomn with derivatives respect to y in inverse Jacobian
    void GetJinvY(Vector<Double> & JinvY);
};

template class Jacobian<Point2D>;
template class Jacobian<Point3D>;

//-----------------------------------------------------------------------------

  //! Base class for description of elements 
  /*! From this class all special finite elements as triangles, quadrilaterals,
    tetraheader, etc. are derived. All memeber functions of this class are pure virtual,
    just the constructor is not virtual.
  */

class BaseElem
{
public:

  //! constructor (does nothing)
  BaseElem();
 
  //! deconstructor 
  virtual ~BaseElem();
 
  //! 
  virtual void Init()=0;

  //!
  virtual void SetIntPoints()=0;

  //!
//  virtual void SetTransFncAtIntPoints()=0;

  //!
//  virtual void SetDerTransFncAtIntPoints()=0;

  //!
  virtual  void GetGradientShFnc(Vector<Double> & ,const Integer i, const Integer j)=0;

   //! Calculation of Jacobian in 2D
  virtual void CalcJacobian(Jacobian<Point2D> & J, const Integer ip,
         const Point2D * const ptCoord, const Boolean NeedJinv=TRUE)=0;

  //! Calculation of Jacobian in 3D
 virtual void CalcJacobian(Jacobian<Point3D> & J, const Integer ip,
              const Point3D * const ptCoord, const Boolean NeedJinv=TRUE)=0;

  //!
  virtual Vector<Double> &  GetShFncAtIP(const Integer iShFnc)=0 ;

  //! Return number of nodes   
  ShortInt GetDim() { return Dim; };
 
  //! Return number of integration pointes
  ShortInt GetNumNodes() {return NumNodes; }; 

  //! Returns number of integration points
  ShortInt GetNumIntPoints() {return NumIntPoints; }; 

  //! Return pointer to coordinate x of integration pointes
  Double GetIntPointsX(Integer i){return IntPoints[i][0];}

  //! Return pointer to coordinate y of integration pointes
  Double GetIntPointsY(Integer i){return IntPoints[i][1];}

  //! Return pointer to integration weights
  Vector<Double> * GetIntWeights(){ return IntWeights;}

protected:

  ShortInt Dim;                //!< space dimension
  ShortInt NumNodes;           //!< number of nodes
  ShortInt NumEdges;           //!< number of edges 
  ShortInt NumFaces;           //!< number of faces
  ShortInt NumIntPoints;       //!< number of integration points
  ShortInt DegreeInteg;        //!< numerical integration order
  Matrix<Double> IntPoints;  //!< integration points
  Vector<Double> * IntWeights; //!< integration weights

private:
 
};
 
}

#endif // FILE_BASEELEM_2001
