#ifndef FILE_BASEELEM_2001
#define FILE_BASEELEM_2001

#include "matrix.hh"
#include "jacobian.hh"

namespace CoupledField
{
  //! Base class for description of elements 
  /*! From this class all special finite elements as triangles, quadrilaterals,
    tetraheader, etc. are derived. All memeber functions of this class are pure virtual,
    just the constructor is not virtual.
  */

class BaseElem
{
public:
  virtual void test(){;} 

  //! constructor (does nothing)
  BaseElem();

  //! deconstructor 
  virtual ~BaseElem();
 
  //! 
  virtual void Init()=0;

  //!
  virtual void SetIntPoints()=0;

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
  enum IntegrationType IntegType;

  enum IntegrationType String2EnumIntegrationType(const Char * inttype);
private:
 
};
 
}

#endif // FILE_BASEELEM_2001
