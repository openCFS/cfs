#ifndef FILE_BASEELEM_2001
#define FILE_BASEELEM_2001

#include <Matrix/matrix.hh>
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

  //! constructor (does nothing)
  BaseElem();

  //! deconstructor 
  virtual ~BaseElem();
 
  //! 
  virtual void Init()=0;

  //!
  virtual void SetIntPoints()=0;



  /// get matrix of shape functions in all integration points
  /*! 
    @param shapeMat: nrIntPts \f$\times \f$ nrSpaceDim
  */
  virtual void GetShapeMatIP(Matrix<Double> shapeMat) //=0;
  { Error("Not implemented ",__FILE__,__LINE__);}



  /// get matrix of shape functions derived after local coordinates in all integration points
  /*! 
    @param derivShapeMat: nrIntPts \f$\times \f$ nrSpaceDim
  */
  virtual void GetLocDerivShapeMatIP(Matrix<Double> derivShapeMat) //=0;
  { Error("Not implemented ",__FILE__,__LINE__);}



  /// get matrix of shape functions derived after global coordinates in all integration points
  /*! 
    @param derivShapeMat: nrIntPts \f$\times \f$ nrSpaceDim
  */
  virtual void GetGlobDerivShapeMatIP(Matrix<Double> derivShapeMat) //=0;
  { Error("Not implemented ",__FILE__,__LINE__);}



  //! get gradient of shape fnc for node i at integration point j
  virtual  void GetGradientShFnc(Vector<Double> & ,const Integer i, const Integer j)=0;

  //! get gradient of shape fnc for node i at center of element
  virtual  void GetGradientShFncAtCenter(Vector<Double> & ,const Integer i)
  { Error("Not implemented ",__FILE__,__LINE__);}

   //! Calculation of Jacobian in 2D
  virtual void CalcJacobian(Jacobian<2> & J, const Integer ip,
          Point<2> *  ptCoord, const Boolean NeedJinv=TRUE)=0;

  //! Calculation of Jacobian in 3D
 virtual void CalcJacobian(Jacobian<3> & J, const Integer ip,
              Point<3> *  ptCoord, const Boolean NeedJinv=TRUE)=0;

  //! Calculation of Jacobian for center point in 2D
  virtual void CalcJacobianAtCenter(Jacobian<2> & J, Point<2> * ptCoord, const Boolean NeedJinv=TRUE)
  { Error(" This function is not implemented ", __FILE__,__LINE__);}

   //! Calculation of Jacobian for center point in 3D
  virtual void CalcJacobianAtCenter(Jacobian<3> & J, Point<3> *  ptCoord, const Boolean NeedJinv=TRUE)
  { Error(" This function is not implemented ", __FILE__,__LINE__);}

  //! get value of shape fnc number iShFnc at integration points
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

  //! return intergral of shape function over line. only for 1d-line
  virtual Double getIntVal(const Point<2> * const ptCoord)
  { 
    Error(" This function is implemented only for 1D elements",__FILE__,__LINE__);
    return Ddummy;
  }	

///////////////////////////////////////////////////////////////////////
  virtual Vector<Double> *  GetDxShFncAtIP(const Integer iShFnc)
  { 
    Error("Not implemented") ;
    return Dvec;
  }

  virtual Vector<Double> *  GetDyShFncAtIP(const Integer iShFnc) 
  { 
    Error("Not implemented") ;
    return Dvec;
  }

  virtual Vector<Double> *  GetDzShFncAtIP(const Integer iShFnc) 
  { 
    Error("Not implemented") ;
    return Dvec;
  }

  //! return FE-Type for LAS++
  virtual Integer feType()=0;

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

  //! Converts the string used for the integration type to an integer
  enum IntegrationType String2EnumIntegrationType(const Char * inttype);

  //just dummies for SUN compiler
  Double Ddummy;
  Vector<Double> *Dvec;

private:
 
};
 
}

#endif // FILE_BASEELEM_2001
