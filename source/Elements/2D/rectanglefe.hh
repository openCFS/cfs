// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_RECTANGLEFE_2003
#define FILE_RECTANGLEFE_2003

#include <Elements/basefe.hh>

namespace CoupledField
{

  //! Class with general procedures for quadrilateral finite elements
  /*! This class is derived from BaseFE. It stores general procedures for each
    type of finite element on quadralateral, such as calculation of Jacobian of
    transformation in standart and information about integration points and
    integration weights */
  class RectangleFE : public BaseFE {

  public:

    //! Constructor with type of integration rule
    RectangleFE(); 
  
    //! Deconstructor
    virtual ~RectangleFE();

    //! return FE-Type for CLA++
    virtual FEType feType() { return ET_UNDEF;}

    //! Calculates a measure for the geometric distortion of an element
    /*!
      \param CornerCoords (input) Corner coordinates of the element
      \param Displacement (input) Displacement of the corner points (same ordering as CornerCoords!!)
    */
    virtual Double CalcDistortion(Matrix<Double> &cornerCoords, Vector<Double> &size, Matrix<Double> &displacements);

    //! Calculates corresponding volume point of neighbouring surface

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
    void GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                   const StdVector<UInt> & volConnect,
                                   const Vector<Double> & surfIntPoint,
                                   Vector<Double> & volIntPoint);
     
    //! compute length of edge with maximal size
    void GetMaxMinEdgeLength( Matrix<Double> &ptCoord, Double &Lmax, Double &Lmin );
    
    /** Set cartesian integration points */                               
    void SetCartesianIntegration(int order1, int order2)
    {
        SetCartesianInteg(order1, order2, 0, false);
    }

    //! Returns wether a given local coordinate is
    //! within this element
    //! \param localCoords (input) local coordinates
    //! \param coordsInside (output) which local coordinates are inside
    virtual void CoordsInsideElem(const Matrix<Double> & localCoords,
                                  const Double tolerance,
                                  StdVector<bool> & coordsInside) const;
                                   
   protected:
      /** the childs fill here the integration points map via AddIntegrationPoints() */    
    virtual void FillIntegrationPoints();
 
  };

} // end of namespace

#endif // FILE_RECTANGLEFE_2003
