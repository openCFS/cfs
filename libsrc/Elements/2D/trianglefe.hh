#ifndef FILE_TRIANGLEFE_2003
#define FILE_TRIANGLEFE_2003

#include <Elements/basefe.hh>

namespace CoupledField
{
  //! Class with general procedures for triangle finite elements
  /*! This class is derived from BaseFE. It stores general procedures for each
    type of finite element on quadralateral, such as calculation of Jacobian of
    transformation in standart and information about integration points and
    integration weights */
  class TriangleFE : public BaseFE {

  public:

    //! Constructor with type of integration rule
    TriangleFE(); 
  
    //! Deconstructor
    virtual ~TriangleFE();

    //! return FE-Type
    virtual FEType feType() { return TRIA;}

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
  protected:

    //! Set integration points
    virtual void SetIntPoints();


  
  };

} // end of namespace

#endif // FILE_RECTANGLEFE_2003
