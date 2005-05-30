#ifndef FILE_HEXAFE_2003
#define FILE_HEXAFE_2003

#include <Elements/basefe.hh>

namespace CoupledField
{

  //! Class with general description of hexahedral element

  //! This class is derived from BaseFE. It stores general procedures for each
  //! type of finite element based on <em>hexahedra</em>, such as information
  //! on integration points and integration weights
  class HexaFE : public BaseFE {

  public:

    //! Constructor with type of integration rule
    HexaFE();

    //! Destructor
    virtual ~HexaFE(); 

    //! return FE-Type
    virtual FEType feType() {
      return HEX;
    }

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

}

#endif
