// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PYRAFE
#define FILE_PYRAFE

#include <Elements/basefe.hh>

namespace CoupledField
{

  //! Class with general description of pyramidal element

  //! This class is derived from BaseFE. It stores general procedures for each
  //! type of finite element based on <em>pyramids</em>, such as information
  //! on integration points and integration weights
  class PyraFE : public BaseFE {

  public:

    //! Constructor with type of integration rule
    PyraFE();

    //! Destructor
    virtual ~PyraFE(); 

    //! return FE-Type for OLAS
    virtual Elem::FEType feType() {
      return Elem::UNDEF;
    }

    //! Returns wether a given local coordinate is
    //! within this element
    //! \param localCoords (input) local coordinates
    //! \param coordsInside (output) which local coordinates are inside
    virtual void CoordsInsideElem(const Matrix<Double> & localCoords,
                                  const Double tolerance,
                                  StdVector<bool> & coordsInside) const;
    
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

    virtual void FillIntegrationPoints();
        
  };

}

#endif
