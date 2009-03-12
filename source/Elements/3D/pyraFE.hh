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

  protected:

    virtual void FillIntegrationPoints();
    
    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        IntegMethod = CLASSICAL;
        IntegOrder  = 2;
    }    

    /** Sets the default reduced integration But we don't have more elements! */ 
    void SetDefaultReducedIntegration()
    {
        IntegMethod = CLASSICAL;
        IntegOrder  = 2;
    }    
    
  };

}

#endif
