#ifndef FILE_WEDGEFE
#define FILE_WEDGEFE

#include <Elements/basefe.hh>

namespace CoupledField
{

  //! Class with general description of six-node wedge element

  //! This class is derived from BaseFE. It stores general procedures for each
  //! type of finite element based on <em>wedges</em>, such as information
  //! on integration points and integration weights
  class WedgeFE : public BaseFE {

  public:

    //! Constructor with type of integration rule
    WedgeFE();

    //! Destructor
    virtual ~WedgeFE(); 

    //! return FE-Type
    virtual FEType feType() {
      return ET_UNDEF;
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
    void GetLocalIntPoints4Surface(const StdVector<Integer> & surfConnect,
                                   const StdVector<Integer> & volConnect,
                                   const Vector<Double> & surfIntPoint,
                                   Vector<Double> & volIntPoint);

    /** Sets the default numerical integration - can be overwritten in XML with integRules */ 
    void SetDefaultIntegration()
    {
        // no alternative :)
        IntegMethod = CLASSICAL;
        IntegOrder  = 2;
    }

    /** Sets the default reduced integration  */ 
    void SetDefaultReducedIntegration()
    {
        // no alternative :)
        IntegMethod = CLASSICAL;
        IntegOrder  = 2;
    }


  protected:

    /** the childs fill here the integration points map via AddIntegrationPoints() */    
    virtual void FillIntegrationPoints();

  private:
     
     /** Creates higher order elements and prints them. This method is only used manually
      * if higher integration is required -> Probable it will never be called again, so kill it if
      * it is in the way - Fabian */
     void CreateHigherOrderElements();
     

  };

}

#endif
