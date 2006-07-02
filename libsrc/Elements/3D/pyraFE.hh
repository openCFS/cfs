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
    virtual FEType feType() {
      return PYR;
    }

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
