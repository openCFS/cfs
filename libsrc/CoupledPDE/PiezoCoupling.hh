#ifndef FILE_PIEZOCOUPLING_HH
#define FILE_PIEZOCOUPLING_HH

#include "BasePairCoupling.hh"


namespace CoupledField
{

  // Forward declarations
  class BaseForm;
  class MaterialData;

  //! Implements the definition of the pairwise piezo coupling
  //! between mechanic and electrostatic field
  class PiezoCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2 );

    //! Destructor
    virtual ~PiezoCoupling();

    //! Trigger calculation of postprocessing results
    void PostProcess();
    
  protected:
    
    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    void ReadStoreResults();
    
    //! Definition of the (bi)linear forms
    void DefineIntegrators();
    
    //! Get correct stiffness integrator
    BaseForm * GetStiffIntegrator( MaterialData * actSDMat,
                                   Boolean reducedInt = FALSE ,
                                   Boolean isdamping = FALSE );

  private:

  };

} // end of namespace

#endif
