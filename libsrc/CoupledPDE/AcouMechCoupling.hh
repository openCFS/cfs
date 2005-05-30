#ifndef FILE_ACOU_MECH_COUPLING_HH
#define FILE_ACOU_MECH_COUPLING_HH

#include "BasePairCoupling.hh"


namespace CoupledField
{

  // Forward declarations
  class BaseForm;
  class MaterialData;

  //! Implements the definition of the pairwise acoustic-mechanic
  //! coupling via surface elements.
  class AcouMechCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    AcouMechCoupling( SinglePDE *pde1, SinglePDE *pde2 );

    //! Destructor
    virtual ~AcouMechCoupling();

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
    
    //! Method for obtaingin material data

    //! This method is overwritten here, since the material of the 
    //! acoustic domain is already written in.
    void ReadMaterialData();
  private:

  };

} // end of namespace

#endif
