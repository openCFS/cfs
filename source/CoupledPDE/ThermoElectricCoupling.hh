#ifndef FILE_THERMO_ELECTRIC_COUPLING_HH
#define FILE_THERMO_ELECTRIC_COUPLING_HH

#include <string>

#include "BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"


namespace CoupledField
{

  // Forward declarations
  class SinglePDE;

  //! Implements the definition of the pairwise thermo-electric coupling
  class ThermoElectricCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    //! \param paramNode pointer to "couplinglist/direct/thermoMechDirect" element
    ThermoElectricCoupling( SinglePDE *pde1, SinglePDE *pde2, 
                      PtrParamNode paramNode  );

    //! Destructor
    virtual ~ThermoElectricCoupling();

    //! Trigger calculation of postprocessing results
    void CalcResults();

    void WriteResultsInFile(const UInt kstep,
                            const Double asteptime,
                            UInt stepOffset,
                            Double timeOffset);    
    
  protected:
    
    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    void ReadStoreResults();
    
    //! Definition of the (bi)linear forms
    void DefineIntegrators();
    
    //! initialize nonlinearities
    virtual void InitNonLin(){;};

	//! Check non linearities of thermoelectric coupling
	bool nonLinThermoElecCoupling_;

	//! SubType of pyroelectric section
    std::string subType_;
 
  private:


  };


#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class ThermoElectricCoupling
  //! 
  //! \purpose This class implements the direct coupling of mechanical and
  //! thermal field via surface elements?.
  //! 
  //! \collab Objects of this class are instantiated by the class Domain.
  //! Afterwards they are passed to an instance of DirectCoupledPDE, which 
  //! can have arbitrarely many of them.
  //! 
  //! \implement 
  //! 
  //! \status In developing
  //! 
  //! \unused
  //! 
  //! \improve 
  //! 

#endif

} // end of namespace

#endif
