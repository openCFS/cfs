// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ACOU_MECH_COUPLING_HH
#define FILE_ACOU_MECH_COUPLING_HH

#include "BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"

namespace CoupledField {
class SinglePDE;
}  // namespace CoupledField


namespace CoupledField
{

  // Forward declarations

  //! Implements the definition of the pairwise acoustic-mechanic
  class AcouMechCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    //! \param paramNode pointer to "couplinglist/direct/acouMechDirect" element
    AcouMechCoupling( SinglePDE *pde1, SinglePDE *pde2, 
                      PtrParamNode paramNode  );

    //! Destructor
    virtual ~AcouMechCoupling();

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
    
    //! Method for obtaingin material data

    //! This method is overwritten here, since the material of the 
    //! acoustic domain is already written in.
    void ReadMaterialData();

    //! initialize nonlinearities
    virtual void InitNonLin(){;};
 
  private:

  };


#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class AcouMechCoupling
  //! 
  //! \purpose This class implements the direct coupling of mechanical and
  //! acoustic field via surface elements.
  //! 
  //! \collab Objects of this class are instantiated by the class Domain.
  //! Afterwards they are passed to an instance of DirectCoupledPDE, which 
  //! can have arbitrarely many of them.
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused
  //! 
  //! \improve 
  //! 

#endif

} // end of namespace

#endif
