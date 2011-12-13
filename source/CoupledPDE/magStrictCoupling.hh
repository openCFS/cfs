// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#ifndef FILE_MAGNETOSTRICTION_COUPLING_HH_
#define FILE_MAGNETOSTRICTION_COUPLING_HH_

#include "BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  // Forward declarations
class BaseResult;
  class SinglePDE;

  //! Class for pairwise magnetostrictive coupling
  
  //! This class implements the direct magnetostrictive coupling 
  //! between the mechanic and magnetic PDE.
  //! \note: Up to now, only the coupling of the magnetic scalar PDE
  //! is defined. The vector-potential case should be integrated as well.
  class MagStrictCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE (Mechanic)
    //! \param pde2 pointer to second coupling PDE (Magnetic)
    //! \param paramNode pointer to "couplinglist/direct/magnetostriction" element(xml file)
    MagStrictCoupling( SinglePDE *pde1, SinglePDE *pde2, 
                       PtrParamNode paramNode );

    //! Destructor
    virtual ~MagStrictCoupling();

    //! Trigger calculation of postprocessing results
    void CalcResults( shared_ptr<BaseResult> result );


  protected:

    //! Definition of the (bi)linear forms
    void DefineIntegrators();

    //! Define available results
    void DefineAvailResults();

  private:

    // *****************
    //  POSTPROCESSING
    // *****************
    
    //! Calculate B Field
    template <class TYPE>
    void CalcBField( shared_ptr<BaseResult> result );    

    //! Calculate coupled stress 
    template <class TYPE>
    void CalcStress( shared_ptr<BaseResult> result );
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class
  // =========================================================================

  //! \class MagStrictCoupling
  //!
  //! \purpose
  //!
  //! \collab
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


#endif /* FILE_MAGNETOSTRICTION_COUPLING_HH_ */
