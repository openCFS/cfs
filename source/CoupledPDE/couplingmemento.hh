// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_COUPLINGMEMENTO
#define FILE_COUPLINGMEMENTO

#include "General/environment.hh"
#include "pdecoupling.hh"

namespace CoupledField
{

  //! Class for saving the internal state of a PDECoupling object
  class CouplingMemento
  {
  public:

    // friend declarations
    friend class PDECoupling;

    //! Constructor
    CouplingMemento();

    //! Copy Constructor
    CouplingMemento( const CouplingMemento & x ) {
      EXCEPTION( "CouplingMemento: Copy constructor not implemented" );
    }

    //! Destructor
    ~CouplingMemento();

    //! Clear all data
    void Clear();
  
    //! Query if information is saved
    bool IsSet() {return isSet_;};
  protected:

    //! true, if information is saved
    bool isSet_;
  
    //! vector containing types of coupling input
    StdVector<CouplingInputType> inputTypes_;        

    //! vector conatining quantities of coupling input
    StdVector<SolutionType> inputQuantities_;         

    //! vector containing pointer to coupling interfaces  
    StdVector<PDECoupling::CouplingInterface> inputInterfaces_;  
  
    // =======================================================================
    // SERIALIZATION FUNCTIONS
    // =======================================================================
    // These functions allow us to write a memento directly
    // into an boost::archive, for saving on a disk or in a 
    // iostream object
    
    //! allow serialization class to access memento entries
    friend class boost::serialization::access;
    
    //! Saving internal state into a boost::archive
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class CouplingMemento
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

 template<class Archive>
  void CouplingMemento::serialize(Archive & ar, 
                                  const unsigned int version) {
   EXCEPTION( "Not imeplemented at the momement" );
    
    // For further details, why this is currently not imeplemented,
    // see PDECoupling::CouplingInterface::serialize
  }

} // end of namespace

#endif
