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
  CouplingMemento(const CouplingMemento & x)
  { Error("CouplingMemento: Copy constructor not implemented");}

  //! Destructor
  ~CouplingMemento();

  //! Clear all data
  void Clear();
  
  //! Query if information is saved
  Boolean IsSet() {return isSet_;};
protected:

  //! TRUE, if information is saved
  Boolean isSet_;
  
  //! vector containing types of coupling input
  StdVector<CouplingInputType> inputTypes_;        

  //! vector conatining quantities of coupling input
  StdVector<SolutionType> inputQuantities_;         

  //! vector containing pointer to coupling interfaces  
  StdVector<PDECoupling::CouplingInterface> inputInterfaces_;  
  
};


} // end of namespace

#endif
