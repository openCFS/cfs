#ifndef FILE_PDEMEMENTO
#define FILE_PDEMEMENTO

#include "General/environment.hh"
#include "CoupledPDE/couplingmemento.hh"
#include "Utils/basenodestoresol.hh"

namespace CoupledField
{

// forward class declaration
class BasePDE;

//! Class for saving the internal state of a PDE
class PDEMemento
{
public:
  // Friend declarations
  friend class BasePDE;

  //! Constructor
  PDEMemento();

  //! Copy constructor
  PDEMemento(const PDEMemento & x)
  {
    Error("PDEMemento: Copy constructor not implemented");
  };

  //! Destructor
  ~PDEMemento();

  //! Clear all data
  void Clear();
  
  //! Query if information is saved
  Boolean IsSet() {return isSet_;};

protected:

  //! TRUE, if information is saved
  Boolean isSet_;

  //! Contains analysistype of PDE
  AnalysisType analysisType_;

  //! Contains solution of PDE
  CFSVector * sol_;
  
  //! Contains first derivative of PDE solution
  Vector<Double> solDeriv1_;
  
  //! Contains second derivative of PDE solution
  Vector<Double> solDeriv2_;

  //! Contains the state of the coupling object
  CouplingMemento couplingMemento_;

};


} // end of namespace

#endif
