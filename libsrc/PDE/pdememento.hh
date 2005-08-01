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
    friend class StdPDE;

    //! Constructor
    PDEMemento();

    //! Copy constructor
    PDEMemento( const PDEMemento &x ) {
      Error( "PDEMemento: Copy constructor not implemented",
             __FILE__, __LINE__ );
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
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class PDEMemento
  //! 
  //! \purpose 
  //! This class encapsulates the internal state of an SinglePDE, i.e. the
  //! solution vector and (in transient case) the time derivative(s). This is
  //! used to transfer the internal state of a PDE from one MultiSequence-step
  //! to the next.
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! This class implements the famous memento concept (ref. Gamma: Design 
  //! patterns)
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! For a generalization this class should also be hierarchical, in order to
  //! also perform multisequence analysis with a DirectCooupledPDE.

#endif

} // end of namespace

#endif
