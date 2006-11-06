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
    friend class SinglePDE;

    //! Type of usage of stored values: usage as start values or
    //! dirichlet values
    typedef enum { NO_TYPE, START_VALUE, DIRICHLET_VALUE } ValueUsageType;

    //@{
    //! Conversion of ValueUsageType
    static std::string Enum2String( ValueUsageType type );
    static ValueUsageType String2Enum( std::string typeString );
    //@}
    
    
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
    bool IsSet() {return isSet_;};

    //! Get restart step
    UInt GetRestartStep() { return stepNum_; }

  protected:

    //! true, if information is saved
    bool isSet_;
    
    //! Contains analysistype of PDE
    AnalysisType analysisType_;

    //! Name of the related grid 
    std::string gridFileName_;

    //! Step number within current analysistype
    UInt stepNum_;

    //! Frequency of solution (harmonic case)
    Double freq_;

    //! Stores the nodal solution for each region
    std::map<std::string,CFSVector*> solution_;
  
    //! Contains first derivative of PDE solution for each region
    std::map<std::string, Vector<Double> > solDeriv1_;
  
    //! Contains second derivative of PDE solution for each region
    std::map<std::string, Vector<Double> > solDeriv2_;

    //! Flag indicating iterative coupling
    bool isIterCoupled_;

    //! Contains the state of the coupling object
    CouplingMemento couplingMemento_;

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
