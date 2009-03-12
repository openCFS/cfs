// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PDEMEMENTO
#define FILE_PDEMEMENTO

#include "General/environment.hh"

// include serialization headers
#include "CoupledPDE/couplingmemento.hh"
#include "Utils/basenodestoresol.hh"
#include "Utils/boost-serialization.hh"

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

    //! Constructor
    PDEMemento();

    //! Copy constructor
    PDEMemento( const PDEMemento &x ) {
      EXCEPTION( "PDEMemento: Copy constructor not implemented" );
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
    BasePDE::AnalysisType analysisType_;

    //! Name of the related grid 
    std::string gridFileName_;

    //! Step number within current analysistype
    UInt stepNum_;

    //! Frequency of solution (harmonic case)
    Double freq_;

    //! Stores the nodal solution for each region
    std::map<std::string, SingleVector*> solution_;
  
    //! Contains first derivative of PDE solution for each region
    std::map<std::string, Vector<Double> > solDeriv1_;
  
    //! Contains second derivative of PDE solution for each region
    std::map<std::string, Vector<Double> > solDeriv2_;

    //! Contains PDE solution of last step for each region
    std::map<std::string, Vector<Double> > sol_tn_1_;

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

  template <class Archive>
  void PDEMemento::serialize(Archive & ar, const unsigned int version) {

#if 0
    ar & BOOST_SERIALIZATION_NVP(isSet_);
    ar & BOOST_SERIALIZATION_NVP(analysisType_); // Problems with HDF5
    ar & BOOST_SERIALIZATION_NVP(gridFileName_);
    ar & BOOST_SERIALIZATION_NVP(stepNum_);
    ar & BOOST_SERIALIZATION_NVP(freq_);
    ar & BOOST_SERIALIZATION_NVP(solution_); // Problems with HDF5
    ar & BOOST_SERIALIZATION_NVP(solDeriv1_); // Problems with HDF5
    ar & BOOST_SERIALIZATION_NVP(solDeriv2_); // Problems with HDF5
    ar & BOOST_SERIALIZATION_NVP(sol_tn_1_); // Problems with HDF5
    ar & BOOST_SERIALIZATION_NVP(isIterCoupled_);
#endif

    ar & isSet_;
    ar & analysisType_; // Problems with HDF5
    ar & gridFileName_;
    ar & stepNum_;
    ar & freq_;
    ar & solution_; // Problems with HDF5
    ar & solDeriv1_; // Problems with HDF5
    ar & solDeriv2_; // Problems with HDF5
    ar & sol_tn_1_; // Problems with HDF5
    ar & isIterCoupled_;

  }

} // end of namespace

#endif
