// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_OUTRST_HH
#define FILE_CFS_SIM_OUTRST_HH

#include <DataInOut/simOutput.hh>

namespace CoupledField
{
  class DynamicLibrary;
  class AnsysBinlibIface;
  
  //! Class for writing results in ANSYS RST postprocessing format
  class SimOutputRST: virtual public SimOutput {
    
  public:


  public:

    //! Constructor
    SimOutputRST( const std::string& fileName, PtrParamNode outputNode );
  
    //! Destructor
    virtual ~SimOutputRST();
  
    //! Initialize class
    void Init( Grid* grid, bool printGridOnly );

    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step,
                                         BasePDE::AnalysisType type,
                                         UInt numSteps);

    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd,
                         bool isHistory );
    
    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    void FinishStep( );

    //! End multisequence step
    virtual void FinishMultiSequenceStep( );

    //! Deallocate last resources
    void Finalize();

  private:
    void GetAnsysRuntimeInfos();
    void RemoveFile(const std::string& fileName,
                    const std::string& exception);
    
    UInt msStep_;
    Double ansysBinlibRev_;
    Double compiledAnsysRev_;
    std::string ansysMachineId_;

    PtrParamNode outputNode_;
    std::string sysPathSep_;

    DynamicLibrary* dynLibrary_;
    AnsysBinlibIface* binlibIface_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SimOutputRST
  //! 
  //! \purpose 
  //! This class provides the interface for writing meshes and elements to
  //! ANSYS postprocessing files .rst. All element types are mapped to four
  //! different ANSYS structural analysis element types PLANE42, SOLID45,
  //! PLANE82 and SOLID95. These are linear/quadratic quads/hexas.
  //! All other element types become degenerated versions of the above types.
  //! At the moment this format is capable of visualizing volume elements, and
  //! surface elements. 
  //! The node solutions get mapped to ANSYS Dofs. There are 32 such Dofs
  //! available (see AnsysNodalDof) and each of them gets written for every
  //! node in the grid in each time step.
  //! This is a tremendous waste of memory and should be redone.
  //! Element results get mapped to one of the 25 ANSYS element result types
  //! (see. AnsysElemDof)
  //! /Note In the moment only transient results will be written. Trying
  //! to write harmonic results will cause an EXCEPTION!
  //! For details about ANSYS binary files see Chapter 1: Format of Binary data
  //! files in the Guide to Interfacing with ANSYS (Programmer's Manual for ANSYS)
  //! The routines used to actually write the files are described in Chapter 2:
  //! Accessing Binary data files. The contents of the files can either be dumped
  //! with the bintst utility compiled with CFS/NACS or in ANSYS with
  //! File->List->Binary Files...
  //! The source code of the unv2rst utility was a good source of info, too.
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

#endif
