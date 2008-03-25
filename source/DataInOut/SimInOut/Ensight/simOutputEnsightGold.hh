// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIMOUTPUTENSIGHT_HH
#define FILE_CFS_SIMOUTPUTENSIGHT_HH

#include <set>

#include <Domain/grid.hh>
#include <Domain/resultInfo.hh>
#include <DataInOut/simOutput.hh>

namespace CoupledField {

  //! Class for reading in a mesh created by the ANSYS mkmesh-extension.

  //! Class, that is derived from class FileType for reading mesh-input data,
  //! which is produced by Ansys mkmesh-interface. 
  class SimOutputEnsightGold: virtual public SimOutput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimOutputEnsightGold(std::string fileName, ParamNode * inputNode);
    
    //! Destructor
    virtual ~SimOutputEnsightGold();

    //! Initialize class 
    virtual void Init( Grid* ptGrid,
                       bool printGridOnly );

    //@}

    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step,
                                         BasePDE::AnalysisType type,
                                         UInt numSteps  );
    
    //! Register result (within one multisequence step)
    virtual void RegisterResult( shared_ptr<BaseResult> sol,
                                 UInt saveBegin, UInt saveInc,
                                 UInt saveEnd,
                                 bool isHistory );

    //! Begin single analysis step
    virtual void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    virtual void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    virtual void FinishStep( );

    //! End multisequence step
    virtual void FinishMultiSequenceStep( );

    //! Finalize the output
    virtual void Finalize();

    //! Initialize 
    virtual void InitModule();

    //! Write grid
    virtual void WriteGrid();
    
  private:

    //! This method maps the type number of an element - as given in the 
    //! mesh file - to a pointer to a reference finite element.
    //! \param itype (input) element type number as read in from the mesh
    UInt RegionElems2EnsightElems(const RegionIdType region,
                                  std::map<FEType, std::vector<UInt> > & elemMap,
                                  std::map<FEType, std::vector<UInt> > & elemIds);
    //@}

    //! Stream for case file
    std::ofstream caseFile_;
    
    //! Stream for geo file
    std::ofstream geoFile_;

    //! Flag indicating if grid is already written
    bool gridWritten_;

    //! Flag indicating if only grid is to be printed
    bool printGridOnly_;

    //! Map internal element types to Ensight names
    std::map<FEType, std::string > elemNames_;
    
    BasePDE::AnalysisType analysisType_;

    struct TimeFileSetInfo
    {
      UInt setNumber;
      UInt actStepIdx;
      StdVector<UInt> stepNums;
      StdVector<Double> stepVals;
      std::string baseFileName;
//      std::ofstream* resultFile;
      StdVector< shared_ptr<BaseResult> > regionResults;
      ResultInfo::EntityUnknownType entityType;
      ResultInfo::EntryType entryType;
    };
    
    
    //! Map available result names to Time / File info set structures.
    std::map<std::string, TimeFileSetInfo > resultSetInfo_;
    
    //! Platform specific path separator
    std::string pathSep_;
    
    //! Maps region ids to parts
    std::map<UInt, UInt> regionIdToPart_;
  };

}

#endif // FILE_ANSYSFILE
