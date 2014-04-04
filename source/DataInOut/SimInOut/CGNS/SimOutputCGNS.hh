// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_OUT_CGNS_HH
#define FILE_CFS_SIM_OUT_CGNS_HH

#include <Domain/Mesh/Grid.hh>
#include <Domain/Results/ResultInfo.hh>
#include <DataInOut/SimOutput.hh>

#include <cgnslib.h>

#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

namespace CoupledField
{

  //! Class for writing output-information(grid,results) in unverg-format file
  
  class SimOutputCGNS : virtual public SimOutput {
    
  public:
    
    //! Constructor
    SimOutputCGNS( std::string& fileName,
                   PtrParamNode outputNode,
                   PtrParamNode infoNode,
                   bool isRestart=false );
    
    //! Destructor
    virtual ~SimOutputCGNS();
  
    //! Initialize class
    void Init( Grid* grid, bool printGridOnly );

    //! Write grid definition in file
    void WriteGrid();

    void BeginMultiSequenceStep( UInt step,
                                 BasePDE::AnalysisType type,
                                 UInt numSteps  );

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
    
    //! Finish multisequence step
    void FinishMultiSequenceStep( );

  private:

    void WriteNodesAndElements();
    
    void WriteMixedSection(const StdVector<Elem*>& elems);
    void WritePureSection(const StdVector<Elem*>& elems);

    void TranslateConnectivity(Elem::FEType feType,
                               cgsize_t* cgnsConn,
                               Elem* elem);

    void InitElemTypeMap();
    
    std::map<Elem::FEType,CGNSLIB_H::ElementType_t> elemTypeMap_;

    int indexFile_, indexBase_;
    StdVector<int> idxZone_;
    StdVector<int> idxNodeSol_;
    StdVector<int> idxElemSol_;
    int cellDim_;
    bool outputFileOK_;

    std::map< RegionIdType, std::map<UInt, UInt> > regionNodeMap_;

    //! Map with result objects for each result type
    ResultMapType resultMap_;
    
    //! Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
            
    //! Offset for step value in case of multisequence analysis
    Double stepValOffset_;

    bool writeQuadElems_;

    typedef std::map<std::string, std::map< RegionIdType, Vector<Double> > > RegionSolsType;
    

    //! for printing nodal results of simulation (static/transient)
    /*!
      \param definedOnNode is data defined on nodes?
      \param title title of the results.
      \param x array with nodal results
      \param step number of the step of the calculation
      \param time time of the calculation
    */
    void NodeElemDataTransient(const bool definedOnNode,
                               RegionSolsType& regionSols,
                               const UInt step, 
                               const Double time);
  
    void FillRegionSols(RegionSolsType& regionSols, 
                        const StdVector<shared_ptr<BaseResult> > & solList,
                        ResultInfo::EntityUnknownType entityType );

  };

} // end of namespace
 
#endif
