// =====================================================================================
// 
//       Filename:  SimInputCGNS.hh (originally FileReader_CGNS.hh in cplreader)
// 
//    Description:  File Reader for the cgns file format
//                  Based on OpenFOAM reader
// 
//        Authors:  Simon Triebenbacher, simon.triebenbacher@tuwien.ac.at
//                  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//                  Jens Grabinger, jens.grabinger@simetris.de
//        Company:  TU Wien, Universitaet Klagenfurt, SIMetris GmbH
// 
// =====================================================================================


#ifndef SIMINPUT_CGNS_HH
#define SIMINPUT_CGNS_HH

#include <string>
#include <set>
#include <unordered_map>

#include <cgnslib.h>

#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

#include "DataInOut/SimInput.hh"

namespace CoupledField{

  //! Class for reading unstructured simulation grids and CFD data from CGNS files.

  //! This class is for reading (C)FD (G)eneral (N)otation (S)ystem file(s).
  //!
  //! Among others, CGNS files can be created in the following pre-processors
  //! and CFD codes:
  //! ANSYS Workbench mesher: File -&gt; Export Mesh -&gt; CGNS Format. ANSYS
  //!          &lt;= 15.0 can only write linear element types. Named components
  //!          are written as CGNS sections using all capital letters.
  //! ANSYS CFX: Native CFX result files can be converted to CGNS by using
  //!          cfx5export -cgns input.res
  //! Pointwise: Uses CGNS as its native file format (www.pointwise.com).
  //!
  //! As the name CGNS states, the format is mostly used for exchange of CFD data.
  //! Since most CFD codes are using linear cell types, most pre- (ANSYS WB) and
  //! post-processors have problems handling quadratic element types often used in
  //! FE codes, even if the CGNS standard properly defines such element types.
  //! Nonetheless, this class is able to read all complete and incomplete element
  //! types supported by CGNS and openCFS. This class handles CGNS sections as regions.
  //! Therefore, it is not yet able to read CGNS files written with the SimOutputCGNS
  //! class, since it maps regions to zones.
  //!
  //! Besides supporting unstructured mesh reading, this class also supports reading
  //! the VelocityX, VelocityY, VelocityZ and Pressure fields for a single time step.
  //!
  //! \note ATTENTION: This class has been implemented for the exchange of CFD data
  //!                  during the E+H harmonic FSI project by Simon Triebenbacher.
  //!                  In this project we only need to transfer results for a single
  //!                  RANS step from CFX. Therefore, no support for more time/freq steps
  //!                  has been implemented and only one multisequence step is supported!
  class SimInputCGNS : public SimInput {
  public:
    
    //! Constructor
    SimInputCGNS( std::string fileName, PtrParamNode inputNode,
                  PtrParamNode infoNode );
    
    //! Destructor
    virtual ~SimInputCGNS();
    
    //! Initializes the reader and reads files' contents
    virtual void InitModule();
    
    //! Fill a Grid object with stored mesh
    virtual void ReadMesh(Grid* ptGrid);
    
    //! @copydoc SimInput::GetDim()
    virtual UInt GetDim() { return dim_; }

    //! @copydoc SimInput::GetNumNodes()
    virtual UInt GetNumNodes() { return numNodes_; }

    //! @copydoc SimInput::GetNumElems()
    virtual UInt GetNumElems(Integer);

    //! @copydoc SimInput::GetNumRegions()
    virtual UInt GetNumRegions() { return numRegions_; }

    //! @copydoc SimInput::GetNumNamedNodes()
    virtual UInt GetNumNamedNodes() { return namedNodes_.size(); }

    //! @copydoc SimInput::GetNumNamedElems()
    virtual UInt GetNumNamedElems() { return namedElems_.size(); }

    //! @copydoc SimInput::GetAllRegionNames()
    virtual void GetAllRegionNames(StdVector<std::string> &regionNames) {
      regionNames = regionNames_;
    }

    //! @copydoc SimInput::GetNodeNames()
    virtual void GetNodeNames(StdVector<std::string>& nodeNames) {
      nodeNames = namedNodes_;
    }

    //! @copydoc SimInput::GetElemNames()
    virtual void GetElemNames(StdVector<std::string>& elemNames) {
      elemNames = namedElems_;
    }

    //! @copydoc SimInput::GetNumNodesElemsByName()
    virtual void GetNumNodesElemsByName(const std::string &name,
                                        UInt &numNodes, UInt &numElems);


    // =========================================================================
    //  GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information

    //! @copydoc SimInput::GetNumMultiSequenceSteps()
    void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                   std::map<UInt, UInt>& numSteps,
                                   bool isHistory = false );

    //! @copydoc SimInput::GetResultTypes()
    void GetResultTypes( UInt sequenceStep, 
                         StdVector<shared_ptr<ResultInfo> >& infos,
                         bool isHistory = false );

    //! @copydoc SimInput::GetStepValues()
    virtual void GetStepValues( UInt sequenceStep,
                                shared_ptr<ResultInfo> info,
                                std::map<UInt, Double>& steps,
                                bool isHistory = false );

    //! @copydoc SimInput::GetResultEntities()
    void GetResultEntities( UInt sequenceStep,
                            shared_ptr<ResultInfo> info,
                            StdVector<shared_ptr<EntityList> >& list,
                            bool isHistory = false );

    //! @copydoc SimInput::GetResultEntityNames()
    void GetResultEntityNames( UInt sequenceStep,
                               shared_ptr<ResultInfo> info,
                               StdVector<std::string> &names,
                               bool isHistory = false);

    //! @copydoc SimInput::GetResult()
    void GetResult( UInt sequenceStep,
                    UInt stepNum,
                    shared_ptr<BaseResult> result,
                    bool isHistory = false );
    //@}

  protected:

    //! Open a CGNS file
    void OpenFile(const std::string &fileName);
    
    //! Read contents of a CGNS file
    //! There is no return value, info is stored in member variables results_,
    //! stepMap_, resultSteps_, resultsEntities_, analysis_
    void ReadFileContents(const std::string &fileName);
    
    //! Read an unstructured zone from the current CGNS file and store its
    //! topology in the Grid object pointed to by this->mi_
    void ReadUnstructuredZone(Integer base, Integer zone);
    
    //! Print the name of an element type on the trace logging stream
    void PrintElementType(ElementType_t eType);
    
    //! Initialize static class variable SinInputCGNS::elemTypeMap_
    void InitElemTypeMap();
    
    //! Re-order CGNS element connectivity to match NACS ordering
    static void TranslateConnectivity(Elem::FEType feType,
                               cgsize_t* cgnsConn,
                               StdVector<UInt>& connect,
                               UInt offset = 0);
    
    //! Convert a polygon into proper element(s)
    static Elem::FEType PolygonToFE(cgsize_t numNodes, cgsize_t *nodes,
                                    StdVector< StdVector<UInt> > &elems,
                                    UInt offset = 0);
    
    //! Triangulate a polygon
    //! IMPORTANT: Extends the connectivity vector, nothing will be overwritten
    template<typename T>
    static void Triangulate(cgsize_t numNodes, T* nodes,
                            StdVector< StdVector<UInt> > &elems,
                            UInt offset = 0);
    
    //! Convert a polyhedron into a proper element
    static Elem::FEType PolyhedronToFE(cgsize_t numFaces, cgsize_t *faces,
        std::unordered_map<cgsize_t, StdVector< StdVector<UInt> > > & polyConn,
        StdVector< StdVector<UInt> > &connect, UInt offset = 0);
    
    //! Convert CGNS field name to NACS solution type
    static SolutionType GetFieldType(const char* field,
                                     ResultInfo::EntryType & entryType);
    
    //! Convert NACS solution type to CGNS field name 
    static void GetFieldNames(SolutionType solType,
                             StdVector<std::string> &fieldNames);
    
    //! Error callback routine called by CGNS library
    static void ErrorCallback(int isError, char *errMsg);
    
    //! Mapping of element types from CGNS to NACS
    static std::map<Integer, Elem::FEType> elemTypeMap_;
 
    //! Mapping from time value to file name
    std::map<Double, std::string>  fileNames_;
 
    //! File handle to be used with calls to the CGNS library
    Integer fileHandle_;

    //! File version of the CGNS file
    float fileVersion_ = -1.0;

    //! Dimension of the cells; 3 for volume cells, 2 for surface cells and 1
    //! for line cells. 
    Integer cellDim_;

    //! Number of bases in the CGNS file
    Integer numBases_;
    
    //! Base index that contains the grid and result data
    Integer baseIndex_;
    
    //! Number of zones in the current base
    Integer numZones_;
    
    //! Number of nodes
    UInt numNodes_;

    //! Number of elements
    UInt numElems_;

    //! Number of regions
    UInt numRegions_;

    //! Remember number of nodes in grid before we add our own ones.
    UInt nodeOffset_;

    //! Remember number of elements in grid before we add our own ones.
    UInt elemOffset_;

    //! Remember the number of polygons, because StarCCM+ apparently does not
    //! count these when exporting to CGNS.
    Integer numPolyFaces_;
    
    //! List of region names
    StdVector<std::string> regionNames_;
    
    // Number of nodes per region or group
    std::map<std::string, UInt> numNodesPerName_;
    
    // Number of elements per region or group
    std::map<std::string, UInt> numElemsPerName_;
    
    // Mapping from region or group name to zone index
    std::map<std::string, Integer> nameToZoneIdx_;
    
    //! Mapping (per zone) from NACS element number to CGNS cell index 
    StdVector< std::unordered_map<UInt, cgsize_t> > elemToCellIdx_;
    
    //! Range of cell indices per region or group
    std::map<std::string, std::pair<cgsize_t, cgsize_t> > cellRange_;
    
    //! Analysis type
    BasePDE::AnalysisType analysis_;
    
    //! List of available results
    StdVector< shared_ptr<ResultInfo> > results_;
    
    //! Mapping from step numbers to step values
    std::map<UInt, Double> stepMap_;

    //! Steps number per result
    std::map< shared_ptr<ResultInfo>, std::set<Double> > resultSteps_;

    //! Entities per result
    std::map< shared_ptr<ResultInfo>, std::set<std::string> > resultEntities_;
  };
}

#endif
