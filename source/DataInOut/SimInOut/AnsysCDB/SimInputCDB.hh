// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SIMINPUTCDB_2014
#define FILE_SIMINPUTCDB_2014

#include <set>
#include <fstream>

#include <DataInOut/SimInput.hh>

namespace CoupledField {

  //! Class for reading in a .cdb mesh file created by the ANSYS CDWRITE command.

  //! Class, that is derived from class FileType for reading mesh-input data,
  //! which is produced by the ANSYS CDWRITE command.
  class SimInputCDB: virtual public SimInput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimInputCDB(std::string fileName, PtrParamNode inputNode, 
                PtrParamNode infoNode );
    
    //! Destructor
    virtual ~SimInputCDB();

    //@}

    virtual void InitModule();

    virtual void ReadMesh(Grid *mi);
  
    // =======================================================================
    // GENERAL MESH INFORMATION
    // =======================================================================
    //@{ \name General Mesh Information

    //! Return dimension of the mesh
    UInt GetDim();
    
    //! Get total number of nodes in mesh
    UInt GetNumNodes();
    
    //! Get total number of elements in mesh
    UInt GetNumElems( const Integer dim = 0 );
    //@}
  
  protected:

    Elem::FEType DegenTypeToNativeType(UInt type, UInt numNodes);

    //! Function  for performing  an  initial scan  of an  APDL  file for  the
    //! locations of nodes, elements, compontents, etc.
    void PreScanCDBFile();

    //! Since APDL is  a very general language  we have to ignore  some of the
    //! blocks during the initial scan.
    void IgnoreBlock(std::string& line, UInt& count);

    void AddIgnoreBlock(const std::string& begin,
                        const std::string& end);

    std::vector< std::pair<std::string, std::string> > ignoreBlks_;

    void SetAnsys2FETypeMap();

    void ReadCoordinatesBlocked();
    void ReadCoordinatesUnBlocked();

    void ReadElementsBlocked();
    void ReadElementsUnBlocked();
    void ReadElementsFromMatBlocks();

    void ReadRegionsAndGroups();

    void StoreSingleNode(UInt fileNodeNum, double x, double y, double z, UInt &nodeNum,
                         UInt &numNodes, UInt &maxNodeNum);

    void GetRegionAndGroupNamesBlocked();
    void ReadBlockedRegionsAndGroups();
    void ReadUnBlockedRegionsAndGroups();
    void ReadUnBlockedElemRegionOrGroup(UInt cmLinePos, std::string regnam);
    void ReadUnBlockedNodeGroup(UInt cmLinePos, std::string regnam);
    void StoreRegion(std::string regname, UInt numdata, int* dataVal);
    void StoreRegion(std::string regname, UInt numdata,
                     const StdVector<UInt>& dataVal);
    void StoreElemGroup(std::string regname, UInt numdata, int* dataVal);
    void StoreElemGroup(std::string regname, UInt numdata, 
                        const StdVector<UInt>& dataVal);
    void StoreNodeGroup(std::string regname, UInt numdata, int* dataVal);
    void StoreNodeGroup(std::string regname, UInt numdata, 
                        const StdVector<UInt>& dataVal);
    void GenElGroupFromPrsSrfElems();
    void GenElGroupFromSurfForceElems();
    void GenerateSurfElGroup(UInt elblock, StdVector<UInt> elemNumbers);

    //! When writing the input file from  Workbench FE Modeler all regions are
    //! written  only as  node  components.  Therefore,  at  least the  volume
    //! regions have to be reconstructed from the previously read node groups.
    void GenerateVolRegionsFromNodeComponents();

    //! Set containing  all ANSYS  element number,  which are  referenced from regions and element groups. 
    // This must not be an unordered_set, as otherwise generated elements, e.g. via GridCFS::GenerateExternalLayer() 
    // have different/undefined ordering. See discussion in https://gitlab.com/openCFS/cfs/-/merge_requests/338
    // TODO: std::set is quite slow with respect to memory layout and complexity. Probably a vector
    // would be faster if one does not need to check uniqueness too often. The read elements do not need to be sorted
    // (which set provides), but just need to be deterministic. So a vector would also work. 
    // But if the cdb file is not sorted by chance, one would again have to generate new reference files.
    std::set<UInt> referencedElems_;

    typedef std::vector< std::vector<UInt> > ElemFacesType;
    typedef std::map<Elem::FEType, ElemFacesType > FEType2ElemFacesType;
    typedef std::map<Elem::FEType, std::vector< Elem::FEType > > FEType2FacesFEType;
    
    typedef std::pair<Elem::FEType, std::vector<UInt> > FaceType;
    typedef std::multimap<std::size_t, FaceType > FacesType;

    FEType2ElemFacesType elemFaces_;
    FEType2FacesFEType elemFaceTypes_;

    typedef enum {NODE_CMP_TO_NAMED_ELEMS,
                  NODE_CMP_TO_SURF_REGIONS,
                  SURF_BCS_TO_NAMED_ELEMS} SurfElemsModeType;
    SurfElemsModeType surfElemsMode_;

    void InitElemFaceTypeMaps();
    void GenerateVolElemFaces(FacesType& faces);    
    void GenerateElemGroupsFromVolElemFacesAndNodeGroups(
      FacesType& faces,
      std::map<std::string, StdVector<UInt> >& elemGroupData,
      std::map<UInt, std::vector<UInt> >& surfTopo,
      std::map<UInt, Elem::FEType >& surfTypes);

#if(defined(WIN32))
    bool GetLine(std::string& line, __int64 pos);
#else
    bool GetLine(std::string& line, std::streampos pos);
#endif
    bool GetNextLine(std::string& line);
    UInt SplitLine(const std::string& line,
                   std::vector< std::string >& tokens,
                   const std::string& addSplitChars = "",
                   std::vector<int>* chunkSizes = NULL,
                   bool trim = false,
                   const std::string& trimChars = "") const;

    UInt DecodeBlockFormatLine(const std::string& line,
                               std::vector<int>& chunkSizes) const;
    
    void OpenCDBFile(std::string fileName);
    void CloseCDBFile();

    void ResortNodes(std::vector<UInt>& elemNodes);
    void DegenerateElement(const Elem::FEType elemTypeIn,
                           Elem::FEType& elemTypeOut,
                           std::vector<UInt>& elemNodes);

    //! Returns true if arr2[] is a subset of arr1[]
    static bool IsSubset(const UInt arr1[], const UInt arr2[],
                         UInt m, UInt n);

    bool strict_;
    bool degen_;
    bool surfElemsFromNodeComps_;    

    std::ifstream inFile_;
    std::streampos fSize_;

    std::vector<std::string> nodeGroupNames_;
    std::vector< StdVector<UInt> > nodeGroupData_;
    UInt numNodeGroups_;

    std::vector<std::string> elemGroupNames_;
    std::vector< StdVector<UInt> > elemGroupData_;
    UInt numElemGroups_;

    UInt numRegions_;
    std::vector<std::string> regionNames_;
    std::vector<StdVector<UInt> > regionData_;

    //! Map from ANSYS element type to FEType
    std::map<UInt, UInt> ans2FEMap_;

    void InitAnsys2FETypes();
    UInt ansysSubTypes_, ansysElTypes_;
    UInt *ans2FEType_;
//    UInt *ansETDim_;

    //! Coordinates of nodes.
    std::vector<Double> nodalCoords_;

    //! Map from ANSYS node number to internal node number
    std::map<UInt, UInt> nodeNumsMap_;

    //! Map from ANSYS element number to FEType
    std::map<UInt, UInt> elemTypes_;

    //! Map from ANSYS element number to ANSYS material number
    std::map<UInt, UInt> elemMaterials_;

    //! Map from ANSYS element number to region id
    std::map<UInt, UInt> elemRegionMap_;

    //! Map from ANSYS element number to connectivity
    std::map<UInt, std::vector<UInt> > topology_;

    //! Maximum number of element nodes
    UInt maxNumElemNodes_;


    // file entries resulting from scan
    std::vector<std::string > lineETCmnds_;
    std::vector<std::string > lineKEYOPTCmnds_;
    UInt numNBlocks_, numEBlocks_, numNCmnds_, numENCmnds_, numMATCmnds_, numTYPECmnds_, numCMBlocks_, numSFECmnds_;
    UInt numEselCmnds_, numNselCmnds_, numCMCmnds_, numPtsPSECmnds_, linePtsPSEStop_;
#if(defined(WIN32))
    std::vector<__int64> linePtsNBlocks_;
    std::vector<__int64> linePtsEBlocks_;
    std::vector<__int64> linePtsNCmnds_;
    std::vector<__int64> linePtsENCmnds_;
    std::vector<__int64> linePtsMATCmnds_;
    std::vector<__int64> linePtsTYPECmnds_;
    std::vector<__int64> linePtsCMBlocks_;
    std::vector<__int64> linePtsEselCmnds_;
    std::vector<__int64> linePtsNselCmnds_;
    std::vector<__int64> linePtsCMCmnds_;
    std::vector<__int64> linePtsPSECmnds_;
    std::vector<__int64> linePtsSFECmnds_;
    FILE *inStream_;
    __int64 inSize_;
    __int64 GetFilePosition();
#else
    std::vector<unsigned long> linePtsNBlocks_;
    std::vector<unsigned long> linePtsEBlocks_;
    std::vector<unsigned long> linePtsNCmnds_;
    std::vector<unsigned long> linePtsENCmnds_;
    std::vector<unsigned long> linePtsMATCmnds_;
    std::vector<unsigned long> linePtsTYPECmnds_;
    std::vector<unsigned long> linePtsCMBlocks_;
    std::vector<unsigned long> linePtsEselCmnds_;
    std::vector<unsigned long> linePtsNselCmnds_;
    std::vector<unsigned long> linePtsCMCmnds_;
    std::vector<unsigned long> linePtsPSECmnds_;
    std::vector<unsigned long> linePtsSFECmnds_;
    unsigned long GetFilePosition();
#endif

  };

}

#endif
