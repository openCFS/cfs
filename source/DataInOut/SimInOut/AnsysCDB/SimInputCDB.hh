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
    
    //! Get total number of regions
    UInt GetNumRegions();

    //! Get total number of named nodes
    UInt GetNumNamedNodes();

    //! Get total number of named elements
    UInt GetNumNamedElems();

    //@}
  
    // =======================================================================
    // ENTITY NAME ACCESS
    // =======================================================================
    //@{ \name Entity Name Access
  
    //! Get vector with all region names in mesh
    
    //! Returns a vector with the names of regions in the mesh of all
    //! dimensions.
    //! \param regionNames (output) vector containing names of regions
    //! \note Since the RegionIdType is guaranteed to be defined by
    //! a number type (UInt, UInt), the regionId of an element can
    //! be directly used as index to the regions-vector
    void GetAllRegionNames( StdVector<std::string> & regionNames );
    
    //! Get vector with region names of given dimension
    
    //! Returns a vector with the names of regions of a given dimension.
    //! This makes it possible to get for example all names of 
    //! 3D, 2D or 1D elements
    //! \param regionNames (output) vector containing names of regions
    //! \param dim (input) dimension of the region (1,2, or 3)
    void GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                              const UInt dim );
   
    //! Get vector with all names of named nodes

    //! Returns a vector which contains all names of named nodes.
    //! \param nodeNames (output) vector with names of named nodes
    void GetNodeNames( StdVector<std::string> & nodeNames );
  
    //! Get vector with all names of named elements

    //! Returns a vector which contains all names of named elements.
    //! \param elemNames (output) vector with names of named elements
    void GetElemNames( StdVector<std::string> & elemNames );

    //@}

    // =======================================================================
    // ENTITY ACCESS
    // =======================================================================
    //@{ \name Entity Access
    
    //! Get all nodal coordinates from 3D grid
    
    //! This method reads all nodal coordinates into a vector of 3D-Points.
    //! \param nodeCoords (output) vector containing nodal coordinates
    void GetCoordinates( std::vector< double > & nodeCoords );
    
    //! Get vector of nodes for each region

    //! This method reads the node numbers of each region into a 
    //! separate vector. 
    //! \param nodes (output) vector containing the node numbers for each
    //!                       region. The access is like \c elems 
    //!                       \c [regionNr] \c [nodeNr]
    //! \param regionId (output) vector containing the region Ids of the
    //!                          nodes corresponding to the outer index in the
    //!                          nodes vector
    void GetNodesOfRegions( std::vector<std::vector<UInt> > &nodes,
                            const std::vector<RegionIdType> & regionId );
    //! Read all elements of given dimension

    //! This method reads all elements of a given dimension (1D, 2D or 3D).
    //! The output is a vector of vectors, where the outer index corresponds
    //! to the different regions and the inner one to the different elements
    //! per region.
    //! \param elems (output) vector containing vectors of pointers to elements
    //!                       per region. The access is like \c elems 
    //!                       \c [regionNr] \c [elemeNr]
    //! \param regionId (output) vector containing the region Ids of the
    //!                          elements corresponding to the outer index in 
    //!                          the elems vector
    //! \param dim (input) dimension of the elements to be read (1,2 or 3)
    void GetElements( std::vector< std::vector<UInt> > & elems,
                      std::vector< std::vector<Elem::Elem::FEType> > & elemTypes,
                      std::vector< std::vector<UInt> > & elemNums,
                      std::vector<RegionIdType> & regionIds,
                      const UInt dim );

    //! Read all named nodes
    
    //! This method reads in all named nodes with their according names.
    //! \param nodes (output) vector containing node numbers for each region.
    //!                       The access is like \c nodes \c [nameNr] 
    //!                       \c [nodeNr]
    //! \param nodeNames (output) vector containing the corresponding
    //!                           node names 
    void GetNamedNodes(StdVector<StdVector<UInt> > & nodes,
                       StdVector<std::string> & nodeNames );

    //! Read all named elements

    //! This method reads in all named elements with their according names.
    //! \param elems (output) vector containing node numbers for each region.
    //!                       The access is like \c elems \c [nameNr] 
    //!                       \c [elemNr]
    //! \param elemNames (output) vector containing the corresponding
    //!                           element names 
    void GetNamedElems(StdVector<StdVector<UInt> > &elems,
                       StdVector<std::string> &elemNames );
    //@}
    

    // =========================================================================
    // GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information
    //! Fill pre-initialized results object with values of specified step
    virtual void GetResult( UInt sequenceStep,
                            UInt stepValue,
                            shared_ptr<BaseResult> result,
                            bool isHistory = false ) {
      
    } 
    //@}


  protected:
    
    Elem::FEType DegenTypeToNativeType(UInt type, UInt numNodes);
    void PreScanCDBFile();
    void SetAnsys2NacsElementTypeMap();

//    void ReadAnsys2NacsElementTypeMap();

    void ReadCoordinatesBlocked();
    void ReadCoordinatesUnBlocked();

    void ReadElementsBlocked();
    void ReadElementsUnBlocked();

    void ReadRegionsAndGroups();

    void StoreSingleNode(UInt fileNodeNum, double x, double y, double z, UInt &nodeNum,
                         UInt &numNodes, UInt &maxNodeNum);

    void GetRegionAndGroupNamesBlocked();
    void ReadBlockedRegionsAndGroups();
    void ReadUnBlockedRegionsAndGroups();
    void ReadUnBlockedElemRegionOrGroup(UInt cmLinePos, std::string regnam);
    void ReadUnBlockedNodeGroup(UInt cmLinePos, std::string regnam);
    void StoreRegion(std::string regname, UInt numdata, int* dataVal);
    void StoreRegion(std::string regname, UInt numdata, StdVector<UInt> dataVal);
    void StoreElemGroup(std::string regname, UInt numdata, int* dataVal);
    void StoreElemGroup(std::string regname, UInt numdata, StdVector<UInt> dataVal);
    void StoreNodeGroup(std::string regname, UInt numdata, int* dataVal);
    void StoreNodeGroup(std::string regname, UInt numdata, StdVector<UInt> dataVal);
    void GenElGroupFromPrsSrfElems();
    void GenElGroupFromSurfForceElems();
    void GenerateSurfElGroup(UInt elblock, StdVector<UInt> elemNumbers);

#ifdef WIN32
    bool GetLine(std::string& line, __int64 pos);
#else
    bool GetLine(std::string& line, std::streampos pos);
#endif
    bool GetNextLine(std::string& line);
    void OpenCDBFile(std::string fileName);
    void CloseCDBFile();

    void ResortNodes(std::vector<UInt>& elemNodes);
    void DegenerateElement(const Elem::FEType elemTypeIn,
    		               Elem::FEType& elemTypeOut,
                           std::vector<UInt>& elemNodes);

    bool strict_;
    bool degen_;

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
    std::map<UInt, UInt> elemRegionMap_;
//    std::vector< std::vector<UInt> > elemNumsOrig_;

    std::map<UInt, UInt> ansNacsMap_;

    void InitAnsys2NacsTypes();
    UInt ansysSubTypes_, ansysElTypes_;
    UInt *ans2NacsET_;
//    UInt *ansETDim_;
    std::vector<UInt> elemMaterials_;

    std::vector<Double> nodalCoords_;
    std::map<UInt, UInt> nodeNumsMap_;

    std::map<UInt, UInt> elemTypes_;
    std::map<UInt, UInt> elemNumsMap_;
    UInt maxNumElemNodes_;
    std::map<UInt, std::vector<UInt> > topology_;

    // file entries resulting from scan
    std::vector<std::string > lineETCmnds_;
    std::vector<std::string > lineKEYOPTCmnds_;
    UInt numNBlocks_, numEBlocks_, numNCmnds_, numENCmnds_, numCMBlocks_, numSFECmnds_;
    UInt numEselCmnds_, numNselCmnds_, numCMCmnds_, numPtsPSECmnds_, linePtsPSEStop_;
#ifdef WIN32
    std::vector<__int64> linePtsNBlocks_;
    std::vector<__int64> linePtsEBlocks_;
    std::vector<__int64> linePtsNCmnds_;
    std::vector<__int64> linePtsENCmnds_;
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
