// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_CFX_2006
#define FILE_FILEREADER_CFX_2006

#include <map>
#include <stack>
#include <string>
#include <vector>

#include "def_cplreader.hh"
#include "cplreader/FileReader.hh"

namespace CoupledField
{

  class FileReader_CFX : public FileReader
  {
  public:

    //! Constructor
    FileReader_CFX(const std::string& name,
                   const UInt dim,
                   const UInt numFiles);

    //! Destructor
    virtual ~FileReader_CFX();

    virtual void Init();

    //! get node coordinates from the corresponding file
    virtual void ReadNodalCoords(std::vector<Double> & NODECOORD);


    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                 std::vector<UInt> & elemTypes);

    //! get nodal values from the corresponding fluid datafile the new way
    virtual void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                 const std::vector<bool>& activeParts,
                                 const UInt timeStepIdx);

    //! return name of region
    virtual std::string GetRegionName(const UInt regionIdx);

    //! get user data from file reader
    virtual void GetUserData(std::map<std::string, std::string>& userData);

    virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                       const UInt regionIdx);

    //! Get node groups
    virtual void GetNodeGroups(std::map<std::string,
                                        std::vector<UInt> >& nodeGroups){};

    //! Get element groups
    virtual void GetElemGroups(std::map<std::string,
                                        std::vector<UInt> >& elemGroups){};
    
    UInt GetNumRegions() { return numRegions_;}
    
  protected:

    void GetInfosFromCommand(std::vector<std::string> commands);
    void ParseCCLLine(const std::string& line);

    PtrParamNode rootNode_;
    PtrParamNode currNode_;
    PtrParamNode latestNode_;
    std::stack<PtrParamNode> parents_;
    bool multiLine_;

    void ParseCommand(std::vector<char>& cmdstr,
                      int& pos,
                      std::string& cmd,
                      std::string& attrib,
                      std::string indent,
                      std::ostream& outFile);

    void ParseOption(std::vector<char>& cmdstr,
                     int& pos,
                     std::string& option,
                     std::string& value,
                     std::string indent,
                     std::ostream& outFile);

    void IOErrorToString(int ioerr, std::string& errStr);

    void CheckTransientFiles();
    
    //! Get element type and connectivity of some face of a volume element
    UInt GetFaceOfElement( UInt elemType, UInt face,
                             std::vector<UInt>::const_iterator &elemConnect,
                             std::vector<UInt> &faceConnect );
    
    //! Extracts result data of a volume region from the whole zone's data
    template<typename A, typename B>
    void MapZoneResultToVolume( std::vector<A> &volResult,
                                const std::vector<B> &zoneResult,
                                const std::set<UInt> &regionNodes,
                                UInt numDOFs );

    //! Opens a CFX file
    bool OpenCFXFile(const std::string &filename, bool throwEx = true);
    
    //! Closes the currently opened CFX file
    bool CloseCFXFile(bool throwEx = true);
    
    //! Reads a scalar (int, float, double) value from a CFX file
    template<typename T>
    T ReadScalar(std::string what,
                 std::string where = "EVERY",
                 Integer when = -1,
                 bool throwEx = true);

    //! Reads a single-line string from a CFX file
    std::string ReadString(std::string what,
                           std::string where = "EVERY",
                           Integer when = -1,
                           bool throwEx = true);
    
    //! Reads a short vector (i.e. stored without compression) from a CFX file
    template<typename T>
    bool ReadShortVector(std::vector<T> &out,
                         Integer size,
                         std::string what,
                         std::string where = "EVERY",
                         Integer when = -1,
                         bool throwEx = true);

    //! Reads a long vector (i.e. stored with compression) from a CFX file
    template<typename T>
    bool ReadLongVector(std::vector<T> &out,
                        Integer size,
                        std::string what,
                        std::string where = "EVERY",
                        Integer when = -1,
                        bool throwEx = true);
    
  private:

    double timeStep_;
    std::string timeStepStr_;
    std::string solTimeUnit_;
    std::string defFile_;
    std::string userDataCFX_COMMANDS_;
    std::map<std::string, std::string> exprMap_;
    std::vector< std::string > trnFilenames_;
    std::vector< int > timeStepNumbers_;
    std::ifstream inFile_;
    bool determineFloatDS_;
    
    //! Number of zones ("Domains" in CFX GUI)
    UInt numZones_;
    
    //! Number of volume regions (numRegions_ includes surface regions)
    UInt numVolumes_;
    
    //! Number of boundary condition patches (surfRegion in CFS++ nomenclature)
    UInt numBCPs_;

    //! Total number of Nodes
    UInt numNodes_;
    
    //! Total number of Elements
    UInt numElems_;
    
    //! Number of volume elements
    UInt numVolElems_;
    
    //! Number of surface elements
    UInt numSurfElems_;
    
    //! Offset per zone for node numbers
    std::vector<UInt> nodeOffsetPerZone_;
    
    //! Offset per zone for element numbers
    std::vector<UInt> elemOffsetPerZone_;
    
    //! Struct for storing info about an element set
    typedef struct {
      UInt          num;
      UInt          numElems;
      Elem::FEType  elemType;
    } ElemSet;
    
    //! Struct for grouping element sets into volumes
    typedef struct {
      UInt                  zone;
      std::vector<ElemSet>  elemSets;
      std::set<UInt>        nodes;
    } Volume;
    
    //! Volume info
    std::vector<Volume> volumes_;
    
    //! Struct for storing info about a face set
    typedef struct {
      UInt          zone;
      UInt          num;
      UInt          numElems;
    } FaceSet;
    
    //! Face sets grouped to BC patches
    std::vector< std::vector<FaceSet> > faceSetsPerBCP_;
    
    //! Vector of element numbers per region
    std::vector< std::vector<UInt> > regionElems_;

    //! Names of regions
    std::vector<std::string> regionNames_;
    
    //! CFX Release number for UserData in .h5 file
    std::string userDataCFXRelease_;
    
    char fn_[4096];
    int whatfile_;
#   define IO_BUFSIZE 80
    char what_[IO_BUFSIZE];
    char where_[IO_BUFSIZE];
    float rarr_[IO_BUFSIZE];
    int iarr_[IO_BUFSIZE];
    char carr_[2*IO_BUFSIZE];
    double darr_[IO_BUFSIZE];
    char sarr_[2*IO_BUFSIZE];

    static std::vector<int> intvec_;
    static std::vector<double> doublevec_;
    static std::vector<float> floatvec_;
    static std::vector<char> charvec_;

  };
} // end of namespace
#endif
