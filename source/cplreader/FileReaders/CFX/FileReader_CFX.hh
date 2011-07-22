// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_CFX_2006
#define FILE_FILEREADER_CFX_2006

#include <map>

#include <def_cplreader.hh>
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

    //! Deconstructor
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
    
  protected:

    void GetInfosFromCommand();

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
    
    UInt GetFaceOfElement( UInt elemType, UInt face,
                             std::vector<UInt>::const_iterator &elemConnect,
                             std::vector<UInt> &faceConnect );

  private:

    double timeStep;
    std::string timeStepStr;
    std::string solTimeUnit;
    std::string defFile;
    std::string userDataCFX_COMMANDS;
    std::map<std::string, std::string> exprMap;
    std::vector< std::string > transientFNs_;
    std::vector< int > timeStepNumbers_;
    std::vector< Elem::FEType > regionElemTypes_;
    std::ifstream inFile_;
    bool determineFloatDS_;
    UInt numElems_;

    // Vector of element numbers per region
    std::map<UInt, std::vector<UInt> > regionElems_;

    // number of volume regions (numRegions_ includes surface regions)
    UInt numVolRegions_;
    
    // names of regions
    std::vector<std::string> regionNames_;

    // CFX Release number for UserData in .h5 file
    std::string userDataCFXRelease;
    
    char fn[4096];
    int nerr;
    int whatfile;
    char what[80];
    char where[80];
    int when;
    int dattyp;
    int length;
    int nsize;
    int iopt;
    int ioptar;
    float rarr[80];
    int iarr[80];
    int nvx;
    char carr[80];
    int larr[80];
    double darr[80];
    char sarr[80];

    static std::vector<int> intvec;
    static std::vector<double> doublevec;
    static std::vector<float> floatvec;
    static std::vector<char> charvec;

  };
} // end of namespace
#endif
