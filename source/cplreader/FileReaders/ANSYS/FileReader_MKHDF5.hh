// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_MKHDF5_2008
#define FILE_FILEREADER_MKHDF5_2008

#include <def_cplreader.hh>
#include "FileReader_ANSYS.hh"

namespace CoupledField
{

  class FileReader_MKHDF5 : public FileReader_ANSYS
  {
  public:

    //! Constructor
    FileReader_MKHDF5(const std::string& name,
                     const UInt dim,
                     const UInt numFiles,
                     const UInt startIndex);

    //! Deconstructor
    virtual ~FileReader_MKHDF5();

    virtual void Init();

    //! Get node groups
    virtual void GetNodeGroups(std::map<std::string,
                               std::vector<UInt> >& nodeGroups);

    //! Get element groups
    virtual void GetElemGroups(std::map<std::string,
                               std::vector<UInt> >& elemGroups);

  protected:
    Elem::FEType ANSYSTypeToFEType(UInt type, UInt numNodes, bool& readAnotherLine);

    void ReadEWRITEFile(std::string fileName,
                        std::map< std::string, std::vector<UInt> >& elemSets);

    void ReadNWRITEFile(std::string fileName,
                        std::map< std::string, std::vector<UInt> >& nodeSets);

    std::map<std::string, std::vector<UInt> > nodeGroups_;
    std::map<std::string, std::vector<UInt> > elemGroups_;
    
    UInt maxOrigElemNum_;
  };


} // end of namespace
#endif
