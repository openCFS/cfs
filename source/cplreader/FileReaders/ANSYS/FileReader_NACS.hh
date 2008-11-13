// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_NACS_2008
#define FILE_FILEREADER_NACS_2008

#include <def_cplreader.hh>
#include "FileReader_ANSYS.hh"

namespace CoupledField
{

  class FileReader_NACS : public FileReader_ANSYS
  {
  public:

    //! Constructor
    FileReader_NACS(const std::string& name,
                     const UInt dim,
                     const UInt numFiles,
                     const UInt startIndex);

    //! Deconstructor
    virtual ~FileReader_NACS();

    virtual void Init();

    //! Get node groups
    virtual void GetNodeGroups(std::map<std::string,
                                        std::vector<UInt> >& nodeGroups);

    //! Get element groups
    virtual void GetElemGroups(std::map<std::string,
                                        std::vector<UInt> >& elemGroups);

  protected:

    FEType DegenTypeToNativeType(UInt type, UInt numNodes);
    void GetRegionAndGroupNames();

    std::vector<std::string> nodeGroupNames_;
    std::vector<std::string> elemGroupNames_;
  };


} // end of namespace
#endif
