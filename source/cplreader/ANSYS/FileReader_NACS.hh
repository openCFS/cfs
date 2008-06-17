#ifndef FILE_FILEREADER_NACS_2008
#define FILE_FILEREADER_NACS_2008

#include <cplreaderdefs.hh>
#include "filereader_ANSYS.hh"

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
