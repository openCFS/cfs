// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_FLUENT_2011
#define FILE_FILEREADER_FLUENT_2011

#include "def_cplreader.hh"
#include "../FileReader_VTKMultiBlock.hh"

namespace CoupledField
{

  class FileReader_FLUENT : public FileReader_VTKMultiBlock
  {
  public:

    //! Constructor
    FileReader_FLUENT(const std::string& name,
                        const UInt dim,
                        const UInt numFiles);

    //! Deconstructor
    virtual ~FileReader_FLUENT();

    /* get nodal values from the corresponding fluid datafile the new way */
    void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                         const std::vector<bool>& activeParts,
                         const UInt timeStepIdx);

    //! get user data from file reader
    virtual void GetUserData(std::map<std::string, std::string>& userData);
  protected:
    void CreateReader();
    void EnableRegions();
    void GetTimeValues();
    void SetTimeValue(Double val);

  private:
    void GUnzipFile(const std::string& fn,
                    const std::string& outFN);

    std::map<UInt, std::string> tsCasFiles_;
    std::map<UInt, std::string> tsDatFiles_;
    std::string workDir_;
    bool oneCasToManyDat_;
  };

} // end of namespace
#endif
