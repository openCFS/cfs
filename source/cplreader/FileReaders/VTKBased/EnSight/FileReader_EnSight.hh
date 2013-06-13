// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_ENSIGHT_2011
#define FILE_FILEREADER_ENSIGHT_2011

#include "def_cplreader.hh"
#include "../FileReader_VTKMultiBlock.hh"

namespace CoupledField
{

  class FileReader_EnSight : public FileReader_VTKMultiBlock
  {
  public:

    //! Constructor
    FileReader_EnSight(const std::string& name,
                        const UInt dim,
                        const UInt numFiles);

    //! Deconstructor
    virtual ~FileReader_EnSight();

    /* get nodal values from the corresponding fluid datafile the new way */
    void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                         const std::vector<bool>& activeParts,
                         const UInt timeStepIdx);

    /* get element values from the corresponding fluid datafile and interpolate to nodes */
    virtual void ReadElemValues(std::vector<FlowDataType>& nodalFlowData,
                                const std::vector<bool>& activeParts,
                                const UInt timeStepIdx);

    //! get user data from file reader
    virtual void GetUserData(std::map<std::string, std::string>& userData);
  protected:
    void CreateReader();
    void EnableRegions();
    void GetTimeValues();
    void SetTimeValue(Double val);

    virtual void InitElemNodeMapping();

  private:
    std::string vx_, vy_, vz_;
    std::string pres_;
    std::string presD2_;
  };

} // end of namespace
#endif
