<<<<<<< .working
#ifndef FILE_OUTPUTWRITER_MPCCI_2008
#define FILE_OUTPUTWRITER_MPCCI_2008

#include "cplreader/OutputWriter.hh"

namespace CoupledField
{
  /*!
    This class handles the subroutines calls concerning MpCCI for coupling the fluid and acoustic computations.
  */

  class OutputWriter_MpCCI : public OutputWriter {
  public:
    OutputWriter_MpCCI();

    virtual ~OutputWriter_MpCCI();

    // Init output writer
    virtual void Init(int argc, char** argv,
                      const std::string& outputPath);

    // Write nodal coordinates
    virtual void WriteNodalCoords(const std::vector<Double>& coords);

    // Write element connectivity/topology
    virtual void WriteTopology(UInt maxNumElemNodes,
                               const std::vector<UInt>& elemTypes,
                               const std::vector<UInt>& elems);

    // Write regions
    virtual void WriteRegion(std::string regionName,
                             const std::vector<UInt>& regionElems,
                             const std::vector<UInt>& regionNodes,
                             UInt regionDim);

    // Write node groups
    virtual void WriteNodeGroups(const std::map< std::string, std::vector<UInt> >&
                                 nodeGroups);

    virtual void WriteElemGroups(const std::map< std::string, std::vector<UInt> >&
                                 elemGroups,
                                 const std::map< std::string, std::vector<UInt> >&
                                 groupNodes,
                                 const std::map< std::string, UInt >&
                                 groupDims);
    

    virtual void WriteFlowData(CouplingHandler* cplHandler,
                               UInt actRegion,
                               const std::vector<std::string>& outputFields);

    virtual void WriteUserData(const std::map<std::string, std::string>& userData);

    virtual void BeginResults(const std::vector<FlowDataType>* flowData);
    
    virtual void BeginStep(UInt stepNum, Double stepValue);
    virtual void EndStep();

    virtual void EndResults();

  private:
    // MpCCI specific member functions
    int ElemTypes2MpCCI(FEType et);

  private:
    //!MpCCI
    int meshId_;
    int GlobalDim_;
    int MpCCIprocess_;
    int quantityId1_;
    int quantityDim1_;
    int quantityId2_;
    int quantityDim2_;
  };

} // end of namespace
#endif
=======
// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_OUTPUTWRITER_MPCCI_2008
#define FILE_OUTPUTWRITER_MPCCI_2008

#include "cplreader/OutputWriter.hh"

namespace CoupledField
{
  /*!
    This class handles the subroutines calls concerning MpCCI for coupling the fluid and acoustic computations.
  */

  class OutputWriter_MpCCI : public OutputWriter {
  public:
    OutputWriter_MpCCI();

    virtual ~OutputWriter_MpCCI();

    // Init output writer
    virtual void Init(int argc, char** argv,
                      const std::string& outputPath);

    // Write nodal coordinates
    virtual void WriteNodalCoords(const std::vector<Double>& coords);

    // Write element connectivity/topology
    virtual void WriteTopology(UInt maxNumElemNodes,
                               const std::vector<UInt>& elemTypes,
                               const std::vector<UInt>& elems);

    // Write regions
    virtual void WriteRegion(std::string regionName,
                             const std::vector<UInt>& regionElems,
                             const std::vector<UInt>& regionNodes,
                             UInt regionDim);

    // Write node groups
    virtual void WriteNodeGroups(const std::map< std::string, std::vector<UInt> >&
                                 nodeGroups);

    virtual void WriteElemGroups(const std::map< std::string, std::vector<UInt> >&
                                 elemGroups,
                                 const std::map< std::string, std::vector<UInt> >&
                                 groupNodes,
                                 const std::map< std::string, UInt >&
                                 groupDims);
    

    virtual void WriteFlowData(CouplingHandler* cplHandler,
                               UInt actRegion,
                               const std::vector<std::string>& outputFields);

    virtual void WriteUserData(const std::map<std::string, std::string>& userData);

    virtual void BeginResults(const std::vector<FlowDataType>* flowData);
    
    virtual void BeginStep(UInt stepNum, Double stepValue);
    virtual void EndStep();

    virtual void EndResults();

  private:
    // MpCCI specific member functions
    int ElemTypes2MpCCI(FEType et);

  private:
    //!MpCCI
    int meshId_;
    int GlobalDim_;
    int MpCCIprocess_;
    int quantityId1_;
    int quantityDim1_;
    int quantityId2_;
    int quantityDim2_;
  };

} // end of namespace
#endif
>>>>>>> .merge-right.r8654
