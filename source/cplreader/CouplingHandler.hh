// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_COUPLINGHANDLER_2006
#define FILE_COUPLINGHANDLER_2006

#include <def_cplreader.hh>
#include "FileReader.hh"
#include "OutputWriter.hh"

namespace CoupledField
{

  class ElemIntegr;

  //! Class for mesh coupling

  class CouplingHandler
  {

    // Main interface methods for main program
  public:

    //!
    CouplingHandler(shared_ptr<FileReader> ptFileReader,
                    std::vector< shared_ptr<OutputWriter> >& outputWriters);

    //!
    virtual ~CouplingHandler();

    // Perform initialization
    void Init(int argc, char *argv[]);

    //! Hand over mesh from reader to writer
    void ConvertMesh();

    //! Performs the coupled computation phase
    void Couple();

    // Perform last operations
    void Finish();

    // Helper methods for handling vectors
  public:
    // Throw away unneccessary entries in vectors
    void ShrinkNodalVector(const UInt partitionIdx,
                           const UInt numDOFs,
                           const std::vector<Double>& input,
                           std::vector<Double>& output);

  private:

    void CalculateAcouSrcs(const int partitionIdx,
                           FlowDataType& flowData);


    void CollectElementNodes(const std::vector<UInt>& elems,
                             std::vector<UInt>& nodes,
                             UInt& dim);

    shared_ptr<FileReader>  ptFileReader_;
    std::vector<Double> nodalCoords_;
    std::vector<UInt> topology_;
    std::vector<UInt> elemTypes_;

    std::map<UInt, std::vector<UInt> > regionElems_;
    std::map<UInt, std::map<UInt, UInt> > regionNodeIndices_;
    std::vector<UInt> numRegionNodes_;

    UInt dim_;

    std::map<std::string, std::vector<UInt> > nodeGroups_;
    std::map<std::string, std::vector<UInt> > elemGroups_;

    std::map<UInt, ElemIntegr *> ptElemIntegr_;
    std::vector<std::string> outputFields_;
    std::vector<bool> activeParts_;

    OutputWriterVectorType outputWriters_;

    /** scales the FLUIDMECH_VELOCITY vector by the factor scaleFac.
     * @param nodalFlowData the data in which the reults are stored and which
     * should be scaled
     * @param scaleFac the factor by which it should be scaled
     * @param dof_idx the dof which should be scaled (x=0, y=1, z=2)
     */
    inline void velScale_(std::vector<FlowDataType>& nodalFlowData,
                   const Double scaleFac, const UInt dof_idx) const
    {
      Settings& settings = Settings::Instance();
      UInt numRegions = nodalFlowData.size();
      if (settings.GetString("type") == "OPENFOAM")
      {
        numRegions = 1;
      }
      for (UInt actRegion=0; actRegion < numRegions; ++actRegion)
      {
        FlowDataType& fd = nodalFlowData[actRegion];
        FlowDataPartStruct* fdps = &fd[FLUIDMECH_VELOCITY];
        const UInt sizeFdps = fdps->data.size();
        const UInt numDofs = fdps->dofNames.size();
        UInt i = dof_idx;
        while (i < sizeFdps)
        {
          fdps->data[i] *= scaleFac;
          i += numDofs;
        }
      }
    }
  };

} // end of namespace
#endif
