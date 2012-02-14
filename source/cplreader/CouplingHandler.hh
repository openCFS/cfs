// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_COUPLINGHANDLER_2006
#define FILE_COUPLINGHANDLER_2006

#include "def_cplreader.hh"
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

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //surface element section START
    //    only used if wanted by the user
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //! Method for calculating flux terms, basically we give the surface region
    //! along with the fow data struct. this should already do the job
    void CalculateSurfaceIntegral(const int surfRegionIdx,
                                      std::vector<FlowDataType> & flowData);

    //!Method for the determination of element neighbours for a given surface region
    //! the resulting datastructure stores for each element of the surface region the
    //! regionidx and the local elemIdx for the volume element
    void GetNeighborElements(std::map<UInt,std::pair<UInt,UInt> >& volNeighbors,UInt surfRegionIdx);

    //!This function calculates the normal of the given surface element number
    //! returns true on success or false if the surface element is not included in the
    //! volume element
    bool CalcSurfaceNormalOutVolume(UInt globalVolElemNum, UInt globalSurfaceElemNum);

    //!Helper method to calculate the surface normal
    void CalcSurfaceNormal( const std::vector<UInt>& nodeNums,
                               Vector<Double>& nVec);

    ///fills the surface region related datastructures if requested
    void PrepareSurfaceRegions();


    ///excplicit storage of volume regions
    std::vector<UInt> volumeRegionIndices_;
    ///explicit storage of surface regions and their active state
    std::map<UInt,bool> surfaceRegionIndices_;

    ///this is the map which stores for each global surface element number
    ///the corresponding volume region index along with the local volume
    ///element number
    ///we can do this, because of there would be more than one volume neighbor
    ///to a global surface neighbor we would have an internal surface which do not need to be considered
    std::map<UInt,std::pair<UInt,UInt> > surfaceNeighbors_;

    ///map storing for each global surface element number the
    ///normal vector out of volume
    std::map<UInt,Vector<Double> > surfaceNormals_;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //surface element section END 
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    shared_ptr<FileReader>  ptFileReader_;
    std::vector<Double> nodalCoords_;
    std::vector<UInt> topology_;
    std::vector<UInt> elemTypes_;

    std::map<UInt, std::vector<UInt> > regionElems_;
    std::map<UInt, std::map<UInt, UInt> > regionNodeIndices_;
    std::vector<UInt> numRegionNodes_;
    std::map<RegionIdType, UInt > regionDims_;
    std::map<std::string, std::vector<UInt> > regionNodes_;

    UInt dim_;
    bool doIntAverageCentre_;
    bool reduce_elementOrder_;

    std::map<std::string, std::vector<UInt> > nodeGroups_;
    std::map<std::string, std::vector<UInt> > elemGroups_;

    std::map<UInt, ElemIntegr *> ptElemIntegr_;
    std::vector<std::string> outputFields_;
    std::vector<bool> activeParts_;

    OutputWriterVectorType outputWriters_;

    /** scales the physical field (FLUIDMECH_VELOCITY, MECH_DISPLACEMT etc.)
     * vector by the factor scaleFac.
     * @param nodalFlowData the data in which the reults are stored and which
     * should be scaled
     * @param scaleFac the factor by which it should be scaled
     * @param dof_idx the dof which should be scaled (x=0, y=1, z=2)
     */
    inline void scale_PhysField_(std::vector<FlowDataType>& nodalFlowData,
        const SolutionType solType, const Double scaleFac, const UInt dof_idx) const
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
        if ( fd.find(solType) != fd.end() ) { // only if the solType exists
          FlowDataPartStruct* fdps = &fd[solType];
          const UInt sizeFdps = fdps->data.size();
          const UInt numDofs = fdps->dofNames.size();
          if ( dof_idx >= numDofs ) {
            EXCEPTION("You are trying to scale a DoF that does not exist!")
          }
          UInt i = dof_idx;
          while (i < sizeFdps)
          {
            fdps->data[i] *= scaleFac;
            i += numDofs;
          }
        }
      }
    }
    /**
     * Find nodes which occur in multiple regions. This is to calculate the
     * acoustic source correctly for nodes which lie in more than one region
     * @param regionMapNodes A std::map with a gathering of nodes on each region
     * @param multiNodes the return value. Return a map showing which node is
     * multiple on which each and which equation number
     */
    inline void findNodeMultiRegion(const std::map<std::string, std::vector<UInt>* >& regionMapNodes,
        std::map<UInt, std::map<std::string, UInt> >& multiNodes)
    {
      std::map<std::string, std::vector<UInt>* >::const_iterator iterRegionMapNodes = regionMapNodes.begin();
      std::map<std::string, std::vector<UInt>* >::const_iterator iterRegionMapNodes2;
      /* go over every region */
      for (; iterRegionMapNodes != regionMapNodes.end(); ++iterRegionMapNodes)
      {
        const std::vector<UInt>& firstRegVec = *iterRegionMapNodes->second;
        iterRegionMapNodes2 = iterRegionMapNodes;
        ++iterRegionMapNodes2;
        /* check every other region */
        for (; iterRegionMapNodes2 != regionMapNodes.end(); ++iterRegionMapNodes2)
        {
          const std::vector<UInt>& secRegVec = *iterRegionMapNodes2->second;
          /* bigger than 50 billion - sorry hardcoded*/
          if (firstRegVec.size() * secRegVec.size() > 50e9)
          {
            EXCEPTION("You are trying to calculate acoustic sources on multiple regions." << std::endl
                << "Cplreader needs to find the common interface of these regions," << std::endl
                << "but since the regions are so big (number of nodes) this will take to long!" << std::endl
                << "If you know what you are doing turn this option of with \"--doCalcMultiNodes 0\"")
          }
          /* nodes of first region */
          for (UInt i = 0; i < firstRegVec.size(); ++i)
          {
            const UInt& firstRegVecElem = firstRegVec[i];
            /* nodes of second region */
            for (UInt j = 0; j < secRegVec.size(); ++j)
            {
              /* find common nodes */
              if (firstRegVecElem == secRegVec[j])
              {
                /* take care of doublicate entrie */
                std::map<std::string, UInt>& regStringVec = multiNodes[firstRegVecElem];
                bool region1Registered = false;
                bool region2Registered = false;
                std::map<std::string, UInt>::iterator iterRegStringVec = regStringVec.begin();
                for (; iterRegStringVec != regStringVec.end(); ++iterRegStringVec)
                {
                  if (iterRegStringVec->first == iterRegionMapNodes->first)
                  {
                    region1Registered = true;
                  }
                  if (iterRegStringVec->first == iterRegionMapNodes2->first)
                  {
                    region2Registered = true;
                  }
                }
                if ( !region1Registered )
                {
                  regStringVec[iterRegionMapNodes->first] = i;
                }
                if ( !region2Registered )
                {
                  regStringVec[iterRegionMapNodes2->first] = j;
                }
              }
            } //end-for second region nodes
          }// end-for first region nodes
        }// end-for second region
      }// end-for first region

    }
  };

} // end of namespace
#endif
