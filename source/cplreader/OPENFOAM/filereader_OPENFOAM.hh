#ifndef FILE_FILEREADER_OPENFOAM_2007
#define FILE_FILEREADER_OPENFOAM_2007

#include <cplreaderdefs.hh>
#include "../filereader.hh"

class vtkOpenFOAMReader;

namespace CoupledField
{

  class FileReader_OPENFOAM : public FileReader
  {
  public:

    //! Constructor
    FileReader_OPENFOAM(const std::string& name,
                        const UInt dim,
                        const UInt numFiles);
    
    //! Deconstructor
    virtual ~FileReader_OPENFOAM();

    virtual void Init();

    //! get node coordinates from the corresponding file
    virtual void ReadNodalCoords(std::vector<Double> & NODECOORD,
                                 const UInt partitionIdx);


    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                              std::vector<UInt> & numNodesPerElem,
                              std::vector<UInt> & elemTypes,
                              const UInt partitionIdx);

    /* get nodal values from the corresponding fluid datafile the new way */
    void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
        const UInt timeStepIdx);

    virtual double GetTimeStep(UInt t);
    virtual std::string GetPartitionName(const UInt partitionIdx);

    //! get user data from file reader
    virtual void GetUserData(std::map<std::string, std::string>& userData);

  protected:
      FEType VTKCellTypeToFEType(UInt cellType);

      std::vector<Integer> dataColumns_;
      /* nodalCoords_ should store the original mesh, which may be needed if we have
       * a moving mesh*/
      std::vector<double> nodalCoords_;
      vtkOpenFOAMReader* reader_;

      /* calculates the mechDisplacement and writes them into newCoords
       * \param origin A vector which has the original position of each node
       * \param newCoords A vector which has the new position of each node and
       * which will store the mechDisplacement after the method call */
      void calcMechDisplacement(const std::vector<double>& origin, \
          std::vector<double>& newCoords) const;
  };

      
} // end of namespace
#endif
