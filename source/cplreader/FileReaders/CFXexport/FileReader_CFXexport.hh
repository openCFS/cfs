#ifndef FILEREADER_CFXEXPORT_H_
#define FILEREADER_CFXEXPORT_H_

#include <map>

#include <def_cplreader.hh>
#include "cplreader/FileReader.hh"

namespace CoupledField {

  class FileReader_CFXexport : public FileReader
  {
    public:

      //! Contructor
      FileReader_CFXexport(const std::string& name,
          const UInt dim,
          const UInt numFiles,
          const UInt startIndex);

      //! Destructor
      virtual ~FileReader_CFXexport();

      virtual void Init();

      //! get node coordinates from the corresponding file
      virtual void ReadNodalCoords(std::vector<Double> & NODECOORD);


      //! get topology information from the corresponding topology file
      virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                std::vector<UInt> & elemTypes);

      //! get nodal values from the corresponding fluid datafile
      virtual void ReadNodalValues(std::vector<FlowDataType > &nodalFlowdata,
                                   const std::vector<bool>& activeParts,
                                   const UInt timeStepIdx);

      virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                     const UInt regionIdx);

    protected:

      enum Quantity { DENSITY, PRESSURE, VELOCITY_X, VELOCITY_Y, VELOCITY_Z };

      UInt startIndex_;
      UInt colX_, colY_, colZ_;
      UInt colDens_, colPres_, colVelX_, colVelY_, colVelZ_;
      std::map<UInt, Quantity> col2Quan_;

      std::ifstream inFile_;

  };

} // namespace CoupledField

#endif /*FILEREADER_CFXEXPORT_H_*/
