#ifndef FILEREADER_CFXEXPORT_H_
#define FILEREADER_CFXEXPORT_H_

#include <map>

#include <cplreaderdefs.hh>
#include "../filereader.hh"

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
      virtual void ReadNodalCoords(std::vector<Double> & NODECOORD,
          const UInt partitionIdx);


      //! get topology information from the corresponding topology file
      virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
          std::vector<UInt> & numNodesPerElem,
          std::vector<UInt> & elemTypes,
          const UInt partitionIdx);

      //! get nodal values from the corresponding fluid datafile
      virtual void ReadNodalValues(std::vector<Double> & flowdata,
          const UInt partitionIdx,
          const UInt timeStepIdx);

    protected:
      
      enum Quantity { DENSITY, PRESSURE, VELOCITY_X, VELOCITY_Y, VELOCITY_Z };

      UInt startIndex_;
      UInt colX_, colY_, colZ_;
      UInt colDens_, colPres_, colVelX_, colVelY_, colVelZ_;
      std::map<UInt, Quantity> col2Quan_;

  };

} // namespace CoupledField

#endif /*FILEREADER_CFXEXPORT_H_*/
