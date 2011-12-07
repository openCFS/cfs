// =====================================================================================
// 
//       Filename:  FileReader_CGNS.hh
// 
//    Description:  File Reader for the cgns file format
//                  Based on OpenFOAM reader
// 
//        Version:  1.0
//        Created:  02/22/2011 12:51:01 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================


#ifndef FILE_FILEREADER_CGNS
#define FILE_FILEREADER_CGNS


#include <def_cplreader.hh>
#include "cgnslib.h"
#include "cplreader/FileReader.hh"

#include <string>
#include <set>

namespace CoupledField{

  class FileReader_CGNS : public FileReader {
    public:

      struct CGNSElem{
         StdVector<UInt> connect;
         //one-based element number
         UInt elemNum;
         Elem::FEType eType;
      };

      //! Constructor
      FileReader_CGNS(  const std::string& name,
                        const UInt dim,
                        const UInt numFiles);

      virtual ~FileReader_CGNS();

      virtual void Init();

      //! get node coordinates from the corresponding file
      virtual void ReadNodalCoords(std::vector<Double> & NODECOORD);


      //! get topology information from the  topology file
      virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                std::vector<UInt> & elemTypes);

      //! get nodal values from the corresponding fluid datafile the new way
      virtual void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                   const std::vector<bool>& activeParts,
                                   const UInt timeStepIdx);

      virtual Double GetTimeStep(UInt stepNumber);

      virtual std::string GetRegionName(const UInt partitionIdx);

      //! get user data from file reader
      virtual void GetUserData(std::map<std::string, std::string>& userData);

      virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx);

    protected:
      UInt numElems_;
      UInt numVertices_;
      std::map<Double, std::string>  fileNames_;
      std::string fileDir_;
      bool gridInitialized_;
      //storing the coordinates
      StdVector< StdVector<Double> > nodeCoords_;
      std::map<CGNSLIB_H::ElementType_t,Elem::FEType> elemTypeMap_;
      std::map<Integer,std::string> regionIndexToNameMap_;
      std::map<Integer,StdVector<CGNSElem> > elemRegionMap_;
      std::map<Integer,Elem::FEType > elemTypeToRegionMap_;
      std::map<Integer,std::set<UInt> > nodesPerRegionMap_;

    private:
      void ReadCGNSDirectory(std::string dirname, std::map<Double, std::string> & fileNames);
      Integer GetFileHandle(std::string fName);
      void CheckFileValidity(Integer fileHandle);
      UInt MapCoordinateIndex(char* coordName);
      void ReadUnstructuredGrid(Integer fileHandle, Integer dim, cgsize_t* size);
      void PrintElementType(ElementType_t eType);
      void InitElemTypeMap();
      void ReadGrid();
      void CalcNumNodesPerRegion();
      UInt MapVelocityIndex(char* coordName);
      UInt MapFrictionIndex(char* coordName);

  };
}

#endif
