// =====================================================================================
// 
//       Filename:  SimInputCGNS.hh (originally FileReader_CGNS.hh in cplreader)
// 
//    Description:  File Reader for the cgns file format
//                  Based on OpenFOAM reader
// 
//        Version:  1.0
//        Created:  02/22/2011 12:51:01 PM
//       Revision:  none
//       Compiler:  g++
// 
//        Authors:  Simon Triebenbacher, simon.triebenbacher@tuwien.ac.at
//                  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  TU Wien, Universitaet Klagenfurt
// 
// =====================================================================================


#ifndef SIMINPUT_CGNS_HH
#define SIMINPUT_CGNS_HH

#include <string>
#include <set>

#include <cgnslib.h>

#include "DataInOut/SimInput.hh"


namespace CoupledField{

  class SimInputCGNS : public SimInput {
  public:
    
    //! Constructor
    SimInputCGNS( std::string fileName, PtrParamNode inputNode,
                  PtrParamNode infoNode );
    
    virtual ~SimInputCGNS();
    
    virtual void InitModule();
    
    virtual void ReadMesh(Grid* ptGrid);
    
    virtual UInt GetDim() { return dim_; }

    virtual UInt GetNumNodes()
    { EXCEPTION("Not implemented!"); }

    virtual UInt GetNumElems(Integer)
    { EXCEPTION("Not implemented!"); }

    virtual UInt GetNumRegions() { return numRegions_; }

    virtual UInt GetNumNamedNodes()
    { EXCEPTION("Not implemented!"); }

    virtual UInt GetNumNamedElems()
    { EXCEPTION("Not implemented!"); }

    virtual void GetAllRegionNames(StdVector<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >&)
    { EXCEPTION("Not implemented!"); }

    virtual void GetRegionNamesOfDim(StdVector<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >&, UInt)
    { EXCEPTION("Not implemented!"); }

    virtual void GetNodeNames(StdVector<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >&)
    { EXCEPTION("Not implemented!"); }

    virtual void GetElemNames(StdVector<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >&)
    { EXCEPTION("Not implemented!"); }

 #if 0
    //! get nodal values from the corresponding fluid datafile the new way
    virtual void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                 const std::vector<bool>& activeParts,
                                 const UInt timeStepIdx);

    virtual Double GetTimeStep(UInt stepNumber);
 #endif
    
  protected:
    UInt numElems_;
    UInt numVertices_;
    std::map<Double, std::string>  fileNames_;
    std::string fileDir_;
    std::map<CGNSLIB_H::ElementType_t,Elem::FEType> elemTypeMap_;
 
    //! Number of regions
    UInt numRegions_;

    //! Number of time steps
    UInt numSteps_;
   
    //! Dimension of the cells; 3 for volume cells, 2 for surface cells and 1
    //! for line cells. 
    Integer dim_;

    //! Number of coordinates required to define a vector in the field.
    Integer physDim_;

    //! Number of vertices, cells, and boundary vertices in each (index)-
    //! dimension.
    Integer vertSize_[9];

  private:
    void ReadCGNSDirectory(std::string dirname, std::map<Double, std::string> & fileNames);
    Integer GetFileHandle(std::string fName);
    void CheckFileValidity(Integer fileHandle);
    UInt MapCoordinateIndex(char* coordName);
    void ReadUnstructuredGrid(Integer fileHandle, Integer dim, cgsize_t* size);
    void PrintElementType(ElementType_t eType);
    void InitElemTypeMap();
    void CalcNumNodesPerRegion();
    UInt MapVelocityIndex(char* coordName);
    UInt MapFrictionIndex(char* coordName);
    UInt MapForceIndex(char* coordName);

    void AddRegionToGrid(RegionIdType& regionId,
                         const Elem::FEType feType,
                         const std::string& sectionName);

    void TranslateConnectivity(Elem::FEType feType,
                               cgsize_t* cgnsConn,
                               StdVector<UInt>& connect);

  };
}

#endif
