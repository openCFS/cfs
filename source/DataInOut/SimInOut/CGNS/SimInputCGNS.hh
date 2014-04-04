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

#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

#include "DataInOut/SimInput.hh"

namespace CoupledField{

  //! Class for reading unstructured simulation grids and CFD data from CGNS files.

  //! This class is for reading (C)FD (G)eneral (N)otation (S)ystem file(s).
  //!
  //! Among others, CGNS files can be created in the following pre-processors
  //! and CFD codes:
  //! ANSYS Workbench mesher: File -&gt; Export Mesh -&gt; CGNS Format. ANSYS
  //!          &lt;= 15.0 can only write linear element types. Named components
  //!          are written as CGNS sections using all capital letters.
  //! ANSYS CFX: Native CFX result files can be converted to CGNS by using
  //!          cfx5export -cgns input.res
  //! Pointwise: Uses CGNS as its native file format (www.pointwise.com).
  //!
  //! As the name CGNS states, the format is mostly used for exchange of CFD data.
  //! Since most CFD codes are using linear cell types, most pre- (ANSYS WB) and
  //! post-processors have problems handling quadratic element types often used in
  //! FE codes, even if the CGNS standard properly defines such element types.
  //! Nonetheless, this class is able to read all complete and incomplete element
  //! types supported by CGNS and CFS++. This class handles CGNS sections as regions.
  //! Therefore, it is not yet able to read CGNS files written with the SimOutputCGNS
  //! class, since it maps regions to zones.
  //!
  //! Besides supporting unstructured mesh reading, this class also supports reading
  //! the VelocityX, VelocityY, VelocityZ and Pressure fields for a single time step.
  //!
  //! \note ATTENTION: This class has been implemented for the exchange of CFD data
  //!                  during the E+H harmonic FSI project by Simon Triebenbacher.
  //!                  In this project we only need to transfer results for a single
  //!                  RANS step from CFX. Therefore, no support for more time/freq steps
  //!                  has been implemented and only one multisequence step is supported!
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

    // =========================================================================
    //  GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information

    //! Return multisequence steps and their analysistypes
    void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                   std::map<UInt, UInt>& numSteps,
                                   bool isHistory = false );

    //! Obtain list with result types in each sequence step
    void GetResultTypes( UInt sequenceStep, 
                         StdVector<shared_ptr<ResultInfo> >& infos,
                         bool isHistory = false );

    //! Return list with time / frequency values and step for a given result
    virtual void GetStepValues( UInt sequenceStep,
                                shared_ptr<ResultInfo> info,
                                std::map<UInt, Double>& steps,
                                bool isHistory = false );

    //! Return entitylist the result is defined on
    void GetResultEntities( UInt sequenceStep,
                            shared_ptr<ResultInfo> info,
                            StdVector<shared_ptr<EntityList> >& list,
                            bool isHistory = false );

    //! Fill pre-initialized results object with values of specified step
    void GetResult( UInt sequenceStep,
                    UInt stepNum,
                    shared_ptr<BaseResult> result,
                    bool isHistory = false );
    //@}


  protected:
    UInt numElems_;
    UInt numVertices_;
    std::map<Double, std::string>  fileNames_;
    std::string fileDir_;
    std::map<int,Elem::FEType> elemTypeMap_;
 
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
    cgsize_t vertSize_[9];

    //! Remember number of nodes in grid before we add our own ones.
    UInt nodeOffset_;

    float cgVersion_;
    
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
