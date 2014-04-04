// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HDF5_COMMON_HH
#define FILE_HDF5_COMMON_HH

#include <exception>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <H5Cpp.h>

using namespace boost;

namespace H5CFS {

  // =======================================================================
  //  ENUM DEFINITIONS
  // =======================================================================
  typedef enum  
  {
    ET_UNDEF      = 0,
    ET_POINT      = 1,
    ET_LINE2      = 2,
    ET_LINE3      = 3,
    ET_TRIA3      = 4,
    ET_TRIA6      = 5,
    ET_QUAD4      = 6,
    ET_QUAD8      = 7,
    ET_QUAD9      = 8,
    ET_TET4       = 9,
    ET_TET10      = 10,
    ET_HEXA8      = 11,
    ET_HEXA20     = 12,
    ET_HEXA27     = 13,
    ET_PYRA5      = 14,
    ET_PYRA13     = 15,
    ET_PYRA14     = 19,
    ET_WEDGE6     = 16,
    ET_WEDGE15    = 17,
    ET_WEDGE18    = 18
  } ElemType;
  
  const unsigned int NUM_ELEM_NODES[] = 
  {
    0,  // ET_UNDEF  
    1,  // ET_POINT
    2,  // ET_LINE2  
    3,  // ET_LINE3  
    3,  // ET_TRIA3  
    6,  // ET_TRIA6  
    4,  // ET_QUAD4  
    8,  // ET_QUAD8  
    9,  // ET_QUAD9  
    4,  // ET_TET4   
    10, // ET_TET10  
    8,  // ET_HEXA8  
    20, // ET_HEXA20 
    27, // ET_HEXA27 
    5,  // ET_PYRA5  
    13, // ET_PYRA13 
    6,  // ET_WEDGE6
    15, // ET_WEDGE15 
    18, // ET_WEDGE18
    14  // ET_PYRA14 
  };
  
  //! Typedef describing the entryType of the result
  typedef enum { 
    UNKNOWN = 0, 
    SCALAR = 1, 
    VECTOR = 2, 
    TENSOR = 3, 
    STRING = 4
  } EntryType;

  //! Entity types a result is related to
  typedef enum{ 
    NODE,
    EDGE, 
    FACE, 
    ELEMENT, 
    SURF_ELEM,
    //PFEM = 6, 
    REGION, 
    SURF_REGION,
    NODELIST,
    //ELEM_GROUP = 10,
    COIL, 
    FREE
  } EntityType;

  //! Type of analysis
  typedef enum {
    STATIC = 1,
    TRANSIENT = 2,
    HARMONIC = 3,
    EIGENFREQUENCY = 4
  } AnalysisType;

  // =======================================================================
  //  HELPER CLASSES
  // =======================================================================
  
  //! Description of a result 
  struct ResultInfo {
    std::string name;
    std::string unit;
    std::vector<std::string> dofNames;
    EntryType entryType;
    EntityType listType;
    std::string listName;
    bool isHistory;
  };
  
  //! Result object itself
  struct Result {
    boost::shared_ptr<ResultInfo> resultInfo;
    bool isComplex;
    std::vector<double> realVals;
    std::vector<double> imagVals;
  };
  
  // =======================================================================
  //  HELPER METHODS 
  // =======================================================================

  // helper function to generate a meaningful std::string exception
#define H5CFS_EXCEPTION(STR){                                          \
  std::ostringstream ostr;                                             \
  ostr << STR << "\nIn file " << __FILE__ << " at line " << __LINE__;  \
  throw ostr.str();                                                    \
}

#define H5CFS_RETHROW_EXCEPTION(REASON, STR){                         \
  std::ostringstream ostr;                                            \
  ostr << REASON << "\n";                                             \
  ostr << STR << "\nIn file " << __FILE__ << " at line " << __LINE__; \
  throw ostr.str();                                                   \
}

  // define commodity method for converting a hdf5 exception
  // to a std::string exception
#define H5_CATCH( STR )                                                 \
  catch (H5::Exception& h5Ex ) {                                        \
    H5CFS_EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );            \
  }  
  
  
  //! Common I/O routines for accessing hdf5 files in CFS format
  
  //! This helper class provides convenience methods and abstractions from
  //! the hdf5-native datatypes to native C++ ones, like e.g. reading and
  //! writing arrays (int, double, string). It also allows the control
  //! of the compression and chuncking settings of the hdf5 format.
  //! Currently this api can treat 1D and 2D arrays, as well as compound
  //! data types with the help of the boost::any-library.
  class Hdf5Common {

  public:

    // =======================================================================
    //  TYPE DEFINITIONS
    // =======================================================================
    typedef std::vector<std::pair<std::string, boost::any> >
    CompoundType;

    typedef std::vector<std::pair<std::string, std::vector<boost::any> >  >
    CompoundArrayType;

    // =======================================================================
    //  WRITE METHODS
    // =======================================================================

    //! Create attribute
    template<typename TYPE>
    static void WriteAttribute( H5::H5Object& obj,
                                const std::string& name,
                                const TYPE& data,
                                const H5::PropList &create_plist
                                = H5::PropList::DEFAULT );

    //! Commodity function for writing a 1D dataset
    template<typename TYPE>
    static void Write1DArray( H5::CommonFG &loc,
                              const std::string& name,
                              unsigned int size,
                              const TYPE * buffer,
                              const H5::DSetCreatPropList &create_plist
                              = H5::DSetCreatPropList::DEFAULT );

    //! Reserve space for a 1D array, but do not write anything to it
    template<typename TYPE>
    static void Reserve1DArray( H5::CommonFG &loc,
                                const std::string& name,
                                unsigned int size,
                                const H5::DSetCreatPropList &create_plist
                                = H5::DSetCreatPropList::DEFAULT );

    //! Set entries entries in 1D array
    template<typename TYPE>
    static void SetEntries1DArray( H5::CommonFG &loc,
                                   const std::string& name,
                                   unsigned int begin, unsigned int end,
                                   const TYPE * buffer  );

    //! Commodity function for writing a 2D dataset
    template<typename TYPE>
    static void Write2DArray( H5::CommonFG &loc,
                              const std::string& name,
                              unsigned int rowSize,
                              unsigned int colSize,
                              const TYPE * buffer,
                              const H5::DSetCreatPropList &create_plist
                              = H5::DSetCreatPropList::DEFAULT );

    //! Reserve space for a 2D array, but do not write anything to it
    template<typename TYPE>
    static void Reserve2DArray( H5::CommonFG &loc,
                                const std::string& name,
                                unsigned int rowSize,
                                unsigned int colSize,
                                const H5::DSetCreatPropList &create_plist
                                = H5::DSetCreatPropList::DEFAULT );

    //! Commodity function for writing a 2D dataset
    template<typename TYPE>
    static void SetEntries2DArray( H5::CommonFG &loc,
                                   const std::string& name,
                                   unsigned int rowBegin, unsigned int rowEnd,
                                   unsigned int colBegin, unsigned int colEnd,
                                   const TYPE * buffer );

    //! Commoditiy function for writing a scalar compound (1x1 rank)
    static void WriteCompound( H5::CommonFG& loc,
                               const std::string& name,
                               const CompoundType comp,
                               const H5::DSetCreatPropList &create_plist
                               = H5::DSetCreatPropList::DEFAULT );

    // =======================================================================
    //  READ METHODS
    // =======================================================================

    //! Get name of an object in a group with given index

    //! This method retrieves the name of an object ob a group with given
    //! index.
    //! \note This methods replaces the buggy method
    //! H5::CommonFG::getObjnameByIdx(), which skips by default the last
    //! character of the groupname
    static std::string GetObjNameByIdx( const H5::CommonFG& loc, hsize_t idx );

    //! Read data from an attribute
    template<typename TYPE>
    static void ReadAttribute( H5::H5Object& obj,
                               const std::string& name,
                               TYPE& data );

    //! Retrieve rank and dimensionality of a dataset and return
    //! total number of entries in the dataset
    static std::vector<unsigned int> GetArrayDims( const H5::CommonFG &loc,
                                                   const std::string& name );

    //! Return number of entries of a dataset / rray
    static unsigned int GetNumEntries( const H5::CommonFG &loc,
                                       const std::string& name );

    //! Retrieve array data from a dataset

    //! Read data from a an dataset of arbitrary dimension into a linear buffer
    //! Note, that the memory has to be allocated from outside
    template<typename TYPE>
    static void ReadArray( H5::CommonFG &loc,
                           const std::string& name,
                           TYPE* data );


    //! Read data from a dataset into a stl vector
    template<typename TYPE>
    static void ReadArray( H5::CommonFG &loc,
                           const std::string& name,
                           std::vector<TYPE>& data );

//     //! Read data from a dataset into a vector
//     template<typename TYPE>
//     static void ReadArray( H5::CommonFG &loc,
//                            const std::string& name,
//                            Vector<TYPE>& data );

//     //! Read data from a dataset into a matrix
//     template<typename TYPE>
//     static void ReadArray( H5::CommonFG& loc,
//                            const std::string& name,
//                            Matrix<TYPE>& buffer );


    // =======================================================================
    //  GENERAL ACCESS METHODS
    // =======================================================================

    //! Obtain grid result group for specified multisequence step
    static H5::Group GetMultiStepGroup( H5::H5File& file, 
                                        unsigned int msStep, 
                                        bool isHistory );

    //! Obtain grid result group for specified step in a given multistep
    static H5::Group GetStepGroup( H5::H5File& file, 
                                   unsigned int msStep, unsigned int 
                                   stepNum );
    
//    // =======================================================================
//    //  CONVERSION METHODS
//    // =======================================================================
//
//    //! Map SimOuput::Capability class to hdf5 type
//    static int MapCapabilityType( SimOutput::Capability c );
//    
//    //! Map hdf5 representation of simOutput::Capability to enum representat
//    static SimOutput::Capability MapCapabilityType( int c );
//    
    //! Map EntityType enum to hdf5 type
    static int MapUnknownType( EntityType t );

    //! Map EntityUnknown hdf5 type to enum
    static EntityType MapUnknownType( int t );
    
    //! Map EntityType enum to string representation
    static std::string MapUnknownTypeAsString( EntityType t );

    //! Map entryType from enum to hdf5 type
    static int MapEntryType( EntryType t );

    //! Map entryType from hdf5 type to enum
    static EntryType MapEntryType( int t );

    // =======================================================================
    //  MISCELANEOUS METHODS
    // =======================================================================

    //! Set chunksize to be used for Array data
    static void SetMaxChunkSize( unsigned int chunkSize );

    //! Check for open objects in the hdf5 file and close them
    static void CheckOpenObjects(H5::H5File& file, bool verbose);
    
  private:

    // =======================================================================
    //  MAXIMUM CHUNKSIZE
    // =======================================================================
    static hsize_t maxChunkSize_;

    // =======================================================================
    //  HELPER CLASSES FOR DEFINING CONVERSION C++ <-> HDF5 DATATYPES
    // =======================================================================

    //! Struct for defining mapping of atom datatypes
    template<typename TYPE>
    struct HdfAtomTypeMap {

      //! Native HDF5 datatype of template parameter
      static H5::PredType HdfNativeType;

      //! Standard HDF5 datatype of template parameter
      static H5::PredType HdfStdType;
    };

    //! Base class for performing conversion
    class BaseHdfTypeConversion {
    public:

      //! Constructor
      BaseHdfTypeConversion() :
        isSet_(false),
        nativeType_(NULL),
        stdType_(NULL),
        size_( 0 ),
        numElems_( 0 )
      {
        if(atomTypeMap_.empty()) {
          InitAtomTypeMap();
        }
      }

      //! Destructor
      virtual ~BaseHdfTypeConversion() {
        CleanUp();
      };

      void AdaptToTypeBeforeRead(const H5::DataType& data) {}
      
      //! Query, if data is set
      bool IsSet() { return isSet_; }

      //! Get platform dependent HDF5 datatype
      H5::DataType* GetNativeType() { return nativeType_; }

      //! Get platform independent HDF5 datatype
      H5::DataType* GetStdType() { return stdType_; }

      //! Get raw pointer to data to be written to hdf5 file
      virtual const void * GetOutBufferPtr() = 0;

      //! Obtain data pointer to internal buffer for
      //! converting standard hdf5 data into native one
      virtual void* GetInBufferPtr( unsigned int numData ) = 0;

      //! Return size of raw data in bytes (used for compounds)
      unsigned int GetRawSize() { return size_; }

      //! Return number of data elements in the buffer
      unsigned int GetNumElems() { return numElems_; }

      //! Clean up conversion method
      virtual void CleanUp() {}

    protected:

      //! Flag indicating if data is set
      bool isSet_;

      //! Native HDF5 DataType
      H5::DataType *nativeType_;

      //! Standard HDF5 DataType
      H5::DataType *stdType_;

      //! Size of the array
      unsigned int size_;

      //! Numer of elements in the array
      unsigned int numElems_;

      static std::map< std::string, std::pair<const H5::PredType*, const H5::PredType*> > atomTypeMap_;

      void InitAtomTypeMap();
    };

    //! Templatized class for type conversion
    template<typename TYPE>
    class HdfTypeConversion : public BaseHdfTypeConversion {
    public:
      HdfTypeConversion() {
        H5CFS_EXCEPTION( "Type conversion not implemented for type '"
                   << typeid(TYPE).name() );
      }
    };


    //! Get conversion object for given boost::any object
    static void GetAnyConversion( const boost::any& anyType,
                                  shared_ptr<BaseHdfTypeConversion>& conv );
  }; // end of class H5IO

} // end of namespace CoupledField

#endif // FILE_CFS_HDF5_IO_HH
