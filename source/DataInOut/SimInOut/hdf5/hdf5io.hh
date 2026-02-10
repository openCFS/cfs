// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_HDF5_IO_HH
#define FILE_CFS_HDF5_IO_HH

#include <set>
#include <map>
#include <any>

#include "General/Environment.hh"
#include "Domain/Results/ResultInfo.hh"
#include "DataInOut/SimOutput.hh"

#include "H5Cpp.h"

namespace CoupledField {

  //! Helper class for providing convenience access to hdf5 data files
  class H5IO {

  public:

    // =======================================================================
    //  TYPE DEFINITIONS
    // =======================================================================
    typedef StdVector<std::pair<std::string, std::any> >
    CompoundType;

    typedef StdVector<std::pair<std::string, StdVector<std::any> >  >
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
                              UInt size,
                              const TYPE * buffer,
                              const H5::DSetCreatPropList &create_plist
                              = H5::DSetCreatPropList::DEFAULT );

    //! Write a 1D array, which can be extended afterwards 
    template<typename TYPE>
    static void Write1DArrayExt( H5::CommonFG &loc,
                                 const std::string& name,
                                 UInt size,
                                 const TYPE * buffer,
                                 const H5::DSetCreatPropList &create_plist
                                 = H5::DSetCreatPropList::DEFAULT );
    
    //! Append entries to extendible array. If not present, it is created.
    template<typename TYPE>
    static void Extend1DArray( H5::CommonFG &loc,
                               const std::string& name,
                               UInt newSize,
                               const TYPE * buffer,
                               const H5::DSetCreatPropList &create_plist
                               = H5::DSetCreatPropList::DEFAULT );
    
    //! Reserve space for a 1D array, but do not write anything to it
    template<typename TYPE>
    static void Reserve1DArray( H5::CommonFG &loc,
                                const std::string& name,
                                UInt size,
                                const H5::DSetCreatPropList &create_plist
                                = H5::DSetCreatPropList::DEFAULT );

    //! Set entries entries in 1D array
    template<typename TYPE>
    static void SetEntries1DArray( H5::CommonFG &loc,
                                   const std::string& name,
                                   UInt begin, UInt end,
                                   const TYPE * buffer  );

    //! Commodity function for writing a 2D dataset
    template<typename TYPE>
    static void Write2DArray( H5::CommonFG &loc,
                              const std::string& name,
                              UInt rowSize,
                              UInt colSize,
                              const TYPE * buffer,
                              const H5::DSetCreatPropList &create_plist
                              = H5::DSetCreatPropList::DEFAULT );


    //! Write a 2D array, which can be extended afterwards 
    template<typename TYPE>
    static void Write2DArrayExt( H5::CommonFG &loc,
                                 const std::string& name,
                                 UInt rowSize,
                                 UInt colSize,
                                 const TYPE * buffer,
                                 const H5::DSetCreatPropList &create_plist
                                 = H5::DSetCreatPropList::DEFAULT );

    //! Append entries to extendible array. If not present, it is created.
    template<typename TYPE>
    static void Extend2DArray( H5::CommonFG &loc,
                               const std::string& name,
                               UInt newRowSize,
                               UInt newColSize,
                               const TYPE * buffer,
                               const H5::DSetCreatPropList &create_plist
                               = H5::DSetCreatPropList::DEFAULT );

    //! Reserve space for a 2D array, but do not write anything to it
    template<typename TYPE>
    static void Reserve2DArray( H5::CommonFG &loc,
                                const std::string& name,
                                UInt rowSize,
                                UInt colSize,
                                const H5::DSetCreatPropList &create_plist
                                = H5::DSetCreatPropList::DEFAULT );

    //! Commodity function for writing a 2D dataset
    template<typename TYPE>
    static void SetEntries2DArray( H5::CommonFG &loc,
                                   const std::string& name,
                                   UInt rowBegin, UInt rowEnd,
                                   UInt colBegin, UInt colEnd,
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

    //! This method retrieves the name of an object of a group with given
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
    static StdVector<UInt> GetArrayDims( const H5::CommonFG &loc,
                                           const std::string& name );

    //! Return number of entries of a dataset / array
    static UInt GetNumEntries( const H5::CommonFG &loc,
                               const std::string& name );

    //! Retrieve array data from a dataset

    //! Read data from a an dataset of arbitrary dimension into a linear buffer
    //! Note, that the memory has to be allocated from outside
    template<typename TYPE>
    static void ReadArray( H5::CommonFG &loc,
                           const std::string& name,
                           TYPE* data );

    //! Read data from a dataset into a openCFS vector
    template<typename TYPE>
    static void ReadArray( H5::CommonFG &loc,
                           const std::string& name,
                           StdVector<TYPE>& data );

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
                                        UInt msStep,
                                        bool isHistory );

    //! Obtain grid result group for specified step in a given multistep
    static H5::Group GetStepGroup( H5::H5File& file,
                                   UInt msStep, UInt
                                   stepNum );

    //! Open group and create it if not yet present
    static H5::Group OpenCreateGroup(H5::CommonFG &curGroup,
                                     const std::string& name );

    //! Check if group exists
    static bool GroupExists( const H5::CommonFG &loc,
                             const std::string& name );

    //! Check if dataset exists
    static bool DatasetExists( const H5::CommonFG &loc,
                               const std::string& name );

    //! Check if attribute exists
    static bool AttrExists( const H5::H5Object& obj,
                            const std::string& name );


    // =======================================================================
    //  CONVERSION METHODS
    // =======================================================================

    //! Map SimOuput::Capability class to hdf5 type
    static Integer MapCapabilityType( SimOutput::Capability c );

    //! Map hdf5 representation of simOutput::Capability to enum representation
    static SimOutput::Capability MapCapabilityType( Integer c );

    //! Map EntityUnknownType enum to hdf5 type
    static Integer MapUnknownType( ResultInfo::EntityUnknownType t );

    //! Map EntityUnknown hdf5 type to enum
    static ResultInfo::EntityUnknownType MapUnknownType( Integer t );

    //! Map EntityUnknownType enum to string representation
    static std::string MapUnknownTypeAsString( ResultInfo::EntityUnknownType t );

    //! Map entryType from enum to hdf5 type
    static Integer MapEntryType( ResultInfo::EntryType t );

    //! Map entryType from hdf5 type to enum
    static ResultInfo::EntryType MapEntryType( Integer t );

    // =======================================================================
    //  MISCELANEOUS METHODS
    // =======================================================================

    //! Set minimum chunksize in case of extensible dataset
    static void SetMinChunkSize( UInt chunkSize );

    //! Set maximum chunksize to be used for Array data
    static void SetMaxChunkSize( UInt chunkSize );
    

    //! Check for open objects in the hdf5 file and close them
    static void CheckOpenObjects(H5::H5File& file, bool verbose);
    
  private:

    // =======================================================================
    //  CHUNKSIZE BOUNDARIES
    // =======================================================================
    //! Minimum chunk size (only used for extensible data arrays)
    static hsize_t minChunkSize_;
    
    //! Maximum chunk size
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

      void AdaptToTypeBeforeRead(const H5::DataType& data){ }

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
      virtual void* GetInBufferPtr( UInt numData ) = 0;

      //! Return size of raw data in bytes (used for compounds)
      UInt GetRawSize() { return size_; }

      //! Return number of data elements in the buffer
      UInt GetNumElems() { return numElems_; }

      //! Clean up conversion method
      virtual void CleanUp() {
        delete nativeType_; nativeType_ = NULL;
        delete stdType_; stdType_ = NULL;
      }

    protected:

      //! Flag indicating if data is set
      bool isSet_;

      //! Native HDF5 DataType
      H5::DataType *nativeType_;

      //! Standard HDF5 DataType
      H5::DataType *stdType_;

      //! Size of the array
      UInt size_;

      //! Number of elements in the array
      UInt numElems_;

      //! Map data type name (key, string representation) to (native,std)-datatype
      static std::map< std::string, std::pair<const H5::PredType*, const H5::PredType*> > atomTypeMap_;

      void InitAtomTypeMap();
    };

    //! Templatized class for type conversion
    template<typename TYPE>
    class HdfTypeConversion : public BaseHdfTypeConversion {
    public:
      HdfTypeConversion() {
        EXCEPTION( "Type conversion not implemented for type '"
                   << typeid(TYPE).name() );
      }
    };


    //! Get conversion object for given std::any object
    static void GetAnyConversion( const std::any& anyType,
                                  shared_ptr<BaseHdfTypeConversion>& conv );
  }; // end of class H5IO

} // end of namespace CoupledField

#endif // FILE_CFS_HDF5_IO_HH
