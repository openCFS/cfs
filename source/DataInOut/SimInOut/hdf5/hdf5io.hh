// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_HDF5_IO_HH
#define FILE_CFS_HDF5_IO_HH

#include <set>
#include <map>
#include <boost/any.hpp>

#include "General/environment.hh"
#include "Domain/resultInfo.hh"
#include "DataInOut/simOutput.hh"


#include "H5Cpp.h"

namespace CoupledField {

  //! Helper class for providing conveniecne access to hdf5 data files
  class H5IO {

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
                                const H5::DSetCreatPropList &create_plist
                                = H5::DSetCreatPropList::DEFAULT );
    
    //! Commodity function for writing a 1D dataset
    template<typename TYPE>
    static void Write1DArray( H5::CommonFG &loc,
                              const std::string& name,
                              UInt size,
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
    static std::vector<UInt> GetArrayDims( const H5::CommonFG &loc,
                                           const std::string& name );

    //! Return number of entries of a dataset / rray
    static UInt GetNumEntries( const H5::CommonFG &loc,
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
                                            UInt msStep, 
                                            bool isHistory );

    //! Obtain grid result group for specified step in a given multistep
    static H5::Group GetStepGroup( H5::H5File& file, 
                                   UInt msStep, UInt 
                                   stepNum );
    
    // =======================================================================
    //  CONVERSION METHODS
    // =======================================================================

    //! Map SimOuput::Capability class to hdf5 type
    static Integer MapCapabilityType( SimOutput::Capability c );
    
    //! Map hdf5 representation of simOutput::Capability to enum representat
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
    
    //! Set chunksize to be used for Array data
    static void SetMaxChunkSize( UInt chunkSize );
    
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
      static H5::DataType HdfNativeType;
      
      //! Standard HDF5 datatype of template parameter
      static H5::DataType HdfStdType;
    };
    
    //! Base class for performing conversion
    class BaseHdfTypeConversion {
    public:
      
      //! Constructor
      BaseHdfTypeConversion() : 
        isSet_(false),
        size_( 0 ),
        numElems_( 0 )
      {}
      
      //! Destructor
      virtual ~BaseHdfTypeConversion() {
        CleanUp();
      };
      
      //! Query, if data is set
      bool IsSet() { return isSet_; }
      
      //! Get platform dependent HDF5 datatype
      H5::DataType GetNativeType() { return nativeType_; }
      
      //! Get platform independent HDF5 datatype
      H5::DataType GetStdType() { return stdType_; }

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
      virtual void CleanUp() {}
     
    protected:
      
      //! Flag indicating if data is set
      bool isSet_;

      //! Native HDF5 DataType
      H5::DataType nativeType_;

      //! Standard HDF5 DataType
      H5::DataType stdType_;

      //! Size of the array
      UInt size_;

      //! Numer of elements in the array
      UInt numElems_;

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
    

    //! Get conversion object for given boost::any object
    static void GetAnyConversion( const boost::any& anyType,
                                  shared_ptr<BaseHdfTypeConversion>& conv );
  }; // end of class H5IO
  
} // end of namespace CoupledField

#endif // FILE_CFS_HDF5_IO_HH
