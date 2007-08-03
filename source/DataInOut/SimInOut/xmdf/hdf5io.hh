// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_HDF5_IO_HH
#define FILE_CFS_HDF5_IO_HH

#include <set>
#include <map>


#include <General/environment.hh>
#include <boost/any.hpp>

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
    
    //! Commodity function for writing a 2D dataset
    template<typename TYPE>
    static void Write2DArray( H5::CommonFG &loc,
                              const std::string& name,
                              UInt rowSize,
                              UInt colSize,
                              const TYPE * buffer,
                              const H5::DSetCreatPropList &create_plist
                              = H5::DSetCreatPropList::DEFAULT );
    
    //! Commoditiy function for writing a scalar compound (1x1 rank)
    static void WriteCompound( H5::CommonFG& loc,
                               const std::string& name,
                               const CompoundType comp,
                               const H5::DSetCreatPropList &create_plist
                               = H5::DSetCreatPropList::DEFAULT );

    // =======================================================================
    //  READ METHODS
    // =======================================================================

    // ... to be implemented ...

    // =======================================================================
    //  CONVERSION METHODS
    // =======================================================================
    
    //! Map XMDF element type to CFS one
    static FEType XMDFElemType2ElemType( const Integer type );

    //! Map CFS element type to XMD one
    static Integer ElemType2XMDFElemType( const FEType type );

    //! Map connectivity CFS <-> XMDF
    static void ReorderConnectivity( const Integer eType,
                                     const bool toXMDF,
                                     const UInt* in,
                                     UInt* out);

  private:

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
        size_( 0 )
      {}
      
      //! Destructor
      virtual ~BaseHdfTypeConversion() {};
      
      //! Query, if data is set
      bool IsSet() { return isSet_; }
      
      //! Get platform ependent HDF5 datatype
      H5::DataType GetNativeType() { return nativeType_; }
      
      //! Get platform independent HDF5 datatype
      H5::DataType GetStdType() { return stdType_; }

      //! Get raw pointer to given object
      virtual const void * GetVoidPtr() = 0;

      //! Return size of raw data in bytes 
      UInt GetSize() { return size_; }

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
