// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "hdf5Common.h"
#include <boost/lexical_cast.hpp>
#include <iostream>


namespace H5CFS {

// define commodity method for converting a hdf5 exception
// to a CFS one
#define H5_CATCH( STR )                                                 \
  catch (H5::Exception& h5Ex ) {                                        \
    H5CFS_EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );                  \
  }


// =================================
//    Initialize Static Variables
// =================================
hsize_t Hdf5Common::maxChunkSize_= 100; 

  std::map< std::string, std::pair<const H5::PredType*, const H5::PredType*> > Hdf5Common::BaseHdfTypeConversion::atomTypeMap_;

  // ====================================
  //    Initialize Atom Data Typemaps
  // ====================================
  void Hdf5Common::BaseHdfTypeConversion::InitAtomTypeMap() 
  {
    atomTypeMap_["bool"].first = &H5::PredType::NATIVE_INT32;
    atomTypeMap_["bool"].second = &H5::PredType::STD_I32LE;

    atomTypeMap_["Integer"].first = &H5::PredType::NATIVE_INT32;
    atomTypeMap_["Integer"].second = &H5::PredType::STD_I32LE;
    atomTypeMap_["int"].first = &H5::PredType::NATIVE_INT32;
    atomTypeMap_["int"].second = &H5::PredType::STD_I32LE;

    atomTypeMap_["UInt"].first = &H5::PredType::NATIVE_UINT32;
    atomTypeMap_["UInt"].second = &H5::PredType::STD_U32LE;
    atomTypeMap_["unsigned int"].first = &H5::PredType::NATIVE_UINT32;
    atomTypeMap_["unsigned int"].second = &H5::PredType::STD_U32LE;

    atomTypeMap_["Double"].first = &H5::PredType::NATIVE_DOUBLE;
    atomTypeMap_["Double"].second = &H5::PredType::IEEE_F64LE;
    atomTypeMap_["double"].first = &H5::PredType::NATIVE_DOUBLE;
    atomTypeMap_["double"].second = &H5::PredType::IEEE_F64LE;

    atomTypeMap_["Float"].first = &H5::PredType::NATIVE_FLOAT;
    atomTypeMap_["Float"].second = &H5::PredType::IEEE_F32LE;
    atomTypeMap_["float"].first = &H5::PredType::NATIVE_FLOAT;
    atomTypeMap_["float"].second = &H5::PredType::IEEE_F32LE;
  }

  // ========================================
  //    Initialize General Datatype mapping
  // ========================================

  // -----------
  //  AtomTypes
  // -----------

#define DECL_HDF_ATOM_TYPE_CONV(TYPE)                           \
  template<>                                                    \
  class Hdf5Common::HdfTypeConversion<TYPE> :                   \
    public Hdf5Common::BaseHdfTypeConversion {                  \
                                                                \
  private:                                                      \
                                                                \
    TYPE* buffer_;                                              \
                                                                \
  public:                                                       \
    HdfTypeConversion()                                         \
      : buffer_( NULL ) {                                       \
                                                                \
      nativeType_ = new H5::PredType(*atomTypeMap_[std::string(#TYPE)].first);  \
      stdType_    = new H5::PredType(*atomTypeMap_[std::string(#TYPE)].second); \
    }                                                           \
                                                                \
    virtual ~HdfTypeConversion() {                              \
      CleanUp();                                                \
    };                                                          \
                                                                \
    const void * GetOutBufferPtr() {                            \
      if( !isSet_ ) {                                           \
        H5CFS_EXCEPTION( "Data buffer is empty" );              \
      }                                                         \
      return (const void*) buffer_;                             \
    }                                                           \
                                                                \
    void* GetInBufferPtr( unsigned int numData ) {                      \
      CleanUp();                                                \
      numElems_ = numData;                                      \
      buffer_ = new TYPE[numData];                              \
      return (void*) buffer_;                                   \
    }                                                           \
                                                                \
    void GetNativeData( TYPE * data) {                          \
      for( unsigned int i = 0; i < numElems_; i++ ) {                   \
        data[i] = buffer_[i];                                   \
      }                                                         \
    }                                                           \
                                                                \
    void SetNativeData( const TYPE& t ) {                       \
      CleanUp();                                                \
      buffer_ = new TYPE[1];                                    \
      buffer_[0] = t;                                           \
      size_ = sizeof( TYPE );                                   \
      numElems_ = 1;                                            \
      isSet_ = true;                                            \
    }                                                           \
                                                                \
    void SetNativeData( const TYPE* t, unsigned int size  ) {           \
      CleanUp();                                                \
      buffer_ = new TYPE[size];                                 \
      for( unsigned int i = 0; i < size; i ++ ) {                       \
        buffer_[i] = t[i];                                      \
      }                                                         \
      size_ = sizeof( TYPE* );                                  \
      numElems_ = size;                                         \
      isSet_ = true;                                            \
    }                                                           \
                                                                \
    void CleanUp() {                                            \
      if( buffer_ ) {                                           \
        delete[] buffer_;                                       \
      }                                                         \
      buffer_ = NULL;                                           \
      size_ = 0;                                                \
      numElems_ = 0;                                            \
      isSet_ = false;                                           \
    }                                                           \
  }

  DECL_HDF_ATOM_TYPE_CONV(int);
  DECL_HDF_ATOM_TYPE_CONV(unsigned int);
  DECL_HDF_ATOM_TYPE_CONV(double);
  DECL_HDF_ATOM_TYPE_CONV(float);

#undef DECL_HDF_ATOM_TYPE_CONV

  // ---------------
  //  bool
  // ---------------
  template<>
  class Hdf5Common::HdfTypeConversion<bool> :
    public Hdf5Common::BaseHdfTypeConversion  {

  private:

    int* buffer_;

  public:

    HdfTypeConversion()
      : buffer_ ( NULL ) {

      nativeType_ = new H5::PredType(*atomTypeMap_["bool"].first);
      stdType_    = new H5::PredType(*atomTypeMap_["bool"].second);
    }

    const void * GetOutBufferPtr() {
      if( !isSet_ ) {
        H5CFS_EXCEPTION( "Data buffer is empty" );
      }
      return buffer_;
    }

    void* GetInBufferPtr( unsigned int numData ) {
      CleanUp();
      numElems_ = numData;
      buffer_ = new int[numElems_];
      return buffer_;
    }

    void GetNativeData( bool * data) {
      for( unsigned int i = 0; i < numElems_; i++ ) {
        data[i] = (buffer_[i] == 1) ? true : false;
      }
    }

    void SetNativeData( const bool& t ) {
      CleanUp();
      buffer_ = new int[1];
      buffer_[0] = t ? 1 : 0;
      size_ = sizeof( int );
      numElems_ = 1;
      isSet_ = true;
    }

    void SetNativeData( const bool* t, unsigned int size  ) {
      CleanUp();
      buffer_ = new int[size];
      for( unsigned int i = 0; i < size; i ++ ) {
        buffer_[i] = t[i] ? 1 : 0;
      }
      size_ = sizeof( int* );
      numElems_ = size;
      isSet_ = true;
    }

    void CleanUp() {
      if( buffer_ ) {
        delete[] buffer_;
      }
      buffer_ = NULL;
      size_ = 0;
      numElems_ = 0;
      isSet_ = false;
    }

  };

  // ---------------
  //  std::string
  // ---------------
  template<>
  class Hdf5Common::HdfTypeConversion<std::string> :
    public Hdf5Common::BaseHdfTypeConversion  {

  private:

    const char** buffer_;

    bool isFixedLength;

    unsigned int fixedLength;

  public:
    HdfTypeConversion()
      : buffer_( NULL ),
        isFixedLength (false),
        fixedLength(0) {
      nativeType_ = new H5::StrType( H5::PredType::C_S1, H5T_VARIABLE);
      stdType_ = new H5::StrType( H5::PredType::C_S1, H5T_VARIABLE);
    }

    void AdaptToTypeBeforeRead(const H5::DataType& data) {
      //only do special things if we have not variable length
      if(!data.isVariableStr()){
        //we need to check for fixed or variable length
        delete nativeType_;
        delete stdType_;
        nativeType_ = NULL;
        stdType_ = NULL;

        unsigned int dim = H5Tget_size(data.getId());
        nativeType_ = new H5::StrType( H5::PredType::C_S1, dim+1);
        stdType_    = new H5::StrType( H5::PredType::C_S1, dim+1);

        isFixedLength = true;
        fixedLength = dim;
        if( buffer_ ) {
          delete[] buffer_;
          buffer_ = NULL;
        }

        if(numElems_) 
        {
          buffer_ = new const char*[numElems_];
          for(unsigned int i = 0;i < numElems_;i++){
            buffer_[i] = new char[fixedLength+1];
          }
        } else 
        {
          isSet_ = false;
        }
      }
    }

    const void * GetOutBufferPtr() {
      if( !isSet_ ) {
        H5CFS_EXCEPTION( "Data buffer is empty" );
      }
      return buffer_;
    }

    void * GetInBufferPtr( unsigned int numData ) {
      CleanUp();
      numElems_ = numData;
      buffer_ = new const char*[numElems_];
      if(isFixedLength){
        for(unsigned int i = 0;i < numElems_;i++){
          buffer_[i] = new char[fixedLength+1];
        }
      }
      //only adjust here was we do not want to write
      // in fixed length string...
      if(isFixedLength){
        return const_cast<char *>(*buffer_);
      } else {
        return buffer_;
      }
      return NULL;
    }

    void GetNativeData( std::string * data) {
      for( unsigned int i = 0; i < numElems_; i++ ) {
        data[i].assign( buffer_[i]);
      }
    }

    void SetNativeData( const std::string& t ) {
      CleanUp();
      buffer_ = new const char*[1];
      buffer_[0] = t.c_str();
      size_ = sizeof( char* );
      numElems_ = 1;
      isSet_ = true;
    }

    void SetNativeData( const std::string* t, unsigned int size  ) {
      CleanUp();
      buffer_ = new const char*[size];
      for( unsigned int i = 0; i < size; i ++ ) {
        buffer_[i] = t[i].c_str();
      }
      size_ = sizeof( char* );
      numElems_ = size;
      isSet_ = true;
    }

    void CleanUp() {
      if( buffer_ ) {
        delete[] buffer_;
      }
      buffer_ = NULL;
      size_ = 0;
      numElems_ = 0;
      isSet_ = false;
    }

    ~HdfTypeConversion() 
    {
      CleanUp();
      delete stdType_;
      delete nativeType_;
    }
    
  };


  // --------------------------
  //  std::vector<std::string>
  // --------------------------
  template<>
  class Hdf5Common::HdfTypeConversion< std::vector<std::string> > :
    public Hdf5Common::BaseHdfTypeConversion {

    private:
      hvl_t * buffer_;

  public:
    HdfTypeConversion()
      : buffer_( NULL ) {
      H5::StrType sType( H5::PredType::C_S1, H5T_VARIABLE);
      nativeType_ = new H5::VarLenType( &sType );
      stdType_ = new H5::VarLenType( &sType );
    }

    void AdaptToTypeBeforeRead(const H5::DataType& data){

      //we need to check for fixed or variable length
      //this should work for fixed string lengths but i am not sure
      //therefore we throw a warning
      if(!data.isVariableStr()){
        delete nativeType_;
        delete stdType_;
        nativeType_ = NULL;
        stdType_ = NULL;

        unsigned int dim = H5Tget_size(data.getId());

        nativeType_ = new H5::StrType( H5::PredType::C_S1, dim+1);
        stdType_    = new H5::StrType( H5::PredType::C_S1, dim+1);
        
        // WARN("HDF5 File input: Converting vectors of fixed length strings is experimental!");
      }
    }


      const void * GetOutBufferPtr() {
        if( !isSet_ ) {
          H5CFS_EXCEPTION( "Data buffer is empty" );
        }
        return buffer_;
      }

      void* GetInBufferPtr( unsigned int numData ) {
        CleanUp();
        numElems_ = numData;
        buffer_ = new hvl_t[numData];
        return buffer_;
      }

      void GetNativeData( std::vector<std::string> * data) {
        for( unsigned int i = 0; i < numElems_; i++ ) {
          data[i].resize( buffer_[i].len );
          for( unsigned int j = 0; j < buffer_[i].len; j++ ) {
            data[i][j].assign( ((const char**)buffer_[i].p)[j] );
          }
        }
      }

      void SetNativeData( const std::vector<std::string>& t ) {
        CleanUp();
        buffer_ = new hvl_t[1];
        buffer_[0].p = (void*) new const char*[t.size()];
        buffer_[0].len = t.size();
        for( unsigned int j = 0; j < t.size(); j++ ) {
          ((const char **)buffer_[0].p)[j] = t[j].c_str();
        }
        size_ = sizeof( hvl_t );
        numElems_ = 1;
        isSet_ = true;
      }


      void SetNativeData( const std::vector<std::string>* t, unsigned int size ) {
        CleanUp();
        buffer_ = new hvl_t[size];
        for( unsigned int i = 0; i < size; i++ ) {
          buffer_[i].p = (void*) new const char*[t[i].size()];
          buffer_[i].len = t[i].size();
          for( unsigned int j = 0; j < t[i].size(); j++ ) {
            ((const char **)buffer_[i].p)[j] = ((t[i])[j]).c_str();
          }
        }
        size_ = sizeof( hvl_t );
        numElems_ = size;
        isSet_ = true;
      }

      void CleanUp() {

        if( buffer_ ) {

          // delete pointers to characters
          for( unsigned int i = 0; i < numElems_; i++ ) {
            delete[] (const char*) buffer_[i].p;
            buffer_[i].p = NULL;
          }

          // delete buffer itself
          delete[] buffer_;
        }
        buffer_ = NULL;
        size_ = 0;
        numElems_ = 0;
        isSet_ = false;
      }

    ~HdfTypeConversion() 
    {
      CleanUp();
      delete stdType_;
      delete nativeType_;
    }

  }; // end of class definition



  // ----------------------
  //  std::vector<TYPE>
  // ----------------------
#define DECL_STL_VECTOR_CONVERSION( TYPE )                              \
  template<>                                                            \
  class Hdf5Common::HdfTypeConversion<std::vector<TYPE> > :             \
    public Hdf5Common::BaseHdfTypeConversion {                          \
  private:                                                              \
    hvl_t * buffer_;                                                    \
                                                                        \
  public:                                                               \
    HdfTypeConversion()                                                 \
      : buffer_( NULL ) {                                               \
      nativeType_ = new H5::VarLenType (atomTypeMap_[std::string(#TYPE)].first); \
      stdType_    = new H5::VarLenType (atomTypeMap_[std::string(#TYPE)].second); \
    }                                                                   \
                                                                        \
    const void * GetOutBufferPtr( ) {                                   \
      if( !isSet_ ) {                                                   \
        H5CFS_EXCEPTION( "Data buffer is empty" );                      \
      }                                                                 \
      return (const void*) buffer_;                                     \
    }                                                                   \
                                                                        \
    void* GetInBufferPtr( unsigned int numData ) {                      \
      CleanUp();                                                        \
      numElems_ = numData;                                              \
      buffer_ = new hvl_t[numData];                                     \
      return buffer_;                                                   \
    }                                                                   \
                                                                        \
    void GetNativeData( std::vector<TYPE> * data) {                     \
      for( unsigned int i = 0; i < numElems_; i++ ) {                   \
        data[i].resize( buffer_[i].len );                               \
        for( unsigned int j = 0; j < buffer_[i].len; j++ ) {            \
          data[i][j] = ((TYPE*)buffer_[i].p)[j];                        \
        }                                                               \
      }                                                                 \
    }                                                                   \
                                                                        \
    void SetNativeData( const std::vector<TYPE>& t ) {                  \
      CleanUp();                                                        \
      buffer_ = new hvl_t[1];                                           \
      buffer_[0].p =  (void*) new TYPE[t.size()];                       \
      buffer_[0].len = t.size();                                        \
      for( unsigned int j = 0; j < t.size(); j++ ) {                    \
        ((TYPE*)buffer_[0].p)[j] = t[j];                                \
      }                                                                 \
      numElems_ = 1;                                                    \
      size_ = sizeof( hvl_t );                                          \
      isSet_ = true;                                                    \
    }                                                                   \
                                                                        \
    void SetNativeData( const std::vector<TYPE>* t, unsigned int size ) {         \
      CleanUp();                                                        \
      buffer_ = new hvl_t[size];                                        \
      for( unsigned int i = 0; i < size; i++ ) {                        \
        buffer_[i].p =  new TYPE[t[i].size()];                          \
        buffer_[i].len = t[i].size();                                   \
        for( unsigned int j = 0; j < t[i].size(); j++ ) {               \
          ((TYPE*)buffer_[i].p)[j] = t[i][j];                           \
        }                                                               \
      }                                                                 \
      numElems_ = size;                                                 \
      size_ = sizeof( hvl_t );                                          \
      isSet_ = true;                                                    \
    }                                                                   \
                                                                        \
  void CleanUp() {                                                      \
    if( buffer_ ) {                                                     \
      for( unsigned int i = 0; i < numElems_; i++ ) {                           \
        delete[] (TYPE*) buffer_[i].p;                                  \
        }                                                               \
        delete[] buffer_;                                               \
      }                                                                 \
                                                                        \
      buffer_ = 0;                                                      \
      size_ = 0;                                                        \
      numElems_ = 0;                                                    \
      isSet_ = false;                                                   \
                                                                        \
    }                                                                   \
                                                                        \
    ~HdfTypeConversion() {                                              \
      CleanUp();                                                        \
      delete nativeType_;                                               \
      delete stdType_;                                                  \
    }                                                                   \
                                                                        \
  }

  DECL_STL_VECTOR_CONVERSION( int );
  DECL_STL_VECTOR_CONVERSION( unsigned int );
  DECL_STL_VECTOR_CONVERSION( double );
  DECL_STL_VECTOR_CONVERSION( float );

#undef DECL_STL_VECTOR_CONVERSION

  template<typename TYPE>
  void Hdf5Common::WriteAttribute( H5::H5Object& obj,
                             const std::string& name,
                             const TYPE& data,
                             const H5::PropList &create_plist ) {

    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType* stdType = conv.GetStdType();
      H5::DataType* nativeType = conv.GetNativeType();

      // create memory data space
      H5::DataSpace space;

      // generate attribute
      conv.SetNativeData( data );
      if( !conv.IsSet() ) {
        H5CFS_EXCEPTION( "Could not convert data for attribute '"
                   << name << "' of type " << typeid(TYPE).name() );
      }
			
			// if an attribute already exists, we have to open it, instead of create
			// so we can now call WriteAttribute multiple times
      H5::Attribute attr;
      try{
        attr = obj.openAttribute( name );
      }catch(H5:: Exception&){
        attr = obj.createAttribute( name, *stdType,
            space, create_plist );
      }

      // write attribute
      attr.write( *nativeType, conv.GetOutBufferPtr() );

      // reset conversion object
      conv.CleanUp();

      // close attribute, dataspace- and types
      space.close();
      attr.close();


    }  catch (H5::Exception& h5ex) {
      H5CFS_EXCEPTION("Could not write attribute '" 
                << name << "':\n" << h5ex.getCDetailMsg());
    } catch( std::string& ex ) {
      H5CFS_RETHROW_EXCEPTION(ex, "Could not write attribute '" << name << "'" );
    }
  }

  template<typename TYPE>
  void Hdf5Common::Write1DArray( H5::CommonFG &loc,
                           const std::string& name,
                           unsigned int size,
                           const TYPE * buffer,
                           const H5::DSetCreatPropList &create_plist ) {

    // check, that size is greate than zero
    if( size == 0 || buffer == NULL ) {
      H5CFS_EXCEPTION( "Attribute data buffer of 1D array '" << name 
                 << "' is NULL or has zero size" );
    }
    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType* stdType = conv.GetStdType();
      H5::DataType* nativeType = conv.GetNativeType();

      // create memory data space
      const hsize_t dims = size;
      H5::DataSpace space( 1, &dims );

      // generate dataset and fill it
      conv.SetNativeData( buffer, size );
      if( !conv.IsSet() ) {
        H5CFS_EXCEPTION( "Could not convert data for 1D array '"
                   << name << "' of type " << typeid(TYPE).name() );
      }

      // set chunking of dataset
      H5::DSetCreatPropList newList(create_plist);
      const hsize_t chunk = std::min( (unsigned int) size, (unsigned int) maxChunkSize_ );
      newList.setChunk( 1, &chunk);
      H5::DataSet dataset = loc.createDataSet( name, *stdType,
                                               space, newList );
      dataset.write( conv.GetOutBufferPtr(), *nativeType  );

      // reset conversion object
      conv.CleanUp();

      // close dataset, dataspace- and types
      space.close();
      dataset.close();


    } catch (H5::Exception& h5ex) {
      H5CFS_EXCEPTION("Could not write 1D-Array '" 
                << name << "':\n" << h5ex.getCDetailMsg());
    } catch( std::string& ex ) {
      H5CFS_RETHROW_EXCEPTION(ex, "Could not write 1D-Array '" << name << "'" );
    }

  }

  template<typename TYPE>
  void Hdf5Common::Reserve1DArray( H5::CommonFG &loc,
                             const std::string& name,
                             unsigned int size,
                             const H5::DSetCreatPropList &create_plist ) {

    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType* stdType = conv.GetStdType();

      // create memory data space
      hsize_t dims[1] = {size};
      const hsize_t maxDims[1] = {H5S_UNLIMITED};
      H5::DataSpace space( 1, dims, maxDims );

      // set chunking of dataset
      H5::DSetCreatPropList newList(create_plist);

      hsize_t chunk = std::min( (unsigned int) size, (unsigned int) maxChunkSize_ );
      newList.setChunk( 1, &chunk);
      H5::DataSet dataset = loc.createDataSet( name, *stdType,
                                               space, newList );

      // reset conversion object
      conv.CleanUp();

      // close dataset, dataspace- and types
      space.close();
      dataset.close();

    } catch (H5::Exception& h5ex) {
      H5CFS_EXCEPTION("Could not reserve 1D-Array '" 
                << name << "':\n" << h5ex.getCDetailMsg());
    } catch( std::string& ex ) {
      H5CFS_RETHROW_EXCEPTION(ex, "Could not reserve 1D-Array '" << name << "'" );
    }
  }


  template<typename TYPE>
  void Hdf5Common::SetEntries1DArray( H5::CommonFG &loc,
                                const std::string& name,
                                unsigned int start, unsigned int end,
                                const TYPE * buffer ) {

    // check, that size is greate than zero
    if( start > end || buffer == NULL ) {
      H5CFS_EXCEPTION( "Attribute data buffer of 1D array '" << name 
                 << "' is NULL or has zero size" );
    }
    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType* nativeType = conv.GetNativeType();

      // create memory data space
      const hsize_t size = ( end - start ) + 1;
      const hsize_t maxDims = H5S_UNLIMITED;
      H5::DataSpace memSpace( 1, &size, &maxDims );

      // fill conversion object
      conv.SetNativeData( buffer, static_cast<unsigned int>(size) );
      if( !conv.IsSet() ) {
        H5CFS_EXCEPTION( "Could not convert data for 1D array '"
                   << name << "' of type " << typeid(TYPE).name() );
      }
      // open dataset and dataspace of file dataset
      H5::DataSet dataset = loc.openDataSet( name );
      H5::DataSpace fileSpace = dataset.getSpace();

      hsize_t Offset[1] = {start};
      hsize_t Mysize[1] = {size};
      fileSpace.selectHyperslab(  H5S_SELECT_SET, Mysize, Offset );

      // write data
      dataset.write( conv.GetOutBufferPtr(), *nativeType, memSpace, fileSpace );

      // reset conversion object
      conv.CleanUp();

      // close dataset, dataspace- and types
      fileSpace.close();
      memSpace.close();
      dataset.close();

    } catch (H5::Exception& h5ex) {
      H5CFS_EXCEPTION("Could not set entries in 1D-Array '" 
                << name << "':\n" << h5ex.getCDetailMsg());
    } catch( std::string& ex ) {
      H5CFS_RETHROW_EXCEPTION(ex, "Could not set entries in 1D-Array '" 
                        << name << "'" );
    }
  }

  template<typename TYPE>
  void Hdf5Common::Write2DArray( H5::CommonFG &loc,
                           const std::string& name,
                           unsigned int rowSize,
                           unsigned int colSize,
                           const TYPE * buffer,
                           const H5::DSetCreatPropList &create_plist
                           ) {

    // check, that size is greate than zero
    if( rowSize == 0 || colSize == 0 || buffer == NULL ) {
      H5CFS_EXCEPTION( "Data buffer of 2D array '" << name 
                 << "' is NULL or has zero size" );
    }

    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType* stdType = conv.GetStdType();
      H5::DataType* nativeType = conv.GetNativeType();

      // create memory data space
      const hsize_t dims[] = {rowSize, colSize};
      const int rank = 2;

      H5::DataSpace space( rank, dims );

      // generate dataset and fill it
      conv.SetNativeData( buffer, rowSize * colSize );
      if( !conv.IsSet() ) {
        H5CFS_EXCEPTION( "Could not convert data for 2D array '"
                   << name << "' of type " << typeid(TYPE).name() );
      }

      H5::DSetCreatPropList newList(create_plist);
      const hsize_t chunk[2] = { std::min( (unsigned int) rowSize, 
                                           (unsigned int) maxChunkSize_ ),
                                 std::min( (unsigned int) colSize, 
                                           (unsigned int) maxChunkSize_ ) };
      newList.setChunk( 2, chunk);
      H5::DataSet dataset = loc.createDataSet( name, *stdType,
                                               space, newList );
      dataset.write( conv.GetOutBufferPtr(), *nativeType  );

      // reset conversion object
      conv.CleanUp();

      // close dataset, dataspace- and types
      space.close();
      dataset.close();

    } catch (H5::Exception& h5ex) {
      H5CFS_EXCEPTION("Could not write 2D-Array '" << name 
                << "':\n" << h5ex.getCDetailMsg());
    } catch( std::string& ex ) {
      H5CFS_RETHROW_EXCEPTION(ex, "Could not write 2D-Array '" << name << "'" );
    }
  }

  template<typename TYPE>
  void Hdf5Common::Reserve2DArray( H5::CommonFG &loc,
                             const std::string& name,
                             unsigned int rowSize,
                             unsigned int colSize,
                             const H5::DSetCreatPropList &create_plist ) {
    try {

       // create conversion helper object and get native / std hdf5 datatype
       HdfTypeConversion<TYPE> conv;
       H5::DataType* stdType = conv.GetStdType();

       // create memory data space
       const hsize_t dims[] = {rowSize, colSize};
       const hsize_t maxDims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
       const int rank = 2;
       H5::DataSpace space( rank, dims, maxDims );

       // set chunking of dataset
       H5::DSetCreatPropList newList(create_plist);
       const hsize_t chunk[2] = { std::min( (unsigned int) rowSize, 
                                            (unsigned int) maxChunkSize_ ),
                                  std::min( (unsigned int) colSize, 
                                            (unsigned int) maxChunkSize_) };
       newList.setChunk( 2, chunk);
       H5::DataSet dataset = loc.createDataSet( name, *stdType,
                                                space, newList );

       // reset conversion object
       conv.CleanUp();

       // close dataset, dataspace- and types
       space.close();
       dataset.close();

     } catch (H5::Exception& h5ex) {
       H5CFS_EXCEPTION("Could not write 2D-Array '" << name
                 << "':\n" << h5ex.getCDetailMsg());
     } catch( std::string& ex ) {
       H5CFS_RETHROW_EXCEPTION(ex, "Could not write 2D-Array '" << name << "'" );
     }

  }

  template<typename TYPE>
  void Hdf5Common::SetEntries2DArray( H5::CommonFG &loc,
                                const std::string& name,
                                unsigned int rowBegin, unsigned int rowEnd,
                                unsigned int colBegin, unsigned int colEnd,
                                const TYPE * buffer ) {
    // check, that size is greate than zero
      if( rowBegin > rowEnd || colBegin > colEnd || buffer == NULL ) {
        H5CFS_EXCEPTION( "Data buffer of 2D array '" << name
                   << "' is NULL or has zero size" );
      }
      try {

        // create conversion helper object and get native / std hdf5 datatype
        HdfTypeConversion<TYPE> conv;
        H5::DataType* nativeType = conv.GetNativeType();

        // create memory data space
        const hsize_t size[2] = { (( rowEnd - rowBegin ) + 1),
                                  (( colEnd - colBegin ) + 1)};
        const hsize_t maxDims[2] = {H5S_UNLIMITED, H5S_UNLIMITED };
        H5::DataSpace memSpace( 2, size, maxDims );


        // fill conversion object
        conv.SetNativeData( buffer, static_cast<unsigned int>(size[0]*size[1]) );
        if( !conv.IsSet() ) {
          H5CFS_EXCEPTION( "Could not convert data for 1D array '"
                     << name << "' of type " << typeid(TYPE).name() );
        }
        // open dataset and dataspace of file dataset
        H5::DataSet dataset = loc.openDataSet( name );
        H5::DataSpace fileSpace = dataset.getSpace();

        hsize_t Offset[2] = {rowBegin, colBegin};
        hsize_t Mysize[2] = {size[0], size[1] };
        fileSpace.selectHyperslab(  H5S_SELECT_SET, Mysize, Offset );

        // write data
        dataset.write( conv.GetOutBufferPtr(), *nativeType, memSpace, fileSpace );

        // reset conversion object
        conv.CleanUp();

        // close dataset, dataspace- and types
        fileSpace.close();
        memSpace.close();
        dataset.close();

      } catch (H5::Exception& h5ex) {
        H5CFS_EXCEPTION("Could not set entries in 1D-Array '" 
                  << name << "':\n" << h5ex.getCDetailMsg());
      } catch( std::string& ex ) {
        H5CFS_RETHROW_EXCEPTION(ex, "Could not set entries in 1D-Array '" 
                          << name << "'" );
      }
  }

  void Hdf5Common::WriteCompound( H5::CommonFG& loc,
                            const std::string& name,
                            const CompoundType comp,
                            const H5::DSetCreatPropList &create_plist ) {

    try {
      // collect datatypes for compound array
      std::vector<shared_ptr<BaseHdfTypeConversion > > conv (comp.size() );  
      std::vector<H5::CompType> memCompTypes( comp.size() );
      unsigned int totalSize = 0;
      for( unsigned int i = 0; i < comp.size(); i++ ) {

        // get name of compound member
        std::string memName = comp[i].first;

        // obtain conversion object for generic compound member
        GetAnyConversion( comp[i].second, conv[i] );
        if( !conv[i]->IsSet() ) {
          H5CFS_EXCEPTION( "Could not map datatype for member '" << memName
                     << "' of compound '" << name << "'" );
        }
        totalSize += conv[i]->GetRawSize();

        // create new compound for memory datatype
        memCompTypes[i] = H5::CompType( (size_t)conv[i]->GetRawSize() );
        memCompTypes[i].insertMember( memName, 0, *conv[i]->GetNativeType() );
      }

      // Create file compound data type
      H5::CompType fileCompType( (size_t) totalSize );
      unsigned int actOffset = 0;
      for( unsigned int i = 0; i < comp.size(); i++ ) {
        fileCompType.insertMember( comp[i].first,
                                   actOffset,
                                   *conv[i]->GetStdType() );
        actOffset += conv[i]->GetRawSize();
      }

      // create  data space
      hsize_t dims = 1;
      H5::DataSpace space( 1, &dims );

      // generate dataset
      H5::DataSet dataset = loc.createDataSet( name, fileCompType,
                                               space, create_plist );

      // iterate again over all entries and fill in values
      for( unsigned int i = 0; i < comp.size(); i++ ) {
        dataset.write( conv[i]->GetOutBufferPtr(), memCompTypes[i] );
        conv[i]->CleanUp();
      }

      // close datatset and dataspace
      dataset.close();
      space.close();

    }  catch (H5::Exception& h5ex) {
      H5CFS_EXCEPTION("Could not write compound '" << name
                       << "':\n" << h5ex.getCDetailMsg());
    } catch( std::string& ex ) {
      H5CFS_RETHROW_EXCEPTION(ex, "Could not write compound '" << name << "'" );
    }
  }

  std::string Hdf5Common::GetObjNameByIdx( const H5::CommonFG& loc, hsize_t idx ) {
    ssize_t name_len = H5Gget_objname_by_idx(loc.getLocId(), idx, NULL, 0);
    if(name_len < 0) {
      H5CFS_EXCEPTION("Was not able to determine name" );
    }

    // now, allocate C buffer to get the name
    // note: obvisously we have read one more byte, as the C API seems not
    // to account for the trailing '\0' for strings, when determining the
    // length of the blank c-array
    char* name_C = new char[name_len+1];

    name_len = H5Gget_objname_by_idx(loc.getLocId(), idx, name_C, name_len+1);

    // clean up and return the string
    std::string name = std::string(name_C);
    delete []name_C;

    return name;
  }

  template<typename TYPE>
  void Hdf5Common::ReadAttribute( H5::H5Object& obj,
                            const std::string& name,
                            TYPE& data ) {
    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;

      // open attribute
      H5::Attribute attribute = obj.openAttribute( name );

      conv.AdaptToTypeBeforeRead(attribute.getDataType());
      H5::DataType* nativeType = conv.GetNativeType();

      // read data in buffer of conversion object
      attribute.read( *nativeType, conv.GetInBufferPtr(1) );

      // obtain data pointer from conversion object and
      // copy it to return buffer
      conv.GetNativeData( &data );

      // reset conversion object
      conv.CleanUp();

      // close attribute
      attribute.close();

    } catch (H5::Exception& h5ex) {
      H5CFS_EXCEPTION("Could not read Attribute '" 
                      << name << "':\n" << h5ex.getCDetailMsg());
    } catch( std::string& ex ) {
      H5CFS_RETHROW_EXCEPTION(ex, "Could not read Attribute '" << name << "'" );
    }

  }


  std::vector<unsigned int> Hdf5Common::GetArrayDims( const H5::CommonFG &loc,
                                      const std::string& name ) {

    H5::DataSet dataset = loc.openDataSet( name );
    H5::DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();

    hsize_t * myDims = new hsize_t[rank];
    rank = dataspace.getSimpleExtentDims( myDims, NULL);

    std::vector<unsigned int> dims( rank );
    for( unsigned int i = 0; i < (unsigned int) rank; i++ ) {
      dims[i] = static_cast<unsigned int>(myDims[i]);
    }

    delete[] myDims;
    return dims;
  }

  unsigned int Hdf5Common::GetNumEntries( const H5::CommonFG &loc,
                            const std::string& name ) {

    H5::DataSet dataset = loc.openDataSet( name );
    H5::DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();

    hsize_t * myDims = new hsize_t[rank];
    rank = dataspace.getSimpleExtentDims( myDims, NULL);

    unsigned int numEntries = 1;
    for( unsigned int i = 0; i < (unsigned int) rank; i++ ) {
      numEntries *= static_cast<unsigned int>(myDims[i]);
    }

    delete[] myDims;
    return numEntries;
  }


  template<typename TYPE>
  void Hdf5Common::ReadArray( H5::CommonFG &loc,
                        const std::string& name,
                        TYPE* data ) {

    // check, that data buffer is not empty
    if( data == NULL ) {
      H5CFS_EXCEPTION( "Data buffer for reading array '" << name 
                       << "' is NULL" );
    }
    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;

      // open dataset
      H5::DataSet dataset = loc.openDataSet( name );

      // get number of elements in the dataset
      H5::DataSpace dataspace = dataset.getSpace();

      // calculate absolute number of entries
      unsigned int size = static_cast<unsigned int>(dataspace.getSelectNpoints());
      // check, that standard hdf5 datatype is the same as the datatype
      // of the stored dataset
      // .. not sure if this can be implemented properly ...
      conv.AdaptToTypeBeforeRead(dataset.getDataType());
      H5::DataType* nativeType = conv.GetNativeType();

      // read data in buffer of conversion object
      dataset.read( conv.GetInBufferPtr(size), *nativeType );

      // obtain data pointer from conversion object and
      // copy it to return buffer
      conv.GetNativeData( data );

      // reset conversion object
      conv.CleanUp();

      // close dataset and space
      dataset.close();
      dataspace.close();

    } catch (H5::Exception& h5ex) {
      H5CFS_EXCEPTION( "Could not read Array '" 
                       << name << "':\n" << h5ex.getCDetailMsg());
    } catch( std::string& ex ) {
      H5CFS_RETHROW_EXCEPTION(ex, "Could not read Array '" << name << "'" );
    }

  }


  template<typename TYPE>
  void Hdf5Common::ReadArray( H5::CommonFG &loc,
                        const std::string& name,
                        std::vector<TYPE>& data ) {

    // clear data
    data.clear();

    // obtain information about dimension of dataset
    unsigned int numEntries = GetNumEntries( loc, name );

    // create temporary buffer
    TYPE * buffer = new TYPE[numEntries];

    // read data into buffer
    ReadArray( loc, name, buffer );

    // copy buffer data to vector
    data.resize( numEntries );
    for( unsigned int i = 0; i < numEntries; i++ ) {
      data[i] = buffer[i];
    }

    // delete buffer
    delete[] buffer;
  }
  
  
  void Hdf5Common::GetAnyConversion( const boost::any& anyType,
                               shared_ptr<BaseHdfTypeConversion>& conv ) {

    // query type of any
#define ANY_CONVERSION( TYPE )                          \
    if( anyType.type() == typeid(TYPE) ){               \
      shared_ptr<HdfTypeConversion< TYPE > >              \
        myConv ( new HdfTypeConversion< TYPE >() );       \
      myConv->SetNativeData( any_cast< TYPE >(anyType) );       \
      conv = myConv;                                    \
      return;                                           \
    }
    // define conversion for standard data types
    ANY_CONVERSION( bool );
    ANY_CONVERSION( unsigned int );
    ANY_CONVERSION( int );
    ANY_CONVERSION( double );
    ANY_CONVERSION( float );
    ANY_CONVERSION( std::string );
    ANY_CONVERSION( std::vector<unsigned int> );
    ANY_CONVERSION( std::vector<int> );
    ANY_CONVERSION( std::vector<double> );
    ANY_CONVERSION( std::vector<float> );
    ANY_CONVERSION( std::vector<std::string> );

#undef ANY_CONVERSION

    H5CFS_EXCEPTION( "Could not convert ANY-type "
                     << anyType.type().name() << "to known value" );
      
  }



  // instantiate all members for different types of values to be written
#define DECL_IO_METHODS( TYPE )                                 \
  template                                                      \
  void Hdf5Common::WriteAttribute<TYPE>( H5::H5Object& obj,           \
                                   const std::string& name,     \
                                   const TYPE& data,            \
                                   const H5::PropList           \
                                   &create_plist );             \
    template                                                    \
    void Hdf5Common::Write1DArray<TYPE>(H5::CommonFG &loc,            \
                                  const std::string& name,      \
                                  unsigned int size,                    \
                                  const TYPE * buffer,          \
                                  const H5::DSetCreatPropList   \
                                  &create_plist);               \
    template                                                    \
    void Hdf5Common::Reserve1DArray<TYPE>(H5::CommonFG &loc,          \
                                    const std::string& name,    \
                                    unsigned int size,                  \
                                    const H5::DSetCreatPropList \
                                    &create_plist);             \
    template                                                    \
    void Hdf5Common::SetEntries1DArray( H5::CommonFG &loc,            \
                                  const std::string& name,      \
                                  unsigned int start, unsigned int end,         \
                                  const TYPE * buffer );        \
    template                                                    \
    void Hdf5Common::Write2DArray<TYPE>( H5::CommonFG &loc,           \
                                   const std::string& name,     \
                                   unsigned int rowSize,                \
                                   unsigned int colSize,                \
                                   const TYPE * buffer,         \
                                   const H5::DSetCreatPropList  \
                                   &create_plist );             \
    template                                                    \
    void Hdf5Common::Reserve2DArray<TYPE>( H5::CommonFG &loc,         \
                                     const std::string& name,   \
                                     unsigned int rowSize,              \
                                     unsigned int colSize,              \
                                     const H5::DSetCreatPropList\
                                     &create_plist );           \
    template                                                    \
    void Hdf5Common::SetEntries2DArray( H5::CommonFG &loc,            \
                                  const std::string& name,      \
                                  unsigned int rowBegin, unsigned int rowEnd,   \
                                  unsigned int colBegin, unsigned int colEnd,   \
                                  const TYPE * buffer );        \
    template                                                    \
    void Hdf5Common::ReadAttribute( H5::H5Object& obj,                \
                              const std::string& name,          \
                              TYPE& data );                     \
                                                                \
    template                                                    \
    void Hdf5Common::ReadArray<TYPE>( H5::CommonFG &loc,              \
                                const std::string& name,        \
                                std::vector<TYPE>& data );        \
                                                                \
    template                                                    \
    void Hdf5Common::ReadArray<TYPE>( H5::CommonFG &loc,              \
                                const std::string& name,        \
                                TYPE* data )
  DECL_IO_METHODS( bool );
  DECL_IO_METHODS( int );
  DECL_IO_METHODS( unsigned int );
  DECL_IO_METHODS( double );
  DECL_IO_METHODS( float );
  DECL_IO_METHODS( std::string );
  DECL_IO_METHODS( std::vector<int> );
  DECL_IO_METHODS( std::vector<unsigned int> );
  DECL_IO_METHODS( std::vector<double> );
  DECL_IO_METHODS( std::vector<float> );
  DECL_IO_METHODS( std::vector<std::string> );

#undef DECL_IO_METHODS


  // =======================================================================
  //  GENERAL ACCESS METHODS
  // =======================================================================
  
  H5::Group Hdf5Common::GetMultiStepGroup( H5::H5File& file, 
                                         unsigned int msStep,
                                         bool isHistory ) {

    // open group with multisteps
    H5::Group resultGroup;
    try {
      if( !isHistory ) {
        resultGroup = file.openGroup("/Results/Mesh");
      } else {
        resultGroup = file.openGroup("/Results/History");
      }
    } H5_CATCH( "Could not open result group" );

    // open specified msgroup
    H5::Group actMsGroup;
    std::string actMsName = "MultiStep_"
      + boost::lexical_cast<std::string>( msStep );
    try {
      actMsGroup = resultGroup.openGroup( actMsName );
    } H5_CATCH( "Could not open group for multistep " << msStep );

    resultGroup.close();
    return actMsGroup;
  }

  H5::Group Hdf5Common::GetStepGroup( H5::H5File& file, 
                                unsigned int msStep, 
                                unsigned int stepNum ) {
    // get multistep group
    H5::Group actMsGroup = GetMultiStepGroup( file, msStep, false );

    std::string groupName = "Step_" +
      boost::lexical_cast<std::string> (stepNum );

    H5::Group stepGroup;
    try {
      stepGroup = actMsGroup.openGroup( groupName );
    } H5_CATCH( "Could not open group for results of step " << stepNum
                << " in multiStep " << msStep );
    actMsGroup.close();
    return stepGroup;
  }
  
  void Hdf5Common::SetMaxChunkSize( unsigned int chunkSize ) {
    Hdf5Common::maxChunkSize_ = chunkSize;
    
  }
//  
//
//  int Hdf5Common::MapCapabilityType( SimOutput::Capability c ) {
//    int ret = 0;
//    switch( c ) {
//    case SimOutput::NONE:
//      ret = 0;
//      break;
//    case SimOutput::MESH:
//      ret = 1;
//      break;
//    case SimOutput::MESH_RESULTS:
//      ret = 2;
//      break;
//    case SimOutput::HISTORY:
//      ret = 3;
//      break;
//    case SimOutput::USERDATA:
//      ret = 4;
//      break;
//    default:
//      H5CFS_EXCEPTION "Could not map capability '" << c 
//                 << "' to hdf5 representation" );
//    }
//    return ret;
//  }
//  
//  SimOutput::Capability Hdf5Common::MapCapabilityType( int c ) {
//    SimOutput::Capability ret = SimOutput::NONE;
//    switch ( c ) {
//    case 0:
//      ret = SimOutput::NONE;
//      break;
//    case 1:
//      ret = SimOutput::MESH;
//      break;
//    case 2:
//      ret = SimOutput::MESH_RESULTS;
//      break;
//    case 3:
//      ret = SimOutput::HISTORY;
//      break;
//    case 4:
//      ret = SimOutput::USERDATA;
//      break;
//    default:
//      break;
//    }
//    return ret;
//  }
//  
//  
//  int Hdf5Common::MapUnknownType( ResultInfo::EntityUnknownType t ) {
//    int definedOn = 0;
//    switch(t) {
//    case ResultInfo::NODE:
//      definedOn = 1;
//      break;
//    case ResultInfo::EDGE:
//      definedOn = 2;
//      break;
//    case ResultInfo::FACE:
//      definedOn = 3;
//      break;
//    case ResultInfo::ELEMENT:
//      definedOn = 4;
//      break;
//    case ResultInfo::SURF_ELEM:
//      definedOn = 5;
//      break;
//    case ResultInfo::PFEM:
//      definedOn = 6;
//      break;
//    case ResultInfo::REGION:
//      definedOn = 7;
//      break;
//    case ResultInfo::SURF_REGION:
//      definedOn = 8;
//      break;
//    case ResultInfo::NODELIST:
//      definedOn = 9;
//      break;
//    case ResultInfo::COIL:
//      definedOn = 10;
//      break;
//    case ResultInfo::FREE:
//      definedOn = 11;
//      break;
//    }
//
//    return definedOn;
//
//  }
//  
//  std::string Hdf5Common::MapUnknownTypeAsString( ResultInfo::EntityUnknownType t ) {
//     std::string definedOn = "";
//     switch(t) {
//     case ResultInfo::NODE:
//       definedOn = "Nodes";
//       break;
//     case ResultInfo::EDGE:
//       definedOn = "Edges";
//       break;
//     case ResultInfo::FACE:
//       definedOn = "Faces";
//       break;
//     case ResultInfo::ELEMENT:
//       definedOn = "Elements";
//       break;
//     case ResultInfo::SURF_ELEM:
//       definedOn = "Elements";
//       break;
//     case ResultInfo::PFEM:
//       definedOn = "Nodes";
//       break;
//     case ResultInfo::REGION:
//       definedOn = "Regions";
//       break;
//     case ResultInfo::SURF_REGION:
//       definedOn = "ElementGroup";
//       break;
//     case ResultInfo::NODELIST:
//       definedOn = "NodeGroup";
//       break;
//     case ResultInfo::COIL:
//       definedOn = "Coils";
//       break;
//     case ResultInfo::FREE:
//       definedOn = "Unknowns";
//       break;
//     }
//
//     return definedOn;
//
//   }
//  
//  ResultInfo::EntityUnknownType Hdf5Common::MapUnknownType( int t ) {
//
//    ResultInfo::EntityUnknownType definedOn;
//    switch(t) {
//    case 1:
//      definedOn = ResultInfo::NODE;
//      break;
//    case 2:
//      definedOn = ResultInfo::EDGE;
//      break;
//    case 3: 
//      definedOn = ResultInfo::FACE;
//      break;
//    case 4:
//      definedOn = ResultInfo::ELEMENT;
//      break;
//    case 5:
//      definedOn = ResultInfo::SURF_ELEM;
//      break;
//    case 6:
//      definedOn = ResultInfo::PFEM;
//      break;
//    case 7:
//      definedOn = ResultInfo::REGION;
//      break;
//    case 8:
//      definedOn = ResultInfo::SURF_REGION;
//      break;
//    case 9:
//      definedOn = ResultInfo::NODELIST;
//      break;
//    case 10:
//      definedOn = ResultInfo::COIL;
//      break;
//    case 11:
//      definedOn = ResultInfo::FREE;
//      break;
//    }
//
//    return definedOn;
//  }
//  
//  
//  int Hdf5Common::MapEntryType( ResultInfo::EntryType t ) {
//    int entryType = 0;
//    
//    switch(t) {
//    case ResultInfo::UNKNOWN:
//      entryType = 0;
//      break;
//    case ResultInfo::SCALAR:
//      entryType = 1;
//      break;
//    case ResultInfo::VECTOR:
//      entryType = 3;
//      break;
//    case ResultInfo::TENSOR:
//      entryType = 6;
//      break;
//    case ResultInfo::STRING:
//      entryType = 32;
//      break;
//    }
//    return entryType;
//  }
//
//  
//  ResultInfo::EntryType Hdf5Common::MapEntryType( int t ) {
//
//    ResultInfo::EntryType entryType;
//    
//    switch(t) {
//    case 0:
//      entryType =ResultInfo::UNKNOWN;
//      break;
//    case 1:
//      entryType = ResultInfo::SCALAR;
//      break;
//    case 3:
//      entryType = ResultInfo::VECTOR;
//      break;
//    case 6:
//      entryType = ResultInfo::TENSOR;
//      break;
//    case 32: 
//      entryType = ResultInfo::STRING;
//      break;
//    }
//    return entryType;
//  }
//
  void Hdf5Common::CheckOpenObjects(H5::H5File& file, bool verbose) {
    std::vector<unsigned int> types;
    std::vector<std::string> typeNames;
    hid_t* ids;

    types.push_back(H5F_OBJ_DATASET); typeNames.push_back("Dataset");
    types.push_back(H5F_OBJ_GROUP); typeNames.push_back("Group");
    types.push_back(H5F_OBJ_DATATYPE); typeNames.push_back("DataType");
    types.push_back(H5F_OBJ_ATTR); typeNames.push_back("Attribute");

    // check for open groups, datasets etc.
    if(verbose)
      std::cerr << "Number of open objects: "
                << "in file " << file.getFileName() << "\n"
                << "--------------------------\n";

    for(unsigned int t=0; t<types.size(); t++)
    {
      unsigned int numObjs = file.getObjCount(types[t]);

      if(verbose)
        std::cerr << typeNames[t] << "s: "<<  numObjs << std::endl;

      ids = new hid_t[numObjs];
      file.getObjIDs(types[t], numObjs, ids);

      for(unsigned int i=0; i<numObjs; i++)
      {
        H5::DataSet ds;
        H5::Group group;
        H5::DataType dt;
        H5::Attribute attr;
        H5::IdComponent *idComp = NULL;

        switch(types[t])
        {
        case H5F_OBJ_DATASET:
          ds.setId((ids[i]));
          idComp = &ds;
          if(verbose)
            std::cerr << "  " << ds.fromClass() << std::endl;
          break;
        case H5F_OBJ_GROUP:
          group.setId((ids[i]));
          idComp = &group;
          if(verbose)
            std::cerr << "  " << group.fromClass() << std::endl;
          break;
        case H5F_OBJ_DATATYPE:
          dt.setId((ids[i]));
          idComp = &dt;
          if(verbose)
            std::cerr << "  " << dt.fromClass() << std::endl;
          break;
        case H5F_OBJ_ATTR:
          attr.setId((ids[i]));
          idComp = &attr;
          if(verbose)
            std::cerr << "  " << attr.fromClass() << std::endl;
          break;
        }

        while ( idComp->getCounter() > 0 )
        {
          if(verbose)
            std::cerr << "  Trying to dec refcounter of HDF5 ID " << ids[i]
                      << "..." << std::endl;            
          idComp->decRefCount();
        }

      }

      delete[] ids;
    }
  }


} // end of namespace CoupledField 
