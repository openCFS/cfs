#include "hdf5io.hh"


namespace CoupledField {


  // ====================================
  //    Initialize Atom Data Typemaps
  // ====================================
#define DECL_HDF_ATOM_TYPE(TYPE, NATIVE_TYPE, STD_TYPE )        \
  template<>                                                    \
  H5::DataType H5IO::HdfAtomTypeMap<TYPE>::HdfNativeType =      \
    H5::IntType( NATIVE_TYPE );                                 \
    template<>                                                  \
    H5::DataType H5IO::HdfAtomTypeMap<TYPE>::HdfStdType =       \
      H5::IntType( STD_TYPE )
  
  
  // Define mapping
  DECL_HDF_ATOM_TYPE( bool,
                      H5::PredType::NATIVE_INT32,
                      H5::PredType::STD_I32LE );
  
  DECL_HDF_ATOM_TYPE( Integer,
                      H5::PredType::NATIVE_INT32,
                      H5::PredType::STD_I32LE );
  
  DECL_HDF_ATOM_TYPE( UInt,
                      H5::PredType::NATIVE_UINT32,
                      H5::PredType::STD_U32LE );
  
  DECL_HDF_ATOM_TYPE( Double,
                      H5::PredType::NATIVE_DOUBLE,
                      H5::PredType::IEEE_F64LE );
  
#undef DECL_HDF_ATOM_TYPE
  
  // ========================================
  //    Initialize General Datatype mapping
  // ========================================

  // -----------
  //  AtomTypes
  // -----------
  
#define DECL_HDF_ATOM_TYPE_CONV(TYPE )                          \
  template<>                                                    \
  class H5IO::HdfTypeConversion<TYPE> :                         \
    public H5IO::BaseHdfTypeConversion {                        \
                                                                \
  private:                                                      \
                                                                \
    TYPE* buffer_;                                              \
                                                                \
  public:                                                       \
    HdfTypeConversion()                                         \
      : buffer_( NULL ) {                                        \
      nativeType_ =                                              \
        HdfAtomTypeMap<TYPE>::HdfNativeType;                     \
      stdType_ =                                                 \
        HdfAtomTypeMap<TYPE>::HdfStdType;                        \
    }                                                           \
                                                                \
    const void * GetOutBufferPtr() {                            \
      if( !isSet_ ) {                                           \
        EXCEPTION( "Data buffer is empty" );                    \
      }                                                         \
      return (const void*) buffer_;                             \
    }                                                           \
                                                                \
    void* GetInBufferPtr( UInt numData ) {                      \
      CleanUp();                                                \
      numElems_ = numData;                                      \
      buffer_ = new TYPE[numData];                              \
      return (void*) buffer_;                                   \
    }                                                           \
                                                                \
    void GetNativeData( TYPE * data) {                          \
      for( UInt i = 0; i < numElems_; i++ ) {                   \
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
    void SetNativeData( const TYPE* t, UInt size  ) {           \
      CleanUp();                                                \
      buffer_ = new TYPE[size];                                 \
      for( UInt i = 0; i < size; i ++ ) {                       \
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

  DECL_HDF_ATOM_TYPE_CONV(Integer);
  DECL_HDF_ATOM_TYPE_CONV(UInt);
  DECL_HDF_ATOM_TYPE_CONV(Double);
  
#undef DECL_HDF_ATOM_TYPE_CONV
  
  // ---------------
  //  bool
  // ---------------
  template<>
  class H5IO::HdfTypeConversion<bool> :
    public H5IO::BaseHdfTypeConversion  {
    
  private:

    Integer* buffer_;
    
  public:
    
    HdfTypeConversion() 
      : buffer_ ( NULL ) {
      nativeType_ =  HdfAtomTypeMap<bool>::HdfNativeType;
      stdType_ = HdfAtomTypeMap<bool>::HdfStdType;
    }

    const void * GetOutBufferPtr() {
      if( !isSet_ ) {
        EXCEPTION( "Data buffer is empty" );
      }
      return buffer_;
    }
    
    void* GetInBufferPtr( UInt numData ) {
      CleanUp();
      numElems_ = numData;
      buffer_ = new Integer[numElems_];
      return buffer_;
    }

    void GetNativeData( bool * data) {
      for( UInt i = 0; i < numElems_; i++ ) {
        data[i] = (buffer_[i] == 1) ? true : false;
      }
    }
    
    void SetNativeData( const bool& t ) {
      CleanUp();
      buffer_ = new Integer[1];
      buffer_[0] = t ? 1 : 0;
      size_ = sizeof( Integer );
      numElems_ = 1;
      isSet_ = true;
    }
    
    void SetNativeData( const bool* t, UInt size  ) {
      CleanUp();
      buffer_ = new Integer[size];
      for( UInt i = 0; i < size; i ++ ) {
        buffer_[i] = t[i] ? 1 : 0;
      }
      size_ = sizeof( Integer* );
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
  class H5IO::HdfTypeConversion<std::string> :
    public H5IO::BaseHdfTypeConversion  {
    
  private:

    const char** buffer_;
    
  public:
    HdfTypeConversion() 
      : buffer_( NULL ) {
      nativeType_ = H5::StrType( H5::PredType::C_S1, H5T_VARIABLE);
      stdType_ = H5::StrType( H5::PredType::C_S1, H5T_VARIABLE);
    }
    
    const void * GetOutBufferPtr() {
      if( !isSet_ ) {
        EXCEPTION( "Data buffer is empty" );
      }
      return buffer_;
    }

    void * GetInBufferPtr( UInt numData ) {
      CleanUp();
      numElems_ = numData;
      buffer_ = new const char*[numElems_];
      return buffer_;
    }

    void GetNativeData( std::string * data) {
      for( UInt i = 0; i < numElems_; i++ ) {
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
    
    void SetNativeData( const std::string* t, UInt size  ) {
      CleanUp();
      buffer_ = new const char*[size];
      for( UInt i = 0; i < size; i ++ ) {
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
  };
  
  
  // --------------------------
  //  std::vector<std::string>
  // --------------------------
  template<>
  class H5IO::HdfTypeConversion<std::vector<std::string> > :
    public H5IO::BaseHdfTypeConversion {

  private:
    hvl_t * buffer_;

  public:
    HdfTypeConversion() 
      : buffer_( NULL ) { 
      H5::StrType sType( H5::PredType::C_S1, H5T_VARIABLE);
      H5::VarLenType vType( &sType );
      nativeType_ = vType;
      stdType_ = vType;
    }
    
    const void * GetOutBufferPtr() {
      if( !isSet_ ) {
        EXCEPTION( "Data buffer is empty" );
      }
      return buffer_;
    }

    void* GetInBufferPtr( UInt numData ) {
      CleanUp();
      numElems_ = numData;
      buffer_ = new hvl_t[numData];
      return buffer_;
    }
    
    void GetNativeData( std::vector<std::string> * data) {
      for( UInt i = 0; i < numElems_; i++ ) {
        data[i].resize( buffer_[i].len );
        for( UInt j = 0; j < buffer_[i].len; j++ ) {
          data[i][j].assign( ((const char**)buffer_[i].p)[j] );
        }
      }
    }

    void SetNativeData( const std::vector<std::string>& t ) {
      CleanUp();
      buffer_ = new hvl_t[1];
      buffer_[0].p = (void*) new const char*[t.size()];
      buffer_[0].len = t.size();
      for( UInt j = 0; j < t.size(); j++ ) {
        ((const char **)buffer_[0].p)[j] = t[j].c_str();
      }
      size_ = sizeof( hvl_t );
      numElems_ = 1;
      isSet_ = true;
    }

    
    void SetNativeData( const std::vector<std::string>* t, UInt size ) {
      CleanUp();
      buffer_ = new hvl_t[size];
      for( UInt i = 0; i < size; i++ ) {
        buffer_[i].p = (void*) new const char*[t[i].size()];
        buffer_[i].len = t[i].size();
        for( UInt j = 0; j < t[i].size(); j++ ) {
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
        for( UInt i = 0; i < numElems_; i++ ) {
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
    
  }; // end of class definition
  
  
  
  // ----------------------
  //  std::vector<TYPE>
  // ----------------------
#define DECL_STL_VECTOR_CONVERSION( TYPE )                              \
  template<>                                                            \
  class H5IO::HdfTypeConversion<std::vector<TYPE> > :                   \
    public H5IO::BaseHdfTypeConversion {                                \
  private:                                                              \
    hvl_t * buffer_;                                                    \
                                                                        \
  public:                                                               \
    HdfTypeConversion()                                                 \
      : buffer_( NULL ) {                                               \
      nativeType_ =                                                     \
        H5::VarLenType ( &HdfAtomTypeMap<TYPE>::HdfNativeType );        \
      stdType_ =                                                        \
        H5::VarLenType ( &HdfAtomTypeMap<TYPE>::HdfStdType );           \
    }                                                                   \
                                                                        \
    const void * GetOutBufferPtr( ) {                                   \
      if( !isSet_ ) {                                                   \
        EXCEPTION( "Data buffer is empty" );                            \
      }                                                                 \
      return (const void*) buffer_;                                     \
    }                                                                   \
                                                                        \
    void* GetInBufferPtr( UInt numData ) {                              \
      CleanUp();                                                        \
      numElems_ = numData;                                              \
      buffer_ = new hvl_t[numData];                                     \
      return buffer_;                                                   \
    }                                                                   \
                                                                        \
    void GetNativeData( std::vector<TYPE> * data) {                     \
      for( UInt i = 0; i < numElems_; i++ ) {                           \
        data[i].resize( buffer_[i].len );                               \
        for( UInt j = 0; j < buffer_[i].len; j++ ) {                    \
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
      for( UInt j = 0; j < t.size(); j++ ) {                            \
        ((TYPE*)buffer_[0].p)[j] = t[j];                                \
      }                                                                 \
      numElems_ = 1;                                                    \
      size_ = sizeof( hvl_t );                                          \
      isSet_ = true;                                                    \
    }                                                                   \
                                                                        \
    void SetNativeData( const std::vector<TYPE>* t, UInt size ) {       \
      CleanUp();                                                        \
      buffer_ = new hvl_t[size];                                        \
      for( UInt i = 0; i < size; i++ ) {                                \
        buffer_[i].p =  new TYPE[t[i].size()];                          \
        buffer_[i].len = t[i].size();                                   \
        for( UInt j = 0; j < t[i].size(); j++ ) {                       \
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
      for( UInt i = 0; i < numElems_; i++ ) {                           \
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
    }                                                                   \
                                                                        \
  }
  
  DECL_STL_VECTOR_CONVERSION( Integer );
  DECL_STL_VECTOR_CONVERSION( UInt );
  DECL_STL_VECTOR_CONVERSION( Double );
    
#undef DECL_ST_VECTOR_CONVERSION
    
  template<typename TYPE>
  void H5IO::WriteAttribute( H5::H5Object& obj,
                             const std::string& name,
                             const TYPE& data,
                             const H5::DSetCreatPropList &create_plist ) {
    
    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType stdType = conv.GetStdType();
      H5::DataType nativeType = conv.GetNativeType();

      // create memory data space
      H5::DataSpace space;
      
      // generate attribute
      conv.SetNativeData( data );
      if( !conv.IsSet() ) {
        EXCEPTION( "Could not convert data for attribute '"
                   << name << "' of type " << typeid(TYPE).name() );
      }
      H5::Attribute attr = obj.createAttribute( name, stdType, 
                                                space, create_plist );
      
      // write attribute
      attr.write( nativeType, conv.GetOutBufferPtr() );
      
      // reset conversion object
      conv.CleanUp();
      
      // close attribute, dataspace- and types
      space.close();
      attr.close();
      

    }  catch (H5::Exception& h5ex) {
      EXCEPTION("Could not write attribute '" 
                << name << "':\n" << h5ex.getCDetailMsg());
    } catch( Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not write attribute '" << name << "'" );
    } 
  }

  template<typename TYPE>
  void H5IO::Write1DArray( H5::CommonFG &loc,
                           const std::string& name,
                           UInt size,
                           const TYPE * buffer,
                           const H5::DSetCreatPropList &create_plist ) {
    
    // check, that size is greate than zero
    if( size == 0 || buffer == NULL ) {
      EXCEPTION( "Attribute data buffer of 1D array '" << name 
                 << "' is NULL or has zero size" );
    }
    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType stdType = conv.GetStdType();
      H5::DataType nativeType = conv.GetNativeType();
  
      // create memory data space
      const hsize_t dims = size;
      H5::DataSpace space( 1, &dims );

      // generate dataset and fill it
      conv.SetNativeData( buffer, size );
      if( !conv.IsSet() ) {
        EXCEPTION( "Could not convert data for 1D array '"
                   << name << "' of type " << typeid(TYPE).name() );
      }
      H5::DataSet dataset = loc.createDataSet( name, stdType, 
                                               space, create_plist );
      dataset.write( conv.GetOutBufferPtr(), nativeType  );

      // reset conversion object
      conv.CleanUp();

      // close dataset, dataspace- and types
      space.close();
      dataset.close();
      
      
    } catch (H5::Exception& h5ex) {
      EXCEPTION("Could not write 1D-Array '" 
                << name << "':\n" << h5ex.getCDetailMsg());
    } catch( Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not write 1D-Array '" << name << "'" );
    }
    
  }

  template<typename TYPE>
  void H5IO::Write2DArray( H5::CommonFG &loc,
                           const std::string& name,
                           UInt rowSize,
                           UInt colSize,
                           const TYPE * buffer,
                           const H5::DSetCreatPropList &create_plist
                           ) {
    
    // check, that size is greate than zero
    if( rowSize == 0 || colSize == 0 || buffer == NULL ) {
      EXCEPTION( "Data buffer of 2D array '" << name 
                 << "' is NULL or has zero size" );
    }
    
    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType stdType = conv.GetStdType();
      H5::DataType nativeType = conv.GetNativeType();
      
      // create memory data space
      const hsize_t dims[] = {rowSize, colSize};
      const Integer rank = 2;

      H5::DataSpace space( rank, dims );
      
      // generate dataset and fill it
      conv.SetNativeData( buffer, rowSize * colSize );
      if( !conv.IsSet() ) {
        EXCEPTION( "Could not convert data for 2D array '"
                   << name << "' of type " << typeid(TYPE).name() );
      }
      H5::DataSet dataset = loc.createDataSet( name, stdType, 
                                               space, create_plist );
      dataset.write( conv.GetOutBufferPtr(), nativeType  );

      // reset conversion object
      conv.CleanUp();

      // close dataset, dataspace- and types
      space.close();
      dataset.close();

    } catch (H5::Exception& h5ex) {
      EXCEPTION("Could not write 2D-Array '" << name 
                << "':\n" << h5ex.getCDetailMsg());
    } catch( Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not write 2D-Array '" << name << "'" );
    }
    

  }

  void H5IO::WriteCompound( H5::CommonFG& loc,
                            const std::string& name,
                            const CompoundType comp,
                            const H5::DSetCreatPropList &create_plist ) {

    try {
      // collect datatypes for compound array
      std::vector<shared_ptr<BaseHdfTypeConversion > > conv (comp.size() );  
      std::vector<H5::CompType> memCompTypes( comp.size() );
      UInt totalSize = 0;
      for( UInt i = 0; i < comp.size(); i++ ) {
        
        // get name of compound member
        std::string memName = comp[i].first;
        
        // obtain conversion object for generic compound member
        GetAnyConversion( comp[i].second, conv[i] );
        if( !conv[i]->IsSet() ) {
          EXCEPTION( "Could not map datatype for member '" << memName
                     << "' of compound '" << name << "'" );
        }
        totalSize += conv[i]->GetRawSize();
      
        // create new compound for memory datatype
        memCompTypes[i] = H5::CompType( (size_t)conv[i]->GetRawSize() );
        memCompTypes[i].insertMember( memName, 0, conv[i]->GetNativeType() );
      }

      // Create file compound data type
      H5::CompType fileCompType( (size_t) totalSize );
      UInt actOffset = 0;
      for( UInt i = 0; i < comp.size(); i++ ) {
        fileCompType.insertMember( comp[i].first,
                                   actOffset,
                                   conv[i]->GetStdType() );
        actOffset += conv[i]->GetRawSize();
      }

      // create  data space
      hsize_t dims = 1;
      H5::DataSpace space( 1, &dims );
    
      // generate dataset
      H5::DataSet dataset = loc.createDataSet( name, fileCompType, 
                                               space, create_plist );

      // iterate again over all entries and fill in values
      for( UInt i = 0; i < comp.size(); i++ ) {
        dataset.write( conv[i]->GetOutBufferPtr(), memCompTypes[i] );
        conv[i]->CleanUp();
      }

      // close datatset and dataspace
      dataset.close();
      space.close();

    }  catch (H5::Exception& h5ex) {
      EXCEPTION("Could not write compound '" << name 
                << "':\n" << h5ex.getCDetailMsg());
    } catch( Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not write compound '" << name << "'" );
    }
  }

  template<typename TYPE>
  void H5IO::ReadAttribute( H5::H5Object& obj,
                            const std::string& name,
                            TYPE& data ) {
    try {
      
      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType stdType = conv.GetStdType();
      H5::DataType nativeType = conv.GetNativeType();
      
      // open attribute
      H5::Attribute attribute = obj.openAttribute( name );
      
      // read data in buffer of conversion object
      attribute.read( nativeType, conv.GetInBufferPtr(1) );
      
      // obtain data pointer from conversion object and
      // copy it to return buffer
      conv.GetNativeData( &data );

      // reset conversion object
      conv.CleanUp();

      // close attribute
      attribute.close();
      
    } catch (H5::Exception& h5ex) {
      EXCEPTION("Could not read Attribute '" 
                << name << "':\n" << h5ex.getCDetailMsg());
    } catch( Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not read Attribute '" << name << "'" );
    }
    
  }


  UInt GetArrayDims( H5::CommonFG &loc,
                     const std::string& name,
                     std::vector<UInt>& dims ) {
    
    H5::DataSet dataset = loc.openDataSet( name );
    H5::DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    
    hsize_t * myDims = new hsize_t[rank];
    rank = dataspace.getSimpleExtentDims( myDims, NULL);

    dims.resize( rank );
    UInt numEntries = 1;
    for( UInt i = 0; i < (UInt) rank; i++ ) {
      dims[i] = myDims[i];
      numEntries *= dims[i];
    }
    delete[] myDims;
    return numEntries;
  }


  template<typename TYPE>
  void H5IO::ReadArray( H5::CommonFG &loc,
                        const std::string& name,
                        TYPE* data ) {
    
    // check, that data buffer is not empty
    if( data == NULL ) {
      EXCEPTION( "Data buffer for reading array '" << name 
                 << "' is NULL" );
    }
    try {

      // create conversion helper object and get native / std hdf5 datatype
      HdfTypeConversion<TYPE> conv;
      H5::DataType stdType = conv.GetStdType();
      H5::DataType nativeType = conv.GetNativeType();
  
      // open dataset
      H5::DataSet dataset = loc.openDataSet( name );

      // get number of elements in the dataset
      H5::DataSpace dataspace = dataset.getSpace();

      // calculate absolute number of entries
      UInt size = dataspace.getSelectNpoints();
      
      // check, that standard hdf5 datatype is the same as the datatype
      // of the stored dataset
      // .. not sure if this can be implemented properly ...

      // read data in buffer of conversion object
      dataset.read( conv.GetInBufferPtr(size), nativeType );
      
      // obtain data pointer from conversion object and
      // copy it to return buffer
      conv.GetNativeData( data );

      // reset conversion object
      conv.CleanUp();

      // close dataset and space
      dataset.close();
      dataspace.close();
      
    } catch (H5::Exception& h5ex) {
      EXCEPTION("Could not read Array '" 
                << name << "':\n" << h5ex.getCDetailMsg());
    } catch( Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not read Array '" << name << "'" );
    }
    
  }
  
  void H5IO::GetAnyConversion( const boost::any& anyType,
                               shared_ptr<BaseHdfTypeConversion>& conv ) {
    
    // query type of any
#define ANY_CONVERSION( TYPE )                          \
    if( anyType.type() == typeid(TYPE) ){               \
      shared_ptr<HdfTypeConversion<TYPE> >              \
        myConv ( new HdfTypeConversion<TYPE>() );       \
      myConv->SetNativeData( any_cast<TYPE>(anyType) );       \
      conv = myConv;                                    \
      return;                                           \
    }
    // define conversion for standard data types
    ANY_CONVERSION( bool );
    ANY_CONVERSION( UInt );
    ANY_CONVERSION( Integer );
    ANY_CONVERSION( Double );
    ANY_CONVERSION( std::string );
    ANY_CONVERSION( std::vector<UInt> );
    ANY_CONVERSION( std::vector<Integer> );
    ANY_CONVERSION( std::vector<Double> );
    ANY_CONVERSION( std::vector<std::string> );

#undef ANY_CONVERSION

    EXCEPTION( "Could not convert ANY-type "
               << anyType.type().name() << "to known value" );
      
  }



  // instantiate all members for different types of values to be written
#define DECL_IO_METHODS( TYPE )                                 \
  template                                                      \
  void H5IO::WriteAttribute<TYPE>( H5::H5Object& obj,           \
                                   const std::string& name,     \
                                   const TYPE& data,            \
                                   const H5::DSetCreatPropList  \
                                   &create_plist );             \
    template                                                    \
    void H5IO::Write1DArray<TYPE>(H5::CommonFG &loc,            \
                                  const std::string& name,      \
                                  UInt size,                    \
                                  const TYPE * buffer,          \
                                  const H5::DSetCreatPropList   \
                                  &create_plist);               \
    template                                                    \
    void H5IO::Write2DArray<TYPE>( H5::CommonFG &loc,           \
                                   const std::string& name,     \
                                   UInt rowSize,                \
                                   UInt colSize,                \
                                   const TYPE * buffer,         \
                                   const H5::DSetCreatPropList  \
                                   &create_plist );             \
    template                                                    \
    void H5IO::ReadAttribute( H5::H5Object& obj,                \
                              const std::string& name,          \
                              TYPE& data );                     \
                                                                \
    template                                                    \
    void H5IO::ReadArray<TYPE>( H5::CommonFG &loc,              \
                                const std::string& name,        \
                              TYPE* data )
  
  DECL_IO_METHODS( bool );
  DECL_IO_METHODS( Integer );
  DECL_IO_METHODS( UInt );
  DECL_IO_METHODS( Double );
  DECL_IO_METHODS( std::string );
  DECL_IO_METHODS( std::vector<Integer> );
  DECL_IO_METHODS( std::vector<UInt> );
  DECL_IO_METHODS( std::vector<Double> );
  DECL_IO_METHODS( std::vector<std::string> );

#undef DECL_IO_METHODS

  Integer H5IO::MapUnknownType( ResultInfo::EntityUnknownType t ) {
    Integer definedOn = 0;
    switch(t) {
    case ResultInfo::NODE:
      definedOn = 1;
      break;
    case ResultInfo::EDGE:
      definedOn = 2;
      break;
    case ResultInfo::FACE:
      definedOn = 3;
      break;
    case ResultInfo::ELEMENT:
      definedOn = 4;
      break;
    case ResultInfo::SURF_ELEM:
      definedOn = 5;
      break;
    case ResultInfo::PFEM:
      definedOn = 6;
      break;
    case ResultInfo::REGION:
      definedOn = 7;
      break;
    case ResultInfo::SURF_REGION:
      definedOn = 8;
      break;
    case ResultInfo::NODELIST:
      definedOn = 9;
      break;
    case ResultInfo::COIL:
      definedOn = 10;
      break;
    case ResultInfo::FREE:
      definedOn = 11;
      break;
    }

    return definedOn;

  }
  
  ResultInfo::EntityUnknownType H5IO::MapUnknownType( Integer t ) {
    return ResultInfo::NODE;
  }
  
  
  Integer H5IO::MapEntryType( ResultInfo::EntryType t ) {
    Integer entryType = 0;
    
    switch(t) {
    case ResultInfo::UNKNOWN:
      entryType = 0;
      break;
    case ResultInfo::SCALAR:
      entryType = 1;
      break;
    case ResultInfo::VECTOR:
      entryType = 3;
      break;
    case ResultInfo::TENSOR:
      entryType = 6;
      break;
    case ResultInfo::STRING:
      entryType = 32;
      break;
    }
    return entryType;
  }

  
  ResultInfo::EntryType H5IO::MapEntryType( Integer t ) {
    return ResultInfo::UNKNOWN;
  }



} // end of namespace CoupledField
