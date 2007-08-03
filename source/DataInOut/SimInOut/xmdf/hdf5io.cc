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
    std::vector<TYPE> buffer_;                                  \
                                                                \
  public:                                                       \
    HdfTypeConversion() {                                       \
      nativeType_ =  HdfAtomTypeMap<TYPE>::HdfNativeType;       \
      stdType_ =  HdfAtomTypeMap<TYPE>::HdfStdType;             \
    }                                                           \
                                                                \
    ~HdfTypeConversion() {};                                    \
                                                                \
    const void * GetVoidPtr() {                                 \
      if( !isSet_ ) {                                           \
        EXCEPTION( "Data buffer is empty" );                    \
      }                                                         \
      return (const void*) &(buffer_[0]);                       \
    }                                                           \
                                                                \
    void SetData( const TYPE& t ) {                             \
      buffer_.resize(1);                                        \
      buffer_[0] = t;                                           \
      size_ = sizeof( TYPE );                                   \
      isSet_ = true;                                            \
    }                                                           \
                                                                \
    void SetData( const TYPE* t, UInt size  ) {                 \
      buffer_.resize( size );                                   \
      for( UInt i = 0; i < size; i ++ ) {                       \
        buffer_[i] = t[i];                                      \
      }                                                         \
      size_ = sizeof( TYPE* );                                  \
      isSet_ = true;                                            \
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
    std::vector<Integer> buffer_;
    
  public:
    HdfTypeConversion()  {
      nativeType_ =  HdfAtomTypeMap<bool>::HdfNativeType;
      stdType_ = HdfAtomTypeMap<bool>::HdfStdType;
    }
    
    const void * GetVoidPtr() {
      if( !isSet_ ) {
        EXCEPTION( "Data buffer is empty" );
      }
      return &buffer_[0];
    }
    
    void SetData( const bool& t ) {
      buffer_.resize(1);
      buffer_[0] = t ? 1 : 0;
      size_ = sizeof( Integer );
      isSet_ = true;
    }
    
    void SetData( const bool* t, UInt size  ) {
      buffer_.resize( size );
      for( UInt i = 0; i < size; i ++ ) {
        buffer_[i] = t[i] ? 1 : 0;
      }
      size_ = sizeof( Integer* );
      isSet_ = true;
    }
  };

  // ---------------
  //  std::string
  // ---------------
  template<>
  class H5IO::HdfTypeConversion<std::string> :
    public H5IO::BaseHdfTypeConversion  {
    
  private:
    std::vector<const char*> buffer_;
    
  public:
    HdfTypeConversion()  {
      nativeType_ = H5::StrType( H5::PredType::C_S1, H5T_VARIABLE);
      stdType_ = H5::StrType( H5::PredType::C_S1, H5T_VARIABLE);
    }
    
    const void * GetVoidPtr() {
      if( !isSet_ ) {
        EXCEPTION( "Data buffer is empty" );
      }
      return &buffer_[0];
    }
    
    void SetData( const std::string& t ) {
      buffer_.resize(1);
      buffer_[0] = t.c_str();
      size_ = sizeof( char* );
      isSet_ = true;
    }
    
    void SetData( const std::string* t, UInt size  ) {
      buffer_.resize( size );
      for( UInt i = 0; i < size; i ++ ) {
        buffer_[i] = t[i].c_str();
      }
      size_ = sizeof( char* );
      isSet_ = true;
    }
  };
  
  
  // --------------------------
  //  std::vector<std::string>
  // --------------------------
  template<>
  class H5IO::HdfTypeConversion<std::vector<std::string> > :
    public H5IO::BaseHdfTypeConversion {

  private:
    std::vector<hvl_t> buffer_;

  public:
    HdfTypeConversion() { 
      H5::StrType sType( H5::PredType::C_S1, H5T_VARIABLE);
      H5::VarLenType vType( &sType );
      nativeType_ = vType;
      stdType_ = vType;
    }
    
    const void * GetVoidPtr() {
      if( !isSet_ ) {
        EXCEPTION( "Data buffer is empty" );
      }
      return &(buffer_[0]);
    }
    
    void SetData( const std::vector<std::string>& t ) {
      CleanUp();
      buffer_.resize( 1 );
      buffer_[0].p = (void*) new const char*[t.size()];
      buffer_[0].len = t.size();
      for( UInt j = 0; j < t.size(); j++ ) {
        ((const char **)buffer_[0].p)[j] = t[j].c_str();
      }
      size_ = sizeof( hvl_t );
      isSet_ = true;
    }

    
    void SetData( const std::vector<std::string>* t, UInt size ) {
      CleanUp();
      buffer_.resize( size );
      for( UInt i = 0; i < size; i++ ) {
        buffer_[i].p = (void*) new const char*[t[i].size()];
        buffer_[i].len = t[i].size();
        for( UInt j = 0; j < t[i].size(); j++ ) {
          ((const char **)buffer_[i].p)[j] = ((t[i])[j]).c_str();
          size_ += t[i][j].size() * sizeof( char );
        }
      }
      isSet_ = true;
    }
 
    void CleanUp() {
      
      // delete buffer for arrays of vectors
      for( UInt i = 0; i < buffer_.size(); i++ ) {
        delete[] (const char*) buffer_[i].p;
      }
      buffer_.clear();
      
      isSet_ = false;
    }
    
    ~HdfTypeConversion() {
      CleanUp();
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
    std::vector<hvl_t> ptBuffer_;                                       \
                                                                        \
  public:                                                               \
    HdfTypeConversion() {                                               \
      nativeType_ =                                                     \
        H5::VarLenType ( &HdfAtomTypeMap<TYPE>::HdfNativeType );        \
      stdType_ =                                                        \
        H5::VarLenType ( &HdfAtomTypeMap<TYPE>::HdfStdType );           \
    }                                                                   \
                                                                        \
    const void * GetVoidPtr( ) {                                        \
      if( !isSet_ ) {                                                   \
        EXCEPTION( "Data buffer is empty" );                            \
      }                                                                 \
      return (const void*) &ptBuffer_[0];                               \
    }                                                                   \
                                                                        \
    void SetData( const std::vector<TYPE>& t ) {                        \
      CleanUp();                                                        \
      ptBuffer_.resize( 1 );                                            \
      ptBuffer_[0].p =  (void*) new TYPE[t.size()];                     \
      ptBuffer_[0].len = t.size();                                      \
      for( UInt j = 0; j < t.size(); j++ ) {                            \
        ((TYPE*)ptBuffer_[0].p)[j] = t[j];                              \
      }                                                                 \
      size_ = sizeof( hvl_t );                                          \
      isSet_ = true;                                                    \
    }                                                                   \
                                                                        \
    void SetData( const std::vector<TYPE>* t, UInt size ) {             \
      CleanUp();                                                        \
      ptBuffer_.resize( size );                                         \
      for( UInt i = 0; i < size; i++ ) {                                \
        ptBuffer_[i].p =  new TYPE[t[i].size()];                        \
        ptBuffer_[i].len = t[i].size();                                 \
        for( UInt j = 0; j < t[i].size(); j++ ) {                       \
          ((TYPE*)ptBuffer_[i].p)[j] = t[i][j];                         \
        }                                                               \
      }                                                                 \
      size_ = sizeof( hvl_t );                                          \
      isSet_ = true;                                                    \
    }                                                                   \
                                                                        \
    void CleanUp() {                                                    \
      for( UInt i = 0; i < ptBuffer_.size(); i++ ) {                    \
        delete[] (TYPE*) ptBuffer_[i].p;                                \
      }                                                                 \
      ptBuffer_.clear();                                                \
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
      conv.SetData( data );
      if( !conv.IsSet() ) {
        EXCEPTION( "Could not convert data for attribute '"
                   << name << "' of type " << typeid(TYPE).name() );
      }
      H5::Attribute attr = obj.createAttribute( name, stdType, 
                                                space, create_plist );
      
      // write attribute
      attr.write( nativeType, conv.GetVoidPtr() );
      
      // reset conversion object
      conv.CleanUp();
      
      // close attribute and dataspace
      attr.close();
      space.close();

    }  catch (H5::Exception& h5ex) {
      EXCEPTION("Could not write attribute '" 
                << name << "': " << h5ex.getCDetailMsg());
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
      conv.SetData( buffer, size );
      if( !conv.IsSet() ) {
        EXCEPTION( "Could not convert data for 1D array '"
                   << name << "' of type " << typeid(TYPE).name() );
      }
      H5::DataSet dataset = loc.createDataSet( name, stdType, 
                                               space, create_plist );
      dataset.write( conv.GetVoidPtr(), nativeType  );

      // reset conversion object
      conv.CleanUp();

      // close dataset and space
      dataset.close();
      space.close();
      
    } catch (H5::Exception& h5ex) {
      EXCEPTION("Could not write 1D-Array '" 
                << name << "': " << h5ex.getCDetailMsg());
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
      conv.SetData( buffer, rowSize * colSize );
      if( !conv.IsSet() ) {
        EXCEPTION( "Could not convert data for 2D array '"
                   << name << "' of type " << typeid(TYPE).name() );
      }
      H5::DataSet dataset = loc.createDataSet( name, stdType, 
                                               space, create_plist );
      dataset.write( conv.GetVoidPtr(), nativeType  );
      
    }  catch (H5::Exception& h5ex) {
      EXCEPTION("Could not write 2D-Array '" << name 
                << "': " << h5ex.getCDetailMsg());
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
        totalSize += conv[i]->GetSize();
      
        // create new compound for memory datatype
        memCompTypes[i] = H5::CompType( (size_t)conv[i]->GetSize() );
        memCompTypes[i].insertMember( memName, 0, conv[i]->GetNativeType() );
      }

      // Create file compound data type
      H5::CompType fileCompType( (size_t) totalSize );
      UInt actOffset = 0;
      for( UInt i = 0; i < comp.size(); i++ ) {
        fileCompType.insertMember( comp[i].first,
                                   actOffset,
                                   conv[i]->GetStdType() );
        actOffset += conv[i]->GetSize();
      }

      // create  data space
      hsize_t dims = 1;
      H5::DataSpace space( 1, &dims );
    
      // generate dataset
      H5::DataSet dataset = loc.createDataSet( name, fileCompType, 
                                               space, create_plist );

      // iterate again over all entries and fill in values
      for( UInt i = 0; i < comp.size(); i++ ) {
        dataset.write( conv[i]->GetVoidPtr(), memCompTypes[i] );
        conv[i]->CleanUp();
      }

      // close datatset and dataspace
      dataset.close();
      space.close();

    }  catch (H5::Exception& h5ex) {
      EXCEPTION("Could not write compound '" << name 
                << "': " << h5ex.getCDetailMsg());
    } catch( Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not write compound '" << name << "'" );
    }
  }
  
  
  void H5IO::GetAnyConversion( const boost::any& anyType,
                               shared_ptr<BaseHdfTypeConversion>& conv ) {
    
    // query type of any
#define ANY_CONVERSION( TYPE )                          \
    if( anyType.type() == typeid(TYPE) ){               \
      shared_ptr<HdfTypeConversion<TYPE> >              \
        myConv ( new HdfTypeConversion<TYPE>() );       \
      myConv->SetData( any_cast<TYPE>(anyType) );       \
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
                                   &create_plist )
  
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


  FEType H5IO::XMDFElemType2ElemType( const Integer type ) {
    switch (type) {
    case 100:   return ET_LINE2;   // 1D linear
    case 101:   return ET_LINE3;   // 1D quadratic
    case 200:   return ET_TRIA3;   // 2D linear triangle
    case 201:   return ET_TRIA6;   // 2D quadratic triangle
    case 210:   return ET_QUAD4;   // 2D linear quadrilateral
    case 211:   return ET_QUAD8;   // 2D quadratic quadrilateral
    case 212:   return ET_QUAD9;   // 2D quadr. quad. with center node
    case 300:   return ET_TET4;    // 3D linear tetrahedron
    case 30010: return ET_TET10;   // 3D quadratic tetrahedron
    case 310:   return ET_WEDGE6;  // 3D linear prism
    case 31010: return ET_WEDGE15; // 3D quadratic prism
    case 320:   return ET_HEXA8;   // 3D linear hexahedron
    case 32010: return ET_HEXA20;  // 3D quadratic hexahedron
    case 32011: return ET_HEXA27;  // 3D quadr. hexa. with center nodes
    case 330:   return ET_PYRA5;   // 3D linear pyramid
    case 33010: return ET_PYRA13;  // 3D quadratic pyramid
    }

    return ET_UNDEF;
    // This place should never be reached!
  }

  Integer H5IO::ElemType2XMDFElemType( const FEType type ) {
    switch (type) {
    case ET_LINE2:   return 100;   // 1D linear
    case ET_LINE3:   return 101;   // 1D quadratic
    case ET_TRIA3:   return 200;   // 2D linear triangle
    case ET_TRIA6:   return 201;   // 2D quadratic triangle
    case ET_QUAD4:   return 210;   // 2D linear quadrilateral
    case ET_QUAD8:   return 211;   // 2D quadratic quadrilateral
    case ET_QUAD9:   return 212;   // 2D quadr. quad. with center node
    case ET_TET4:    return 300;   // 3D linear tetrahedron
    case ET_TET10:   return 30010; // 3D quadratic tetrahedron
    case ET_WEDGE6:  return 310;   // 3D linear prism
    case ET_WEDGE15: return 31010; // 3D quadratic prism
    case ET_HEXA8:   return 320;   // 3D linear hexahedron
    case ET_HEXA20:  return 32010; // 3D quadratic hexahedron
    case ET_HEXA27:  return 32011; // 3D quadr. hexa. with center nodes
    case ET_PYRA5:   return 330;   // 3D linear pyramid
    case ET_PYRA13:  return 33010; // 3D quadratic pyramid
    default: return ET_UNDEF;
    }

    // This place should never be reached!
    return -1;
  }

  void H5IO::ReorderConnectivity( const Integer eType,
                                  const bool toXMDF,
                                  const UInt* in,
                                  UInt* out) {
    static Integer toXMDFIdxsLine[] = {0, 2, 1};
    static Integer fromXMDFIdxsLine[] = {0, 2, 1};
    static Integer toXMDFIdxsTria[] = {0, 3, 1, 4, 2, 5};
    static Integer fromXMDFIdxsTria[] = {0, 2, 4, 1, 3, 5};
    static Integer toXMDFIdxsQuad[] = {0, 4, 1, 5, 2, 6, 3, 7, 8};
    static Integer fromXMDFIdxsQuad[] = {0, 2, 4, 6, 1, 3, 5, 7, 8};
    
    std::vector<UInt> tmp;
    Integer* it;

    UInt numNodes = NUM_ELEM_NODES[eType];
    tmp.resize(numNodes);
    memcpy(&tmp[0], in, numNodes*sizeof(UInt));

    switch (eType) {
    case ET_LINE3: // 1D quadratic line
      if(toXMDF)
        it = toXMDFIdxsLine;
      else
        it = fromXMDFIdxsLine;
      break;
    case ET_TRIA6: // 2D quadratic triangle
      if(toXMDF)
        it = toXMDFIdxsTria;
      else
        it = fromXMDFIdxsTria;
      break;
    case ET_QUAD8: // 2D quadratic quadrilateral
    case ET_QUAD9: // 2D quadratic quadrilateral with center node
      if(toXMDF)
        it = toXMDFIdxsQuad;
      else
        it = fromXMDFIdxsQuad;
      break;
    }

    for(UInt i=0; i<numNodes; i++, it++)
    {
      out[i] = tmp[*it];
    }
  }




} // end of namespace CoupledField
