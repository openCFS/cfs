// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/lexical_cast.hpp>

#include "entityList.hh"
#include "Domain/grid.hh"



namespace CoupledField {


  EntityList::EntityList( Grid* grid ) {
    grid_ = grid;
    type_ = NO_LIST;
    defineType_ = NO_TYPE;
    size_ = 0;
  }

  EntityList::~EntityList() {
  }

  void EntityList::Enum2String( ListType type, std::string& out ) {
    
    std::string ret;
    
    switch( type ) {

    case ELEM_LIST:
      ret = "elemList";
      break;
      
    case SURF_ELEM_LIST:
      ret = "surfElemList";
      break;
    case NC_ELEM_LIST:
      ret = "ncElemList";
      break;

    case NODE_LIST:
      ret = "nodeList";
      break;
      
    case REGION_LIST:
      ret = "regionList";
      break;
      
    case NUMBER_LIST:
      ret = "numberList";
      break;
      
    default:
      EXCEPTION( "Type '" << type << "' describes no entityList!" );
    }

    out = ret;
  }

  void EntityList::String2Enum( const std::string& type,
                                EntityList::ListType & out ) {

    if( type == "elemList" ) {
      out = ELEM_LIST;
    } else if( type == "surfElemList" ||
               type == "surfRegion") {
      out = SURF_ELEM_LIST;
    } else if( type == "ncElemList" ) {
      out = NC_ELEM_LIST;
    } else if( type == "nodeList" ) {
      out =  NODE_LIST;
    } else if( type == "regionList" ||
               type == "region") {
      out = REGION_LIST;
    } else if( type == "numberList" ) {
      out = NUMBER_LIST;
    } else {
      EXCEPTION( "Type '" << type << "' describes no EntityList." );
    }
  }
  
  void EntityList::Enum2String( DefineType type, std::string& out ) {
  }
  

  void EntityList::String2Enum( const std::string& type,
                                EntityList::DefineType & out ) {
    
  }

  // --- Elem List ---
  ElemList::ElemList( Grid* grid ) 
    : EntityList( grid ) {

    type_ = ELEM_LIST;
    region_ = NO_REGION_ID;
  }
  ElemList::~ElemList() {
  }

  std::string ElemList::GetName() const {
    
    if( defineType_ == NO_TYPE ) {
      return std::string("*anoymous list*");
    } else {
      return name_;
    }
  }

  void ElemList::SetNamedElems( const std::string& name ) {
    defineType_ = NAMED_ELEMS;
    name_ = name;
    StdVector<Elem*> elems;
    grid_->GetElemsByName( elems, name );
    for ( UInt i = 0; i<elems.GetSize(); i++ ) {
      list_.Push_back( elems[i]->elemNum);
    }
    size_ = list_.GetSize();
    name_ = name;
  }
      

  void ElemList::SetRegion( RegionIdType region ) {
    
    defineType_ = REGION;
    region_ = region;
    StdVector<Elem*> elems;
    grid_->GetElems( elems, region );
    for ( UInt i = 0; i<elems.GetSize(); i++ ) {
      list_.Push_back( elems[i]->elemNum);
    }
    size_ = list_.GetSize();
    name_ = grid_->RegionIdToName( region );
  }

  void ElemList::SetElement( const Elem* elem) {

    defineType_ = NO_TYPE;
    region_ = NO_REGION_ID;
    list_.Clear();
    list_.Push_back( elem->elemNum );
  }
  
  RegionIdType ElemList::GetRegion() const {
    return region_;
  }

  const Elem * ElemList::GetElem( UInt nr ) const {
    return grid_->GetElem( list_[nr] );
  }

  EntityIterator ElemList::GetIterator() const{
    EntityIterator it;
    it.type_ = ELEM_LIST;
    it.elemList_ = this;
    it.pos_ = 0;
    it.size_ = list_.GetSize();
    return it;
  }
 


  // --- SurfElem List ---
  SurfElemList::SurfElemList( Grid* grid )
    : EntityList( grid ) {

    type_ = SURF_ELEM_LIST;
    region_ = NO_REGION_ID;
  }

  std::string SurfElemList::GetName() const {
    return name_;
  }

  void SurfElemList::SetNamedElems( const std::string& name ) {
    defineType_ = NAMED_ELEMS;
    name_ = name;
    StdVector<Elem*> elems;
    grid_->GetElemsByName( elems, name );
    SurfElem * actElem = NULL;
    for ( UInt i = 0; i<elems.GetSize(); i++ ) {
      actElem = dynamic_cast<SurfElem*>(elems[i]);
      if( actElem == NULL ) {
        EXCEPTION( "Element #" << elems[i]->elemNum 
                   << " in element list '" << name 
                   << "' is no surface element!" );
      }
      list_.Push_back( actElem );
    }
    size_ = list_.GetSize();
    name_ = name;
  }

  void SurfElemList::SetRegion( RegionIdType  region ) { 

    defineType_ = REGION;
    region_ = region;
    StdVector<SurfElem*> elems;
    grid_->GetSurfElems( elems, region );
    for ( UInt i = 0; i<elems.GetSize(); i++ ) {
      list_.Push_back( elems[i]);
    }
    size_ = list_.GetSize();
    name_ = grid_->RegionIdToName( region );
  }

    //! Get RegionId
  RegionIdType SurfElemList::GetRegion() const {
    return region_;
  }

  const SurfElem* SurfElemList::GetSurfElem( UInt nr ) const{
    return list_[nr];
  }

  EntityIterator SurfElemList::GetIterator() const {
    EntityIterator it;
    it.type_ = SURF_ELEM_LIST;
    it.surfElemList_ = this;
    it.pos_ = 0;
    it.size_ = list_.GetSize();
    return it;
  }
 
  // --- Node List ---
  NodeList::NodeList( Grid* grid ) 
    : EntityList( grid) {
    type_ = NODE_LIST;
  }

  std::string NodeList::GetName() const {
    if( name_ == std::string() ) {
      std::stringstream out;
      for (UInt i = 0; i < list_.GetSize()-1; i++ ) {
        out << list_[i] << ", ";
      }
      if ( list_.GetSize() > 0 ) {
        out << list_.Last();
      }
      return out.str();
    } else {
      return name_;
    }
  }
  
  void NodeList::SetNamedNodes( const std::string & name) {
    defineType_ = NAMED_NODES;
    name_ = name;
    grid_->GetNodesByName( list_, name );
    size_ = list_.GetSize();
  }

  void NodeList::SetNodesOfRegion( RegionIdType regionId ) {
    defineType_ = REGION;
    name_ = grid_->RegionIdToName( regionId );
    grid_->GetNodesByRegion( list_, regionId );
    size_ = list_.GetSize();
  }
  
  void NodeList::SetNodes( const StdVector<UInt>& nodeList ) {
    defineType_ = NO_TYPE;
    name_ = "";
    list_ = nodeList;
    size_ = list_.GetSize();

  }


  EntityIterator NodeList::GetIterator() const {
    EntityIterator it;
    it.type_ = NODE_LIST;
    it.nodeList_ = this;
    it.pos_ = 0;
    it.size_ = list_.GetSize();
    return it;
  }

  // --- Region List ---
  RegionList::RegionList( Grid* grid )
    : EntityList( grid ) {
    type_ = REGION_LIST;
  }

  std::string RegionList::GetName() const {
//    std::string ret = "regionList( ";
//    for( UInt i=0; i < list_.GetSize()-1; i++ ) { 
//      ret += grid_->RegionIdToName( list_[i] ) + ", ";
//    }
//    if( list_.GetSize() > 0 ) {
//      ret+= grid_->RegionIdToName( list_.Last() );
//    }
//    return ret + ")";
    
    return grid_->RegionIdToName( list_.Last() );
  }

  void RegionList::SetRegionId( RegionIdType region ) {
    defineType_ = REGION;
    list_.Resize(1);
    list_[0] = region;
    size_ = 1;
  }

  void RegionList::SetRegionNames( const StdVector<RegionIdType>& regions ) {
    defineType_ = REGION;
    list_ = regions;
    size_ = list_.GetSize();
  }

  
  EntityIterator RegionList::GetIterator() const {
    EntityIterator it;
    it.type_ = REGION_LIST;
    it.regionList_ = this;
    it.pos_ = 0;
    it.size_ = list_.GetSize();
    return it;
  }

  // --- Number List ---
  NumberList::NumberList( Grid* grid )
    : EntityList( grid ) {
    type_ = NUMBER_LIST;
  }

  std::string NumberList::GetName() const {
    return std::string( "" );
  }

  EntityIterator NumberList::GetIterator() const {
    EntityIterator it;
    EXCEPTION( "Not implemented yet");
    return it;
  }

  // =================================================
  // E N T I T Y   I T E R A T O R
  // =================================================
  
  EntityIterator::EntityIterator() {

  }
  

  void EntityIterator::Begin() {
    pos_ = 0;
  }

  bool EntityIterator::IsEnd() const {
    return (pos_ == size_);
  }
  
  EntityIterator& EntityIterator::operator++(int) {
    pos_++;
    return *this;
  }
  
  const Elem* EntityIterator::GetElem() const{
    switch(type_)
    {
    case EntityList::ELEM_LIST:
      return elemList_->GetElem( pos_ );
      
    case EntityList::SURF_ELEM_LIST:
      return surfElemList_->GetSurfElem( pos_ );
      
    default: EXCEPTION("type " << type_ << " not implemented");
    }
    return NULL;
  }
  
  const SurfElem* EntityIterator::GetSurfElem() const {
    return surfElemList_->list_[ pos_ ];
  }
  
  RegionIdType EntityIterator::GetRegion() const {
    return regionList_->list_[ pos_ ];
  }
  
  UInt EntityIterator::GetNode() const {
    return nodeList_->list_[ pos_ ];
  }
  
  UInt EntityIterator::GetUnknown() const {
    EXCEPTION( "Not Implemented" );
    return 0;
  }
  
  std::string EntityIterator::GetIdString() const {
    
    std::string id = "";
    switch( type_ ) {
    case EntityList::ELEM_LIST:
    case EntityList::SURF_ELEM_LIST:
    case EntityList::NC_ELEM_LIST:
      id = lexical_cast<std::string>(GetElem()->elemNum);
      break;
    case EntityList::NODE_LIST:
      id = lexical_cast<std::string>(GetNode());
      break;
    case EntityList::REGION_LIST:
      id = regionList_->grid_->RegionIdToName( GetRegion() );
      break;
    default:
        EXCEPTION( "Not implemented" );
    }
      
     return id; 
  }


} // end of namespace
