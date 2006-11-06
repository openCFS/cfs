#include "entityList.hh"
#include "Domain/grid.hh"


namespace CoupledField {


  EntityList::EntityList( Grid* grid ) {
    grid_ = grid;
    type_ = NO_TYPE;
    size_ = 0;
  }

  EntityList::~EntityList() {
  }

  std::string EntityList::TypeToString( Type type ) {
    
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
      *error << "Type '" << type << "' describes no entityList!";
      Error( __FILE__, __LINE__ );
    }

    return ret;
  }

  EntityList::Type EntityList::StringToType( const std::string& type ) {

    Type ret;
    
    if( type == "elemList" ) {
      ret = ELEM_LIST;
    } else if( type == "surfElemList" ) {
      ret = SURF_ELEM_LIST;
    } else if( type == "ncElemList" ) {
      ret = NC_ELEM_LIST;
    } else if( type == "nodeList" ) {
      ret =  NODE_LIST;
    } else if( type == "regionList" ) {
      ret = REGION_LIST;
    } else if( type == "numberList" ) {
      ret = NUMBER_LIST;
    } else {
      *error << "Type '" << type << "' describes no EntityList.";
      Error( __FILE__, __LINE__ );
    }

    return ret;
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
    
    return grid_->RegionIdToName( region_ );
  }

  void ElemList::SetRegion( RegionIdType region ) {
    ENTER_FCN( "ElemList::SetRegions" );
    
    region_ = region;
    StdVector<Elem*> elems;
    grid_->GetElems( elems, region );
    for ( UInt i = 0; i<elems.GetSize(); i++ ) {
      list_.Push_back( elems[i]->elemNum);
    }
    size_ = list_.GetSize();
  }

  void ElemList::SetElement( const Elem* elem) {
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
    return grid_->RegionIdToName( region_ );
  }

  void SurfElemList::SetRegion( RegionIdType  region ) { 
    region_ = region;
    StdVector<SurfElem*> elems;
    grid_->GetSurfElems( elems, region );
    for ( UInt i = 0; i<elems.GetSize(); i++ ) {
      list_.Push_back( elems[i]);
    }
    size_ = list_.GetSize();
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
    name_ = name;
    grid_->GetNodesByName( list_, name );
    size_ = list_.GetSize();
  }

  void NodeList::SetNodesOfRegion( RegionIdType regionId ) {
    name_ = grid_->RegionIdToName( regionId );
    grid_->GetNodesByRegion( list_, regionId );
    size_ = list_.GetSize();
  }
  
  void NodeList::SetNodes( const StdVector<UInt>& nodeList ) {
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
    return std::string( "" );
  }

  EntityIterator RegionList::GetIterator() const {
    EntityIterator it;
    Error( "Not implemented yet", __FILE__, __LINE__ );
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
    Error( "Not implemented yet", __FILE__, __LINE__ );
    return it;
  }

  // =================================================
  // E N T I T Y   I T E R A T O R
  // =================================================
  
  EntityIterator::EntityIterator() {

  }
  

  void EntityIterator::Begin() {
    ENTER_FCN( "EntityIterator::Begin" );
    pos_ = 0;
  }

  bool EntityIterator::IsEnd() const {
    ENTER_FCN( "EntityIterator::IsEnd" );
    return (pos_ == size_);
  }
  
  EntityIterator& EntityIterator::operator++(int) {
    ENTER_FCN( "EntityIterator::operator++" );
    pos_++;
    return *this;
  }
  
  const Elem* EntityIterator::GetElem() const{
    ENTER_FCN( "EntityIterator:: GetElem" );
    //std::cerr << "elemList_ = " << elemList_ << std::endl;
    if( type_ == EntityList::ELEM_LIST ) {
      return elemList_->GetElem( pos_ );
    } else if( type_ == EntityList::SURF_ELEM_LIST ) {
      return surfElemList_->GetSurfElem( pos_ );
    }
  }
  
  const SurfElem* EntityIterator::GetSurfElem() const {
    ENTER_FCN( "EntityIterator::GetSurfElem" );
    return surfElemList_->list_[ pos_ ];
  }
  
  RegionIdType EntityIterator::GetRegion() const {
    ENTER_FCN( "EntityIterator::GetRegion" );
    Error( "Not implemented", __FILE__, __LINE__ );
    return NO_REGION_ID;
  }
  
  UInt EntityIterator::GetNode() const {
    ENTER_FCN( "EntityIterator::GetNode" );
    return nodeList_->list_[ pos_ ];
  }
  
  UInt EntityIterator::GetUnknown() const {
    ENTER_FCN( "EntityIterator::GetUnknown" );
    Error( "Not Implemented", __FILE__, __LINE__ );
    return 0;
  }


}
