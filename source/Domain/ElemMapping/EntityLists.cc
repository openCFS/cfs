
#include "EntityLists.hh"
#include "Domain/Mesh/Grid.hh"


namespace CoupledField {

  Enum<EntityList::ListType>   EntityList::listType;
  Enum<EntityList::DefineType> EntityList::defineType;

  EntityList::EntityList( Grid* grid ) {
    grid_ = grid;
    type_ = NO_LIST;
    defineType_ = NO_TYPE;
    size_ = 0;
  }

  EntityList::~EntityList() {
  }

  void EntityList::SetEnums()
  {
    EntityList::listType.SetName("EntityList::ListType");
    EntityList::listType.Add(EntityList::ELEM_LIST, "elemList");
    EntityList::listType.Add(EntityList::SURF_ELEM_LIST, "surfElemList", false);
    EntityList::listType.Add(EntityList::SURF_ELEM_LIST, "surfRegion", false);
    EntityList::listType.Add(EntityList::NC_ELEM_LIST, "ncElemList");
    EntityList::listType.Add(EntityList::NODE_LIST, "nodeList");
    EntityList::listType.Add(EntityList::NAME_LIST, "nameList" );
    EntityList::listType.Add(EntityList::REGION_LIST, "regionList", false);
    EntityList::listType.Add(EntityList::REGION_LIST, "region", false);
    EntityList::listType.Add(EntityList::NUMBER_LIST, "numberList");
    EntityList::listType.Add(EntityList::COIL_LIST, "coilList");

    EntityList::defineType.SetName("EntityList::DefineType");
    EntityList::defineType.Add(EntityList::NO_TYPE, "no_type");
    EntityList::defineType.Add(EntityList::REGION, "region");
    EntityList::defineType.Add(EntityList::NAMED_NODES, "named_nodes");
    EntityList::defineType.Add(EntityList::NAMED_ELEMS, "named_elems");


  }
  
  void EntityList::Intersect(const StdVector<shared_ptr<EntityList> >& set1,
                             const StdVector<shared_ptr<EntityList> >& set2,
                             StdVector<shared_ptr<EntityList> >& intersect ) {

    for( UInt i = 0; i < set1.GetSize(); ++i ) {
      std::string name1 = set1[i]->GetName();
      for( UInt j = 0; j < set2.GetSize(); ++j ) {
        if( set2[j]->GetName() == name1 ) {
          intersect.Push_back( set1[i]);
          break;
        }
      }
    }
    intersect.Trim();
  }

  void EntityList::Union(const StdVector<shared_ptr<EntityList> >& set1,
                         const StdVector<shared_ptr<EntityList> >& set2,
                         StdVector<shared_ptr<EntityList> >& unionSet) {
    // grab entries of set1
    for( UInt i = 0; i < set1.GetSize(); ++i ) {
      bool found = false;
      std::string name1 = set1[i]->GetName();
      for( UInt k = 0; k < unionSet.GetSize(); ++k ) {
        if( unionSet[k]->GetName() == name1 ) {
          found = true;
          break;
        }
      }
      if(!found)
        unionSet.Push_back(set1[i]);
    }
    // grab entries of set2
    for( UInt j = 0; j < set2.GetSize(); ++j ) {
      bool found = false;
      std::string name2 = set1[j]->GetName();
      for( UInt k = 0; k < unionSet.GetSize(); ++k ) {
        if( unionSet[k]->GetName() == name2 ) {
          found = true;
          break;
        }
      }
      if(!found)
        unionSet.Push_back(set2[j]);
    }
    unionSet.Trim();
  }

  std::string EntityList::ToString() const
  {
    std::stringstream ss;
    ss << "type=" << listType.ToString(GetType());
    ss << " dt=" << defineType.ToString(GetDefineType());
    ss << " name=" << GetName();
    ss << " reg=" << GetRegion();
    ss << " size=" << GetSize();
    return ss.str();
  }



  // --- Elem List ---
  ElemList::ElemList( Grid* grid ) 
    : EntityList( grid ) {

    type_ = ELEM_LIST;
    region_ = NO_REGION_ID;
  }
  
  ElemList::ElemList(  const Elem* elem, Grid* grid ) 
      : EntityList( grid ) {

      type_ = ELEM_LIST;
      region_ = NO_REGION_ID;
      defineType_ = NO_TYPE;
      list_.Push_back( elem->elemNum );
  }
  
  
  ElemList::~ElemList() {
  }

  std::string ElemList::GetName() const {
    
    if( defineType_ == NO_TYPE ) {
      return std::string("*anonymous list*");
    } else {
      return name_;
    }
  }

  void ElemList::SetNamedElems( const std::string& name ) {
    StdVector<Elem*> elems;
    grid_->GetElemsByName( elems, name );
    UInt numElems = elems.GetSize();
    
    list_.Clear();
    list_.Reserve(numElems);
    
    for ( UInt i=0; i<numElems; ++i ) {
      list_.Push_back(elems[i]->elemNum);
    }
    
    defineType_ = NAMED_ELEMS;
    size_ = list_.GetSize();
    region_ = NO_REGION_ID;
    name_ = name;
  }
      

  void ElemList::SetRegion( RegionIdType region ) {
    StdVector<Elem*> elems;
    grid_->GetElems( elems, region );
    UInt numElems=elems.GetSize();
    
    list_.Clear();
    list_.Reserve(numElems);
    
    for ( UInt i=0; i<numElems; ++i ) {
      list_.Push_back( elems[i]->elemNum);
    }

    defineType_ = REGION;
    size_ = list_.GetSize();
    region_ = region;
    name_ = grid_->GetRegion().ToString( region );
  }

  void ElemList::SetElement( const Elem* elem) {
    defineType_ = NO_TYPE;
    size_ = 1;
    region_ = NO_REGION_ID;
    name_ = "";
    list_.Resize(1, elem->elemNum);
  }
  
  RegionIdType ElemList::GetRegion() const {
    return region_;
  }

  const Elem* ElemList::GetElem( UInt nr ) const {
    return grid_->GetElem( list_[nr] );
  }

  EntityIterator ElemList::GetIterator() const{
    EntityIterator it;
    it.type_ = ELEM_LIST;
    it.ptGrid_ = this->grid_;
    it.elemList_ = this;
    it.pos_ = 0;
    it.size_ = list_.GetSize();
    return it;
  }
 
  //! Add an element to the list
  void ElemList::AddElement( const Elem* elem ) {
#pragma omp critical (ENTLISTS_ELEMLIST)
{
    list_.Push_back(elem->elemNum);
    ++size_;
}
  }

  void ElemList::AddElements( const StdVector<Elem*>& elems ){
    list_.Reserve( list_.GetSize() + elems.GetSize() );
    for( unsigned int k = 0; k < elems.GetSize(); ++k ){
      list_.Push_back(elems[k]->elemNum);
    }
    size_ = list_.GetSize();
  }



  // --- SurfElem List ---
  SurfElemList::SurfElemList( Grid* grid )
    : ElemList( grid ) {

    type_ = SURF_ELEM_LIST;
    region_ = NO_REGION_ID;
  }

  void SurfElemList::SetNamedElems( const std::string& name ) {
    StdVector<Elem*> elems;
    grid_->GetElemsByName( elems, name );
    
    surfElemList_.Clear();
    
    SurfElem * actElem = nullptr;
    for ( UInt i=0, numElems=elems.GetSize(); i<numElems; ++i ) {
      actElem = dynamic_cast<SurfElem*>(elems[i]);
      if( actElem == nullptr) {
        EXCEPTION( "Element #" << elems[i]->elemNum 
                   << " in element list '" << name 
                   << "' is no surface element!" );
      }
      surfElemList_.Push_back( actElem );
    }
    
    defineType_ = NAMED_ELEMS;
    size_ = surfElemList_.GetSize();
    region_ = NO_REGION_ID;
    name_ = name;
  }

  void SurfElemList::SetRegion( RegionIdType  region ) { 
    StdVector<SurfElem*> elems;
    grid_->GetSurfElems( elems, region );
    
    surfElemList_.Clear();
    
    for ( UInt i=0, numElems=elems.GetSize(); i<numElems; ++i ) {
      surfElemList_.Push_back( elems[i]);
    }

    defineType_ = REGION;
    size_ = surfElemList_.GetSize();
    region_ = region;
    name_ = grid_->GetRegion().ToString( region );
  }

  void SurfElemList::SetElement(const Elem * elem) {
    const SurfElem * sElem = dynamic_cast<const SurfElem*>(elem);
    if ( !sElem ) {
      EXCEPTION("Cannot use SurfElemList with volume elements");
    }
    SetSurfElem(sElem);
  }

  void SurfElemList::SetSurfElem(const SurfElem * elem) {
    defineType_ = NO_TYPE;
    size_ = 1;
    region_ = NO_REGION_ID;
    name_ = "";
    surfElemList_.Resize(1, elem);
  }
  
  void SurfElemList::AddElement(const SurfElem* elem) {
#pragma omp critical (ENTLISTS_SURFELEMLIST)
{
    surfElemList_.Push_back(elem);
    ++size_;
}
  }
  
  const Elem* SurfElemList::GetElem(UInt nr) const {
    return surfElemList_[nr];
  }
  
  const SurfElem* SurfElemList::GetSurfElem( UInt nr ) const {
    return surfElemList_[nr];
  }

  EntityIterator SurfElemList::GetIterator() const {
    EntityIterator it;
    it.type_ = SURF_ELEM_LIST;
    it.ptGrid_ = this->grid_;
    it.surfElemList_ = this;
    it.pos_ = 0;
    it.size_ = surfElemList_.GetSize();
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
    name_ = grid_->GetRegion().ToString( regionId );
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
    it.ptGrid_ = this->grid_;
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
//      ret += grid_->GetRegion().ToString( list_[i] ) + ", ";
//    }
//    if( list_.GetSize() > 0 ) {
//      ret+= grid_->GetRegion().ToString( list_.Last() );
//    }
//    return ret + ")";
    
    return grid_->GetRegion().ToString( list_.Last() );
  }

  void RegionList::SetRegion( RegionIdType region ) {
    defineType_ = REGION;
    list_.Resize(1);
    list_[0] = region;
    size_ = 1;
  }

  void RegionList::SetRegions( const StdVector<RegionIdType>& regions ) {
    defineType_ = REGION;
    list_ = regions;
    size_ = list_.GetSize();
  }

  
  EntityIterator RegionList::GetIterator() const {
    EntityIterator it;
    it.type_ = REGION_LIST;
    it.ptGrid_ = this->grid_;
    it.regionList_ = this;
    it.pos_ = 0;
    it.size_ = list_.GetSize();
    return it;
  }
  
  // --- Name List ---
   NameList::NameList( Grid* grid )
   : EntityList( grid ) {
     type_ = NAME_LIST;
     size_ = 0;
   }

   std::string NameList::GetName() const {
     if( list_.GetSize() == 0  ) {
       return "";
     } else {
       return list_[0];
     }
   }

   void NameList::SetName( const std::string& name ) {
     list_.Clear();
     list_.Resize( 1 );
     list_[0] = name;
     size_ = 1;
   }
   void NameList::SetNames( const StdVector<std::string>& names ) {
     list_ = names;
     size_ = list_.GetSize();
   }

   EntityIterator NameList::GetIterator() const {
     EntityIterator it;
     it.type_ = NAME_LIST;
     it.ptGrid_ = this->grid_;
     it.nameList_ = this;
     it.pos_ = 0;
     it.size_ = list_.GetSize();
     return it;
   }
  
   
   // --- CoilList ---
   CoilList::CoilList( Grid* grid )
   : EntityList( grid ) {
     type_ = COIL_LIST;
   }

   std::string CoilList::GetName() const {
     std::string id = "";
     if( size_ > 0 ) {
       id = list_[0]->coilId_;
     }
     return id;
   }

   void CoilList::AddCoil( shared_ptr<Coil> aCoil ) {
     defineType_ = NO_TYPE;
     list_.Push_back( aCoil );
     size_ = list_.GetSize();
   }

   EntityIterator CoilList::GetIterator() const {
     EntityIterator it;
     it.type_ = COIL_LIST;
     it.coilList_ = this;
     it.pos_ = 0;
     it.size_ = list_.GetSize();
     it.ptGrid_ = grid_;
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

  // --- NC Element List ---
  // Default constructor
  NcSurfElemList::NcSurfElemList( Grid * grid )
    : SurfElemList(grid)
  {
  }
  
  // Constructor
  NcSurfElemList::NcSurfElemList( Grid * grid, std::string name )
    : SurfElemList(grid)
  {
    this->type_ = NC_ELEM_LIST;
    name_ = name;
  }

  std::string NcSurfElemList::GetName() const {
    if (name_.length() > 0) {
      return name_;
    } else {
      return "*anonymous list*";
    }
  }
  
  void NcSurfElemList::SetName(const std::string & name) {
    name_ = name;
    defineType_ = NO_TYPE;
  }
  
  void NcSurfElemList::SetNamedElems( const std::string& name ) {
    StdVector<Elem*> elems;
    grid_->GetElemsByName( elems, name );
    
    ncElems_.Clear();
    
    NcSurfElem * actElem = nullptr;
    for ( UInt i=0, numElems=elems.GetSize(); i<numElems; ++i ) {
      actElem = dynamic_cast<NcSurfElem*>(elems[i]);
      if( actElem == nullptr ) {
        EXCEPTION( "Element #" << elems[i]->elemNum 
                   << " in element list '" << name 
                   << "' is no NcSurfElem!" );
      }
      ncElems_.Push_back( shared_ptr<NcSurfElem>(actElem) );
    }
    
    defineType_ = NAMED_ELEMS;
    size_ = ncElems_.GetSize();
    region_ = NO_REGION_ID;
    name_ = name;
  }

  void NcSurfElemList::SetRegion( RegionIdType region ) { 
    StdVector<SurfElem*> elems;
    grid_->GetSurfElems( elems, region );
    
    name_ = grid_->GetRegion().ToString( region );
    ncElems_.Clear();
    
    NcSurfElem * curElem = nullptr;
    for ( UInt i=0, numElems=elems.GetSize(); i<numElems; ++i ) {
      curElem = dynamic_cast<NcSurfElem*>(elems[i]);
      if ( ! curElem ) {
        EXCEPTION( "Element #" << elems[i]->elemNum 
                   << " in element list '" << name_
                   << "' is no NcSurfElem!" );
      }
      ncElems_.Push_back( shared_ptr<NcSurfElem>(curElem) );
    }

    defineType_ = REGION;
    size_ = ncElems_.GetSize();
    region_ = region;
  }

  void NcSurfElemList::SetElement( const Elem * elem ) {
    const NcSurfElem * ncElem = dynamic_cast<const NcSurfElem*>(elem);
    if ( !ncElem ) {
      EXCEPTION("Cannot use NcSurfElemList with anything other than NcSurfElems");
    }
    SetNcSurfElem( shared_ptr<NcSurfElem>( const_cast<NcSurfElem*>(ncElem) ) );
  }

  void NcSurfElemList::SetSurfElem( const SurfElem * elem ) {
    const NcSurfElem * ncElem = dynamic_cast<const NcSurfElem*>(elem);
    if ( !ncElem ) {
      EXCEPTION("Cannot use NcSurfElemList with anything other than NcSurfElems");
    }
    SetNcSurfElem( shared_ptr<NcSurfElem>( const_cast<NcSurfElem*>(ncElem) ) );
  }

  void NcSurfElemList::SetNcSurfElem( const shared_ptr<NcSurfElem> elem ) {
    defineType_ = NO_TYPE;
    size_ = 1;
    region_ = NO_REGION_ID;
    ncElems_.Resize(1, elem);
  }
  
  const Elem * NcSurfElemList::GetElem( UInt nr ) const
  {
    return ncElems_[nr].get();
  }

  const SurfElem * NcSurfElemList::GetSurfElem( UInt nr ) const
  {
    return ncElems_[nr].get();
  }

  //! returns const reference to NcElem
  NcSurfElem* NcSurfElemList::GetNcSurfElem( UInt nr ) const
  {
    return ncElems_[nr].get();
  }

  //! Adds an element using a shared pointer which is better suited here
  void NcSurfElemList::AddElement( const shared_ptr<NcSurfElem> elem ) {
#pragma omp critical (ENTLISTS_NCELEMLIST)
{
    ncElems_.Push_back(elem);
    ++size_;
}
  }

  //! Get iterator
  EntityIterator NcSurfElemList::GetIterator() const {
    EntityIterator it;
    it.type_ = NC_ELEM_LIST;
    it.ptGrid_ = this->grid_;
    it.ncElemList_ = this;
    it.pos_ = 0;
    it.size_ = ncElems_.GetSize();
    return it;
  }

  // =================================================
  // E N T I T Y   I T E R A T O R
  // =================================================
  
  EntityIterator::EntityIterator() {
    type_          = EntityList::NO_LIST;
    ptGrid_        = NULL;
    elemList_      = NULL; 
    surfElemList_  = NULL;
    nodeList_      = NULL;
    regionList_    = NULL;
    nameList_      = NULL;
    coilList_      = NULL;
    numberList_    = NULL;
    ncElemList_    = NULL;
    pos_ = 0;
    size_ = 0;
  }
  
  std::string EntityIterator::ToString() const
  {
     return EntityList::listType.ToString(type_) + " -> " + std::to_string(pos_) + "/" + std::to_string(size_);
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

  EntityIterator& EntityIterator::operator+=(int val){
    assert((pos_+val) < size_);
    pos_ += val;
    return *this;
  }
  
  const Elem* EntityIterator::GetElem() const{
    switch(type_)
    {
    case EntityList::ELEM_LIST:
      return elemList_->GetElem( pos_ );
      
    case EntityList::SURF_ELEM_LIST:
      return surfElemList_->GetSurfElem( pos_ );
      
    case EntityList::NC_ELEM_LIST:
       return ncElemList_->GetNcSurfElem( pos_ );

    default:
      EXCEPTION("type " << type_ << " not implemented");
      break;
    }
    return nullptr;
  }
  

  bool EntityIterator::IsElemType() const
  {
    switch(type_)
    {
    case EntityList::ELEM_LIST:
    case EntityList::SURF_ELEM_LIST:
    case EntityList::NC_ELEM_LIST:
      return true;

    default:
      return false;
    }
  }




  const SurfElem* EntityIterator::GetSurfElem() const {
    switch(type_)
    {
    case EntityList::SURF_ELEM_LIST:
      return surfElemList_->GetSurfElem( pos_ );
      break;
    case EntityList::NC_ELEM_LIST:
      return ncElemList_->GetNcSurfElem( pos_ );
      break;
    default:
      EXCEPTION("type " << type_ << " not implemented for GetSurfElem");
      return nullptr;
    }

  }
  
  const NcSurfElem* EntityIterator::GetNcSurfElem() const {
    if ( type_ == EntityList::NC_ELEM_LIST ) {
      return ncElemList_->GetNcSurfElem( pos_ );
    } else {
      EXCEPTION("Type " << type_ << " not implemented for GetNcSurfElem");
      return nullptr;
    }
  }
  
  RegionIdType EntityIterator::GetRegion() const {
    return regionList_->list_[ pos_ ];
  }
  
  std::string  EntityIterator::GetName() const {
     return nameList_->list_[ pos_ ];
   }
  
  UInt EntityIterator::GetNode() const {
    return nodeList_->list_[ pos_ ];
  }
  
  UInt EntityIterator::GetUnknown() const {
    EXCEPTION( "Not Implemented" );
    return 0;
  }
  
  shared_ptr<Coil> EntityIterator::GetCoil() const {
    return coilList_->list_[ pos_ ];
  }

  std::string EntityIterator::GetIdString() const {
    
    std::string id = "";
    switch( type_ ) {
    case EntityList::ELEM_LIST:
    case EntityList::SURF_ELEM_LIST:
    case EntityList::NC_ELEM_LIST:
      id = std::to_string(GetElem()->elemNum);
      break;
    case EntityList::NODE_LIST:
      id = std::to_string(GetNode());
      break;
    case EntityList::NAME_LIST:
      id = nameList_->GetName();
      break;
    case EntityList::REGION_LIST:
      id = regionList_->grid_->GetRegion().ToString( GetRegion() );
      break;
    case EntityList::COIL_LIST:
      id = coilList_->list_[ pos_ ]->coilId_;
      break;
    default:
      EXCEPTION( "Not implemented" );
      break;
    }
      
     return id; 
  }


} // end of namespace
