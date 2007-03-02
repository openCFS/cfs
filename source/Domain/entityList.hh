#ifndef FILE_CFS_ENTITYLIST_HH
#define FILE_CFS_ENTITYLIST_HH

#include "General/environment.hh"
#include "Utils/StdVector.hh"


namespace CoupledField {

  // Forward class declaration
  class EntityIterator;
  class Grid;
  class SurfElem;
  class Elem;
  class Coil;


  //! Base class for providing lists of geometric entities
  class EntityList {

  public:
    
    //! Typedef for all geometric entities
    typedef enum {NO_LIST, ELEM_LIST, SURF_ELEM_LIST, NC_ELEM_LIST, 
                  NODE_LIST, REGION_LIST, COIL_LIST, 
                  NUMBER_LIST} ListType;

    //! Typedef describing the type the list is defined by
    typedef enum {NO_TYPE, REGION, NAMED_NODES, NAMED_ELEMS}  DefineType;

    //! Constructor
    EntityList( Grid *grid);

    //! Destructor
    virtual ~EntityList();

    //! Return type of EntityList
    ListType GetType() const { return type_;}

    //! Return DefineType of EntityList
    DefineType GetDefineType() const { return defineType_; }
    
    //! Return name of entityList

    //! This method returns the name of the entitylist,
    //! which can be either the name of the Region, the elements or the
    //! nodes contained.
    virtual std::string GetName() const = 0;
    
    
    //! Return size of list
    UInt GetSize() const { return size_; }

    //! Get iterator
    virtual EntityIterator GetIterator() const = 0;

    static std::string TypeToString( ListType type );

    static ListType StringToType( const std::string& type );

  protected:
    
    //! ListType of entities
    ListType type_;

    //! DefineType of entities
    DefineType defineType_;
    
    //! Pointer to grid
    Grid* grid_;
    
    //! Number of entities
    UInt size_;

  };

  //! List for standard elements
  class ElemList : public EntityList {

  public:
    friend class EntityIterator;
    
    //! Constructor
    ElemList( Grid * grid);

    virtual ~ElemList();
    
    //! Returns the name of the region contained in the entitylist
    std::string GetName() const;

    //! Set name of named element list
    void SetNamedElems( const std::string& name);

    //! Set RegionId
    void SetRegion( RegionIdType  region );

    //! Set single element
    void SetElement( const Elem* elem );

    //! Get RegionId
    RegionIdType GetRegion() const;
    
    //! Get one element of the list
    const Elem * GetElem( UInt nr ) const;

    //! Get iterator
    EntityIterator GetIterator() const;

  private:
    StdVector<UInt> list_;

    RegionIdType region_;
    
    std::string name_;
  };
  
  
  //! List for surface elements
  class SurfElemList : public EntityList {
 
  public:
    
    friend class EntityIterator;

    //! Constructor
    SurfElemList( Grid * grid);

    //! Returns the name of the region contained in the entitylist
    std::string GetName() const;

    //! Set name of named element list
    void SetNamedElems( const std::string& name);

    //! Set RegionId
    void SetRegion( RegionIdType  region );

    //! Get RegionId
    RegionIdType GetRegion() const;
    
    //! Get one surface element of the list
    const SurfElem * GetSurfElem( UInt nr ) const;

    //! Get iterator
    EntityIterator GetIterator() const;
    
  private:
    StdVector<SurfElem*> list_;

    RegionIdType region_;

    std::string name_;
  };

  //! List for nodes
  class NodeList : public EntityList {

  public:
    friend class EntityIterator;

    //! Constructor
    NodeList( Grid* grid);

    //! Returns the name of the named nodes contained in the entitylist
    std::string GetName() const;

    //! Set name of named node list
    void SetNamedNodes( const std::string& name);

    //! Set nodes by region ID
    void SetNodesOfRegion( RegionIdType region );

    //! Set explicitly the nodal vector
    void SetNodes( const StdVector<UInt>& nodeList );
    
    //! Get Name of nodes
    const std::string& GetNodesName() { return name_;}

    //! Get list of nodes
    const StdVector<UInt> & GetNodes() { return list_; }
    
    //! Get iterator
    EntityIterator GetIterator() const;
    
  private:
    StdVector<UInt> list_;
    std::string name_;
  };

  
  //! List for regions
  class RegionList : public EntityList {

  public:
    friend class EntityIterator;
    
    //! Constructor
    RegionList( Grid* grid );
    
    //! Returns the name of the region contained in the entitylist
    std::string GetName() const;

    //! Set Region id
    void SetRegionId( RegionIdType name );

    //! Set Region ids
    void SetRegionNames( const StdVector<RegionIdType>& names );

    //! Get Region anmes
    const StdVector<RegionIdType>& GetRegionIds() const {return list_; }
    

    //! Get iterator
    EntityIterator GetIterator() const;
    
  private:
    StdVector<RegionIdType> list_;
  };

  //! List for magnetic coils
  class CoilList : public EntityList {
  public:
    friend class EntityIterator;
    
    //! Constructor
    CoilList( Grid* grid );
    
    //! Set coil
    void AddCoil( shared_ptr<Coil> aCoil );
    
    //! Returns the name of the coil
    std::string GetName() const;

    //! Get iterator
    EntityIterator GetIterator() const;
    
  private:
    StdVector<shared_ptr<Coil> > list_;

  };

  //! List for non-geometric entities
  class NumberList : public EntityList {
  public:
    friend class EntityIterator;
    
    //! Constructor
    NumberList( Grid* grid );
    
    //! Returns the name of the free entities
    std::string GetName() const;

    //! Get iterator
    EntityIterator GetIterator() const;
    
  private:
    StdVector<UInt> list_;

  };

  //! Iterator class for all kinds of iterable entities
  class EntityIterator {
    
  public:
    
    friend class ElemList;
    friend class SurfElemList;
    friend class NodeList;
    friend class RegionList;
    friend class CoilList;
    friend class NumberList;

    EntityIterator();
    
    //! Get type of iterator
    EntityList::ListType GetType() const { return type_; }
    
    //! Set iterator to begin of list
    void Begin();
    
    //! True, if iterator is at end of list
    bool IsEnd() const;
    
    //! Increment iterator
    EntityIterator& operator++(int);
    
    //! Return current element
    const Elem* GetElem() const;
    
    //! Return current surface element
    const SurfElem* GetSurfElem() const;

    //! Return current region id
    RegionIdType GetRegion() const;

    //! Return current node number
    UInt GetNode() const;
    
    //! Return current coil
    shared_ptr<Coil> GetCoil() const;

    //! Return unknown number
    UInt GetUnknown() const;

    //! Return position
    UInt GetPos() const { return pos_; }
    
  protected:
    EntityList::ListType type_;
    const ElemList* elemList_;
    const SurfElemList* surfElemList_;
    const NodeList* nodeList_;
    const RegionList* regionList_;
    const CoilList * coilList_;
    const NumberList* numberList_;
    UInt pos_;
    UInt size_;
  };

}

#endif
