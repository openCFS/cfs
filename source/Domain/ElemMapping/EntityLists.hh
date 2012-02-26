// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_ENTITYLIST_HH
#define FILE_CFS_ENTITYLIST_HH

#include "General/Environment.hh"
#include "Utils/StdVector.hh"


namespace CoupledField {

  // Forward class declaration
  class EntityIterator;
  class Grid;
  struct SurfElem;
  struct NcSurfElem;
  struct Elem;
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

    static Enum<ListType> listType;

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
    
    virtual RegionIdType GetRegion() const { return -1; }

    //! Return pointer to grid.
    Grid* GetGrid() const { return grid_; }
    
    //! Return size of list
    UInt GetSize() const { return size_; }

    //! Get iterator
    virtual EntityIterator GetIterator() const = 0;

    // =======================================================================
    //  CONVERSION METHODS
    // =======================================================================

    static void SetEnums();

    //! Conversion from ListType to string
    static void Enum2String( ListType in, std::string& out );

    //! Conversion from string to ListType
    static void String2Enum( const std::string& in, ListType& out );

    //! Conversion from DefineType to string
    static void Enum2String( DefineType in, std::string& out );

    //! Conversion from string to DefineType
    static void String2Enum( const std::string& in, DefineType& out );

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
    
    //! Create explicitly element list from single element
    explicit ElemList( const Elem* elem, Grid* grid );

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
    virtual ~SurfElemList() {}

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
    virtual ~NodeList() {}

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
    virtual ~RegionList() {}
    
    //! Returns the name of the region contained in the entitylist
    std::string GetName() const;

    //! Set Region id
    void SetRegion( RegionIdType name );

    //! Set Region ids
    void SetRegions( const StdVector<RegionIdType>& names );

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
    virtual ~CoilList() {}
    
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
    virtual ~NumberList() {}
    
    //! Returns the name of the free entities
    std::string GetName() const;

    //! Get iterator
    EntityIterator GetIterator() const;
    
  private:
    StdVector<UInt> list_;

  };

  //! Class for storing NC elems
  //! IMPORTANT: in difference to the other entitylists,
  //!   this list stroes directly the pointers to the elements
  class NcElemList : public EntityList{
  public:

    //! Constructor
    NcElemList( Grid * grid,std::string name);

    virtual ~NcElemList(){

    }

    //! returns the name of the list
    std::string GetName() const;

    //! returns const reference to NcElem
    const NcSurfElem * GetNcSurfElem( UInt nr ) const;

    //! Adds an element using a shared pointer which is better suited here
    void SetElement( const shared_ptr<NcSurfElem> elem );

    //! Get iterator
    EntityIterator GetIterator() const;

  private:
    StdVector< shared_ptr<NcSurfElem> > ncElems_;

    std::string name_;
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
    friend class NcElemList;

    //! Constructor
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

    //! Return current surface element
    const NcSurfElem* GetNcSurfElem() const;

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
    
    /** The total size */
    UInt GetSize() const { return size_; }

    //! This method returns for each type of unknown a characteristic string 
    
    //! This method resturn the following information (as string) for the 
    //! following entity types:
    //! ELEM_LIST, SURF_ELEM_LIST, NC_ELEM_LIST: number of element 
    //! NODE_LIST: node number
    //! REGION_LIST: region name
    //! COIL_LIST: coild id
    //! NUMBER_LIST: unknown number
    std::string GetIdString() const;
    
    //! Returns the name of the list, i.e. region name, etc.
    std::string GetName() const;
    
    
  protected:
    EntityList::ListType type_;
    const ElemList* elemList_;
    const SurfElemList* surfElemList_;
    const NodeList* nodeList_;
    const RegionList* regionList_;
    const CoilList * coilList_;
    const NumberList* numberList_;
    const NcElemList* ncElemList_;
    UInt pos_;
    UInt size_;
    std::string name_;
  };

}

#endif
