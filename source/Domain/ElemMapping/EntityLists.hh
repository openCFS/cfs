// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_ENTITYLIST_HH
#define FILE_CFS_ENTITYLIST_HH

#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "Utils/Coil.hh"

namespace CoupledField {

  // Forward class declaration
  class EntityIterator;
  class Grid;
  struct SurfElem;
  struct NcSurfElem;
  struct Elem;
  class Coil;


  // -----------------------
  //  Base Class EntityList
  // -----------------------
  //! Base class for providing lists of geometric entities
  class EntityList {

  public:
    
    //! Typedef for all geometric entities
    typedef enum {NO_LIST, ELEM_LIST, SURF_ELEM_LIST, NC_ELEM_LIST, 
                  NODE_LIST, REGION_LIST, COIL_LIST, 
                  NUMBER_LIST, NAME_LIST} ListType;

    //! Typedef describing the type the list is defined by
    typedef enum {NO_TYPE, REGION, NAMED_NODES, NAMED_ELEMS}  DefineType;

    static Enum<ListType> listType;

    static Enum<DefineType> defineType;

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
    
    virtual RegionIdType GetRegion() const { return NO_REGION_ID; }

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


    // =======================================================================
    //  UTITLY METHODS
    // =======================================================================
    
    //! Get intersection of two sets of EntityLists
    static void Intersect(const StdVector<shared_ptr<EntityList> >& set1,
                          const StdVector<shared_ptr<EntityList> >& set2,
                          StdVector<shared_ptr<EntityList> >& intersect );
    
    //! Get union of two sets of tEntityLists
    static void Union(const StdVector<shared_ptr<EntityList> >& set1,
                      const StdVector<shared_ptr<EntityList> >& set2,
                      StdVector<shared_ptr<EntityList> >& intersect );

    /** dumps some data for logging and debugging purpose */
    std::string ToString() const;
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
  
  // --------------
  //  Element List
  // --------------

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
    virtual void SetNamedElems( const std::string& name );

    //! Set RegionId
    virtual void SetRegion( RegionIdType region );

    //! Set single element
    virtual void SetElement( const Elem* elem );

    //! Add an element to the list
    virtual void AddElement( const Elem* elem );
    
    //! Add vector of elements to the list
    void AddElements( const StdVector<Elem*>& elems );

    //! Get RegionId
    RegionIdType GetRegion() const;
    
    //! Get one element of the list
    virtual const Elem * GetElem( UInt nr ) const;

    //! Get iterator
    virtual EntityIterator GetIterator() const;

  protected:
    StdVector<UInt> list_;

    RegionIdType region_;
    
    std::string name_;
  };
  
  // ----------------------
  //  Surface Element List
  // ----------------------
  
  //! List for surface elements
  class SurfElemList : public ElemList {
 
  public:
    
    friend class EntityIterator;

    //! Constructor
    SurfElemList( Grid * grid);
    
    //! Destructor
    virtual ~SurfElemList() {}

    //! Set name of named element list
    virtual void SetNamedElems( const std::string& name );

    //! Set RegionId
    virtual void SetRegion( RegionIdType region );

    //! Set single element (if it can be casted to a SurfElem)
    virtual void SetElement( const Elem* elem );

    //! Set single surface element
    virtual void SetSurfElem( const SurfElem* elem );

    //! Add an element to the list
    virtual void AddElement( const SurfElem* elem );

    //! Get one element of the list
    virtual const Elem * GetElem( UInt nr ) const;

    //! Get one surface element of the list
    virtual const SurfElem * GetSurfElem( UInt nr ) const;

    //! Get iterator
    virtual EntityIterator GetIterator() const;
    
  private:
    StdVector<const SurfElem*> surfElemList_;
  };
  
  // ------------------------------
  //  Non Conforming Elements List
  // ------------------------------
  
  //! Class for storing NC elems
  //! IMPORTANT: in difference to the other entitylists,
  //!   this list stores directly the pointers to the elements
  class NcSurfElemList : public SurfElemList {
  public:

    //! Default constructor
    NcSurfElemList( Grid * grid);
    
    //! Constructor
    NcSurfElemList( Grid * grid, std::string name );

    virtual ~NcSurfElemList(){};

    //! Returns the name of the NcSurfElemList
    std::string GetName() const;

    //! Sets the name of the NcSurfElemList
    virtual void SetName(const std::string &name);
    
    //! Set name of named element list
    virtual void SetNamedElems( const std::string& name );

    //! Set RegionId
    virtual void SetRegion( RegionIdType region );

    //! Set single element (if it can be casted to a NcSurfElem)
    virtual void SetElement( const Elem* elem );

    //! Set single surface element (if it can be casted to a NcSurfElem)
    virtual void SetSurfElem( const SurfElem* elem );

    //! Set single NcSurfElem
    virtual void SetNcSurfElem( const shared_ptr<NcSurfElem> elem );

    //! Get one element of the list
    virtual const Elem * GetElem( UInt nr ) const;

    //! Get one surface element of the list
    virtual const SurfElem * GetSurfElem( UInt nr ) const;

    //! returns const reference to NcElem
    virtual NcSurfElem * GetNcSurfElem( UInt nr ) const;

    //! Adds an element using a shared pointer which is better suited here
    virtual void AddElement( const shared_ptr<NcSurfElem> elem );

    //! Deletes all elements from the list
    virtual void Clear(bool keepCapacity=false) { ncElems_.Clear(keepCapacity); size_ = 0; }
    
    //! Get iterator
    virtual EntityIterator GetIterator() const;

  private:
    StdVector< shared_ptr<NcSurfElem> > ncElems_;
    
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

  // -------------
  //  Region List
  // -------------
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

  
  // -----------
  //  Coil List
  // -----------
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

  // -------------
  //  Number List
  // -------------
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
  
  
  // -------------
  //  Name List
  // -------------

  //! List for names (regions or element/node groups)
  class NameList : public EntityList {

  public:
    friend class EntityIterator;

    //! Constructor
    NameList( Grid* grid );

    //! Returns the name of the region / group contained in the entitylist
    std::string GetName() const;

    //! Set Name
    void SetName( const std::string& name );

    //! Set list of names
    void SetNames( const StdVector<std::string>& names );


    //! Get iterator
    EntityIterator GetIterator() const;

  private:
    StdVector<std::string> list_;
  };

  // -----------------
  //  Entity Iterator
  // -----------------

  //! Iterator class for all kinds of iterable entities
  class EntityIterator {
    
  public:
    
    friend class ElemList;
    friend class SurfElemList;
    friend class NcSurfElemList;
    friend class NodeList;
    friend class RegionList;
    friend class NameList;
    friend class CoilList;
    friend class NumberList;

    //! Constructor
    EntityIterator();
    
    //! Get type of iterator
    EntityList::ListType GetType() const { return type_; }
    
    //! Set iterator to begin of list
    void Begin();
    
    //! True, if iterator is at end of list
    bool IsEnd() const;
    
    
    //! Get pointer to grid
    Grid* GetGrid() {
      return ptGrid_;
    }
    
    //! Increment iterator
    EntityIterator& operator++(int);
    
    //! Increment iterator by fixed value
    EntityIterator& operator+=(int val);

    //! Return current element
    const Elem* GetElem() const;
    
    //! Return current surface element
    const SurfElem* GetSurfElem() const;

    //! Return current surface element
    const NcSurfElem* GetNcSurfElem() const;

    //! Return current region id
    RegionIdType GetRegion() const;
    
    //! Return current name
    std::string GetName() const;

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
    //! COIL_LIST: coil id
    //! NUMBER_LIST: unknown number
    std::string GetIdString() const;
    
    
  protected:
    EntityList::ListType type_;
    Grid* ptGrid_;
    const ElemList* elemList_;
    const SurfElemList* surfElemList_;
    const NodeList* nodeList_;
    const RegionList* regionList_;
    const NameList* nameList_;
    const CoilList * coilList_;
    const NumberList* numberList_;
    const NcSurfElemList* ncElemList_;
    UInt pos_;
    UInt size_;
  };
  
}

#endif
