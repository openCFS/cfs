// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_RESULT_HH
#define FILE_CFS_RESULT_HH

#include "Domain/entityList.hh"
#include "vector.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Domain/entityList.hh"


namespace CoupledField {

  //! Base class representing a result object
  class BaseResult {
    
  public:
    
    //! constructor
    BaseResult();
    
    //! destructor
    virtual ~BaseResult();
    
    //! Get entry type of data
    virtual EntryType::ScalarType GetEntryType() const = 0;

    //! Return result dof
    shared_ptr<ResultInfo> GetResultInfo() { return resultDof_; }

    //! Set result dof object
    void SetResultInfo( shared_ptr<ResultInfo> resultDof ) {
      resultDof_ = resultDof;
    }

    //! Set entitiylist
    void SetEntityList( shared_ptr<EntityList> list ) {
      entities_ = list;
    }

    //! Return entitylist
    shared_ptr<EntityList> GetEntityList() const {return entities_;}
    
    //! Return vector containing data
    virtual CFSVector* GetCFSVector() = 0 ;
    
  protected:
    
    //! Object describing the type of result
    shared_ptr<ResultInfo> resultDof_;

    //! Entitylist the result is associated with
    shared_ptr<EntityList> entities_;

  };


  //! Class containing the simulation results
  template <class TYPE>
  class Result : public BaseResult {

  public:

    //! Constructor
    Result();
    
    //! Destructor
    virtual ~Result();

    //! Get entry type of data
    EntryType::ScalarType GetEntryType() const {
      return EntryTypeMap<TYPE>::S_TYPE;
    }

    //! Return data vector
    CFSVector* GetCFSVector() { return &values_; }

    //! Return specific data vector
    Vector<TYPE>& GetVector() {return values_; }

  protected:
    
    //! data vector
    Vector<TYPE> values_;

  };

}

#endif
