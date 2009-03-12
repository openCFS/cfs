// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_RESULT_HH
#define FILE_CFS_RESULT_HH

#include "Domain/entityList.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"

#include "MatVec/basematrix.hh"
#include "Domain/resultInfo.hh"

namespace CoupledField {

  //! Base class representing a result object
  class BaseResult {

  public:

    BaseResult();

    //! destructor
    virtual ~BaseResult();

    //! Get entry type of data
    virtual BaseMatrix::EntryType GetEntryType() const = 0;

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
    virtual SingleVector* GetSingleVector() = 0 ;

    /** Set all result values to the null value. Used, if one cannot compute */
    virtual void Init() = 0;

    /** Gives back some information for debug output */
    std::string ToString() const;

    /** Only some special region (integral) results have a InfoNode to
     * write the data also to info.xml
     * @return often NULL */
    InfoNode* GetInfoNode() { return infoNode_; }

    /** To leave the constructor clean, set an InfoNode from external if
     * the result is of type region integral or such like */
    void SetInfoNode(InfoNode* in) { this->infoNode_ = in; }

    /** Dumps a result list */
    static void Dump(StdVector<shared_ptr<BaseResult> >& resultList);

  protected:

    //! Object describing the type of result
    shared_ptr<ResultInfo> resultDof_;

    //! Entitylist the result is associated with
    shared_ptr<EntityList> entities_;

    /** Some results, like region results with a single scalar, get here the InfoNode
     * to write the data at the right position in info.xml */
    InfoNode* infoNode_;
  };


  //! Class containing the simulation results
  template <class TYPE>
  class Result : public BaseResult
  {
  public:

    //! Constructor
    Result();

    //! Destructor
    virtual ~Result();

    //! Get entry type of data
    BaseMatrix::EntryType GetEntryType() const {
      return  EntryType<TYPE>::M_EntryType;
    }

    //! Return data vector
    SingleVector* GetSingleVector() { return &values_; }

    //! Return specific data vector
    Vector<TYPE>& GetVector() {return values_; }

    /** Initialize data vector */
    void Init();

  protected:

    //! data vector
    Vector<TYPE> values_;

  };

}

#endif
