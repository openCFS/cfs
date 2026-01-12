// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ResultManager.hh
 *       \brief    <Description>
 *
 *       \date     Nov 24, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_UTILS_RESULTMANAGER_HH_
#define SOURCE_CFSDAT_UTILS_RESULTMANAGER_HH_



#include <utility>
#include <vector>

//boost
#define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <unordered_map>

//CFS
#include "General/Environment.hh"
#include "Domain/Results/BaseResults.hh"

//Dat
#include "cfsdat/DatUtils/Defines.hh"
#include "cfsdat/DatUtils/DataStructs.hh"
#include "cfsdat/DatUtils/ResultCache.hh"
#include "EqnNumberingSimple.hh"

#define DEBUG

namespace CFSDat{

typedef CF::StdVector<uuids::uuid> ResultIdList;

//Basically a openCFS Resulthandler adapted to our special needs
//its main purpose is to coordinate result access and to give the
//possibility to optimize result storage independent of the
//filter pipeline defined
class ResultManager{

public:

  typedef str1::shared_ptr<ExtendedResultInfo> InfoPtr;

  typedef str1::shared_ptr<const ExtendedResultInfo> ConstInfoPtr;

  typedef str1::shared_ptr<GenericResultCache> ResPtr;

  typedef str1::shared_ptr<const GenericResultCache> ConstResPtr;

  ResultManager();

  ~ResultManager();

  uuids::uuid AddResult(std::string name, uuids::uuid filterTag);
  
  uuids::uuid AddResult(std::string name, uuids::uuid filterTag, 
         Integer minStepOffset, Integer maxStepOffset);
  
  void CombineResults(uuids::uuid resultIdA, uuids::uuid resultIdB);
  
  uuids::uuid GetMasterResult(uuids::uuid resultId);

  // =======================================================================
  // Results Access
  // =======================================================================
  //@{ \name Result(Adapter/Vector) Access

  //! Sets the index of the actual step value for the given id
  //! Results get marked as invalid if the new value does not match any cached result.
  //! \param (in) requestedId unique id of result
  //! \param (in) value new step number
  void SetStepIndex(uuids::uuid requestedId, Integer value);
  
  //! Obtain index of the actual step value
  //! \param (in) requestedId  unique id of result
  //! \return UInt index of the step
  UInt GetStepIndex(uuids::uuid requestedId);

  //! Sets the step value of the given id
  //! Results get marked as invalid if the new value does not match any cached result.
  //! \param (in) requestedId unique id of result
  //! \param (in) value new step value
  void SetStepValue(uuids::uuid requestedId, Double value);
  
  //! Obtain step value of given result
  //! \param (in) requestedId  unique id of result
  //! \return Double value of the step
  Double GetStepValue(uuids::uuid requestedId);
  
  //! Check if the given result is already updated
  //! \param (in) requestedId  unique id of result
  //! \param (in) offsetStep time/frequency step offset for future/past results
  //! \return Bool flag, true if result is valid
  bool IsResultVecUpToDate(uuids::uuid requestedId);

  void SetResultVecUpToDate(uuids::uuid requestedId, bool upToDate=true);
  
  //! Direct access to solution vector, real valued version
  //! \param (in) requestedId  unique id of result
  //! \param (out) eqnNumbers vector of equation numbers in case the result vector is defined globally
  //! \return  reference to result vector
  template<typename T>
  Vector<T>&  GetResultVector(uuids::uuid requestedId, CF::StdVector<UInt>& eqnNumbers);

  //! Direct access to solution vector, real valued version
  //! \param (in) requestedId  unique id of result
  //! \param (out) eqnNumbers vector of equation numbers in case the result vector is defined globally
  //! \return  reference to result vector
  template<typename T>
  Vector<T>&  GetResultVector(uuids::uuid requestedId);

//  //! Direct access to solution vector, complex valued version
//  //! \param (in) requestedId  unique id of result
//  //! \param (out) eqnNumbers vector of equation numbers in case the result vector is defined globally
//  //! \return reference to result vector
//  Vector<Complex>& GetResultVector(uuids::uuid requestedId, CF::StdVector<UInt>& eqnNumbers);

  //! Return the equation map of a result
  //! \param (in) requestedId  unique id of result
  //! \return shared pointer to equation map
  str1::shared_ptr<EqnMapSimple> GetEqnMap(uuids::uuid requestedId);

  //! Return reference to BaseResult vector of result adapter. Used e.g. in output filters
  //! to copy the values in SimOut conformal matter
  //! \param (in) requestedId  unique id of result
  //! \return reference to Result objects (e.g. in case multiple regions)
  CF::StdVector<shared_ptr<BaseResult> >& GetBaseResultVector(uuids::uuid requestedId);

  // =======================================================================
  // Extended Result Info Access
  // =======================================================================
  //@{ \name Result Info Access

  //! Returns const pointer to result Info for quering information
  //! \param (in) resId unique id of Result
  //! \return shared pointer to constant element info
  ConstInfoPtr GetExtInfo(uuids::uuid resId);

  //! Copies everything from the given result ID except the result name and unit
  //!   - Memory will be reallocated in this method.
  //!   - All fields will be overwritten
  //! \param (in) srcId  unique id of result to be copied from
  //! \param (in) trgId unique id of result to be copied to
  //! \return Nothing
  void CopyResultData(uuids::uuid srcId, uuids::uuid trgId);

  //! Returns a copy of the timeline of the current result.
  //!   - Memory will be copied in this method.
  //! \param (in) resId unique id of Result
  //! \return Vector of time step values (copied)
  CF::StdVector<Double> GetTimeLine(uuids::uuid resId);

  //! Returns true if a result is constant in time
  //! \param (in) resId unique id of Result
  //! \return true if a result is constant in time
  bool IsStatic(uuids::uuid resId);
  
  //! Set if a quantity is static in time
  //! \param (in) resId unique id of Result
  //! \param (in) isStatic if the quantity is static
  //! \return Nothing
  void SetStatic(uuids::uuid resId, bool isStatic);

  //bool IsConstant(uuids::uuid resId); Using IsStatic() due to merge from tuWien

  //! Set the time step value vector of the given result
  //!   - Replace if vector has already been set
  //!   - Only valid during initialization phase.
  //! \param (in) resId unique id of Result
  //! \param (in) tVec new time step values
  //! \return Nothing
  void SetTimeLine(uuids::uuid resId, CF::StdVector<Double> tVec);


  //! Returns a copy of the setpnumbers of the current result.
  //!   - Memory will be copied in this method.
  //! \param (in) resId unique id of Result
  //! \return Vector of time step values (copied)
  CF::StdVector<UInt> GetStepNumbers(uuids::uuid resId);

  //! Set the time step number vector of the given result
  //!   - Replace if vector has already been set
  //!   - Only valid during initialization phase.
  //! \param (in) resId unique id of Result
  //! \param (in) tVec new time step vector
  //! \return Nothing
  void SetStepNumbers(uuids::uuid resId, CF::StdVector<UInt> tVec);

  //! Returns the current data type of the result
  //!   - Memory will be copied in this method.
  //! \param (in) resId unique id of Result
  //! \return Current datatype of result (DOUBLE/COMPLEX/INT)
  ExtendedResultInfo::ResDType GetDType(uuids::uuid resId);

  //! Set the datatype of the result.
  //!   - Only valid during initialization phase.
  //!   - Replace if type has already been set
  //! \param (in) resId unique id of Result
  //! \param (in) dType new datatype request
  //! \return Nothing
  void SetDType(uuids::uuid resId, ExtendedResultInfo::ResDType dType);

  //! Returns the current data type of the result
  //!   - Memory will be copied in this method.
  //! \param (in) resId unique id of Result
  //! \return Current datatype of result (DOUBLE/COMPLEX/INT)
  ExtendedResultInfo::EntityUnknownType GetDefOn(uuids::uuid resId);

  //! Set the entity type for the result
  //!   - Only valid during initialization phase.
  //!   - Replace if type has already been set
  //! \param (in) resId unique id of Result
  //! \param (in) eType new entity type
  //! \return Nothing
  void SetDefOn(uuids::uuid resId, ExtendedResultInfo::EntityUnknownType eType);

  //! Some Results are defined on complete regions (e.g. output Results)
  //! Here, we return a copy of the region names
  //!   - Memory will be copied in this method.
  //! \param (in) resId unique id of Result
  //! \return Set of currently defined region names for the result
  std::set<std::string>  GetRegionNames(uuids::uuid resId);

  //! Set the vector of region names.
  //!   - Replaces if already filled
  //!   - Only valid during initialization phase.
  //! \param (in) resId unique id of Result
  //! \param (in) regions names to be added to the result
  //! \return Nothing
  void SetRegionNames(uuids::uuid resId, std::set<std::string>  regions);

  //! Set the info to indicate a result define on volume or surface meshes
  //!   - Replaces if already set
  //!   - Only valid during initialization phase.
  //! \param (in) resId unique id of Result
  //! \param (in) isMesh value is true if it is a mesh result
  //! \return Nothing
  void SetMeshResult(uuids::uuid resId, bool isMesh);

  //! Obtain the vector of entity numbers. Useful if the result is defined
  //! only on some elements/nodes within the given region(s).
  //!   - Memory will be copied in this method.
  //! \param (in) resId unique id of Result
  //! \return Vector of entity numbers referring to nodes or elements
  CF::StdVector<UInt> GetEntityNumbers(uuids::uuid resId);

  //! Inserts specific entity numbers into the set. Make sure that the
  //! numbers are valid w.r.t. the given definedOn Entity type and region names!
  //! \note{If this set is empty at the end of initialization, it will be
  //! filled from region name vector}
  //!   - Inserts if set has already been filled
  //!   - Only valid during initialization phase.
  //! \param (in) resId unique id of Result
  //! \param (in) numbers to be inserted into the set
  //! \return Nothing
  void SetEntityNumbers(uuids::uuid resId, CF::StdVector<UInt> numbers);

  //! Obtain the vector of equation numbers of current result.
  //! This is important as the result vector can be bigger than
  //! the number of equations relevant for the filter (i.e. entity numbers
  //! is subset).
  //!   - Memory will be copied in this method.
  //!   - Only valid after initialization phase.
  //! \param (in) resId unique id of Result
  //! \return Vector of equation numbers
  StdVector<UInt> GetEqnNumbers(uuids::uuid resId);

  //! Returns the complex data format
  //!   - Memory will be copied in this method.
  //! \param (in) resId unique id of Result
  //! \return Current format of complex values (REAL_IMAG or AMPLITUDE_PHASE)
  CF::ComplexFormat GetCFormat(uuids::uuid resId);

  //! Set the desired complex format for the result.
  //!   - Only valid during initialization phase.
  //!   - Replace if type has already been set
  //! \param (in) resId unique id of Result
  //! \param (in) cType new complex format request
  //! \return Nothing
  void SetCFormat(uuids::uuid resId, CF::ComplexFormat cType);

  //! Returns the the spatial type of the result (SALAR/VECTOR/...) data format.
  //!   - Memory will be copied in this method.
  //! \param (in) resId unique id of Result.
  //! \return Current result entry type ( SCALAR/VECTOR/TENSOR/STRING )
  ExtendedResultInfo::EntryType GetEntryType(uuids::uuid resId);

  //! Set the entry type of the requested result.
  //!   - Only valid during initialization phase.
  //!   - Replace if type has already been set.
  //! \param (in) resId unique id of Result
  //! \param (in) rType new complex format request
  //! \return Nothing
  void SetEntryType(uuids::uuid resId, ExtendedResultInfo::EntryType rType);

  //! Returns the physical unit of the result (e.g. m/s).
  //! \param (in) resId unique id of Result
  //! \return Physical unit in string representation
  std::string GetUnit(uuids::uuid resId);

  //! Set the entry type of the requested result.
  //!   - Only valid during initialization phase.
  //!   - Replace if unit has already been set.
  //! \param (in) resId unique id of Result
  //! \param (in) newUnit new result unit
  //! \return Nothing
  void SetUnit(uuids::uuid resId, std::string newUnit);

  //! Returns result name as string
  //! \param (in) resId unique id of Result
  //! \return Result name string
  std::string GetResultName(uuids::uuid resId);

  //! Get the names of the degrees of freedom.
  //! \param (in) resId unique id of Result
  //! \return Vector of DOF names (e.g. "x","y","z")
  CF::StdVector<std::string> GetDofNames(uuids::uuid resId);

  //! Set the names of the result DOFs.
  //!   - Only valid during initialization phase.
  //!   - Replace if names were already set.
  //! \param (in) resId unique id of Result
  //! \param (in) vector of DOF names
  //! \return Nothing
  void SetDofNames(uuids::uuid resId, CF::StdVector<std::string> newDofNames);

  //! Mark the given result as valid i.e. the result is up to date
  //!   - Replace if result was invalid.
  //! \param (in) resId unique id of Result
  //! \return Nothing
  void SetValid(uuids::uuid resId);

  //! Set the grid pointer of the given Result
  //!   - Only valid during initialization
  //!   - field will be overwritten
  //!   - Replace if result was invalid.
  //! \param (in) resId unique id of Result
  //! \param (in) ptGrid pointer to grid object
  //! \return Nothing
  void SetGrid(uuids::uuid resId, Grid* ptGrid);

  //! Set the given result as ouput and give flag if
  //! mesh result. in case of mesh results, the BaseResult struct will
  //! be initialized with the specified regions. Otherwise only a
  //! single baseResult will be created
  //! \params (in) requestedId unique id of result
  //! \params (in) isMeshResult for volume results (e.g. to hdf5)
  void SetAsOutputResult(uuids::uuid requestedId, bool isMeshResult){
    checkNotFinalized();
    resultMap_[requestedId].first->isOutput = true;
    resultMap_[requestedId].first->isMeshResult = isMeshResult;
  }
  //@}


  void CompressResults();

  void Finalize();


  bool IsFinalized(){
    return isFinalized_;
  }

  std::set<uuids::uuid> GetActiveResults(){
    return activeIds_;
  }

  bool isActiveResult(uuids::uuid isActive){
    return activeIds_.find(isActive)!= activeIds_.end();
  }

  void ActivateResult(uuids::uuid toActive){
    if(resultMap_.find(toActive) == resultMap_.end()){
      EXCEPTION("Trying to activate a result which is not generated! Make sure every call "
                      "to ResultManager provides valid resultIds.");
    }
    activeIds_.insert(toActive);
  }

  void DeactivateResult(uuids::uuid toSleep){
    if(resultMap_.find(toSleep) == resultMap_.end()){
      EXCEPTION("Trying to deactivate a result which is not generated! Make sure every call " <<
                      "to ResultManager provides valid resultIds.");
    }
    activeIds_.erase(toSleep);
  }

private:

  inline void checkNotFinalized(){
    if(isFinalized_){
      EXCEPTION("Call not valid if ResultManager is finalized.");
    }
  }

  inline void checkFinalized(){
    if(!isFinalized_){
      EXCEPTION("Call not valid if ResultManager is not finalized.");
    }
  }

  inline bool checkEqualStageStrict(InfoPtr resOne,InfoPtr resTwo){
    bool isSimilar = true;
    isSimilar &= resOne->dType                          == resTwo->dType;
    isSimilar &= (*resOne->regNames.get())              == (*resTwo->regNames.get());
    isSimilar &= (*resOne->entityNumbers.get())         == (*resTwo->entityNumbers.get());
    isSimilar &= (*resOne->eqnNumbers.get())            == (*resTwo->eqnNumbers.get());
    isSimilar &= (*resOne->timeLine.get())              == (*resTwo->timeLine.get());
    isSimilar &= resOne->definedOn                      == resTwo->definedOn;
    isSimilar &= resOne->dofNames                       == resTwo->dofNames;
    isSimilar &= resOne->entryType                      == resTwo->entryType;
    isSimilar &= resOne->resultName                     == resTwo->resultName;
    isSimilar &= resOne->unit                           == resTwo->unit;
    isSimilar &= resOne->ptGrid                         == resTwo->ptGrid;
    return isSimilar;
  }

  inline bool checkEqualEqnMap(std::pair<InfoPtr,ResPtr> resOne,std::pair<InfoPtr,ResPtr> resTwo){
    bool isSimilar = true;
    isSimilar &= resOne.first->dType                          == resTwo.first->dType;
    isSimilar &= (*resOne.first->regNames.get())              == (*resTwo.first->regNames.get());
    isSimilar &= resOne.first->definedOn                      == resTwo.first->definedOn;
    isSimilar &= resOne.first->entryType                      == resTwo.first->entryType;
    isSimilar &= resOne.first->ptGrid                         == resTwo.first->ptGrid;
    isSimilar &= resOne.first->dofNames.GetSize()             == resTwo.first->dofNames.GetSize();
    return isSimilar;
  }

  inline bool checkEqualEqnVec(std::pair<InfoPtr,ResPtr> resOne,std::pair<InfoPtr,ResPtr> resTwo){
    bool isSimilar = true;
    isSimilar &= checkEqualEqnMap(resOne,resTwo);
    isSimilar &= (*resOne.first->entityNumbers.get())         == (*resTwo.first->entityNumbers.get());
    return isSimilar;
  }

  void SetResultVectorSize(InfoPtr cInfo, ResPtr cRes);

  void CreateEqnMapping(InfoPtr cInfo, ResPtr cRes);


  bool isFinalized_;

  typedef std::unordered_map<uuids::uuid,std::pair<InfoPtr,ResPtr> > UuidMap;

  //map storing for each result id a pair of infos and the actual result
  //will get changed during result compression
  UuidMap resultMap_;

  //set of Ids to be activated
  std::set<uuids::uuid> activeIds_;

};

}



#endif /* SOURCE_CFSDAT_UTILS_RESULTMANAGER_HH_ */
