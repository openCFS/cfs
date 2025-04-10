// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SimInputEnsight.hh
 *       \brief    <Description>
 *
 *       \date     21.02.2014
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SIMINPUTENSIGHT_HH_
#define SIMINPUTENSIGHT_HH_

#include "DataInOut/SimInOut/VTKBased/SimInputVTKBased.hh"
#include <boost/unordered_map.hpp>

namespace CoupledField{

class SimInputEnsight : public SimInputVTKBased {

public:

  //! definition of a result with additional information for ensight
  struct ResultDef{
    ResultDef(){
      isComplex = false;
      isStatic = false;
    }
    //! CFS result deifintion
    shared_ptr<ResultInfo> res;

    //! map for assigning each dof of the result
    //! a string for the ensight file
    StdVector<std::string> dofNameMap;

    //!flag indicating if this result is valid i.e. definied by user and given in
    //! ensight file
    StdVector<bool> isValid_;

    bool isComplex;
    
    bool isStatic;
  };

  SimInputEnsight(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode );

  virtual ~SimInputEnsight();

  //! \copydoc SimInput::GetResultTypes
  virtual void GetResultTypes( UInt sequenceStep,
                                  StdVector<shared_ptr<ResultInfo> >& infos,
                                  bool isHistory = false );

  //! \copydoc SimInput::GetResultEntities
  virtual void GetResultEntities( UInt sequenceStep,
                                      shared_ptr<ResultInfo> info,
                                      StdVector<shared_ptr<EntityList> >& list,
                                      bool isHistory = false);

  //! \copydoc SimInput::GetResult
  virtual void GetResult( UInt sequenceStep,
                             UInt stepValue,
                             shared_ptr<BaseResult> result,
                             bool isHistory = false );
protected:
  //! \copydoc SimInputVTKBased::GetTimeValues
  virtual void GetTimeValues();

  //! \copydoc SimInputVTKBased::CreateReader
  virtual void CreateReader();

  //! \copydoc SimInputVTKBased::SetTimeValue
  virtual void SetTimeValue(Double val);

  //! \copydoc SimInputVTKBased::EnableRegions
  virtual void EnableRegions(){}

  //! Check validity of user input and Results by element or by node
  virtual void ValidateResultDefinition();

private:

  //! specialization for nodal results
  template<typename T>
  void GetNodeResult( UInt sequenceStep,
      UInt stepValue,
      shared_ptr<BaseResult> result);

  template<typename T>
  void GetElemResult( UInt sequenceStep,
      UInt stepValue,
      shared_ptr<BaseResult> result);
  
  void FillResultMap(SolutionType sType, ResultInfo::EntryType eType);

  //! creates a map for CFS<->Ensight results
  virtual void FillResultMap();

  virtual void DisableAllArrays();

  //! function for rounding all values to 6 digits
  Double dround(Double a, int ndigits);

  std::map<SolutionType,ResultDef> cfsEnsightResMap_;

  //! store the number of available timesteps
  UInt numAvailSteps_;

  boost::unordered_map<UInt,Double> elemVolumes_;

  //! set of available results
  std::set<SolutionType> availResults_;

  //! neutral TimeValue
  Double neutralTime_;

};

}


#endif /* SIMINPUTENSIGHT_HH_ */
