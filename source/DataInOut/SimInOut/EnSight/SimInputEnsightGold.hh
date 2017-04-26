// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SimInputEnsightGold.hh
 *       \brief    <Description>
 *
 *       \date     Sep 17, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SIMINPUTENSIGHTGOLD_HH_
#define SIMINPUTENSIGHTGOLD_HH_

#include "DataInOut/SimInput.hh"

namespace CoupledField{

//! Class for reading unstructured simulation grids and CFD data from EnSight Gold (Binary and ASCII) files.

//! Some Comments on the Ensight Gold Format
//! - it is part based, meaning each part has its own coordinate and node array. Thereby the meshes are discontinuous by definition
//!   No FEM computation can be conducted on
class SimInputEnSightGold : public SimInput {
 public:

  //! Constructor
  SimInputEnSightGold( std::string fileName, PtrParamNode inputNode,
                       PtrParamNode infoNode );

  virtual ~SimInputEnSightGold();

  virtual void InitModule();

  virtual void ReadMesh(Grid* ptGrid);

  // ========================================================================
  // GENERAL MESH INFORMATION
  // ========================================================================
  //@{ \name General Mesh Information

  //! \copydoc SimInput::GetDim
  virtual UInt GetDim();

  //! \copydoc SimInput::GetNumNodes
  virtual UInt GetNumNodes();

  //! \copydoc SimInput::GetNumElems
  virtual UInt GetNumElems( const Integer dim = -1 );

  //! \copydoc SimInput::GetNumRegions
  virtual UInt GetNumRegions();

  //! \copydoc SimInput::GetNumNamedNodes
  virtual UInt GetNumNamedNodes();

  //! \copydoc SimInput::GetNumNamedElems
  virtual UInt GetNumNamedElems();

  //@}

  // =========================================================================
  // ENTITY NAME ACCESS
  // =========================================================================
  //@{ \name Entity Name Access

  //! \copydoc SimInput::GetAllRegionNames
  virtual void GetAllRegionNames( StdVector<std::string> & regionNames );

  //! \copydoc SimInput::GetRegionNamesOfDim
  virtual void GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                    const UInt dim );

  //! \copydoc SimInput::GetNodeNames
  virtual void GetNodeNames( StdVector<std::string> & nodeNames );

  //! \copydoc SimInput::GetElemNames
  virtual void GetElemNames( StdVector<std::string> & elemNames );

  //@}

  // =========================================================================
  // GENERAL SOLUTION INFORMATION
  // =========================================================================
  //@{ \name General Solution Information

  //! \copydoc SimInput::GetNumMultiSequenceSteps
  virtual void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                         std::map<UInt, UInt>& numSteps,
                                         bool isHistory = false );

  //! \copydoc SimInput::GetResultTypes
  virtual void GetResultTypes( UInt sequenceStep,
                               StdVector<shared_ptr<ResultInfo> >& infos,
                               bool isHistory = false );

  //! \copydoc SimInput::GetStepValues
  virtual void GetStepValues( UInt sequenceStep,
                              shared_ptr<ResultInfo> info,
                              std::map<UInt, Double>& steps,
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
  //@}

};

}



#endif /* SIMINPUTENSIGHTGOLD_HH_ */
