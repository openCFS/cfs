// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SimInputVTKBased.hh
 *       \brief    <Description>
 *
 *       \date     07.02.2014
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SIMINPUTVTKBASED_HH_
#define SIMINPUTVTKBASED_HH_

#include "DataInOut/SimInput.hh"
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkUnstructuredGrid.h>
#include "vtkCompositeDataIterator.h"

// FV stuff
#include <unordered_set>
#include <unordered_map>

namespace CoupledField{

//! This is a special class used to identify a finite volume face by its nodes.
//! The master and slave elements are equal, if it is a boundary element.
//! The two faces created from adjacent cells shoudl have the same nodes in reverse order (see EqualTest)
class VTKFVFace {

public:

  //! Contructor using master Element number and nodes
  VTKFVFace(UInt masterElement, StdVector<UInt>& points);
 
  VTKFVFace();
  
  virtual ~VTKFVFace();
  
  //! returns true if faces use same nodes in same or reverse order
  bool operator==(const VTKFVFace &other) const;
  
  //! Returns 1, if the other face contains the same nodes in same order
  //! -1 if the other face contains teh same nodes in reverse order
  //! 0 else
  Integer EqualTest(const VTKFVFace &other) const;
  
  UInt GetMasterElement();
  
  UInt GetSlaveElement();
  
  void SetSlaveElement(UInt slaveElement);

  bool HasSlaveElement();
  
  UInt GetPointCount();
  
  void GetPoints(StdVector<UInt>& points);
  
  bool HasPoint(UInt point);
  
  //! returrns a hash code unique for faces with the same node numers
  size_t GetHash() const;
  
private: 

  UInt masterElement_;

  UInt slaveElement_;
  
  StdVector<UInt> points_;
  
  std::size_t hash_;

};

}

// for use of VTKFVFace in unordered_set or unordered_map
namespace std {

template<>
struct hash<CoupledField::VTKFVFace>
{
  size_t operator()(const CoupledField::VTKFVFace & obj) const {
        return obj.GetHash();
  }
};
}
// FV stuff end

namespace CoupledField{

class SimInputVTKBased : public SimInput {


public:


  SimInputVTKBased(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode )
  : SimInput(fileName, inputNode,  infoNode){
    mi_ = NULL;

    readFVMesh_ = inputNode->Get("readFVMesh")->As<bool>();
    fixFVPyramids_ = inputNode->Get("fixFVPyramids")->As<bool>();
//    readFVMesh_ = true;
//    fixFVPyramids_ = true;
    numFVTestedFaces_ = 0;
    numFVInternalFaces_ = 0;
    numFVBoundaryFaces_ = 0;
  }

  virtual ~SimInputVTKBased(){

  }

  virtual void InitModule();

  virtual void ReadMesh( Grid *mi );

  // ========================================================================
  // GENERAL MESH INFORMATION
  // ========================================================================
  //@{ \name General Mesh Information

  //! \copydoc SimInput::GetDim
  virtual UInt GetDim();

  //! \copydoc SimInput::GetFileName
  const std::string& GetFileName() const { return fileName_; }

  //! \copydoc SimInput::GetNumNodes
  virtual UInt GetNumNodes(){
    return numNodes_;
  };

  //! \copydoc SimInput::GetNumElems
  virtual UInt GetNumElems( const Integer dim = -1 ){
    return numElems_;
  }

  //! \copydoc SimInput::GetNumRegions
  virtual UInt GetNumRegions(){
    return numRegions_;
  }

  //! \copydoc SimInput::GetNumNamedNodes
  virtual UInt GetNumNamedNodes(){
    EXCEPTION("GetNumNamedNodes not supported by VTK Readers")
    return 0;
  }

  //! \copydoc SimInput::GetNumNamedElems
  virtual UInt GetNumNamedElems(){
    EXCEPTION("GetNumNamedElems not supported by VTK Readers")
    return 0;
  }

  //@}

  // =========================================================================
  // ENTITY NAME ACCESS
  // =========================================================================
  //@{ \name Entity Name Access

  //! \copydoc SimInput::GetAllRegionNames
  virtual void GetAllRegionNames( StdVector<std::string> & regionNames ){
    std::map<UInt, StdVector<std::string> >::iterator mapiter = regionNamesOfDim_.begin();
    while(mapiter != regionNamesOfDim_.end()){
      for(UInt i=0;i<mapiter->second.GetSize();i++)
      regionNames.Push_back(mapiter->second[i]);
    }
  }

  //! \copydoc SimInput::GetRegionNamesOfDim
  virtual void GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                    const UInt dim ){
    regionNames = regionNamesOfDim_[dim];
  }
  //! \copydoc SimInput::GetNodeNames
  virtual void GetNodeNames( StdVector<std::string> & nodeNames ){
    EXCEPTION("GetNodeNames not supported by VTK Readers")
  }

  //! \copydoc SimInput::GetElemNames
  virtual void GetElemNames( StdVector<std::string> & elemNames ){
    EXCEPTION("GetElemNames not supported by VTK Readers")
  }

  //@}

  // =========================================================================
  // GENERAL SOLUTION INFORMATION
  // =========================================================================
  //@{ \name General Solution Information

  //! \copydoc SimInput::GetNumMultiSequenceSteps
  virtual void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                             std::map<UInt, UInt>& numSteps,
                                             bool isHistory = false ){
    //VTK Files can only store one multisequence step?!
    //right now we hard-code to transient results
    analysis[1] = BasePDE::TRANSIENT;
    numSteps[1] = timeValues_.size();
  }

  //! \copydoc SimInput::GetResultTypes
  virtual void GetResultTypes( UInt sequenceStep,
                                   StdVector<shared_ptr<ResultInfo> >& infos,
                                   bool isHistory = false )=0;

  //! \copydoc SimInput::GetElemNames
  virtual void GetStepValues( UInt sequenceStep,
                                  shared_ptr<ResultInfo> info,
                                  std::map<UInt, Double>& steps,
                                  bool isHistory = false ){
    for(UInt i=0;i<timeValues_.size();i++){
      steps[i] = timeValues_[i];
    }
  }

  //! \copydoc SimInput::GetResultEntities
  virtual void GetResultEntities( UInt sequenceStep,
                                      shared_ptr<ResultInfo> info,
                                      StdVector<shared_ptr<EntityList> >& list,
                                      bool isHistory = false)=0;

  //! \copydoc SimInput::GetResult
  virtual void GetResult( UInt sequenceStep,
                             UInt stepValue,
                             shared_ptr<BaseResult> result,
                             bool isHistory = false )=0;

  //@}

protected:
  //! Read available timesteps
  virtual void GetTimeValues()=0;

  //! Create the file Reader
  virtual void CreateReader()=0;

  //! Set the specified time value to the reader
  virtual void SetTimeValue(Double val)=0;

  //! Enable reading of Regions
  virtual void EnableRegions()=0;

  virtual void GetPointIDsByRegion(std::map<std::string, StdVector<UInt> >& pointMap);

  virtual void ReadElemData(std::map<std::string, StdVector<UInt> >& pointMap);

  //! determine the name of the current region, pass an index as fallback to create a region name
  //! The iterator is NOT incremented here
  virtual void GetRegionName(std::string& name,vtkCompositeDataIterator* iter, UInt idx);


  //! ConvertCells which are not FE Standard
  virtual UInt ComputeNumElems(vtkCompositeDataIterator* iter);


  Elem::FEType VTKCellTypeToFEType(UInt cellType);

  // Pointer for storing the reader
  vtkMultiBlockDataSetAlgorithm* reader_;

  //store available timesteps
  std::vector<Double> timeValues_;

  //! store number of nodes
  UInt numNodes_;

  //! store total number of elements
  UInt numElems_;

  //! store number of regions
  UInt numRegions_;

  //! associate RegionIdTypes to block indices of vtk dataset
  std::map< RegionIdType, UInt> regionAssoc_;

  //! associate global element number (index) to region local element number (value)
  std::unordered_map<UInt,UInt> globElemLocElem_;

  //! associate global node number (index) to region local node number (value)
  StdVector<UInt> globNodeLocElem_;


  // =======================================================================
  // FINITE VOLUME REPRESENTATION SECTION
  // =======================================================================

  //! true if the finite volume mesh should be read
  bool readFVMesh_;
  
  //! true if some pyramids need to be fixed
  bool fixFVPyramids_;

  //! Used for faces
  std::unordered_set<VTKFVFace> FVFaces_;

  //! Some counter
  UInt numFVTestedFaces_;
  //! Some counter
  Integer numFVInternalFaces_;
  //! Some counter
  Integer numFVBoundaryFaces_;

  void AddFVEdge(UInt masterElement, UInt pA, UInt pB);
  
  void AddFVTriFace(UInt masterElement, UInt pA, UInt pB, UInt pC);
  
  void AddFVQuadFace(UInt masterElement, UInt pA, UInt pB, UInt pC, UInt pD);
  
  void AddFVFace(UInt masterElement, StdVector<UInt>& points);
  
};


}


#endif /* SIMINPUTVTKBASED_HH_ */
