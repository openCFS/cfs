// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FEBasedInterpolator.hh
 *       \brief    This is a finte element based interpolator. It works similar as the interpolation
 * in CoefFunctionGridNodel but with extended functionality:
 * 
 *  - the input values are taken either from the nodes or from the elements (centroids)
 * 
 *  - for conservative interpolation the elements in the target grid are found in which the source
 *    points (node coordinates of element centroids) lay in, then the shape function of the element nodes are evaluated at the input point to 
 *    get the interpolation coefficients
 *     - In case of 100% interpolation, also target elements/nodes are searched for, if the sources are not located inside a target element
 *           - if connectivity is on, the neighbours of sources which are evaluated successfully are evaluated and the non-interpolated sources get targets accordingly. The neighbour search is done either based on finite elements (element->node->element) or finite volumes (element->face->element)
 *               - if the distance from source to target exceeds maxDistance only numNeighbours neighbours are taken into account to avoid too big matrices. At such places there should be no relavant data.
 *           - if connectivity is off, numNeigbhours defines how many nearest neighbours are searched for to make the interpolation conservative
 *     - if 100% interpolation and no connectivity is selectet is selected, the nearest nodes to the source points are set for the 
 * 
 *  - for non-conservative interpolation the search is done in reversed order. The target points are
 *    located in the source elements and the shape functions are evaluated from these elements.
 *     - In case of 100% interpolation sorces elements/nodes are searched for, if no source element was found for a target. This can also be done either by connectivity or nearest neighbours
 * 
 * The option exponent regulates the weightings in case of 100& option is needed to find a corresponding source or target. The higher "exponent" is , the higher is the weighting for the near replacements in contrast to the far away ones
 *      
 *       \date     Jun 18, 2018
 *       \author   mtautz
 */
//================================================================================================

#pragma once //instead of the #ifndef #def ...


#include "DataInOut/SimInput.hh"
#include <Filters/MeshFilter.hh>
#include "DatUtils/InterpolationMatrix.hh"
#include "AbstractInterpolator.hh"



namespace CFSDat{

class FEBasedInterpolator : public AbstractInterpolator {

public:

  FEBasedInterpolator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~FEBasedInterpolator();

protected:
  
  virtual void FillMatrix(StdVector<CF::UInt>& globSrcEntity, StdVector<CF::UInt>& globTrgEntity);

private:

  //! Stores if interpolation is done conservatively
  bool conservative_;

  //! Stores if the interpolation is done 100 percently
  bool hundredPercent_;
  
  // use connectivity to achieve 100 percently, elsse use nearest neighbour
  bool connectivity_;
  
  // use finite volume connectivity instead of vertices and elements
  bool useFV_;
  
  //! number of nearest neighbours
  UInt numNeighbours_;
  
  //! exponent to use when computing the factors for 100 percent interpolation
  Double exponent_;
  
  //! maximum distance to use when applying connectivity to switch to nearest neighbour
  Double maxDistance_;
  
  UInt nSteps_;
  
  UInt iSteps_;
  
  //! find some search entitites (search***) inside the elements of a target grid (find***), find replacements if not found and 100 %
  void FindTargetElements(Grid* findGrid, std::set<std::string>& findRegions, 
                          StdVector<CF::UInt>& findEntitites, bool findElems, shared_ptr<EqnMapSimple> findMap,
                          Grid* searchGrid, std::set<std::string>& searchRegions, 
                          StdVector<CF::UInt>& searchEntitites, bool searchElems, shared_ptr<EqnMapSimple> searchMap,
                          StdVector<CF::UInt>& foundElements, StdVector< LocPoint >& locCoords, std::vector<bool>& hasFoundElement,
                          StdVector< StdVector<CF::UInt>* >& foundReplacementEntities, StdVector< StdVector<CF::Double>* >& foundReplacementDistances);
  
  //! copy found elements of one entity to another entity
  void CopyElemeVectorEntries(UInt source, UInt target, StdVector< StdVector<UInt>* >& elements, 
                             std::vector<bool>& hasElementsVector, std::vector<bool>& newAddedElements,
                             UInt& countEls, bool& hasPutSomeNewElements);
  
};

}
