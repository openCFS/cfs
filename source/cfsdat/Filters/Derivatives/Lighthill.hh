// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Lighthill.hh
 *       \brief    <Description>
 *
 *       \date     Oct 4, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include <Filters/Derivatives/AeroacousticBase.hh>
#include "DataInOut/SimInput.hh"


namespace CFSDat{

/*********************************************************************************
 * BEFORE IMPLEMENTING PLEASE READ INFORMATION AT AdaptFilterResults() down below
 *********************************************************************************/

//! Class which provides methods to compute the Lamb- and LighthillSource-vector as well as
//! the whole LighthillSource-term
class Lighthill : public AeroacousticBase{

public:

  Lighthill(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~Lighthill();

  virtual bool Run();



protected:

  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();


  //! -) AeroacousticSource_LambVector: output is vector-valued and has the dimension of the grid (2 for 2D, 3 for 3D)
  //!    -) externVorticity == true: NODE-RESULTS
  //!    -) externVorticity == false: ELEMENT-RESULTS
  //! -) AeroacousticSource_LighthillSourceVector: output is vector-valued and has the dimension of the grid (2 for 2D, 3 for 3D)
  //!    -) ELEMENT-RESULTS no matter if externVorticity is provided or not
  //! -) AeroacousticSource_LighthillSourceTerm (only 2D yet): physically it's a scalar but due to the
  //!                                                          result-managing it is vector valued with the scalar
  //!                                                          quantity in x-direction. to extract it and get a real
  //!                                                          scalar value, use the filter PostLighthillSourceTerm
  //!    -) ELEMENT-RESULTS no matter if externVorticity is provided or not
  virtual void AdaptFilterResults();

  std::string res1Name;
  uuids::uuid res1Id;

  std::string res2Name;
  uuids::uuid res2Id;


private:

  void LambVector(Vector<Double>& tempRetVec);

  void LighthillSourceVector(Vector<Double>& tempRetVec);

  void LighthillSourceTerm(Vector<Double>& tempRetVec, const str1::shared_ptr<EqnMapSimple>& upMap);

  //! data struct for node-values
  std::vector<QuantityStruct> derivDataNode_;

  //! data struct for element-values
  std::vector<QuantityStruct> derivDataElem_;


  //! Coordinates of input data
  CF::StdVector< CF::Vector<double> > sourceCoords_;

  //! Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;

  //! Density, if not specified in xml-scheme it is automatically set to one
  Double density_;

  //! String if the full Lighthill or only the Lamb-vector is computed
  std::string Form_;

  //! Boolean if an extern vorticity-input is provided of if we have to compute it
  bool externVorticity_;

  PhysicalEntity data_;

};



}
