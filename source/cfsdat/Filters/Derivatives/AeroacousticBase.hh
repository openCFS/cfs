// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AeroacousticBase.hh
 *       \brief    <Description>
 *
 *       \date     Jan 7, 2017
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include <Filters/MeshFilter.hh>
#include "DataInOut/SimInput.hh"


namespace CFSDat{

//! Base class for aeroacoustic computations, e.g. source-term computations.
//! Regarding the physical quantities, we have some standards, which must be
//! fulfilled:
//! -) Velocity input (e.g. from CFD) must have the same dimension as the grid
//! -) Velocity input must be defined on nodes
//! -) Vorticity input (if provided) for a 2D grid must be scalar-valued, altough it physically
//!    points in z-direction (if grid is x-y). For a 3D grid, the vorticity is also 3D
//! It has some similar methods as MeshFilter but due to the fact that we carry out several
//! operations in one filter, there are some adaptations to these, especially for the
//! derivatives
class AeroacousticBase : public MeshFilter{

public:

  struct PhysicalEntity{
    Vector<Double> U; //Velocity on nodes
    Vector<Double> Om; //Vorticity on nodes
    Vector<Double> Uelem; //Velocity on elements
    CF::StdVector< Vector<Double> >  scatteredData;
    UInt dim; //1 = scalar, 2 = 2D, 3 = 3D
    bool externVorticity; // true if an extern vorticity is provided

    PhysicalEntity() : dim(0), externVorticity(false){
      U.Resize(0);
      Om.Resize(0);
      Uelem.Resize(0);
      scatteredData.Resize(0);
    }

  };


  AeroacousticBase(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~AeroacousticBase();

  virtual bool Run() = 0;



protected:

  virtual void PrepareCalculation() = 0;

  virtual ResultIdList SetUpstreamResults() = 0;

  virtual void AdaptFilterResults() = 0;


  //***************************************************************************
  //********************** Operations with physical quantities    *************
  //***************************************************************************

  //TODO Make the arguments more uniform, e.g. one method computes the scatteredData
  //internally, for the other one it has to be provided

  //! Computes the vector product of vorticity and velocity
  //! The problem is that the input-results, provided by result-manager, is
  //! only a vector containing the results in the form (u1,v1,w1,u2,v2,w2,...)
  //! so we have to take special care of the dimensions

  //! \param data (in) Struct containing the velocity and vorticity
  //!                 (Velocity must be 2D and Vorticity scalar!)
  //! \param derivData (in) Struct containing the information of the entities,
  //!                   for which the computation is performed
  //! \param retVec (out) Vector, containing the vector-product-results (2D vector xy)
  void OmegaVectorProductU_2D(const PhysicalEntity& data,
                              const std::vector<QuantityStruct>& derivData,
                              Vector<Double>& retVec);

  //! Computes the vector product of vorticity and velocity in 3D, similar
  //! to OmegaVectorProductU_2D above

  //! \param data (in) Struct containing the velocity and vorticity
  //!                 (Velocity and vorticity must be 3D!)
  //! \param derivData (in) Struct containing the information of the entities,
  //!                   for which the computation is performed
  //! \param retVec (out) Vector, containing the vector-product-results (3D vector)
  void OmegaVectorProductU_3D(const PhysicalEntity& data,
                              const std::vector<QuantityStruct>& derivData,
                              Vector<Double>& retVec);

  //! Computes the divergence in 2D

  //! \param retVec (out) Vector containing the divergence-values (scalar-valued)
  //! \param inVec (in) Input-vector, must have twice the size as srcCoords (2D-vector)
  //! \param srcCoords (in) Vector of coordinates (node-coordinates)
  //! \param derivData (in)  Struct containing the information of the entities,
  //!                   for which the computation is performed
  void CalcDivergence2D(Vector<Double>& returnVec,
                        const Vector<Double>& inVec,
                        const CF::StdVector< CF::Vector<double> >& srcCoords,
                        const std::vector<QuantityStruct>& derivData);

  //! Computes the divergence in 3D

  //! \param retVec (out) Vector containing the divergence-values (3D vector-valued)
  //! \param inVec (in) Input-vector, must have three times the size of srcCoords (3D-vector)
  //! \param srcCoords (in) Vector of coordinates (node-coordinates)
  //! \param derivData (in)  Struct containing the information of the entities,
  //!                   for which the computation is performed
  void CalcDivergence3D(Vector<Double>& returnVec,
                        const Vector<Double>& inVec,
                        const CF::StdVector< CF::Vector<double> >& srcCoords,
                        const std::vector<QuantityStruct>& derivData){
    EXCEPTION("CalcDivergence3D not implemented yet");
  }



  //! Computes the curl of the given scatteredData

  //! \param retVec (out) Vector containing the curl-values (scalar-valued)
  //! \param scatteredData (in) Vector of scattered data, needed for nearest neighbour search.
  //!                           In the divergence-methods above it's computed internally
  //TODO Calculation of scatteredData should only be carried out, if it's really necessary
  //! \param inVec (in) Input-vector, must have twice the size as srcCoords (2D-vector)
  //! \param srcCoords (in) Vector of coordinates (node-coordinates)
  //! \param derivData (in)  Struct containing the information of the entities,
  //!                   for which the computation is performed
  void CalcCurlU_2D(Vector<Double>& retVec,
                  const CF::StdVector< Vector<Double> >&  scatteredData,
                  const std::vector<QuantityStruct>& derivData,
                  const CF::StdVector< CF::Vector<double> >& coords,
                  const boost::uuids::uuid& resId,
                  const UInt& dim);


  //! NOT TESTED!
  void CalcCurlU_3D(Vector<Double>& retVec,
                  const CF::StdVector< Vector<Double> >&  scatteredData,
                  const std::vector<QuantityStruct>& derivData,
                  const CF::StdVector< CF::Vector<double> >& coords,
                  const boost::uuids::uuid& resId,
                  const UInt& dim);


  //! Computes the gradient of the dot-product u*u

  //! \param retVec (out) Vector containing the gradient (vector-valued)
  //! \param data (in) Struct containing the velocity and vorticity
  //! \param derivData (in)  Struct containing the information of the entities,
  //!                   for which the computation is performed
  //! \param coords (in) Vector of coordinates (node-coordinates)
  //! \param resId (in) Result id of the current filter, needed to get the
  //!                   correct equations (for elements or nodes, as defined in
  //!                   AdaptFilterResults() )
  void CalcGradUScalarU(Vector<Double>& retVec,
                            const PhysicalEntity& data,
                            const std::vector<QuantityStruct>& derivData,
                            const CF::StdVector< CF::Vector<double> >& coords,
                            const boost::uuids::uuid& resId);


  Vector<Double> ScalarToTwoD(const Vector<Double>& inVec);

  Vector<Double> TwoDToScalar(const Vector<Double>& inVec);

};



}
