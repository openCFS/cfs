// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionNodalGrid.hh
 *       \brief    This Coefficient Function is specially suited for the reading of the standard
 *                 Nodal results we mostly have up to now.
 *                 Right now we choose the "little" hack of creating always this function as
 *                 our result definition is not capable of storing higher order results anyway
 *
 *       \date     Jun 27, 2012
 *       \author   ahueppe
 */
//================================================================================================

#ifndef COEFFUNCTIONNODALGRID_HH_
#define COEFFUNCTIONNODALGRID_HH_

#include "CoefFunctionGrid.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include <set>
#include <string>

namespace CoupledField{
/*! \class CoefFunctionGridBase
 *     \brief Coefficient function which obtains and interpolates values from another Grid nSpecialized for nodal results
 *     \tparam DATA_TYPE Can be Complex or Double
 *     @author A. Hueppe
 *     @date 06/2012
 *
 *     This function handles nodal results from other grids the constructor takes for now only a pointer
 *     to an grid xml node which is evaluated. Based on this information it triggers the appropriate
 *     ResultHandler to provide an FeFunction defined on the external Grid in order to have easy access
 *     to the field values.
 *     Basic Functionality:
 *      - In the GetVector/Scalar function we ask the external grid to provide the containing element of the
 *          given local point
 *      - Based on this source element we evaluate an indentity operator to interpolate to the requested point
 *     To get a better computational speed at the cost of having more memory consumption
 *     we store right no a sort of fake equation map to refer to results on the source grid
 */
template<class DATA_TYPE>
class CoefFunctionNodalGrid : public CoefFunction{

public:

  CoefFunctionNodalGrid(PtrParamNode configNode);


  virtual ~CoefFunctionNodalGrid();


  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  virtual void GetTensor(Matrix<DATA_TYPE>& CoefMat,
                         const LocPointMapped& lpm );

  virtual void GetVector(Vector<DATA_TYPE>& CoefMat,
                         const LocPointMapped& lpm );

  virtual void GetScalar(DATA_TYPE& CoefMat,
                         const LocPointMapped& lpm );
  //@}

  //! Dump coefficient function to string
  virtual std::string ToString() const;

  virtual void Recalculate();


  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;

  virtual bool IsComplex();
protected:
   //! reads information from configNode
   void ReadXmlNode(PtrParamNode configNode);

   //! input reader Id for this function
   std::string inputId_;

   //! input grid Id for this function
   std::string gridId_;

   //! pointer to src grid
   Grid* srcGrid_;

   //! list of regions on given grid the coefficients are defined on
   std::set<RegionIdType> srcRegions_;

   //! the solution type given by this function
   SolutionType solType_;

   //! Pointer to math parser instance
   MathParser* mp_;

   //! Handle for expression e.g. to accomplish automatic updates of coefficients
   MathParser::HandleType mHandleStep_;

   //! BOperator to map solutions to arbitrary points
   //! Right now, hardcoded identity operator
   shared_ptr<BaseBOperator > myOperator_;

   //! pointer to FeFunction defined on external Grid
   shared_ptr<FeFunction<DATA_TYPE> > feFunct_;

   //! pointer to resultInfo
   shared_ptr<ResultInfo> resultInfo_;

   //! external step to t map
   std::map<UInt, Double> stepValueMap_;

   //! stores the current stepnumber of external result
   UInt curStep_;

   //! stores current timestep
   Double curTStep_;

   shared_ptr<CoefFunction> factorFnc_;

   StdVector< Vector<Double> > locCoords_;

   StdVector<UInt> elemNums_;

   //! the local eqn map
   std::map<UInt,UInt> nodeIdxMap_;

   //! the eqnNumbers
   StdVector< StdVector<UInt> > eqnNumbers_;

   UInt numEqns_;

   UInt dofDim_;

   Vector<DATA_TYPE> solVec_;

   Vector<DATA_TYPE> solVecFuture_;

private:
   //! list of regions on given grid the coefficients are defined on
   std::set<std::string> srcRegionNames_;

   void MapEqns();

   void ReadSolution(UInt step,Vector<DATA_TYPE> & sol);

   void GetElemSolution(Vector<DATA_TYPE> & sol, UInt eNum);

   virtual UInt GetStepNum(bool & interpolateT,Double & iFactor);

   virtual DATA_TYPE GetFactor( const LocPointMapped& lpm ){
     DATA_TYPE res = 0;
     this->factorFnc_->GetScalar(res,lpm);
     return res;
   }
};
}

#endif /* COEFFUNCTIONNODALGRID_HH_ */
