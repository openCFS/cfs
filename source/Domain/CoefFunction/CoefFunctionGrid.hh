// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     CoefFunctionGrid.hh 
 *       \brief    Coefficient function which obtains and interpolates values from another Grid
 *
 *       \date     01/23/2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef COEFFUNCTIONGRID_HH
#define COEFFUNCTIONGRID_HH

#include "CoefFunction.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "FeBasis/FeFunctions.hh"
#include <boost/tr1/type_traits.hpp>
#include <string>

namespace CoupledField{

/*! \class CoefFunctionGridBase
 *     \brief Coefficient function which obtains and interpolates values from another Grid
 *     \tparam DATA_TYPE Can be Complex or Double
 *     @author A. Hueppe
 *     @date 01/2012
 *
 *     This function handles results from other grids the constructor takes for now only a pointer
 *     to an grid xml node which is evaluated. Based on this information it triggers the appropriate
 *     ResultHandler to provide an FeFunction defined on the external Grid in order to have easy access
 *     to the field values.
 *     Basic Functionality:
 *      - In the GetVector/Scalar function we ask the external grid to provide the containing element of the
 *          given local point
 *      - Based on this source element we evaluate an indentity operator to interpolate to the requested point
 *
 *     Future plans:
 *      \li Make the Operator variable to directly compute e.g. gradients during field interpolation
 *      \li Provide the possibility to map the coefficients to an entire FeFunction using either normal
 *          or conservative interpolation
 *      \li Perform temporal interpolation of the given field values
 */
template<class DATA_TYPE>
class CoefFunctionGridBase : public CoefFunction{
  public:

    CoefFunctionGridBase(PtrParamNode configNode);

    virtual ~CoefFunctionGridBase();

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

    virtual void AddEntities(shared_ptr<EntityList> ent){
      if(!this->entities_.Contains(ent)){
        entities_.Push_back(ent);
        //domain->GetGrid(gridId_)->ComputeContainingElements(ent,srcRegions_);
      }else{
        WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
      }
    }

    //! \copydoc CoefFunction::GetVecSize
    virtual UInt GetVecSize() const;

    //! \copydoc CoefFunction::GetTensorSize
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;

    virtual void Recalculate();

    virtual UInt GetStepNum();

  protected:
    //! reads information from configNode
    virtual void ReadXmlNode(PtrParamNode configNode);

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

    bool isStatic_;

};

}
#endif
