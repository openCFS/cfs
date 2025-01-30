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
#include <string>

namespace CoupledField{

/*! \class CoefFunctionGridBase
 *     \brief Coefficient function which obtains and interpolates values from another Grid
 *     \tparam DATA_TYPE Can be Complex or Double
 *     @author A. Hueppe
 *     @date 01/2012
 *
 *     This function handles results from other grids the constructor takes for now only a pointer
 *     to an grid xml node which is evaluated. Child classes perform the actual task, wrt the fact
 *     which type is needed.
 *     NOTE: Implemented right now are only external Grids providing nodal results!
 */
class CoefFunctionGrid : public CoefFunction{
  public:

    // ========================
    //  ENUMS / TYPEDEFS
    // ========================
    //@{ \name Public type definitions

    //! Enumeration for the Interpolation type
    typedef enum{
      NO_INTERPOLATION, /*!< The source Grid is the default Grid */
      STANDARD,         /*!< Interpolate the Values with given order */
      CONSERVATIVE     /*!< Perform conservative Interpolation */
    } InterpType;
    static Enum<InterpType> InterpType_;
    //@}


    /*! This function generates a Grid Based Coefficient function
     * based on the parameter configuration supplied by the user
     */
    static PtrCoefFct
    Generate( Domain* ptDomain,
              Global::ComplexPart format,
              PtrParamNode infoNode,
              PtrParamNode configNode,
              shared_ptr<RegionList> regions,
              ResultInfo::EntryType type,
              const std::string& resName = "");

    ///Constructor which sets every field to default values
    /// carefully check for each field in an overloaded class
    CoefFunctionGrid(Domain* ptDomain, PtrParamNode configNode, shared_ptr<RegionList> regions);

    ///Destructor
    virtual ~CoefFunctionGrid();

    virtual string GetName() const { return "CoefFunctionGrid"; }

    // ========================
    //  ACCESS METHODS
    // ========================
    //@{ \name Access Methods

    virtual void GetTensor(Matrix<Double>& CoefMat,
                           const LocPointMapped& lpm );

    virtual void GetVector(Vector<Double>& CoefMat,
                           const LocPointMapped& lpm );

    virtual void GetScalar(Double& CoefMat,
                           const LocPointMapped& lpm );

    virtual void GetTensor(Matrix<Complex>& CoefMat,
                           const LocPointMapped& lpm );

    virtual void GetVector(Vector<Complex>& CoefMat,
                           const LocPointMapped& lpm );

    virtual void GetScalar(Complex& CoefMat,
                           const LocPointMapped& lpm );
    //@}


    virtual void AddEntityList(shared_ptr<EntityList> ent){
    	if(!this->entities_.Contains(ent)){
    		entities_.Push_back(ent);
    	}
    	else {
    		WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
    	}
    }

    //! \copydoc CoefFunction::GetVecSize
    virtual UInt GetVecSize() const;

    //! \copydoc CoefFunction::GetTensorSize
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;

    //! Currently only implemented in CoefFunctionGridNodalInterp
    virtual void PrepareInterpolation(){
      EXCEPTION("PrepareInterpolation not implemented here");
    }

    //! Sets the variables for conservative or standard interpolation
    virtual void SetConservative(bool value){
      bool isReset = false;
      std::string interpTyStr;
      if ((this->curInterpType_ == CoefFunctionGrid::CONSERVATIVE) != value
          && this->curInterpType_ != CoefFunctionGrid::NO_INTERPOLATION)
      {
        isReset = true;
      }
      if(value){
        this->curInterpType_ = CoefFunctionGrid::CONSERVATIVE;
      }else{
        this->curInterpType_ = CoefFunctionGrid::STANDARD;
      }
      interpTyStr = CoefFunctionGrid::InterpType_.ToString(this->curInterpType_);
      this->extDataInfo_->Get("interpolation")->Get("type")->SetValue(interpTyStr);
      if(isReset){
        this->extDataInfo_->Get("interpolation")->Get("setBy")->SetValue("CFS++, Overwritten");
      }else{
        this->extDataInfo_->Get("interpolation")->Get("setBy")->SetValue("CFS++");
      }
    }
  protected:

    ///Stores the current type of interpolation
    CoefFunctionGrid::InterpType curInterpType_;

    ///Obtains Result from the domain class according to inputID and
    ///Sequence step
    virtual void DetermineResult(std::string inputID,UInt seqStep);

    void SetEntitiesByRegions(shared_ptr<RegionList> regions){
      Grid* grid = regions->GetGrid();
      StdVector<RegionIdType> regIDs = regions->GetRegionIds();
      for (UInt i = 0; i < regIDs.GetSize(); i++) {
        std::string name = grid->GetRegionName(regIDs[i]);
        
        // entity list
        EntityList::ListType listType = EntityList::ELEM_LIST; 
        if( grid->GetEntityDim( name ) == grid->GetDim() - 1) {
          listType = EntityList::SURF_ELEM_LIST;
        }
        shared_ptr<EntityList> eList = grid->GetEntityList( listType, name );
        
        // add entity list
        if(!this->entities_.Contains(eList)){
          entities_.Push_back(eList);
        }else{
          WARN("entity list " << eList->GetName() << " already contained in CoefFunction")
        }
      }
    }
    
    //! input reader Id for this function
    std::string inputId_;

    //! input grid Id for this function
    std::string gridId_;

    //! pointer to src grid
    Grid* srcGrid_;
    
    //! Pointer to domain
    Domain* domain_;

    //! list of regions on given grid the coefficients are defined on
    std::set<std::string> srcRegions_;
    
    //! string containing the names of all source regions
    std::string allSrcRegionNames_;

    //! name of the solution type given in this function
    std::string solName_;
    
    //! the solution type given by this function
    SolutionType solType_;

    //! pointer to resultInfo
    shared_ptr<ResultInfo> resultInfo_;

    //! Store the dimension of (Vectorial) Result
    UInt dimDof_;

    //! Store the dimension of tensor
    UInt numRowTensor_;
    UInt numColTensor_;

    //! Support of the CoefFunction. Only needed for grid/solution results
    StdVector<shared_ptr<EntityList> > entities_;
    
    //! external step to t map
    std::map<UInt, Double> stepValueMap_;

    //! stores the current stepnumber of external result
    UInt curStep_;

    //! stores the current step value of external result
    Double curTStep_;

    //! stores the current mutlisequence step for external result
    UInt aSeqStep_;

    //!Stores the conguration parameter node
    PtrParamNode myConfigNode_;

    //! Pointer to info node of grid coefFunction
    PtrParamNode extDataInfo_;

    //! Flag indicating verbose level of the coefFunction
    bool verbose_;

};

}
#endif
