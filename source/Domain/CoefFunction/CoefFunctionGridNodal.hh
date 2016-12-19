// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     CoefFunctionGridNodal.hh 
 *       \brief    Coefficient function which obtains and interpolates values from another Grid
 *                 Specialized for nodal Results from the other grid
 *
 *       \date     11/21/2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef COEFFUNCTIONGRIDNODAL_HH
#define COEFFUNCTIONGRIDNODAL_HH

#include "CoefFunctionGrid.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "DataInOut/ResultHandler.hh"

namespace CoupledField{

/*! \class CoefFunctionGridBase
 *     \brief Coefficient function which obtains and interpolates values from another Grid
 *     \tparam DATA_TYPE Can be Complex or Double
 *     @author A. Hueppe
 *     @date 11/2012
 *     
 *     This CoefFunction handels nodal resutls from another Grid. 
 *     The Implementation is specializied thus beeing quite efficient but inflexible wrt 
 *     to higher order results. The class CoefFunctionGridHigher will cope with this. 
 */
template<class DATA_TYPE>
class CoefFunctionGridNodal : public CoefFunctionGrid{
  public:

    ///Constructor with the configuration xml node reading the type of interpolation
    CoefFunctionGridNodal(Domain* ptDomain, PtrParamNode configNode, shared_ptr<RegionList> regions);

    ///Destructor freeing al used data strutures
    virtual ~CoefFunctionGridNodal();

    // ========================
    //  ACCESS METHODS
    // ========================
    //@{ \name Access Methods

    ///Getting a Tensor Value at the requested lpm, INvalid call for Conservative interpolation
    virtual void GetTensor(Matrix<DATA_TYPE>& CoefMat,
                           const LocPointMapped& lpm );

    ///Getting a vector Value at the requested lpm, Invalid call for Conservative interpolation
    virtual void GetVector(Vector<DATA_TYPE>& CoefMat,
                           const LocPointMapped& lpm );

    ///Getting a scalar value at the requested lpm, INvalid call for Conservative interpolation
    virtual void GetScalar(DATA_TYPE& CoefMat,
                           const LocPointMapped& lpm );
    //@}

    //! Dump coefficient function to string
    virtual std::string ToString() const =0 ;

    //! \copydoc CoefFunction::SetDerivativeOperation
    virtual void SetDerivativeOperation(CoefDerivativeType type){
      this->derivType_ = type;

      //make some checks here!
      switch(dimType_){
      case SCALAR:
        //only NONE is valid right now
        //if extended to gradient, this would be fine too
        if(type==VECTOR_DIVERGENCE){
          EXCEPTION("CoefFunctionExpression: VECTOR_DIVERGENCE is not a valid operator for scalar coefFunction");
        }
        break;
      case VECTOR:
        //this is fine in all cases right now
        if(type==VECTOR_DIVERGENCE){
          //change dim type to scalar
          this->dimType_ = SCALAR;
          UInt gDim = this->domain_->GetGrid()->GetDim();
          UInt dDim = this->resultInfo_->dofNames.GetSize();
          this->CreateDivOperator(gDim,dDim);
        }
        break;
      case TENSOR:
        if(type==VECTOR_DIVERGENCE){
          EXCEPTION("CoefFunctionExpression: VECTOR_DIVERGENCE is not a valid operator for tensor coefFunction");
        }
        break;
      default:
        break;
      }
      return;
    }

  protected:

    //=======================================================================================
    // Interpolation operator for nodal results
    //=======================================================================================
    //!Create the interpolation operator
    //!and sets the corresponding class variable
    void CreateOperator(UInt spaceDim, UInt dofDim);

    //!Creates a divergence operator in case we want the divergence of an
    //! external vector field
    void CreateDivOperator(UInt spaceDim, UInt dofDim);

    //! BOperator to map solutions to arbitrary points
    //! Right now, hardcoded identity operator
    shared_ptr<BaseBOperator > myOperator_;

    //=======================================================================================
    // Functions and structure for nodal grid result handling 
    //=======================================================================================

    //! initalize region nodes
    void InitRegionNodes();
    
    //! stores if region nodes have been initialized
    bool initializedRegionNodes_;
    
    //! cached data of nodes per region
    StdVector< StdVector<UInt> > regionNodes_;
    
    //! initialize used region nodes for evaluting the sum
    void InitUsedRegionNodesForSum();
    
    //! stores if the nodes needes for evaluating the sum for each source are initialized
    bool initializedUsedRegionNodesForSum_;
    
    //! stores, if all region nodes should be used for evaluating the sum
    bool useAllRegionNodesForSum_;
    
    //! stores the nodes to take into account for evaluating the sum of sources
    StdVector< std::vector<bool> > usedRegionNodesForSum_;
    
    //! Read solution from sourceFile according to the given stepnumber
    void ReadSolution(UInt step,Vector<DATA_TYPE> & sol);
    
    //! Updates the solution vector
    bool UpdateSolution();

    //! Perform a simple equation mapping for nodal grids
    //! Initialize the solution vector solVec_;
    void InitSolVec();

    //! Extract the solution for a source element in order to apply the interpolation operator
    void GetElemSolution(Vector<DATA_TYPE> & sol, UInt eNum);

    //! Returns the step to be read in along with temporal interpolation factors
    //! If those are necessary
    UInt GetStepNum(bool &interpolateT,Double & iFactor1, Double & iFactor2);

    //! Stores the current solution vector
    Vector<DATA_TYPE> solVec_;

    //! stores the next solution vector for temporal interpolation
    Vector<DATA_TYPE> solVecFuture_;

    //! stores the solution vector at previuous timesteps
    Vector<DATA_TYPE> solVecOld_;

    //! Total number of nodes
    UInt numNodes_;
    
    //! Total number of equations
    UInt numEqns_;

    //! Pointer to math parser instance
    MathParser* mp_;

    //! Handle for expression determines current time/freq value 
    MathParser::HandleType mHandleStep_;
    
    //! stores the stepnumber of the last read solution process
    UInt lastStepRead_;

    //! flag indicating if the timevalue map of input should be ignored
    bool snapToCFSStep_;
    
    //! flag to verbose the sum of the sources
    bool verboseSum_;

    //=======================================================================================
    // Functions and structure for global factor handling
    //=======================================================================================

    //! initalize region node coordinates
    void InitRegionNodeCoordinates();
    
    //! clear region node coordinates
    void ClearRegionNodeCoordinates();
    
    //! stores if region nodes have been initalized
    bool initializedRegionNodeCoordinates_;
    
    //! cached region node coordinates
    StdVector<StdVector<Vector<Double> > > regionNodeCoordinates_;

    //! create one factor function from one factor string 
    shared_ptr<CoefFunction> CreateFactorFunction(std::string factorString);
    
    //! get the dependency type of a given factor string
    CoefFunction::CoefDependType ReadFactorFunctionStringDependency(std::string factorString);

    //! Initialize global factor functions
    void InitGlobalFactorFunctions(PtrParamNode configNode);
    
    //! Write the global factor settings to info XML file
    void WriteGlobalFactorsToXML(PtrParamNode configNode);
    
    
    //! stores if any global factor has beend set
    bool hasFactor_;
    
    
    //! stores if one constant factor is present
    bool hasConstantFactor_;
    
    //! stores the constant global factor (not dependend on time and space) to be evalated
    DATA_TYPE constantFactor_;
    
    
    //! stores if a factor solely depending on space is present
    bool hasSpaceFactor_;
    
    //! factor function depending solely on space
    shared_ptr<CoefFunction> spaceFactorFunction_;
    
    //! initalize space factor
    void InitSpaceFactor();
    
    //! stores if the space factors have already been initialized
    bool initializedSpaceFactors_;
    
    //! cached space factor values for each region
    StdVector< StdVector< DATA_TYPE > > spaceFactor_;
    
    
    //! stores if a factor solely depending on time/frequency is present
    bool hasTimeFreqFactor_;
    
    //! factor function depending solely on time/frequency
    shared_ptr<CoefFunction> timeFreqFactorFunction_;
    
    //! stores if factor values depending solely on time/frequency is present should be verbosed in each time step
    bool verboseTimeFreqFactor_;
    
    
    //! stores if a general factor depending on time/frequency and space is present
    bool hasGeneralFactor_;
    
    //! factor function depending on time/frequency and space
    shared_ptr<CoefFunction> generalFactorFunction_;

    //=======================================================================================
    // Functions and structure for copying the nodal grid result
    //=======================================================================================
    
    //! struct needed for creation of typedef to template function for copying the result vector
    struct CopyResultFunction {
      //! type definition for pointer to CopyResult function
      //typedef void (CoefFunctionGridNodal<DATA_TYPE>::*blafunction)(Vector<DATA_TYPE>&, Vector<DATA_TYPE>&, StdVector<UInt>&,
      typedef void (*Ptr)(Vector<DATA_TYPE>&, Vector<DATA_TYPE>&, StdVector<UInt>&, 
                                StdVector<DATA_TYPE>&, StdVector<DATA_TYPE>&, const DATA_TYPE,
                                StdVector<DATA_TYPE>&, StdVector<DATA_TYPE>&, std::vector<bool>&);
    };
    
    //! Copy the data from the result vector to the solution vector, apply factors and calculate sum
    //!
    //! This function uses templates, so unneeded branches will be optimized easily.
    //! Otherwise one would have to trust the compiler, that all boolean options will be moved to switches 
    //! outside the copy loop, which may not appear due to many loops used. 
    //! For GNU compiler this should usually be done by funswitch-loop options
    //! The setting of dimDof to a template value furthermore allows the compiler the unrolling of 
    //! small inner loops
    template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor, 
              bool countSum, bool countAllValues, UInt dimDof>
    static void CopyResult(Vector<DATA_TYPE>& sol,
                          Vector<DATA_TYPE>& res,
                          StdVector<UInt>& nodeNums, 
                          StdVector<DATA_TYPE>& spaceFactorA, 
                          StdVector<DATA_TYPE>& spaceFactorB,
                          const DATA_TYPE constFactor,
                          StdVector<DATA_TYPE>& sum, 
                          StdVector<DATA_TYPE>& factorSum,
                          std::vector<bool>& countNodes);
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor, 
        bool countSum, bool countAllValues, UInt dimDof>
    static typename CoefFunctionGridNodal<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction() {
      return &CoefFunctionGridNodal<DATA_TYPE>::CopyResult<
          useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,dimDof>;
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor, 
        bool countSum, bool countAllValues>
    static typename CoefFunctionGridNodal<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(const UInt dimDof) {
      if (dimDof == 1) { // scalar
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,1>();
      } else if (dimDof == 2) { // 2D vector
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,2>();
      } else if (dimDof == 3) { // 3D vector and 2D symmtric tensor
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,3>();
      } else if (dimDof == 4) { // 2D tensor
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,4>();
      } else if (dimDof == 6) { // 3D symmetric tensor
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,6>();
      } else if (dimDof == 8) { // 2D tesnor 3rd order 
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,8>();
      } else if (dimDof == 9) { // 3D tensor
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,9>();
      } else if (dimDof == 27) { // 3D tesnor 3rd order 
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,27>();
      }
      EXCEPTION("Cannot find CopyResult function for dimension " << dimDof << " in CoefFunctionGridNodal::GetCopyResultFunction");
      return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,1>();
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor>
    static typename CoefFunctionGridNodal<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(
        const bool countSum, const bool countAllValues, const UInt dimDof) {
      if (countSum) {
        if (countAllValues) {
          return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,true,true>(dimDof);
        } else {
          return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,true,false>(dimDof);
        }  
      }
      return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,false,false>(dimDof);
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA, bool useSpaceFactorB>
    static typename CoefFunctionGridNodal<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(
        const bool useConstFactor, const bool countSum, const bool countAllValues, const UInt dimDof) {
      if (useConstFactor) {
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,true>(countSum, countAllValues, dimDof);
      }
      return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,false>(countSum, countAllValues, dimDof);
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA>
    static typename CoefFunctionGridNodal<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(
        const bool useSpaceFactorB, const bool useConstFactor, 
        const bool countSum, const bool countAllValues, const UInt dimDof) {
      if (useSpaceFactorB) {
        return GetCopyResultFunction<useSpaceFactorA,true>(useConstFactor, countSum, countAllValues, dimDof);
      }
      return GetCopyResultFunction<useSpaceFactorA,false>(useConstFactor, countSum, countAllValues, dimDof);
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    static typename CoefFunctionGridNodal<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(
        const bool useSpaceFactorA, const bool useSpaceFactorB, const bool useConstFactor, 
        const bool countSum, const bool countAllValues, const UInt dimDof) {
      if (useSpaceFactorA) {
        return GetCopyResultFunction<true>(useSpaceFactorB, useConstFactor, countSum, countAllValues, dimDof);
      }
      return GetCopyResultFunction<false>(useSpaceFactorB, useConstFactor, countSum, countAllValues, dimDof);
    }

  private:

};
}

#endif
