// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     CoefFunctionGridElem.hh 
 *       \brief    Coefficient function which obtains and interpolates values from another grid
 *                 Specialized for element Results from the other grid
 *
 *       \date     20/09/2024
 *       \author   Dominik Mayrhofer
 */
//================================================================================================

#ifndef COEFFUNCTIONGRIDELEM_HH
#define COEFFUNCTIONGRIDELEM_HH

#include "CoefFunctionGrid.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "DataInOut/ResultHandler.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif


namespace CoupledField{

/*! \class CoefFunctionGridBase
 *     \brief Coefficient function which obtains and interpolates values from another Grid
 *     \tparam DATA_TYPE Can be Complex or Double
 *     @author Dominik Mayrhofer
 *     @date 09/2024
 *     
 *     This CoefFunction handels element resutls from another Grid. 
 *     Currently, only the default grid is supported. 
 *     Nevertheless, the implementation is designed similarly to the nodal one to make future implementations easier.
 */
template<class DATA_TYPE>
class CoefFunctionGridElem : public CoefFunctionGrid{
  public:

    ///Constructor with the configuration xml node reading the type of interpolation
    CoefFunctionGridElem(Domain* ptDomain, PtrParamNode configNode, shared_ptr<RegionList> regions);

    ///Destructor freeing al used data strutures
    virtual ~CoefFunctionGridElem();

    virtual string GetName() const { return "CoefFunctionGridElem"; }

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

    virtual void AddEntityList(shared_ptr<EntityList> ent){
    	if(!this->entities_.Contains(ent)){
    		entities_.Push_back(ent);
    	}
    	else {
    		WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
    	}
    }


  protected:

    //=======================================================================================
    // Functions and structure for element grid result handling 
    //=======================================================================================

    //! initalize region elements
    void InitRegionElements();
    
    //! stores if region elements have been initialized
    bool initializedRegionElements_;
    
    //! cached data of elements per region
    StdVector< StdVector<UInt> > regionElements_;

    //! initialize used region elements for evaluting the sum
    void InitUsedRegionElementsForSum();

    //! stores if the elements needes for evaluating the sum for each source are initialized
    bool initializedUsedRegionElementsForSum_;
    
    //! stores, if all region elements should be used for evaluating the sum
    bool useAllRegionElementsForSum_;
    
    //! stores the elements to take into account for evaluating the sum of sources
    StdVector< std::vector<bool> > usedRegionElementsForSum_;
    
    //! Read solution from sourceFile according to the given stepnumber and stepvalue
    void ReadSolution(UInt step, Double stepValue, Vector<DATA_TYPE> & sol);
    
    //! Read solution from sourceFile according to the given stepnumber
    void ReadSolution(UInt step,Vector<DATA_TYPE> & sol);
    
    //! Updates the solution vector
    bool UpdateSolution();

#ifdef USE_OPENMP
    //! thread locking for UpdateSolution function
    omp_lock_t updateSolutionLock_;
#endif
    
    //! Extract the solution for a source element
    void GetElemSolution(Vector<DATA_TYPE> & sol, UInt eNum);
    
    //! Initialize the solution vector solVec_;
    void InitSolVec();

    //! Stores the current solution vector
    Vector<DATA_TYPE> solVec_;

    //! stores the step number of a solution vector for temporal interpolation
    UInt stepNumberInterpolationA_;
    
    //! stores the a solution vector for temporal interpolation
    Vector<DATA_TYPE> solVecInterpolationA_;

    //! stores the step number of a second solution vector for temporal interpolation
    UInt stepNumberInterpolationB_;
    
    //! stores the a second solution vector for temporal interpolation
    Vector<DATA_TYPE> solVecInterpolationB_;

    //! Total number of elements
    UInt numElements_;
    
    //! Total number of nodes
    UInt numNodes_;
    
    //! Total number of equations
    UInt numEqns_;

    //! Equation mapping completed
    bool eqnMapComplete_;

    //! Pointer to math parser instance
    MathParser* mp_;

    //! Handle for expression determines current time/freq value 
    unsigned int mHandleStep_;
    
    //! shared pointer to coefFunction representing multiplicative factor
    shared_ptr<CoefFunction> factorFnc_;

    //! stores the stepnumber of the last update solution process
    UInt lastStepUpdate_;
    
    //! flag indicating if the timevalue map of input should be ignored
    bool snapToCFSStep_;
    
    //! flag to verbose the sum of the sources
    bool verboseSum_;

    //=======================================================================================
    // Functions and structure for global factor handling
    //=======================================================================================

    //! initalize region element coordinates
    void InitRegionElementCoordinates();
    
    //! clear region element coordinates
    void ClearRegionElementCoordinates();
    
    //! stores if region elements have been initalized
    bool initializedRegionElementCoordinates_;
    
    //! cached region element coordinates
    StdVector<StdVector<Vector<Double> > > regionElementCoordinates_;
    
    //! function to create on string of global factors from multiple strings
    std::string OneFactorStringFromMultipleFactors(StdVector<std::string>& factors);
    
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
    
    
    //! stores if the special case is present, that constant input data is given but time dependent global factors
    bool hasSpaceInputButTimeFactor_;
    
    //! stores if the constant input data is initialized
    bool initializedConstantInput_;
    
    //! cached constant input data
    StdVector< Vector< DATA_TYPE > > constantInput_;

    //=======================================================================================
    // Functions and structure for copying the element-based grid result
    //=======================================================================================
    
    //! struct needed for creation of typedef to template function for copying the result vector
    struct CopyResultFunction {
      //! type definition for pointer to CopyResult function
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
                          StdVector<UInt>& elemNums, 
                          StdVector<DATA_TYPE>& spaceFactorA, 
                          StdVector<DATA_TYPE>& spaceFactorB,
                          const DATA_TYPE constFactor,
                          StdVector<DATA_TYPE>& sum, 
                          StdVector<DATA_TYPE>& factorSum,
                          std::vector<bool>& countElems);
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor, 
        bool countSum, bool countAllValues, UInt dimDof>
    static typename CoefFunctionGridElem<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction() {
      return &CoefFunctionGridElem<DATA_TYPE>::CopyResult<
          useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,dimDof>;
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor, 
        bool countSum, bool countAllValues>
    static typename CoefFunctionGridElem<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(const UInt dimDof) {
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
      EXCEPTION("Cannot find CopyResult function for dimension " << dimDof << " in CoefFunctionGridElem::GetCopyResultFunction");
      return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,useConstFactor,countSum,countAllValues,1>();
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor>
    static typename CoefFunctionGridElem<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(
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
    static typename CoefFunctionGridElem<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(
        const bool useConstFactor, const bool countSum, const bool countAllValues, const UInt dimDof) {
      if (useConstFactor) {
        return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,true>(countSum, countAllValues, dimDof);
      }
      return GetCopyResultFunction<useSpaceFactorA,useSpaceFactorB,false>(countSum, countAllValues, dimDof);
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    template<bool useSpaceFactorA>
    static typename CoefFunctionGridElem<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(
        const bool useSpaceFactorB, const bool useConstFactor, 
        const bool countSum, const bool countAllValues, const UInt dimDof) {
      if (useSpaceFactorB) {
        return GetCopyResultFunction<useSpaceFactorA,true>(useConstFactor, countSum, countAllValues, dimDof);
      }
      return GetCopyResultFunction<useSpaceFactorA,false>(useConstFactor, countSum, countAllValues, dimDof);
    }
    
    //! Function for getting the right implementation of CopyResult function for the given problem
    //! recursively resolving all parameters into template parameters
    static typename CoefFunctionGridElem<DATA_TYPE>::CopyResultFunction::Ptr GetCopyResultFunction(
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
