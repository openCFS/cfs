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
    CoefFunctionGridNodal(PtrParamNode configNode);

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

  protected:

    //=======================================================================================
    // Interpolation operator for nodal results
    //=======================================================================================
    //!Create the interpolation operator
    //!and sets the corresponding class variable
    void CreateOperator(UInt spaceDim, UInt dofDim);

    //! BOperator to map solutions to arbitrary points
    //! Right now, hardcoded identity operator
    shared_ptr<BaseBOperator > myOperator_;

    //=======================================================================================
    // Functions and structure for nodal grid result handling 
    //=======================================================================================

    //! Read solution from sourceFile according to the given stepnumber
    void ReadSolution(UInt step,Vector<DATA_TYPE> & sol);
    
    //! Updates the solution vector
    void UpdateSolution();

    //! Perform a simple equation mapping for nodal grids
    //! to make solution access simpler
    void MapEqns();

    //! Extract the solution for a source element in order to apply the interpolation operator
    void GetElemSolution(Vector<DATA_TYPE> & sol, UInt eNum);

    //! Returns the step to be read in along with temporal interpolation factors
    //! If those are necessary
    UInt GetStepNum(bool &interpolateT,Double & iFactor1, Double & iFactor2);

     //! Associate node numbers to equations
    std::map<UInt,UInt> nodeIdxMap_;

    //! the equation Numbers, for each node, dimDOF equation numbers
    StdVector< StdVector<UInt> > eqnNumbers_;

    //! Stores the current solution vector
    Vector<DATA_TYPE> solVec_;

    //! stores the next solution vector for temporal interpolation
    Vector<DATA_TYPE> solVecFuture_;

    //! stores the solution vector at previuous timesteps
    Vector<DATA_TYPE> solVecOld_;

    //! Total number of equations
    UInt numEqns_;

    //! Equation mapping completed
    bool eqnMapComplete_;

    //! Pointer to math parser instance
    MathParser* mp_;

    //! Handle for expression determines current time/freq value 
    MathParser::HandleType mHandleStep_;

    //! shared pointer to coefFunction representing multiplicative factor
    shared_ptr<CoefFunction> factorFnc_;

    //! stroes the stepnumber of the las read solution process
    UInt lastStepRead_;

  private:

};
}

#endif
