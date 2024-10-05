#ifndef FILE_CFS_POSTPROC_HH
#define FILE_CFS_POSTPROC_HH

#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "Domain/Results/ResultInfo.hh"
#include "MatVec/Vector.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {


  // forward class declarations
  class Grid;
  class BaseResult;
  class MathParser;

  //! Base class for defining interface for operations on result objects
  class PostProc {

  public :

    //! Type defining the spatial reduction of results
    typedef enum { SPACE, TIME_FREQ, NONE } ReductionType;
    
    //! Constructor
    PostProc( Grid* ptGrid, PtrParamNode postProcNode );
    
    //! Destructor
    virtual ~PostProc();

    //! Return name of postProc
    const std::string& GetName() const { return name_; }

    //! Add new result
    virtual void SetResult( shared_ptr<BaseResult> res ) = 0;

    //! Get output result
    shared_ptr<BaseResult> GetOutputResult() { return output_; }
    
    //! Get input result
    shared_ptr<BaseResult> GetInputResult() { return input_; }

    //! Get reduction type
    ReductionType GetReductionType() { return reducType_; }

    //! Get output destination names
    void GetOutDestNames( StdVector<std::string> & outNames );

    //! Get flag, if result should be written
    bool IsWriteResult();

    //! Get flag if result of postprocessing is history result
    bool IsHistory();
    
    //! Get next postProc for this result
    const std::string& GetNextPostProcName() { return next_;}

    //! Apply procedure
    virtual void Apply( ) = 0;

    //! Finalize 
    virtual void Finalize( ) {};

    //! Create postProcessing routine by reading related xml-specification
    static void CreatePostProc( PtrParamNode procNode,
                                Grid * ptGrid, 
                                StdVector<shared_ptr<PostProc> >& postProcs );
    
  protected:

    //! Get for a given entitylist and reduction type a new entitylist
    void GetReducedList( shared_ptr<EntityList>& outList,
                         ResultInfo::EntityUnknownType& outType,
                         const shared_ptr<EntityList> inList,
                         const ResultInfo::EntityUnknownType inType,
                         const ReductionType reducType );
    
    //! Name of postProc
    std::string name_;

    //! Reduction type
    ReductionType reducType_;

    //! Pointer to grid
    Grid * ptGrid_;

    //! Pointer to parameter node of current postProc
    PtrParamNode myParam_;

    //! Input result
    shared_ptr<BaseResult> input_;

    //! Output result
    shared_ptr<BaseResult> output_;

    //! Flag for writing output
    bool writeResult_;
    
    //! Flag indicatin if result is history
    bool isHistory_;

    //! Name of output destinations for output result
    StdVector<std::string> outputNames_;

    //! Reference to next postprocessing routine
    std::string next_;
  };

  // -------------------------------------------------------------------------

  //! Calculates the sum (sptatial/time) for one result type
  class PostProcSum : public PostProc {
    
  public:
    
    //! Constructor
    PostProcSum( Grid * ptGrid, ReductionType type, 
                 PtrParamNode postProcNode );

    //! Destructor
    virtual ~PostProcSum();
    
    
    void SetResult( shared_ptr<BaseResult> res );

    //! Apply procedure (calculate sum)
    void Apply( );
    

    //! Finalize (return sum of accumulated values)
    void Finalize( );
    

  protected:

    //! Templatized sub-routine
    template<class TYPE>
    void CalcSum();
    
  };

  // -------------------------------------------------------------------------

  //! Calculates the maximum (spatial/time) of a result
  class PostProcMax : public PostProc {

  public:
  
    //! Constructor
    PostProcMax( Grid * ptGrid, ReductionType type,
                 PtrParamNode postProcNode );
    
    //! Destructor
    virtual ~PostProcMax();
    
    //! Initialize postprocessing operator
    void SetResult( shared_ptr<BaseResult> res );
    
    //! Apply procedure (calculate maximum)
    void Apply( );
    
    //! Finalize 
    void Finalize( );

  protected:
    
    //! Templatized sub-routine
    template<class TYPE>
    void CalcMax();

    //! Calculates the operation '<' for each singleDof

    //! Calculates the operation '<' for a single number.
    //! For complex numbers, tha absolute value is compared
    template<class TYPE> 
    static bool opLtSingleDof(TYPE a, TYPE b);
  };


  // -------------------------------------------------------------------------

  //! Creates a new result by applying an arbitrary function to each result
  class PostProcFunc : public PostProc {

  public:

    //! Constructor
    PostProcFunc( Grid * ptGrid, PtrParamNode postProcNode );

    //! Destructor
    virtual ~PostProcFunc();

    //! Initialize postprocessing operator
    void SetResult( shared_ptr<BaseResult> res );

    //! Set dofNames and functions for the new result
    void Initialize( const std::string& resultName,
                     const std::string& unit,
                     const StdVector<std::string>& dofNames,
                     const StdVector<std::string>& rFunctions,
                     const StdVector<std::string>& iFunctions,
                     const StdVector<std::string>& variableNames,
                     const PtrParamNode& postProcNode ); 
    
    //! Apply procedure
    void Apply( );

    //! Finalize
    void Finalize( );
    
  private:

    //! Apply function (real valued)
    void ApplyReal();
    
    //! Apply function (complex values);
    void ApplyComplex();
    
    //! name of new result to be created
    std::string resultName_;

    //! unit of new result to be created
    std::string unit_;
    
    //! dof names of new result
    StdVector<std::string> dofNames_;

    //! function string of new result (real part)
    StdVector<std::string> rFuncs_;

    //! function string of new result (imag part)
    StdVector<std::string> iFuncs_;

    //! variable name for the result handling (will be used to evaluate the mathParser)
    StdVector<std::string> variableNames_;

    //! function string of new result (real part), matched to grid DOFs
    StdVector<std::string> rFuncsSorted_;

    //! function string of new result (real part), matched to grid DOFs, adapted to mathParser naming convetion of variables
    StdVector<std::string> rFuncsSortedMp_;

    //! vector with related handles for math parser (real part)
    StdVector<unsigned int> rHandles_;

    //! vector with related handles for math parser (imag part)
    StdVector<unsigned int> iHandles_;

    //! vector containing values for registered variables (real part)
    Vector<Double> rVals_;

    //! vector containing values for registered variables (imag part)
    Vector<Double> iVals_;

    //! Pointer to MathParser object
    MathParser * mParser_;

    //! ParamNode
    PtrParamNode postProcNode_;
  };

// -------------------------------------------------------------------------

  //! Creates a new result by applying an arbitrary function to each result and calculating the L2 norm
  class PostProcL2Norm : public PostProc {

  public:

    //! Constructor
    PostProcL2Norm( Grid * ptGrid, PtrParamNode postProcNode );

    //! Destructor
    virtual ~PostProcL2Norm();

    //! Initialize postprocessing operator
    void SetResult( shared_ptr<BaseResult> res );

    //! Set dofNames and functions for the new result
    void Initialize( const std::string& resultName,
                     const UInt& integrationOrder,
                     const std::string& mode,
                     const std::string& computeSolutionSteps,
                     const StdVector<std::string>& dofNames,
                     const StdVector<std::string>& rFunctions,
                     const StdVector<std::string>& iFunctions ); 
    
    //! Apply procedure
    void Apply( );

    //! Finalize
    void Finalize( );

  private:

    //! Templatized sub-routine
    template<class TYPE>
    void CalcNorm();
    
    //! name of new result to be created
    std::string resultName_;

    //! integration order for evaluation, needs to be set high enough to accurately get the differenece between analytical function and FE approximation
    UInt integrationOrder_;

    //! mode how the L2 norm function shall operate (absolute/relative)
    std::string mode_;

    //! defines which time/frequency steps should be considered for computing the postProc result
    //! currently, this is only used in the L2Norm and can be "all" or "last"
    std::string computeSolutionSteps_;
    
    //! dof names of new result
    StdVector<std::string> dofNames_;

    //! function string of new result (real part)
    StdVector<std::string> rFuncs_;

    //! function string of new result (imag part)
    StdVector<std::string> iFuncs_;

    //! function string of new result (real part), matched to grid DOFs
    StdVector<std::string> rFuncsSorted_;

    //! function string of new result (imag part), matched to grid DOFs
    StdVector<std::string> iFuncsSorted_;

    //! vector with related handles for math parser (real part)
    StdVector<unsigned int> rHandles_;

    //! vector with related handles for math parser (imag part)
    StdVector<unsigned int> iHandles_;

    //! vector containing values for registered variables (real part)
    Vector<Double> rVals_;

    //! vector containing values for registered variables (imag part)
    Vector<Double> iVals_;

    //! Pointer to MathParser object
    MathParser * mParser_;

    //! coefFunction of the difference (FE solution minus reference solution) for norm evaluation
    PtrCoefFct diffCoef_;

    //! coefFunnction of the reference solution
    PtrCoefFct coef_;

  };

  // -------------------------------------------------------------------------

  //! Cut off values beneath or above a given limit
  class PostProcLimit : public PostProc {

  public:
  
    //! Constructor
    PostProcLimit( Grid * ptGrid, 
                   PtrParamNode postProcNode );
    
    //! Destructor
    virtual ~PostProcLimit();
    
    //! Initialize postprocessing operator
    void SetResult( shared_ptr<BaseResult> res );
    
    //! Set lower limit
    void SetLowerLimit( Double loLim );
    
    //! Set upper limit
    void SetUpperLimit( Double upLim );
    
    //! Apply procedure (calculate maximum)
    void Apply( );

  protected:
    
    //! Templatized sub-routine
    template<class TYPE>
    void CalcLimit();

    //! Calculates the operation '<' for a single number.
    //! For complex numbers, the absolute value is compared
    template<class TYPE1, class TYPE2> 
    static bool opLtSingleDof(TYPE1 a, TYPE2 b);

    //! Calculates the operation '>' for a single number.
    //! For complex numbers, the absolute value is compared
    template<class TYPE1, class TYPE2> 
    static bool opGtSingleDof(TYPE1 a, TYPE2 b);
    
    //! flag indicating if lower limit is set
    bool loLimSet_;

    //! flag indicating if upper limit is set
    bool upLimSet_;

    //! lower limit value
    Double loLim_;

    //! upper limit value
    Double upLim_;

  };

}

    
#endif
