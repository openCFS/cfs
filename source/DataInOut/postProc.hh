#ifndef FILE_CFS_POSTPROC_HH
#define FILE_CFS_POSTPROC_HH

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Domain/resultInfo.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/vector.hh"

namespace CoupledField {


  // forward class declarations
  class Grid;
  class BaseResult;
  class ParamNode;

  //! Base class for defining interface for operations on result objects
  class PostProc {

  public :

    //! Type defining the spatial reduction of results
    typedef enum { SPACE, TIME_FREQ, NONE } ReductionType;
    
    //! Constructor
    PostProc( Grid* ptGrid, ParamNode * postProcNode );
    
    //! Destructor
    virtual ~PostProc();

    //! Return name of postProc
    std::string GetName() { return name_; } 

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

    //! Get next postProc for this result
    const std::string& GetNextPostProcName() { return next_;}

    //! Apply procedure
    virtual void Apply( ) = 0;

    //! Finalize 
    virtual void Finalize( ) {};

    //! Create postProcessing routine by reading related xml-specification
    static void CreatePostProc( ParamNode * procNode,
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
    ParamNode * myParam_;

    //! Input result
    shared_ptr<BaseResult> input_;

    //! Output result
    shared_ptr<BaseResult> output_;
    
    //! Flag for writing output
    bool writeResult_;

    //! Name of output destinations for output result
    StdVector<std::string> outputNames_;

    //! Reference to next postprocessing routine
    std::string next_;
  };


  //! Calculates the sum (sptatial/time) for one result type
  class PostProcSum : public PostProc {
    
  public:
    
    //! Constructor
    PostProcSum( Grid * ptGrid, ReductionType type, 
                 ParamNode * postProcNode );

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


  //! Calculates the maximum (spatial/time) of a result
  class PostProcMax : public PostProc {

  public:
  
    //! Constructor
    PostProcMax( Grid * ptGrid, ReductionType type,
                 ParamNode * postProcNode );
    
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

  //! Creates a new result by applying an arbitrary function to each result
  class PostProcFunc : public PostProc {

  public:

    //! Constructor
    PostProcFunc( Grid * ptGrid, ParamNode * postProcNode );

    //! Destructor
    virtual ~PostProcFunc();

    //! Initialize postprocessing operator
    void SetResult( shared_ptr<BaseResult> res );

    //! Set dofNames and functions for the new result
    void Initialize( const std::string& resultName,
                     const std::string& unit,
                     const StdVector<std::string>& dofNames,
                     const StdVector<std::string>& rFunctions,
                     const StdVector<std::string>& iFunctions ); 
    
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

    //! vector with related handles for math parser (real part)
    StdVector<MathParser::HandleType> rHandles_;

    //! vector with related handles for math parser (imag part)
    StdVector<MathParser::HandleType> iHandles_;

    //! vector containing values for registered variables (real part)
    Vector<Double> rVals_;

    //! vector containing values for registered variables (imag part)
    Vector<Double> iVals_;

    //! Pointer to MathParser object
    MathParser * mParser_;
  };

}

    
#endif
