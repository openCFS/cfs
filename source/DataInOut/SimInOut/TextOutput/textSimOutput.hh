#ifndef FILE_CFS_TEXT_SIM_OUTPUT_HH
#define FILE_CFS_TEXT_SIM_OUTPUT_HH

#include "Domain/grid.hh"
#include "DataInOut/simOutput.hh"

namespace CoupledField {


  class SimOutputText : public SimOutput {

  public:

    //! Constructor
    SimOutputText( const std::string& fileName, ParamNode * outputNode );

    //! Destructor
    ~SimOutputText();

    //! Initialize class
    void Init( Grid * ptGrid );

    //! Write grid (in this case nothing)
    void WriteGrid() {};
    
    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd );

    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    void FinishStep( );

  private:

    //! Create new output stream
    void CreateFiles( shared_ptr<BaseResult> res,
                      UInt step,
                      Double stepVal );
    
    //! Map with result objects for each result type
    ResultMapType resultMap_;
    
    //! map of result and vector of ofstreams
    std::map<shared_ptr<BaseResult>,
             StdVector<std::ofstream*> > outFiles_;
    
    
    //! comment char
    char cmChar_;
  };
}



#endif
