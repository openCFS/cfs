#ifndef SIMOUTPUTINFO_HH_
#define SIMOUTPUTINFO_HH_

#include "DataInOut/simOutput.hh"
#include "Domain/entityList.hh"

namespace CoupledField 
{
  class InfoNode;

  /** This class collects the information, the basic text writer writes to files 
   * in to info.xml output.
   * This writer is intialized by default */
  class SimOutputInfo : public SimOutput 
  {
  public:

    SimOutputInfo( ParamNode * outputNode );

    ~SimOutputInfo();

    //! Initialize class
    void Init( Grid * ptGrid, bool printGridOnly );    

    //! Register result (within one multisequence step)
    void RegisterResult(shared_ptr<BaseResult> sol, UInt saveBegin, UInt saveInc,
                                 UInt saveEnd, bool isHistory );    
    
    //! Begin multisequence step
    void BeginMultiSequenceStep(UInt step, BasePDE::AnalysisType type, UInt numSteps);
    
    
    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );
  private:
    /** to be set via "info" file entry in xml */
    bool output;
    
    /** root in info.xml */
    InfoNode* info_root;
  };
}

#endif /* SIMOUTPUTINFO_HH_ */
