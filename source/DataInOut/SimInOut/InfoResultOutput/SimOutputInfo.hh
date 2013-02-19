#ifndef SIMOUTPUTINFO_HH_
#define SIMOUTPUTINFO_HH_

#include "DataInOut/SimOutput.hh"
#include "Domain/ElemMapping/EntityLists.hh"

namespace CoupledField 
{
  class InfoNode;

  /** This class collects the information, the basic text writer writes to files 
   * in to info.xml output.
   * This writer is intialized by default */
  class SimOutputInfo : public SimOutput 
  {
  public:

    SimOutputInfo( PtrParamNode outputNode );

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
    
    /** root in info.xml */
    PtrParamNode info_root;
  };
}

#endif /* SIMOUTPUTINFO_HH_ */
