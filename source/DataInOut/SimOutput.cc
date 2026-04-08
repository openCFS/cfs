#include "SimOutput.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/Mesh/Grid.hh"
#include "Utils/Timer.hh"

namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simOutput, "simOutput")
    
  
    SimOutput::SimOutput( const std::string& fileName, 
                          PtrParamNode outputNode,
                          PtrParamNode infoNode,
                          bool isRestart ) {
    fileName_ = fileName;
    isRestart_ = isRestart;
    myParam_ = outputNode;
    myInfo_ = infoNode;
    
    timer = make_shared<Timer>();
    myInfo_->Get("timer")->SetValue(timer);
  }

  SimOutput::~SimOutput() {
  }

  template<class TYPE>
  void SimOutput::FillGlobalVec(Vector<TYPE>& gSol, 
                                const StdVector<shared_ptr<BaseResult> > & solList,
                                ResultInfo::EntityUnknownType entityType ) {
    Grid* ptGrid = solList[0]->GetEntityList()->GetGrid();
    ResultInfo & actDof = *(solList[0]->GetResultInfo());
    UInt numDofs = actDof.dofNames.GetSize();
    LOG_DBG(simOutput) << "Filling global vector for result '" 
                       << actDof.resultName << "' on ";
    for( UInt i = 0; i < solList.GetSize(); i++ ) {
      LOG_DBG(simOutput) << solList[i]->GetEntityList()->GetName();
    }

    if( entityType == ResultInfo::ELEMENT ||
        entityType == ResultInfo::SURF_ELEM ) {

      // === Element Results ===
      gSol.Resize( ptGrid->GetNumElems()*numDofs );
      gSol.Init();
      for( UInt i = 0; i < solList.GetSize(); i++ ){
        EntityIterator it = solList[i]->GetEntityList()->GetIterator();
        Vector<TYPE> & actSol = dynamic_cast<Result<TYPE>&>
          (*(solList[i])).GetVector();
        for( it.Begin(); !it.IsEnd(); it++ ) {
          UInt elemNum = it.GetElem()->elemNum;
          for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
            gSol[(elemNum-1)*numDofs+iDof] = actSol[it.GetPos()*numDofs+iDof];
          }
        }
      }
    } else if ( entityType == ResultInfo::NODE ) {

      // === Nodal Results ===
      gSol.Resize( ptGrid->GetNumNodes()*numDofs );
      gSol.Init();
      for( UInt i = 0; i < solList.GetSize(); i++ ){
        EntityIterator it = solList[i]->GetEntityList()->GetIterator();
        Vector<TYPE> & actSol = dynamic_cast<Result<TYPE>&>
          (*(solList[i])).GetVector();
        assert( (UInt) (actSol.GetSize()/numDofs) 
                == solList[i]->GetEntityList()->GetSize());
        for( it.Begin(); !it.IsEnd(); it++ ) {
          UInt nodeNum = it.GetNode();
          for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
            gSol[(nodeNum-1)*numDofs+iDof] = actSol[it.GetPos()*numDofs+iDof];
          }
        }
      }
    } else {
      EXCEPTION( "Can only map nodal / element results to grid" );
    }
      
  }

  bool SimOutput::ValidateNodesAndElements(ResultInfo& actInfo)
  {
    if(actInfo.definedOn != ResultInfo::NODE &&
       actInfo.definedOn != ResultInfo::ELEMENT &&
       actInfo.definedOn != ResultInfo::SURF_ELEM ) 
    {
      std::string msg = formatName_ + " can only write results on element and nodes";
      WARN(msg.c_str());
      return false;
    }
    else 
    {
      return true;
    }
  }

  // explicit template instantiation
  template void SimOutput::
  FillGlobalVec<Double>(Vector<Double>& gSol, 
                        const StdVector<shared_ptr<BaseResult> > & solList,
                        ResultInfo::EntityUnknownType entityType );
  template void SimOutput::
  FillGlobalVec<Complex>(Vector<Complex>& gSol, 
                         const StdVector<shared_ptr<BaseResult> > & solList,
                         ResultInfo::EntityUnknownType entityType );
}
