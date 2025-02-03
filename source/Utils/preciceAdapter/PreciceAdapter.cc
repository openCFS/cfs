
#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


namespace CoupledField {

  // Define / declare logging stream
  DEFINE_LOG(preciceAdapter, "preciceAdapter")

  PreciceAdapter::PreciceAdapter(boost::shared_ptr<ParamNode> paramNode)
        : paramNode_(paramNode),
#ifdef USE_PRECICE
          participant_(nullptr),
#endif
          rank_(0),
          size_(1) 
      {}

  PreciceAdapter::~PreciceAdapter() {
      finalize();
  }


  void PreciceAdapter::initialize() {

#ifdef USE_PRECICE
      if (!paramNode_) {
          EXCEPTION("Parameter node is not set.");
      }

      // Retrieve participantName and configFileName from paramNode_
      std::string participantName;
      std::string configFileName;

      // Assuming Get returns a pointer to a ParamNode or similar
      paramNode_->Get("fileFormats/preciceCoupling/participantName")->GetValue("name", participantName);
      paramNode_->Get("fileFormats/preciceCoupling/precice_configFile")->GetValue("file", configFileName);

      LOG_DBG(preciceAdapter) << "Initializing PreCICE participant: " << participantName
                << " with config file: " << configFileName
                << " at rank: " << rank_ << " out of size: " << size_;

      // Create the PreCICE participant
      participant_ = std::make_unique<precice::Participant>(
          participantName,
          configFileName,
          rank_,
          size_
      );

      participant_->initialize();
    
#endif
  }


    void PreciceAdapter::finalize() {
#ifdef USE_PRECICE
      if (participant_) {
            participant_->finalize();
            participant_.reset();
            LOG_DBG(preciceAdapter) << "PreCICE participant finalized.";
          
      }
#endif
    }

} // end of namespace
