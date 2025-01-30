#include "NodeResult.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include <cmath>
#include <limits>
#include <sstream>

namespace CoupledField {
  DEFINE_LOG(preciceAdapterNodeRes, "preciceAdapterNodeRes")

  NodeResult::NodeResult(const ResultConfig& config)
      : ResultBase(config)
  {}

  void NodeResult::readData(const std::vector<int>& nodeIndices, double deltaT) {
    nodeResultMap_.clear();
    size_t numNodes = nodeIndices.size();
    for (size_t i = 0; i < numNodes; ++i) {
      Vector<double> result;
      result.Resize(config_.quantitydim);
      for (int k = 0; k < config_.quantitydim; ++k) {
        result[k] = flatData_[i * config_.quantitydim + k];
      }
      nodeResultMap_[nodeIndices[i]] = result;
    }
  }

  void NodeResult::writeData(const std::vector<int>& nodeIndices) {
    flatData_.resize(nodeIndices.size() * config_.quantitydim);
    for (size_t i = 0; i < nodeIndices.size(); ++i) {
      auto it = nodeResultMap_.find(nodeIndices[i]);
      if (it != nodeResultMap_.end()) {
        const Vector<double>& result = it->second;
        for (int k = 0; k < config_.quantitydim; ++k) {
          flatData_[i * config_.quantitydim + k] = result[k];
        }
      }
    }
  }

  ResultType NodeResult::getResultType() const {
    return ResultType::NODE;
  }

  void NodeResult::updateFromOpenCFS(
      const std::vector<std::pair<shared_ptr<BaseResult>, shared_ptr<ResultHandler::ResultContext>>>& contexts,
                            const std::vector<int>& nodeIndices)
  {
    nodeResultMap_.clear();

    // Loop over the provided contexts
    for (const auto &pair : contexts) {
      shared_ptr<BaseResult> baseResult = pair.first;
      shared_ptr<ResultHandler::ResultContext> resultContext = pair.second;
      if (baseResult->GetResultInfo()->resultName == config_.cfsname) {
        // Ensure we read the current time-step values even if this result was
        // not scheduled for regular output in ResultHandler::isNeeded_.
        if (resultContext->functor) {
          resultContext->functor->EvalResult(resultContext->result);
        }

        // Assume the result context contains a node-based vector.
        Vector<double>* res_vec = dynamic_cast<Vector<Double>*>(resultContext->result->GetSingleVector());
        if (!res_vec) {
          EXCEPTION("Unable to retrieve node result vector for opencfs result " << config_.cfsname);
        }

        const int vecSize = res_vec->GetSize();
        if (config_.quantitydim <= 0) {
          EXCEPTION("Invalid quantity dimension " << config_.quantitydim
                    << " for opencfs result " << config_.cfsname);
        }
        if ((vecSize % config_.quantitydim) != 0) {
          EXCEPTION("Result vector size " << vecSize
                    << " is not divisible by quantity dimension " << config_.quantitydim
                    << " for opencfs result " << config_.cfsname);
        }

        int entityCount = 0;
        {
          EntityIterator countIt = resultContext->result->GetEntityList()->GetIterator();
          for (countIt.Begin(); !countIt.IsEnd(); countIt++) {
            ++entityCount;
          }
        }
        const int expectedSize = entityCount * config_.quantitydim;
        if (vecSize < expectedSize) {
          EXCEPTION("Result vector for opencfs result " << config_.cfsname
                    << " is too small: size=" << vecSize
                    << ", expected at least " << expectedSize
                    << " (entities=" << entityCount
                    << ", dim=" << config_.quantitydim << ")");
        }

        int pos = 0;
        EntityIterator it = resultContext->result->GetEntityList()->GetIterator();
        for (it.Begin(); !it.IsEnd(); it++) {
          int node = it.GetNode();
          Vector<double> nodeVal;
          nodeVal.Resize(config_.quantitydim);
          const int base = pos * config_.quantitydim;
          for (int k = 0; k < config_.quantitydim; ++k) {
            const int idx = base + k;
            if (idx >= vecSize) {
              EXCEPTION("Out-of-bounds access while reading opencfs result " << config_.cfsname
                        << ": idx=" << idx << ", size=" << vecSize
                        << ", pos=" << pos << ", dim=" << config_.quantitydim << ")");
            }
            nodeVal[k] = (*res_vec)[idx];
          }
          setNodeResult(node, nodeVal);
          ++pos;
        }

        // Ensure that all coupling nodes requested by preCICE are available.
        double minVal = std::numeric_limits<double>::infinity();
        double maxVal = -std::numeric_limits<double>::infinity();
        double maxAbs = 0.0;
        int nonFiniteCount = 0;
        int hugeCount = 0;
        int firstHugeNode = -1;
        std::vector<double> firstHugeVals(config_.quantitydim, 0.0);
        const double hugeThreshold = 1.0e20;

        for (int nodeId : nodeIndices) {
          auto mapIt = nodeResultMap_.find(nodeId);
          if (mapIt == nodeResultMap_.end()) {
            EXCEPTION("Missing node result for coupling node " << nodeId
                      << " in opencfs result " << config_.cfsname);
          }

          const Vector<double>& nodeVal = mapIt->second;
          for (int k = 0; k < config_.quantitydim; ++k) {
            const double v = nodeVal[k];
            if (!std::isfinite(v)) {
              ++nonFiniteCount;
              continue;
            }

            if (v < minVal) {
              minVal = v;
            }
            if (v > maxVal) {
              maxVal = v;
            }

            const double absV = std::fabs(v);
            if (absV > maxAbs) {
              maxAbs = absV;
            }
            if (absV > hugeThreshold) {
              ++hugeCount;
              if (firstHugeNode < 0) {
                firstHugeNode = nodeId;
                for (int kk = 0; kk < config_.quantitydim; ++kk) {
                  firstHugeVals[kk] = nodeVal[kk];
                }
              }
            }
          }
        }

        if (!std::isfinite(minVal)) {
          minVal = 0.0;
        }
        if (!std::isfinite(maxVal)) {
          maxVal = 0.0;
        }

        LOG_DBG(preciceAdapterNodeRes) << "pre-write stats: mesh='" << config_.meshName
                  << "', data='" << config_.precicename
                  << "', cfs='" << config_.cfsname
                  << "', nodes=" << nodeIndices.size()
                  << ", vecSize=" << vecSize
                  << ", min=" << minVal
                  << ", max=" << maxVal
                  << ", maxAbs=" << maxAbs
                  << ", nonFinite=" << nonFiniteCount
                  << ", huge(>|" << hugeThreshold << "|)=" << hugeCount;

        if (nonFiniteCount > 0 || hugeCount > 0) {
          std::ostringstream vals;
          for (int k = 0; firstHugeNode >= 0 && k < config_.quantitydim; ++k) {
            vals << " " << firstHugeVals[k];
          }
          WARN("PreciceAdapter(NodeResult) writes " << nonFiniteCount << " non-finite and "
               << hugeCount << " huge (>|" << hugeThreshold << "|) values for '" << config_.cfsname
               << "'; first huge node " << firstHugeNode << " values:" << vals.str());
        }
      }
    }
  }

  void NodeResult::setNodeResult(int nodeId, const Vector<double>& result) {
    nodeResultMap_[nodeId] = result;
  }

  const std::map<int, Vector<double>>& NodeResult::getNodeResultMap() const {
    return nodeResultMap_;
  }

} // namespace CoupledField
