#include "NodeResult.hh"

namespace CoupledField {

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
    const std::vector<std::pair<boost::shared_ptr<BaseResult>, boost::shared_ptr<ResultHandler::ResultContext>>>& contexts,
                          const std::vector<int>& nodeIndices)
{
    // Loop over the provided contexts
    for (const auto &pair : contexts) {
        shared_ptr<BaseResult> baseResult = pair.first;
        shared_ptr<ResultHandler::ResultContext> resultContext = pair.second;
        if (baseResult->GetResultInfo()->resultName == config_.cfsname) {
            // Assume the result context contains a node-based vector.
            Vector<double>* res_vec = dynamic_cast<Vector<Double>*>(resultContext->result->GetSingleVector());
            if (!res_vec) {
                EXCEPTION("Unable to retrieve node result vector for opencfs result " << config_.cfsname);
            }
            int pos = 0;
            EntityIterator it = resultContext->result->GetEntityList()->GetIterator();
            for (it.Begin(); !it.IsEnd(); it++) {
                int node = it.GetNode();
                Vector<double> nodeVal;
                nodeVal.Resize(config_.quantitydim);
                for (int k = 0; k < config_.quantitydim; ++k) {
                    nodeVal[k] = (*res_vec)[pos * config_.quantitydim + k];
                }
                setNodeResult(node, nodeVal);
                ++pos;
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
