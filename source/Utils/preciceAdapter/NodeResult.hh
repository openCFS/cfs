#ifndef NODE_RESULT_HH
#define NODE_RESULT_HH

#include "DataInOut/ResultHandler.hh"
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/Domain.hh"

#include "ResultBase.hh"
#include "MatVec/Vector.hh"
#include <map>
#include <vector>

namespace CoupledField {

class NodeResult : public ResultBase {
public:
    explicit NodeResult(const ResultConfig& config);
    virtual ~NodeResult() = default;
    
    // Read flatData_ into the node result map
    void readData(const std::vector<int>& nodeIndices, double deltaT) override;
    
    // Populate flatData_ from the node result map
    void writeData(const std::vector<int>& nodeIndices) override;
    
    ResultType getResultType() const override;

    // Update the internal map from OpenCFS result contexts
    // The contexts vector is a collection of (BaseResult, ResultContext) pairs
    void updateFromOpenCFS(const std::vector<std::pair<boost::shared_ptr<BaseResult>, boost::shared_ptr<ResultHandler::ResultContext>>>& contexts,
                          const std::vector<int>& nodeIndices);


    void setNodeResult(int nodeId, const Vector<double>& result);
    const std::map<int, Vector<double>>& getNodeResultMap() const;

private:
    std::map<int, Vector<double>> nodeResultMap_;
};

} // namespace CoupledField

#endif // NODE_RESULT_HH
