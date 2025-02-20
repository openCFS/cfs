#ifndef ELEMENT_RESULT_HH
#define ELEMENT_RESULT_HH

#include "DataInOut/ResultHandler.hh"
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/Domain.hh"

#include "ResultBase.hh"
#include "MatVec/Vector.hh"
#include <map>
#include <vector>
#include <memory>

namespace CoupledField {

class ElementResult : public ResultBase {
public:
    explicit ElementResult(const ResultConfig& config);
    virtual ~ElementResult() = default;
    
    // Read flatData_ into the element result map
    void readData(const std::vector<int>& elemIndices, double deltaT) override;
    
    // Populate flatData_ from the element result map
    void writeData(const std::vector<int>& elemIndices) override;
    
    ResultType getResultType() const override;

    // Update the internal element map from OpenCFS result contexts
    void updateFromOpenCFS(const std::vector<std::pair<boost::shared_ptr<BaseResult>, boost::shared_ptr<ResultHandler::ResultContext>>>& contexts,
                        const std::vector<int>& elemIndices);


    void setElementResult(int elemId, const Vector<double>& result);
    const std::map<int, Vector<double>>& getElementResultMap() const;

private:
    std::map<int, Vector<double>> elementResultMap_;
};

} // namespace CoupledField

#endif // ELEMENT_RESULT_HH
