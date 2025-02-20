#include "ElementResult.hh"


namespace CoupledField {

ElementResult::ElementResult(const ResultConfig& config)
    : ResultBase(config)
{}

void ElementResult::readData(const std::vector<int>& elemIndices, double deltaT) {
    elementResultMap_.clear();
    size_t numElems = elemIndices.size();
    for (size_t i = 0; i < numElems; ++i) {
        Vector<double> result;
        result.Resize(config_.quantitydim);
        for (int k = 0; k < config_.quantitydim; ++k) {
            result[k] = flatData_[i * config_.quantitydim + k];
        }
        elementResultMap_[elemIndices[i]] = result;
    }
}

void ElementResult::writeData(const std::vector<int>& elemIndices) {
    flatData_.resize(elemIndices.size() * config_.quantitydim);
    for (size_t i = 0; i < elemIndices.size(); ++i) {
        auto it = elementResultMap_.find(elemIndices[i]);
        if (it != elementResultMap_.end()) {
            const Vector<double>& result = it->second;
            for (int k = 0; k < config_.quantitydim; ++k) {
                flatData_[i * config_.quantitydim + k] = result[k];
            }
        }
    }
}

ResultType ElementResult::getResultType() const {
    return ResultType::ELEMENT;
}

void ElementResult::updateFromOpenCFS(
    const std::vector<std::pair<boost::shared_ptr<BaseResult>, boost::shared_ptr<ResultHandler::ResultContext>>>& contexts,
                        const std::vector<int>& elemIndices)
{
    for (const auto &pair : contexts) {
        shared_ptr<BaseResult> baseResult = pair.first;
        shared_ptr<ResultHandler::ResultContext> resultContext = pair.second;
        if (baseResult->GetResultInfo()->resultName == config_.cfsname) {
            // Downcast the functor.
            FieldCoefFunctor<double>* fcfunctor = dynamic_cast<FieldCoefFunctor<double>*>(resultContext->functor.get());
            if (!fcfunctor) {
                EXCEPTION("PreciceAdapter: no downcast to FieldCoefFunctor possible");
            }
            resultContext->functor->EvalResult(resultContext->result);
            EntityIterator it = resultContext->result->GetEntityList()->GetIterator();
            for (it.Begin(); !it.IsEnd(); it++) {
                const Elem* el = it.GetElem();
                LocPoint lp = Elem::shapes[el->type].midPointCoord;
                LocPointMapped lpm;
                Vector<double> tempField;
                tempField.Resize(config_.quantitydim);
                shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
                lpm.Set(lp, esm, 0.0);
                fcfunctor->GetVector(tempField, lpm);
                setElementResult(el->GetElemNum(), tempField);
            }
        }
    }
}

void ElementResult::setElementResult(int elemId, const Vector<double>& result) {
    elementResultMap_[elemId] = result;
}

const std::map<int, Vector<double>>& ElementResult::getElementResultMap() const {
    return elementResultMap_;
}

} // namespace CoupledField
