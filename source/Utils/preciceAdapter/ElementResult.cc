#include "ElementResult.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include <cmath>
#include <limits>
#include <sstream>

namespace CoupledField {
  DEFINE_LOG(preciceAdapterElemRes, "preciceAdapterElemRes")

  ElementResult::ElementResult(const ResultConfig& config)
      : ResultBase(config)
  {}

  void ElementResult::readData(const std::vector<int>& elemIndices, double deltaT) {
    elementResultMap_.clear();
    size_t numElems = elemIndices.size();
    double minVal = std::numeric_limits<double>::infinity();
    double maxVal = -std::numeric_limits<double>::infinity();
    double maxAbs = 0.0;
    int nonFiniteCount = 0;
    int hugeCount = 0;
    const double hugeThreshold = 1.0e20;
    int firstHugeElem = -1;
    std::vector<double> firstHugeVals(config_.quantitydim, 0.0);

    for (size_t i = 0; i < numElems; ++i) {
      Vector<double> result;
      result.Resize(config_.quantitydim);
      for (int k = 0; k < config_.quantitydim; ++k) {
        const double v = flatData_[i * config_.quantitydim + k];
        result[k] = v;

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
          if (firstHugeElem < 0) {
            firstHugeElem = elemIndices[i];
            for (int kk = 0; kk < config_.quantitydim; ++kk) {
              firstHugeVals[kk] = result[kk];
            }
          }
        }
      }
      elementResultMap_[elemIndices[i]] = result;
      LOG_DBG(preciceAdapterElemRes) << "Read result "<< config_.cfsname<< " for element "<< elemIndices[i]<< ": " << result.ToString();
    }

    if (!std::isfinite(minVal)) {
      minVal = 0.0;
    }
    if (!std::isfinite(maxVal)) {
      maxVal = 0.0;
    }

    LOG_DBG(preciceAdapterElemRes) << "read stats: mesh='" << config_.meshName
              << "', data='" << config_.precicename
              << "', cfs='" << config_.cfsname
              << "', elems=" << elemIndices.size()
              << ", min=" << minVal
              << ", max=" << maxVal
              << ", maxAbs=" << maxAbs
              << ", nonFinite=" << nonFiniteCount
              << ", huge(>|" << hugeThreshold << "|)=" << hugeCount;

    if (nonFiniteCount > 0 || hugeCount > 0) {
      std::ostringstream vals;
      for (int k = 0; firstHugeElem >= 0 && k < config_.quantitydim; ++k) {
        vals << " " << firstHugeVals[k];
      }
      WARN("PreciceAdapter(ElementResult) read " << nonFiniteCount << " non-finite and "
           << hugeCount << " huge (>|" << hugeThreshold << "|) values for '" << config_.cfsname
           << "'; first huge element " << firstHugeElem << " values:" << vals.str());
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
      const std::vector<std::pair<shared_ptr<BaseResult>, shared_ptr<ResultHandler::ResultContext>>>& contexts,
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
