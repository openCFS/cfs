#ifndef RESULT_BASE_HH
#define RESULT_BASE_HH

#include "General/Environment.hh"
#include "boost/shared_ptr.hpp"
#include <vector>
#include <map>
#include <vector>


namespace CoupledField {

// Define an enum to distinguish result types
enum class ResultType {
    NODE,
    ELEMENT
};

struct ResultConfig {
    std::string precicename;   // e.g., "Temperature"
    std::string cfsname;       // e.g., "heatTemperature"
    SolutionType solutiontype; // e.g., SolutionType::HEAT_TEMPERATURE
    std::string meshName;      // Mesh name where the result is defined
    int quantitydim;           // Dimension of the quantity (scalar, vector, etc.)
};


class ResultBase {
public:
    ResultBase(const ResultConfig& config) : config_(config) {}
    virtual ~ResultBase() = default;

    // Reads flat data into internal storage (e.g., mapping into internal maps)
    virtual void readData(const std::vector<int>& indexMapping, double deltaT) = 0;
    
    // Writes the internal result data into flat data (which can be passed to preCICE)
    virtual void writeData(const std::vector<int>& indexMapping) = 0;

    // Allocate flatData_ vector based on the number of indices and quantity dimension.
    void allocateData(size_t numIndices) {
        flatData_.resize(numIndices * config_.quantitydim);
    }

    // Returns the flat data vector (read-only)
    const std::vector<double>& getFlatData() const { return flatData_; }

    // Returns whether this is a node- or element-based result.
    virtual ResultType getResultType() const = 0;

    // Provides access to the configuration
    const ResultConfig& getConfig() const { return config_; }

protected:
    ResultConfig config_;
    std::vector<double> flatData_; // This is the flatData vector that was in PreciceRuntimeQuantity
};

} // namespace CoupledField

#endif