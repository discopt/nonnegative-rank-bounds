#include "separation_oracle.h"

#include <stdexcept>
#include <iostream>

namespace cpm
{
 
  SeparationOracle::SeparationOracle(std::size_t ambientDimension, int priority)
    : _ambientDimension(ambientDimension), _priority(priority)
  {

  }

  SeparationOracle::~SeparationOracle()
  {

  }

  void SeparationOracle::separate(bool separatePoint, const double* vector, std::vector<double>& lhs, std::vector<double>& rhs,
    std::vector<std::size_t>& begin, std::vector<std::size_t>& indices, std::vector<double>& values, double& violationLowerBound,
    double& violationUpperBound)
  {
    std::cerr << "SeparationOracle::separate() called for " << (separatePoint ? "point" : "ray") << ":";
    for (std::size_t v = 0; v < ambientDimension(); ++v)
    {
      std::cerr << " " << vector[v];
    }
    std::cerr << ".\n" << std::flush;

    violationUpperBound = 0.0;
  }

  void SeparationOracle::getFeasiblePointProperties(bool& scaledDown, bool& scaledUp, std::vector<std::size_t>& modifiableVariables) const
  {
    scaledDown = false;
    scaledUp = false;

    // Mark all variables as modifiable by default.

    modifiableVariables.clear();
    for (std::size_t v = 0; v < ambientDimension(); ++v)
      modifiableVariables.push_back(v);
  }

  std::size_t SeparationOracle::numFeasiblePoints() const
  {
    return 0;
  }

  void SeparationOracle::getFeasiblePoint(std::size_t id, std::vector< double >& point) const
  {
    if (id >= numFeasiblePoints())
      throw std::runtime_error("Invalid index while querying a feasible point.");
  }

} /* namespace cpm */
