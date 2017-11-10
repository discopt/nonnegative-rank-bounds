#ifndef _SEPARATION_ORACLE_H_
#define _SEPARATION_ORACLE_H_

#include <numeric>
#include <vector>

// TODO: Implement variable names array
// TODO: Implement variable index mapping from core to oracle (having a subset).

namespace cpm
{

  class SeparationOracle
  {
  protected:
    std::size_t _ambientDimension;
    const int _priority;

  public:
    SeparationOracle(std::size_t ambientDimension, int priority);

    virtual ~SeparationOracle();

    inline std::size_t ambientDimension() const
    {
      return _ambientDimension;
    }

    inline int priority() const
    {
      return _priority;
    }

    virtual void separate(bool separatePoint, const double* vector, std::vector<double>& lhs, std::vector<double>& rhs,
      std::vector<std::size_t>& begin, std::vector<std::size_t>& indices, std::vector<double>& values, double& violationLowerBound,
      double& violationUpperBound);

    virtual void getFeasiblePointProperties(bool& scaledDown, bool& scaledUp, std::vector<std::size_t>& modifiableVariables) const;

    virtual std::size_t numFeasiblePoints() const;

    virtual void getFeasiblePoint(std::size_t id, std::vector<double>& point) const;
  };

} /* namespace cpm */

#endif /* _SEPARATION_ORACLE_H_ */