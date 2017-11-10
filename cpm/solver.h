#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <vector>
#include <string>

#include <soplex.h>

#include "separation_oracle.h"

namespace cpm
{

  class Solver
  {
  public:
    typedef int Status;
    const static Status ERROR = -4;
    const static Status TIME_LIMIT = -3;
    const static Status ITERATION_LIMIT = -2;
    const static Status NUMERICAL_TROUBLE = -1;
    const static Status OPTIMAL =  1;
    const static Status UNBOUNDED =  2;
    const static Status INFEASIBLE =  3;

    Solver();

    virtual ~Solver();

    virtual std::size_t addVariable(const std::string& name, double objective, double lowerBound, double upperBound) = 0;

    virtual void addInequalities(std::vector<double>& lhs, std::vector<double>& rhs, std::vector<std::size_t>& begin,
      std::vector<std::size_t>& indices, std::vector<double>& values) = 0;

    virtual Status run() = 0;

    inline const std::string& variableName(std::size_t variable) const
    {
      return _variableNames[variable];
    }

    virtual double objectiveCoefficient(std::size_t variable) const = 0;

    inline std::size_t numVariables() const
    {
      return _variableNames.size();
    }

    inline const std::vector<double>& point() const
    {
      return _point;
    }

    inline const std::vector<double>& ray() const
    {
      return _ray;
    }

  protected:
    std::size_t addVariable(const std::string& name);

    std::vector<std::string> _variableNames;
    std::vector<double> _point;
    std::vector<double> _ray;
  };

} /* namespace cpm */

#endif /* _SOLVER_H_ */