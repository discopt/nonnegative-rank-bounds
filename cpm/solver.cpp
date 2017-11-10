#include "solver.h"

namespace cpm {

  Solver::Solver()
  {

  }

  Solver::~Solver()
  {

  }

  std::size_t Solver::addVariable(const std::string& name)
  {
    _variableNames.push_back(name);
    _point.push_back(0.0);
    _ray.push_back(0.0);
    return _variableNames.size() - 1;
  }

}