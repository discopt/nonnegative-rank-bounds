#ifndef _CORE_H_
#define _CORE_H_

#include <memory>
#include <chrono>

#include "solver.h"
#include "separation_oracle.h"

namespace cpm
{
  struct SolutionData
  {
    std::vector<double> values;
    double objectiveValue;

    SolutionData(const std::vector<double>& values, double objectiveValue);
  };
  typedef std::shared_ptr<SolutionData> Solution;
  
  class Core
  {
  public:
    Core(Solver* solver);

    ~Core();

    std::size_t addVariable(const std::string& name, double objective, double lowerBound = -std::numeric_limits<double>::infinity(),
      double upperBound = std::numeric_limits<double>::infinity());

    void addOracle(SeparationOracle* oracle);

    void run();

  protected:
    Solver* _solver;
    std::vector<SeparationOracle*> _oracles;
    std::vector<Solution> _solutions;
    Solution _bestSolution;
    std::chrono::steady_clock::time_point _timeStart; 
  };

}

#endif /* _CORE_H_ */