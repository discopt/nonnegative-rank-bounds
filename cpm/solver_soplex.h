#ifndef _SOLVER_SOPLEX_H_
#define _SOLVER_SOPLEX_H_

#include <soplex.h>

#include "solver.h"
#include "separation_oracle.h"

namespace cpm
{

  class SolverSoPlex : public Solver
  {
  public:
    SolverSoPlex();

    virtual ~SolverSoPlex();

    virtual std::size_t addVariable(const std::string& name, double objective, double lowerBound, double upperBound) override;

    virtual double objectiveCoefficient(std::size_t variable) const override;

    virtual void addInequalities(std::vector<double>& lhs, std::vector<double>& rhs, std::vector<std::size_t>& begin,
      std::vector<std::size_t>& indices, std::vector<double>& values) override;

    virtual Status run() override;

  protected:
    soplex::SoPlex _spx;
    soplex::DVectorReal _vector;
  };

} /* namespace cpm */

#endif /* _SOLVER_SOPLEX_H_ */