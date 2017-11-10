#include "core.h"

namespace cpm
{
  SolutionData::SolutionData(const std::vector<double>& vals, double objVal)
    : values(vals), objectiveValue(objVal)
  {

  }

  Core::Core(Solver* solver)
    : _solver(solver), _bestSolution(nullptr)
  {

  }

  Core::~Core()
  {
    delete _solver;
  }

  std::size_t Core::addVariable(const std::string& name, double objective, double lowerBound, double upperBound)
  {
    return _solver->addVariable(name, objective, lowerBound, upperBound);
  }

  void Core::addOracle(SeparationOracle* oracle)
  {
    _oracles.push_back(oracle);
  }

  void Core::run()
  {
    _timeStart = std::chrono::steady_clock::now();
    
    std::vector<double> lhs;
    std::vector<double> rhs;
    std::vector<std::size_t> begin;
    std::vector<std::size_t> indices;
    std::vector<double> values;

    std::vector<double> vector;
    bool abort = false;

    double primalBound = -std::numeric_limits<double>::min();
    double dualBound = std::numeric_limits<double>::max();
    while (!abort)
    {
      lhs.clear();
      rhs.clear();
      begin.clear();
      indices.clear();
      values.clear();

      abort = true;
      Solver::Status status = _solver->run();
      if (status == Solver::OPTIMAL)
      {
        double violationLowerBound = 0.0;
        double violationUpperBound = std::numeric_limits<double>::max();

        vector = _solver->point();
        
        dualBound = 0.0;
        for (std::size_t v = 0; v < _solver->numVariables(); ++v)
          dualBound += _solver->objectiveCoefficient(v) * vector[v];

        std::cerr << std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::steady_clock::now() - _timeStart).count()
          << ": Dual bound is " << dualBound << ". Primal bound is " << primalBound << ".\n" << std::flush;

        // TODO: Call stabilizer to change the vector or change LP and resolve.

        for (std::size_t o = 0; o < _oracles.size(); ++o)
        {
          std::size_t oldNumInequalities = lhs.size();
          
          _oracles[o]->separate(true, &vector[0], lhs, rhs, begin, indices, values, violationLowerBound, violationUpperBound);
          if (violationLowerBound < 1.0e-3)
            violationLowerBound = 0.0;

          std::cout << "Oracle " << o << " returned " << (lhs.size() - oldNumInequalities) << " cuts and proved " << violationLowerBound
            << " <= maximum cut violation";
          if (violationUpperBound <= 0.5 * std::numeric_limits<double>::max())
            std::cout << " <= " << violationUpperBound;
          std::cout << ".\n" << std::flush;

          if (_oracles[o]->numFeasiblePoints() > 0)
          {
            // TODO: Ensure that this point is also valid w.r.t. the other oracles. For this, dominance information should be there. This helps to
            //       avoid unnecessary calls to separation oracles. Using a sequence of the latter we might generate a sequence of infeasible points,
            //       each being infeasible except for the last. This works if the modified-variable sets of the oracles induce an implication graph
            //       without cycles.

            std::vector<double> vals;
            for (std::size_t i = 0; i < _oracles[o]->numFeasiblePoints(); ++i)
            {
              _oracles[o]->getFeasiblePoint(i, vals);
              double objVal = 0.0;
              for (std::size_t v = 0; v < _solver->numVariables(); ++v)
                objVal += _solver->objectiveCoefficient(v) * vals[v];
              std::cerr << "Oracle " << o << " returned a solution with objective value " <<  objVal << std::endl;
              _solutions.push_back(std::make_shared<SolutionData>(vals, objVal));
              if (objVal > primalBound)
              {
                primalBound = objVal;
              }
            }
          }

          if (!lhs.empty())
          {
            _solver->addInequalities(lhs, rhs, begin, indices, values);
            abort = false;
            break;
          }
          if (violationUpperBound <= 0.0)
            break;
        }

        // TODO: If abort is true, then vector is feasible. Depending on stabilization it might not be optimal, though
      }
      else if (status == Solver::UNBOUNDED)
      {
        double violationLowerBound = 0.0;
        double violationUpperBound = std::numeric_limits<double>::max();

        vector = _solver->ray();

        // TODO: Call stabilizer to change the vector or change LP and resolve.

        for (std::size_t o = 0; o < _oracles.size(); ++o)
        {
          _oracles[o]->separate(false, &vector[0], lhs, rhs, begin, indices, values, violationLowerBound, violationUpperBound);
          if (!lhs.empty())
          {
            _solver->addInequalities(lhs, rhs, begin, indices, values);
            abort = false;
            break;
          }
          if (violationUpperBound <= 0.0)
            break;
        }

        // TODO: If abort is true, then vector is feasible. Depending on stabilization it might not be optimal, though
      }
      else
      {
        throw std::runtime_error("Unhandled Solver::Status.");
      }
    }
    
    std::cerr << "Total time: " << std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::steady_clock::now() - _timeStart).count()
      << std::endl;
  }

} /* namespace cpm */