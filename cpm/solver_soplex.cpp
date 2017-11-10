#include "solver_soplex.h"

namespace cpm
{
  SolverSoPlex::SolverSoPlex()
  {
    _spx.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
    _spx.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_WARNING);
  }

  SolverSoPlex::~SolverSoPlex()
  {

  }

  std::size_t SolverSoPlex::addVariable(const std::string& name, double objective, double lowerBound, double upperBound)
  {
    assert(_spx.numColsReal() == (int) numVariables());

    soplex::LPColReal col;
    col.setLower(lowerBound);
    col.setUpper(upperBound);
    col.setObj(objective);
    _spx.addColReal(col);

    _point.push_back(0.0);
    _ray.push_back(0.0);

    return Solver::addVariable(name);
  }

  double SolverSoPlex::objectiveCoefficient(std::size_t variable) const
  {
    return _spx.objReal(variable);
  }

  void SolverSoPlex::addInequalities(std::vector<double>& lhs, std::vector<double>& rhs, std::vector<std::size_t>& begin,
    std::vector<std::size_t>& indices, std::vector<double>& values)
  {
    assert(lhs.size() == rhs.size());
    assert(lhs.size() == begin.size() || lhs.size() + 1 == begin.size());
    assert(indices.size() == values.size());

    soplex::LPRowSetReal rows(lhs.size());
    soplex::DSVectorReal vector;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
      vector.clear();
      std::size_t first = begin[i];
      std::size_t beyond = (i+1 < lhs.size()) ? begin[i+1] : indices.size();
      for (std::size_t p = first; p < beyond; ++p)
        vector.add(indices[p], values[p]);
      double l = std::max(lhs[i], -soplex::infinity);
      double r = std::min(rhs[i], soplex::infinity);
      rows.add(l, vector, r);
    }
    _spx.addRowsReal(rows);
  }

  Solver::Status SolverSoPlex::run()
  {
    _vector.reDim(numVariables(), false);

    std::cerr << "SolverSoPlex: Solving LP with " << _spx.numRowsReal() << " rows and " << _spx.numColsReal() << " cols..." << std::endl;

    soplex::SPxSolver::Status status = _spx.solve();

    std::cerr << "SolverSoPlex: Status = " << status << ", #iters = " << _spx.numIterations() << ", obj.val = " << _spx.objValueReal() << std::endl;
    
    if (_spx.hasPrimal())
    {
      _spx.getPrimalReal(_vector);
          
      for (std::size_t v = 0; v < numVariables(); ++v)
        _point[v] = _vector[v];
    }
    if (_spx.hasPrimalRay())
    {
      _spx.getPrimalRayReal(_vector);
      for (std::size_t v = 0; v < numVariables(); ++v)
        _ray[v] = _vector[v];
    }
    switch (status)
    {
      case soplex::SPxSolver::NO_RATIOTESTER:
        throw std::runtime_error("SoPlex.solve() did not find a ratio tester.");
      case soplex::SPxSolver::NO_PRICER:
        throw std::runtime_error("SoPlex.solve() did not find a pricer.");
      case soplex::SPxSolver::NO_SOLVER:
        throw std::runtime_error("SoPlex.solve() did not find a solver.");
      case soplex::SPxSolver::NOT_INIT:
        throw std::runtime_error("SoPlex.solve() not initialized.");
      case soplex::SPxSolver::ABORT_EXDECOMP:
        throw std::runtime_error("SoPlex.solve() aborted to exit decomposition simplex.");
      case soplex::SPxSolver::ABORT_DECOMP:
        throw std::runtime_error("SoPlex.solve() aborted due to commence decomposition simplex.");
      case soplex::SPxSolver::ABORT_CYCLING:
        throw std::runtime_error("SoPlex.solve() aborted due to detection of cycling.");
      case soplex::SPxSolver::ABORT_TIME:
        return Solver::TIME_LIMIT;
      case soplex::SPxSolver::ABORT_ITER:
        return Solver::ITERATION_LIMIT;
      case soplex::SPxSolver::ABORT_VALUE:
        throw std::runtime_error("SoPlex.solve() aborted due to objective limit.");
      case soplex::SPxSolver::SINGULAR:
        return Solver::NUMERICAL_TROUBLE;
      case soplex::SPxSolver::NO_PROBLEM:
        throw std::runtime_error("SoPlex.solve() has no problem.");
      case soplex::SPxSolver::REGULAR:
        throw std::runtime_error("SoPlex.solve() has a regular basis.");
      case soplex::SPxSolver::RUNNING:
        throw std::runtime_error("SoPlex.solve() is running.");
      case soplex::SPxSolver::UNKNOWN:
      case soplex::SPxSolver::ERROR:
        return Solver::ERROR;
      case soplex::SPxSolver::OPTIMAL:
        return Solver::OPTIMAL;
      case soplex::SPxSolver::UNBOUNDED:
        return Solver::UNBOUNDED;
      case soplex::SPxSolver::INFEASIBLE:
        return Solver::INFEASIBLE;
      case soplex::SPxSolver::INForUNBD:
        _spx.setIntParam(soplex::SoPlex::SIMPLIFIER, soplex::SoPlex::SIMPLIFIER_OFF);
    }
    return Solver::ERROR;
  }
  
  
/*
  void SoPlexKelley::run()
  {
    std::vector<double> lhs;
    std::vector<double> rhs;
    std::vector<std::size_t> begin;
    std::vector<std::size_t> indices;
    std::vector<double> values;
    std::vector<double> feasiblePoint;

    soplex::DVectorReal vector;
    vector.reDim(_spx.numColsReal());
    bool done = false;
    while (!done)
    {
      soplex::SPxSolver::Status status = _spx.solve();
      if (status == soplex::SPxSolver::OPTIMAL)
      {
        _spx.getPrimalReal(vector);
        std::cerr << "Current point = " << vector << std::endl;
        
        double violationLowerBound = 0.0;
        double violationUpperBound = std::numeric_limits<double>::max();
        done = true;
        for (std::size_t o = 0; o < _oracles.size(); ++o)
        {
          _oracles[o]->separate(true, vector.get_const_ptr(), lhs, rhs, begin, indices, values, violationLowerBound, violationUpperBound, 
            feasiblePoint);
          if (!lhs.empty())
          {
            addInequalities(lhs, rhs, begin, indices, values);
            done = false;
            break;
          }
          if (violationUpperBound <= 0.0)
            break;
        }
      }
      else if (status == soplex::SPxSolver::UNBOUNDED)
      {
        _spx.getPrimalRayReal(vector);
        std::cerr << "Current ray = " << vector << std::endl;

        double violationLowerBound = 0.0;
        double violationUpperBound = std::numeric_limits<double>::max();
        done = true;
        for (std::size_t o = 0; o < _oracles.size(); ++o)
        {
          _oracles[o]->separate(false, vector.get_const_ptr(), lhs, rhs, begin, indices, values, violationLowerBound, violationUpperBound, 
            feasiblePoint);
          if (!lhs.empty())
          {
            addInequalities(lhs, rhs, begin, indices, values);
            done = false;
            break;
          }
          if (violationUpperBound <= 0.0)
            break;
        }
      }
    }
  }*/

}