#include "enum_oracle.h"

#include <cassert>
#include <limits>
#include <iostream>

MaximumWeightRectangleEnumOracle::MaximumWeightRectangleEnumOracle(const Slackmatrix& slackmatrix, int priority)
  : SeparationOracle(slackmatrix.nonzeros.size(), priority)
{
  _slackmatrix = new Slackmatrix(slackmatrix.numRows, slackmatrix.numColumns, slackmatrix.nonzeros);
  _sumPositiveByRow.resize(_slackmatrix->numRows);
  _sumPositiveByColumn.resize(_slackmatrix->numColumns);
  _variables.reserve(2 * std::max(_slackmatrix->numRows, _slackmatrix->numColumns));
}

MaximumWeightRectangleEnumOracle::~MaximumWeightRectangleEnumOracle()
{
  delete _slackmatrix;
}


void MaximumWeightRectangleEnumOracle::separate(bool separatePoint, const double* vector, std::vector<double>& lhs, std::vector<double>& rhs, 
  std::vector<std::size_t>& begin, std::vector<std::size_t>& indices, std::vector<double>& values, double& violationLowerBound,
  double& violationUpperBound)
{
  assert(separatePoint);

  for (std::size_t r = 0; r < _slackmatrix->numRows; ++r)
    _sumPositiveByRow[r] = 0.0;
  for (std::size_t c = 0; c < _slackmatrix->numColumns; ++c)
    _sumPositiveByColumn[c] = 0.0;

  // Sum up positive entries row- and column-wise.

  for (std::size_t i = 0; i < _slackmatrix->nonzeros.size(); ++i)
  {
    double value = vector[i];
    if (value > 0.0)
    {
      _sumPositiveByRow[_slackmatrix->nonzeros[i].row] += value;
      _sumPositiveByColumn[_slackmatrix->nonzeros[i].column] += value;
    }
  }

  // Find all rectangles consisting of 1 row.

  for (std::size_t r = 0; r < _slackmatrix->numRows; ++r)
  {
    double violation = _sumPositiveByRow[r] - 1.0;
    if (violation > 1.0e-3)
    {
      lhs.push_back(-std::numeric_limits<double>::infinity());
      rhs.push_back(1.0);
      begin.push_back(indices.size());
      for (std::size_t c = 0; c < _slackmatrix->numColumns; ++c)
      {
        std::size_t v = _slackmatrix->denseIndices[r][c];
        if (v == std::numeric_limits<std::size_t>::max())
          continue;

        if (vector[v] > 0.0)
        {
          indices.push_back(v);
          values.push_back(1.0);
        }
      }
      if (violation > violationLowerBound)
        violationLowerBound = violation;
    }
  }

  // Find all rectangles consisting of 1 column.

  for (std::size_t c = 0; c < _slackmatrix->numColumns; ++c)
  {
    double violation = _sumPositiveByColumn[c] - 1.0;
    if (violation > 1.0e-3)
    {
      lhs.push_back(-std::numeric_limits<double>::infinity());
      rhs.push_back(1.0);
      begin.push_back(indices.size());
      for (std::size_t r = 0; r < _slackmatrix->numRows; ++r)
      {
        std::size_t v = _slackmatrix->denseIndices[r][c];
        if (v == std::numeric_limits<std::size_t>::max())
          continue;

        if (vector[v] > 0.0)
        {
          indices.push_back(v);
          values.push_back(1.0);
        }
      }
      if (violation > violationLowerBound)
        violationLowerBound = violation;
    }
  }

  // Check pairs of rows.

  if (lhs.empty())
  {
    for (std::size_t r1 = 0; r1 < _slackmatrix->numRows; ++r1)
    {
      for (std::size_t r2 = r1 + 1; r2 < _slackmatrix->numRows; ++r2)
      {
        if (_sumPositiveByRow[r1] + _sumPositiveByRow[r2] - 1.0 <= 1.0e-3)
          continue;

//         std::cerr << "Candidate row pair " << r1 << " " << r2 << std::endl;
        _variables.clear();
        double violation = -1.0;
        for (std::size_t c = 0; c < _slackmatrix->numColumns; ++c)
        {
          std::size_t v1 = _slackmatrix->denseIndices[r1][c];
          std::size_t v2 = _slackmatrix->denseIndices[r2][c];
          
          if (v1 == std::numeric_limits<std::size_t>::max())
            continue;
          if (v2 == std::numeric_limits<std::size_t>::max())
            continue;
          double contribution = vector[v1] + vector[v2];
          if (contribution > 0.0)
          {
            violation += contribution;
            _variables.push_back(v1);
            _variables.push_back(v2);
          }
        }
        
        if (violation > 1.0e-3)
        {
          lhs.push_back(-std::numeric_limits<double>::infinity());
          rhs.push_back(1.0);
          begin.push_back(indices.size());
          for (std::size_t i = 0; i < _variables.size(); ++i)
          {
            indices.push_back(_variables[i]);
            values.push_back(1.0);
          }
        }
      }
    }
  }

  // Check pairs of columns.

  if (lhs.empty())
  {
    for (std::size_t c1 = 0; c1 < _slackmatrix->numColumns; ++c1)
    {
      for (std::size_t c2 = c1 + 1; c2 < _slackmatrix->numColumns; ++c2)
      {
        if (_sumPositiveByColumn[c1] + _sumPositiveByColumn[c2] - 1.0 <= 1.0e-3)
          continue;

//         std::cerr << "Candidate column pair " << c1 << " " << c2 << std::endl;
        _variables.clear();
        double violation = -1.0;
        for (std::size_t r = 0; r < _slackmatrix->numRows; ++r)
        {
          std::size_t v1 = _slackmatrix->denseIndices[r][c1];
          std::size_t v2 = _slackmatrix->denseIndices[r][c2];
          
          if (v1 == std::numeric_limits<std::size_t>::max())
            continue;
          if (v2 == std::numeric_limits<std::size_t>::max())
            continue;
          double contribution = vector[v1] + vector[v2];
          if (contribution > 0.0)
          {
            violation += contribution;
            _variables.push_back(v1);
            _variables.push_back(v2);
          }
        }
        
        if (violation > 1.0e-3)
        {
          lhs.push_back(-std::numeric_limits<double>::infinity());
          rhs.push_back(1.0);
          begin.push_back(indices.size());
          for (std::size_t i = 0; i < _variables.size(); ++i)
          {
            indices.push_back(_variables[i]);
            values.push_back(1.0);
          }
        }
      }
    }
  }
}

