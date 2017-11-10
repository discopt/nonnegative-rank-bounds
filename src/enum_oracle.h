#ifndef _ENUM_ORACLE_H_
#define _ENUM_ORACLE_H_

#include <cpm/separation_oracle.h>

#include "slackmatrix.h"

class MaximumWeightRectangleEnumOracle : public cpm::SeparationOracle
{
public:
  MaximumWeightRectangleEnumOracle(const Slackmatrix& slackmatrix, int priority);

  virtual ~MaximumWeightRectangleEnumOracle();

  virtual void separate(bool separatePoint, const double* vector, std::vector<double>& lhs, std::vector<double>& rhs,
    std::vector<std::size_t>& begin, std::vector<std::size_t>& indices, std::vector<double>& values, double& violationLowerBound,
    double& violationUpperBound);

protected:
  const Slackmatrix* _slackmatrix;
  std::vector<double> _sumPositiveByRow;
  std::vector<double> _sumPositiveByColumn;
  std::vector<std::size_t> _variables;
};

#endif /* _ENUM_ORACLE_H_ */