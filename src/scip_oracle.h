#ifndef _SCIP_ORACLE_H_
#define _SCIP_ORACLE_H_

#include <scip/scip.h>

#include <cpm/separation_oracle.h>

#include "scip_exception.h"
#include "slackmatrix.h"

class MaximumWeightRectangleIPOracle : public cpm::SeparationOracle
{
public:
  MaximumWeightRectangleIPOracle(const Slackmatrix& slackmatrix, int priority);

  virtual ~MaximumWeightRectangleIPOracle();

  virtual void separate(bool separatePoint, const double* vector, std::vector<double>& lhs, std::vector<double>& rhs,
    std::vector<std::size_t>& begin, std::vector<std::size_t>& indices, std::vector<double>& values, double& violationLowerBound,
    double& violationUpperBound);

  virtual void getFeasiblePointProperties(bool& scaledDown, bool& scaledUp, std::vector< std::size_t >& modifiableVariables) const;

  virtual std::size_t numFeasiblePoints() const;

  virtual void getFeasiblePoint(std::size_t id, std::vector<double>& point) const;

  void setDoubleParam(const std::string& param, double value);
  
  void setIntParam(const std::string& param, int value);

protected:
  const Slackmatrix* _slackmatrix;
  SCIP* _scip;
  std::vector<SCIP_VAR*> _rowVariables;
  std::vector<SCIP_VAR*> _columnVariables;
  std::vector<SCIP_VAR*> _nonzeroVariables;
  std::vector<double> _feasiblePoint;
};

#endif /* _SCIP_ORACLE_H_ */