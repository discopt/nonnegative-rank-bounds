#include "scip_oracle.h"

#include <iostream>
#include <sstream>
#include <limits>

#include <scip/scipdefplugins.h>

MaximumWeightRectangleIPOracle::MaximumWeightRectangleIPOracle(const Slackmatrix& slackmatrix, int priority)
  : SeparationOracle(slackmatrix.nonzeros.size(), priority)
{
  _slackmatrix = new Slackmatrix(slackmatrix.numRows, slackmatrix.numColumns, slackmatrix.nonzeros);

  SCIP_CALL_EXC(SCIPcreate(&_scip));
  SCIP_CALL_EXC(SCIPcreateProbBasic(_scip, "max-weight-rectangle"));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));
  SCIP_CALL_EXC(SCIPsetIntParam(_scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPsetObjlimit(_scip, 1.0));
  SCIP_CALL_EXC(SCIPsetBoolParam(_scip, "misc/catchctrlc", false));
  SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));
  
  // Disable trivial cuts. TODO: Check how worthy cuts look like.
  
  SCIP_CALL_EXC(SCIPsetIntParam(_scip, "separating/clique/freq", -1));
  SCIP_CALL_EXC(SCIPsetIntParam(_scip, "separating/impliedbounds/freq", -1));

  _rowVariables.resize(_slackmatrix->numRows);
  for (std::size_t row = 0; row < _slackmatrix->numRows; ++row)
  {
    std::stringstream ss;
    ss << "row#" << row;
    SCIP_CALL_EXC(SCIPcreateVarBasic(_scip, &_rowVariables[row], ss.str().c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
    SCIP_CALL_EXC(SCIPaddVar(_scip, _rowVariables[row]));
  }

  _columnVariables.resize(_slackmatrix->numColumns);
  for (std::size_t column = 0; column < _slackmatrix->numColumns; ++column)
  {
    std::stringstream ss;
    ss << "column#" << column;
    SCIP_CALL_EXC(SCIPcreateVarBasic(_scip, &_columnVariables[column], ss.str().c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
    SCIP_CALL_EXC(SCIPaddVar(_scip, _columnVariables[column]));
  }

  _nonzeroVariables.resize(_slackmatrix->nonzeros.size());
  for (std::size_t i = 0; i < _nonzeroVariables.size(); ++i)
  {
    std::stringstream ss;
    ss << "nonzero#" << i << "#row#" << _slackmatrix->nonzeros[i].row << "#column#" << _slackmatrix->nonzeros[i].column;
    SCIP_CALL_EXC(SCIPcreateVarBasic(_scip, &_nonzeroVariables[i], ss.str().c_str(), 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY));
    SCIP_CALL_EXC(SCIPaddVar(_scip, _nonzeroVariables[i]));
  }

  // x_row + y_column <= 1 for zero entry (row,column).

  for (std::size_t row = 0; row < _slackmatrix->numRows; ++row)
  {
    for (std::size_t column = 0; column < _slackmatrix->numColumns; ++column)
    {
      if (_slackmatrix->denseIndices[row][column] < std::numeric_limits<std::size_t>::max())
        continue;

      SCIP_CONS* cons = NULL;
      std::stringstream ss;
      ss << "zero#" << row << "#" << column;
      SCIP_CALL_EXC(SCIPcreateConsBasicLinear(_scip, &cons, ss.str().c_str(), 0, NULL, NULL, -SCIPinfinity(_scip), 1.0));
      SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _rowVariables[row], 1.0));
      SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _columnVariables[column], 1.0));
      SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
      SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
    }
  }

  // x_row - z_i >= 0 for i'th nonzero (row,column).

  for (std::size_t i = 0; i < _slackmatrix->nonzeros.size(); ++i)
  {
    SCIP_CONS* cons = NULL;
    std::stringstream ss;
    ss << "nonzero#" << i << "#row#" << _slackmatrix->nonzeros[i].row;
    SCIP_CALL_EXC(SCIPcreateConsBasicLinear(_scip, &cons, ss.str().c_str(), 0, NULL, NULL, 0.0, 1.0));
    SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _rowVariables[_slackmatrix->nonzeros[i].row], 1.0));
    SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _nonzeroVariables[i], -1.0));
    SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
    SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
  }

  // y_column - z_i >= 0 for i'th nonzero (row,column).

  for (std::size_t i = 0; i < _slackmatrix->nonzeros.size(); ++i)
  {
    SCIP_CONS* cons = NULL;
    std::stringstream ss;
    ss << "nonzero#" << i << "#column#" << _slackmatrix->nonzeros[i].column;
    SCIP_CALL_EXC(SCIPcreateConsBasicLinear(_scip, &cons, ss.str().c_str(), 0, NULL, NULL, 0.0, 1.0));
    SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _columnVariables[_slackmatrix->nonzeros[i].column], 1.0));
    SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _nonzeroVariables[i], -1.0));
    SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
    SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
  }

  // x_row + y_column - z_i <= 1 for i'th nonzero (row,column).

  for (std::size_t i = 0; i < _slackmatrix->nonzeros.size(); ++i)
  {
    SCIP_CONS* cons = NULL;
    std::stringstream ss;
    ss << "nonzero#" << i << "#row#" << _slackmatrix->nonzeros[i].row << "#column#" << _slackmatrix->nonzeros[i].column;
    SCIP_CALL_EXC(SCIPcreateConsBasicLinear(_scip, &cons, ss.str().c_str(), 0, NULL, NULL, -SCIPinfinity(_scip), 1.0));
    SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _rowVariables[_slackmatrix->nonzeros[i].row], 1.0));
    SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _columnVariables[_slackmatrix->nonzeros[i].column], 1.0));
    SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _nonzeroVariables[i], -1.0));
    SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
    SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
  }
}

MaximumWeightRectangleIPOracle::~MaximumWeightRectangleIPOracle()
{
  for (std::size_t row = 0; row < _slackmatrix->numRows; ++row)
    SCIP_CALL_EXC( SCIPreleaseVar(_scip, &_rowVariables[row]) );
  for (std::size_t column = 0; column < _slackmatrix->numColumns; ++column)
    SCIP_CALL_EXC( SCIPreleaseVar(_scip, &_columnVariables[column]) );
  for (std::size_t i = 0; i < _nonzeroVariables.size(); ++i)
  {
    SCIP_CALL_EXC( SCIPreleaseVar(_scip, &_nonzeroVariables[i]) );
  }
  
  SCIP_CALL_EXC(SCIPfreeProb(_scip));
  SCIP_CALL_EXC(SCIPfree(&_scip));
  delete _slackmatrix;
}

void MaximumWeightRectangleIPOracle::separate(bool separatePoint, const double* vector, std::vector<double>& lhs, std::vector<double>& rhs,
  std::vector<std::size_t>& begin, std::vector<std::size_t>& indices, std::vector<double>& values, double& violationLowerBound,
  double& violationUpperBound)
{
//     double dummy;
//     SeparationOracle::separate(separatePoint, vector, lhs, rhs, begin, indices, values, violationLowerBound, dummy);

  for (std::size_t i = 0; i < _slackmatrix->nonzeros.size(); ++i)
  {
    SCIP_CALL_EXC( SCIPchgVarObj(_scip, _nonzeroVariables[i], vector[i]) );
  }

  SCIP_CALL_EXC( SCIPsolve(_scip) );

  SCIP_SOL* bestSol = SCIPgetBestSol(_scip);
  if (bestSol != NULL)
  {
    double bestValue = SCIPgetSolOrigObj(_scip, SCIPgetBestSol(_scip));
    violationLowerBound = bestValue - 1.0;
    violationUpperBound = SCIPgetDualbound(_scip) - 1.0;
  }
  
//   SCIP_CUT** cuts = NULL;
//   cuts = SCIPgetPoolCuts(_scip);
//   for (int i = 0; i < SCIPgetNPoolCuts(_scip); ++i)
//   {
//     SCIP_ROW* row = SCIPcutGetRow(cuts[i]);
//     SCIPprintRow(_scip, row, stderr);
//     std::cerr << std::endl;
//   }
//   SCIPprintStatistics(_scip, stderr);
//   SCIP_CALL_EXC( SCIPprintBestSol(_scip, stdout, false) );

  int numSols = SCIPgetNSols(_scip);
  SCIP_SOL** sols = SCIPgetSols(_scip);
  for (int sol = 0; sol < numSols; ++sol)
  {
    double objectiveValue = SCIPgetSolOrigObj(_scip, sols[sol]);
    if (objectiveValue - 1.0 <= 1.0e-3) // TODO: epsilon
      continue;

    lhs.push_back(-std::numeric_limits<double>::infinity());
    rhs.push_back(1.0);
    begin.push_back(indices.size());

    for (std::size_t v = 0; v < ambientDimension(); ++v)
    {
      if (SCIPgetSolVal(_scip, sols[sol], _nonzeroVariables[v]) > 0.5)
      {
        indices.push_back(v);
        values.push_back(1.0);
      }
    }
//     std::size_t numRectangleRows = 0;
//     for (std::size_t r = 0; r < _rowVariables.size(); ++r)
//     {
//       if (SCIPgetSolVal(_scip, sols[sol], _rowVariables[r]) > 0.5)
//         ++numRectangleRows;
//     }
//     std::size_t numRectangleColumns = 0;
//     for (std::size_t c = 0; c < _columnVariables.size(); ++c)
//     {
//       if (SCIPgetSolVal(_scip, sols[sol], _columnVariables[c]) > 0.5)
//         ++numRectangleColumns;
//     }
//     std::cerr << " [" << numRectangleRows << "x" << numRectangleColumns << " rect]" << std::flush;
  }

  if (violationUpperBound < 1000)
  {
    _feasiblePoint.resize(ambientDimension());
    for (std::size_t v = 0; v < ambientDimension(); ++v)
      _feasiblePoint[v] = vector[v] / SCIPgetDualbound(_scip);
  }

  SCIP_CALL_EXC( SCIPfreeSolve(_scip, true) );
  SCIP_CALL_EXC( SCIPfreeTransform(_scip) );
}

void MaximumWeightRectangleIPOracle::getFeasiblePointProperties(bool& scaledDown, bool& scaledUp, std::vector<std::size_t>& modifiableVariables) const
{
  cpm::SeparationOracle::getFeasiblePointProperties(scaledDown, scaledUp, modifiableVariables);
  scaledDown = true;
}

std::size_t MaximumWeightRectangleIPOracle::numFeasiblePoints() const
{
  return _feasiblePoint.empty() ? 0 : 1;
}

void MaximumWeightRectangleIPOracle::getFeasiblePoint(std::size_t id, std::vector<double>& point) const
{
  assert(id == 0);
  assert(!_feasiblePoint.empty());
  point = _feasiblePoint;
}

void MaximumWeightRectangleIPOracle::setDoubleParam(const std::string& param, double value)
{
  SCIP_CALL_EXC(SCIPsetRealParam(_scip, param.c_str(), value));
}


void MaximumWeightRectangleIPOracle::setIntParam(const std::string& param, int value)
{
  SCIP_CALL_EXC(SCIPsetIntParam(_scip, param.c_str(), value));
}


