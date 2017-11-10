#include "slackmatrix.h"

#include <limits>

Slackmatrix::Slackmatrix(std::size_t nRows, std::size_t nColumns, const std::vector<Nonzero>& nzs)
  : numRows(nRows), numColumns(nColumns), nonzeros(nzs)
{
  denseIndices.resize(nRows);
  for (std::size_t row = 0; row < nRows; ++row)
  {
    denseIndices[row].resize(nColumns, std::numeric_limits<std::size_t>::max());
  }
  for (std::size_t i = 0; i < nzs.size(); ++i)
  {
    denseIndices[nonzeros[i].row][nonzeros[i].column] = i;
  }
}

Slackmatrix::~Slackmatrix()
{
  
}