#ifndef _SLACKMATRIX_H_
#define _SLACKMATRIX_H_

#include <vector>

class Slackmatrix
{
public:
  struct Nonzero
  {
    std::size_t row;
    std::size_t column;
    std::size_t slack;
  };

  std::size_t numRows;
  std::size_t numColumns;
  std::vector<Nonzero> nonzeros;
  std::vector<std::vector<std::size_t> > denseIndices;

  Slackmatrix(std::size_t numRows, std::size_t numColumns, const std::vector<Nonzero>& nonzeros);
  ~Slackmatrix();
};

#endif /* _SLACKMATRIX_H_ */