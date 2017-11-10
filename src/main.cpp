
#include <cpm/core.h>
#include <cpm/solver_soplex.h>
#include <cpm/separation_oracle.h>

#include "slackmatrix.h"
#include "enum_oracle.h"
#include "scip_oracle.h"

int main(int argc, char** argv)
{
  std::string fileName = argv[1];
  
  std::size_t numRows, numColumns;
  std::vector<Slackmatrix::Nonzero> nonzeros;
  std::ifstream file(fileName.c_str());
  file >> numRows >> numColumns;
  unsigned int maxEntry = 0.0;
  for (std::size_t row = 0; row < numRows; ++row)
  {
    for (std::size_t column = 0; column < numColumns; ++column)
    {
      unsigned int slack;
      file >> slack;
      assert(slack >= 0);
      if (slack > 0)
      {
        Slackmatrix::Nonzero nz = { row, column, slack };
        nonzeros.push_back(nz);
        maxEntry = std::max(maxEntry, slack);
      }
    }
  }
  Slackmatrix slackmatrix(numRows, numColumns, nonzeros);
  file.close();

  std::cout << "Read " << numRows << "x" << numColumns << " matrix with " << nonzeros.size() << " nonzeros." << std::endl;

  if (maxEntry == 0)
  {
    std::cerr << "Error: matrix is the zero matrix!" << std::endl;
    return SCIP_ERROR;
  }

  double scalingFactor = 1.0 / maxEntry;
  cpm::Core* core = new cpm::Core(new cpm::SolverSoPlex());

  // Variables

  for (std::size_t i = 0; i < slackmatrix.nonzeros.size(); ++i)
  {
    std::stringstream ss;
    ss << "nonzero#" << i << "#" << slackmatrix.nonzeros[i].row << "#" << slackmatrix.nonzeros[i].column;
    core->addVariable(ss.str(), scalingFactor * slackmatrix.nonzeros[i].slack, -std::numeric_limits<double>::infinity(), 1.0);
  }

  // Oracles

  MaximumWeightRectangleEnumOracle* enumOracle = new MaximumWeightRectangleEnumOracle(slackmatrix, 2);
  core->addOracle(enumOracle);
  
//   MaximumWeightRectangleIPOracle* heuristicSCIPOracle = new MaximumWeightRectangleIPOracle(slackmatrix, 1);
//   heuristicSCIPOracle->setDoubleParam("limits/time", 1.0);
//   heuristicSCIPOracle->setIntParam("limits/bestsol", 2);
//   core->addOracle(heuristicSCIPOracle);
  
  MaximumWeightRectangleIPOracle* exactSCIPOracle = new MaximumWeightRectangleIPOracle(slackmatrix, -1);
  core->addOracle(exactSCIPOracle);

  // Run
  
  core->run();

  delete exactSCIPOracle;
//   delete heuristicSCIPOracle;
  delete enumOracle;
  delete core;
}
