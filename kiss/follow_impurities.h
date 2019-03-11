#include <iostream>
using namespace std;

#ifndef INPUT_FILE
#define INPUT_FILE
#include "input_file.h"
#endif

#ifndef PLASMA_GRID
#define PLASMA_GRID
#include "plasma_grid.h"
#endif

#ifndef FOLLOW_IMPURITIES
#define FOLLOW_IMPURITIES

class Impurity{
  public:
    int number;
    double r_pos, p_pos, y_pos;

  Impurity(int number, double rad, double pol, double par);
  Impurity();
};

Impurity launch_impurity(int imp_number, InputFile* input_file_ptr,
                         vector<vector<vector<Cube>>> &grid_ptr);

#endif
