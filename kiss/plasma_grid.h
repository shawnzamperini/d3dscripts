#include <iostream>
using namespace std;

#ifndef INPUT_FILE
#define INPUT_FILE
#include "input_file.h"
#endif

#ifndef PLASMA_GRID
#define PLASMA_GRID

class Cube{
  public:
    double rad_mid, pol_mid, par_mid;
    double te, ne;

    Cube();
};

class PlasmaGrid{
  public:
    InputFile input_file;

    PlasmaGrid(InputFile* input_file_ptr);
    vector<vector<vector<Cube>>> construct_grid();
    void assign_background(vector<vector<vector<Cube>>> &grid);
};

#endif
