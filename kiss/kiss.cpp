#include "read_input.h"
#include "plasma_grid.h"
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#define VERSION "0.1"


int main(int argc, char **argv){

  cout << "KISS Version " << VERSION << endl;

  // Get input filename.
  string input_filename = argv[1];

  // The grid to hold the background plasma. This is simply a 3D array, just
  // constructed with vectors (a vector of vectors of vectors).
  vector<vector<vector<Cube>>> grid;

  // Create the InputFile object, loading in all the parameters.
  InputFile input_file(input_filename);

  // Pass in a POINTER to the InputFile into the PlasmaGrid (this is the way
  // you're supposed to do it I think).
  PlasmaGrid plasma_grid(&input_file);

  // Create the grid, assigning each Cube an (R, P, Y) coordinate.
  grid = plasma_grid.construct_grid();

  // Fill in Te, ne, etc. data into the grid. The depends on the input file options,
  // such as what kind of SOL or if the collector probe creates its own plasma.
  plasma_grid.assign_background(grid);

  return 0;
}
