#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <tuple>
#include "read_input.h"
#include "plasma_grid.h"
#include "read_input.h"

#define BOLTZ_CONST 8.617e-5 // eV K-1
#define AMU 931.49e6 / (3e8)^2 // eV s2 m-2

using namespace std;

class Impurity{
  public:
    int number;
    double r_pos, p_pos, y_pos, vel;

  Impurity();
  Impurity(int number, double rad, double pol, double par);
  //Impurity launch_impurity(InputFile* input_file_ptr);
};

// Constructer for an impurity ion.
Impurity::Impurity(int number, double rad, double pol, double par){
  number  = number;
  r_pos = rad; p_pos = pol; y_pos = par;
}

Impurity::Impurity(){}

// Function to see what Cube an Impurity is in, and returns that Cube.
Cube get_my_cube(Impurity &imp, InputFile* input_file_ptr,
                  vector<vector<vector<Cube>>> &grid){

  int rad_bin, pol_bin, par_bin;
  double tmp_dist, tmp_dist_old;
  tuple <int, int, int> bins;

  InputFile* input_file = input_file_ptr;


  tmp_dist_old = 999999;
  for(int i=0; i<input_file->num_rad_bins; i++){
    tmp_dist = abs(grid[i][0][0].rad_mid - imp.r_pos);
    //cout << "Radial tmp_dist = " << tmp_dist << endl;
    if(tmp_dist < tmp_dist_old){
      rad_bin = i;
      tmp_dist_old = tmp_dist;
    }
  }
  tmp_dist_old = 999999;
  for(int i=0; i<input_file->num_pol_bins; i++){
    tmp_dist = abs(grid[0][i][0].pol_mid - imp.p_pos);
    if(tmp_dist < tmp_dist_old){
      pol_bin = i;
      tmp_dist_old = tmp_dist;
    }
  }
  tmp_dist_old = 999999;
  for(int i=0; i<input_file->num_par_bins; i++){
    //cout << "grid[0][0][" << i << "] = " << grid[0][0][i].par_mid << endl;
    //cout << "imp.y_pos = " << imp.y_pos << endl;
    tmp_dist = abs(grid[0][0][i].par_mid - imp.y_pos);
    //cout << "Parallel tmp_dist = " << tmp_dist << endl;
    //cout << "Parallel tmp_dist = " << tmp_dist << endl;
    //cout << "Parallel tmp_dist_old = " << tmp_dist_old << endl;
    if(tmp_dist < tmp_dist_old){
      //cout << "Parallel tmp_dist = " << tmp_dist << endl;
      par_bin = i;
      tmp_dist_old = tmp_dist;
    }
  }

  //cout << "Checkpoint #5" << endl;

  // Can't figure out how to do this correctly. Do tuple instead.
  Cube cube = grid[rad_bin][pol_bin][par_bin];
  return cube;
  //return grid[rad_bin][pol_bin][par_bin];

  // Return as tuple.
  //return make_tuple(rad_bin, pol_bin, par_bin);
}


// Monte Carlo approach to choose initial position, velocity of impurity.
Impurity launch_impurity(int imp_number, InputFile* input_file_ptr,
                         vector<vector<vector<Cube>>> &grid_ptr){

  //cout << "Checkpoint #3" << endl;
  double launch_r, launch_p, launch_y, v_therm;
  //tuple <int, int, int> bins;

  // Load in input file to access input parameters.
  InputFile* input_file = input_file_ptr;

  //cout << "Creating random injection locations..." << endl;;
  // Find initial launch position. Note: All options are probably going to
  // use this approach for the initial position, so no need to rewrite the
  // code for each case. If this changes, may need to move these lines into
  // each case, or just overwrite these variables for the special cases.
  launch_r = input_file->launch_rmin + double(rand()) / RAND_MAX
             * (input_file->launch_rmax - input_file->launch_rmin);
  launch_p = input_file->launch_pmin + double(rand()) / RAND_MAX
             * (input_file->launch_pmax - input_file->launch_pmin);
  launch_y = input_file->launch_ymin + double(rand()) / RAND_MAX
             * (input_file->launch_ymax - input_file->launch_ymin);

  Impurity imp = Impurity(imp_number, launch_r, launch_p, launch_y);

  Cube cube = get_my_cube(imp, input_file_ptr, grid_ptr);

  switch(input_file->inj_opt){

    // Opt 0: Launch with random velocity between -v_therm and +v_therm.
    case 0:
      v_therm = sqrt(BOLTZ_CONST * cube.te / input_file->imp_amu);
      imp.vel = -v_therm + 2 * double(rand()) / RAND_MAX * v_therm;
      break;

    case 1:
      v_therm = sqrt(BOLTZ_CONST * cube.te / input_file->imp_amu);
      imp.vel = - double(rand()) / RAND_MAX * v_therm;
      break;

  cout << "Launch velocity = " << imp.vel << endl;;

  }

  // Create the Impurity particle object and return it.
  return imp;
}
