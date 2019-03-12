#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <tuple>
#include <fstream>
#include "read_input.h"
#include "plasma_grid.h"
#include "read_input.h"

#define BOLTZ_CONST 8.617e-5 // eV K-1
#define AMU         931.49e6 / (3e8 * 3e8) // eV s2 m-2
#define AMU_KG      1.661e-27 // kg
#define ELEC        1.609e-19 // C

using namespace std;

class Impurity{
  public:
    int number, Z;
    double r_pos, p_pos, y_pos, vel, mass;

  Impurity();
  Impurity(int number, double rad, double pol, double par);
  //Impurity launch_impurity(InputFile* input_file_ptr);
};

// Constructer for an impurity ion.
Impurity::Impurity(int imp_number, double rad, double pol, double par){
  number  = imp_number;
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

    // Opt 0: Launch with random velocity between -v_therm and +v_therm in Y direction.
    case 0:
      //v_therm = sqrt(BOLTZ_CONST * cube.te / input_file->imp_amu);
      v_therm = sqrt(cube.te / (input_file->imp_amu * AMU));
      imp.vel = -v_therm + 2 * double(rand()) / RAND_MAX * v_therm;
      break;

    case 1:
      //v_therm = sqrt(BOLTZ_CONST * cube.te / input_file->imp_amu);
      v_therm = sqrt(cube.te / (input_file->imp_amu * AMU));
      imp.vel = - double(rand()) / RAND_MAX * v_therm * 1e3;
      break;
  }
  //cout << "Launch velocity = " << imp.vel << endl;

  // Create the Impurity particle object and return it.
  imp.mass = input_file->imp_amu * AMU_KG;
  imp.Z    = input_file->imp_Z;
  return imp;
}

// Function to determine if particle has hit a surface. Criteria is if the
// particle's new position passes through a surface.
bool not_hit_surface(double old_r, double old_p, double old_y,
                     double new_r, double new_p, double new_y,
                     InputFile* input_file){

  // Check if it hit the right or left wall.
  if(new_y > input_file->par_halfwidth){
    cout << "Hit right wall." << endl;
    return false;
  }
  else if(new_y < -input_file->par_halfwidth){
    cout << "Hit left wall." << endl;
    return false;
  }
  else if(new_r > input_file->rad_fullwidth){
    cout << "Hit bottom wall." << endl;
    return false;
  }

  // Check if it hit the collector probe.

  // If you're in the right radial range of the probe...
  if(old_r < input_file->probe_tip || new_r < input_file->probe_tip){

    // ...and also in the right poloidal range...
    if(old_p < input_file->probe_width/2.0 || new_p < input_file->probe_width/2.0){

      // ...and if you crossed the parallel location of the probe...
      if(old_y < input_file->probe_par_loc && input_file->probe_par_loc < new_y){ // From the left.
        cout << "Hit left side of probe." << endl;
        return false;
      }
      else if(new_y < input_file->probe_par_loc && input_file->probe_par_loc < old_y){  // From the right.
        cout << "Hit right side of probe." << endl;
        return false;
      }
    }
  }
  return true;
}

// Function to follow an impurity ion to its death.
Impurity follow_to_end(Impurity &imp, InputFile* input_file_ptr,
                       vector<vector<vector<Cube>>> &grid){

    // Save the input file into a new pointer (just get rid of the ptr label).
    InputFile* input_file = input_file_ptr;

    // Keep track of how long the particle took to reach the end.
    int time_start = time(nullptr);
    //double time_step = input_file->time_step;

    double old_r = imp.r_pos, old_p = imp.p_pos, old_y = imp.y_pos;
    double new_r, new_p, new_y;
    //double ff, fig, feg, fe;
    int num_iter = 1;
    double time_spent = 0;

    // Open file for debug_track if it's on.
    ofstream track_file;
    if(input_file->debug_track){
      track_file.open("track.txt", fstream::app);
      track_file << "Particle #" << imp.number << endl;
      track_file << old_r << "," << old_p << "," << old_y << endl;
    }

    // Continue following particle until not_hit_surface returns false.
    do{
      //cout << "Start iteration #" << num_iter << endl;
      // Compute diffusion contributions to position change.
      //double delta_rad = 0, delta_pol = 0, delta_par = 0;
      double delta_rad = sqrt(2 * input_file->rad_diff * input_file->time_step);

      // Poloidal diffusion could go either way. Randomly choose.
      double delta_pol = (double(rand()) / RAND_MAX -1) * sqrt(2 * input_file->pol_diff * input_file->time_step);
      //cout << "Delta pol. = " << delta_pol << endl;

      // Calculate the forces. First get the Cube it's in.
      Cube cube = get_my_cube(imp, input_file, grid);

      // Friction force, FF. Assumes deuterium background.
      double vel_i = cube.mach * sqrt(2.0 * cube.te / (2.014 * AMU));
      double vel_z = imp.vel;
      double tau_s = (1.47e13 * input_file->imp_amu * cube.te * sqrt(cube.te / 2.014)) /
                     ((1 + 2.014 / input_file->imp_amu) * cube.ne * double(imp.Z) * double(imp.Z) * 15);
      double ff    = imp.mass * (vel_i - vel_z) / tau_s;

      // Electric field force, FE.
      double fe = imp.Z * ELEC * cube.elec;

      // Electron/ion temperature gradient force, FeG/FiG.
      double te_grad;
        if(cube.par_bin == 0){
        // If we're all the way at the left, use the bin in front.
        te_grad = (grid[cube.rad_bin][cube.pol_bin][cube.par_bin+1].te - cube.te) /
                  (grid[cube.rad_bin][cube.pol_bin][cube.par_bin+1].par_mid - cube.par_mid);
      }
      else{
        // If we're all the way to the right, use the bin behind. Or if we're
        // anywhere else in the volume.
        te_grad = (grid[cube.rad_bin][cube.pol_bin][cube.par_bin-1].te - cube.te) /
                  (grid[cube.rad_bin][cube.pol_bin][cube.par_bin-1].par_mid - cube.par_mid);
      }

      double alpha_e = 0.71 * pow(imp.Z,2);
      double mu = imp.mass / (imp.mass + 2.014 * AMU);
      double beta_i = 3 * (mu+5.0*sqrt(2.0)*pow(imp.Z,2)*(1.1*pow(mu,2.5)-0.35*pow(mu,1.5))-1.0) /
                      (2.6-2*mu+5.4*pow(mu,2));

      double feg = alpha_e * te_grad;
      double fig = beta_i  * te_grad; // Assumes Ti=Te;

      double ftot = ff + fe + feg + fig;

      //cout << "  FF =  " << ff << endl;
      //cout << "  FE =  " << fe << endl;
      //cout << "  FeG = " << feg << endl;
      //cout << "  FiG = " << fig << endl;
      //cout << "  FTOT = " << ftot << endl;

      // Find distance travelled due to force acceleration. a=F/m.
      double delta_par = imp.vel * input_file->time_step
                         + 0.5 * ftot/imp.mass * pow(input_file->time_step, 2);

      // Update impurity velocity.
      imp.vel = imp.vel + ftot/imp.mass * input_file->time_step;

      // Save old positions for checking to see if crossed boundary and then
      // update impurity position.
      old_r = imp.r_pos;
      old_p = imp.p_pos;
      old_y = imp.y_pos;
      imp.r_pos += delta_rad; new_r = imp.r_pos;
      imp.p_pos += delta_pol; new_p = imp.p_pos;
      imp.y_pos += delta_par; new_y = imp.y_pos;

      //cout << "delta_rad = " << delta_rad << endl;
      //cout << "delta_par = " << delta_par << endl;
      cout << new_r << ", " << new_p << ", " << new_y << endl;

      if(input_file->debug_track){
        track_file << new_r << "," << new_p << "," << new_y << endl;
      }

      time_spent = time(nullptr) - time_start;
      ++num_iter;
    }
    while(not_hit_surface(old_r, old_p, old_y, new_r, new_p, new_y, input_file)
          && time_spent < input_file->max_time_per);

    int time_end = time(nullptr);
    int time_following = time_end - time_start;
    cout << "Time following impurity #" << imp.number << "/" << input_file->num_imps << " = "
         << time_following << " ms" << endl;

    if(input_file->debug_track){
      track_file.close();
    }

    return imp;
}
