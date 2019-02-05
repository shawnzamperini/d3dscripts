/*

Will eventually fill this in with documentation. For now an important
explanation. The coordinate system is (Radial, Poloidal, Parallel) or (R,P,Y).
The origin is located at the "top" of the simulation box, which can be considered
the separatrix for now. Further implementation would maybe maybe R=0cm the
separatrix, and then have an absorbing surface kick in at around R>=7cm, i.e.
where the windowed region begins. This would require a different plasma solution
in the 0 > R > 7cm region, but this shouldn't be prohibited to implement.

*/


#include <iostream>
#include <math.h>
#include <fstream>
#include "read_input.h"
using namespace std;

const double epsilon = 1e-6;

// Simple utility function to test if zero, or close enough to it.
bool greater_than_zero(double number){
  if(number >= epsilon){
    return true;
  }
  else{
    return false;
  }
}

bool is_zero(double number){
  if(abs(number) < epsilon){
    return true;
  }
  else{
    return false;
  }
}

bool is_equal(double number1, double number2){
  if(abs(number1 - number2) < epsilon){
    return true;
  }
  else{
    return false;
  }
}

// A Cube object to represent a single 3D cell in the simulation grid.
class Cube{
  public:
    double rad_mid, pol_mid, par_mid;
    double te, ne, mach, elec;

    Cube();
};
Cube::Cube(){}


// Create a plasma grid using a given InputFile.
class PlasmaGrid{
  public:
    InputFile input_file;
    int sol_option;

    PlasmaGrid(InputFile* input_file);
    vector<vector<vector<Cube>>> construct_grid();
    void assign_background(vector<vector<vector<Cube>>> &grid);
};

// Constructor to accept an InputFile passed in.
PlasmaGrid::PlasmaGrid(InputFile* input_file_ptr){

  // Deconstruct the pointer and load in the InputFile.
  input_file = *input_file_ptr;
}

// Create a 3D grid of "Cubes", where each cube o the grid will hold its postion,
// Te, ne, etc.
vector<vector<vector<Cube>>> PlasmaGrid::construct_grid(){

  // Create the 3D grid.
  cout << "Creating grid... ";
  //Cube grid[input_file.num_rad_bins]   // i
  //         [input_file.num_pol_bins]   // j
  //         [input_file.num_par_bins];  // k

  // This is how you create a 3D array with vectors. It doesn't look very intuitive,
  // but it does the same thing as the above commented out code. Vectors are
  // better because you can easily pass them between functions.
  vector<vector<vector<Cube>>> grid(input_file.num_rad_bins,
                                    vector<vector<Cube>>(input_file.num_pol_bins,
                                    vector<Cube>(input_file.num_par_bins, Cube())));
  cout << "done." << endl;

  // Need to assign each grid a position. So let's just make a 1-1 mapping between
  // each index and its radial/poloidal/parallel location. For example, if the R
  // range is [0,0.015] and there are 100 bins:
  // R index   R location
  // 0         0/100   * 0.015
  // 1         1/100   * 0.015
  // 2         2/100   * 0.015
  // ...       ...............
  // 100       100/100 * 0.015
  // Same idea for poloidal and parallel directions, except you need to account
  // for the fact they can be negative values.
  for(int j=0; j<input_file.num_pol_bins; j++){
    for(int k=0; k<input_file.num_par_bins; k++){
      for(int i=0; i<input_file.num_rad_bins; i++){
        grid[i][j][k].rad_mid = double(i) / double(input_file.num_rad_bins - 1) * input_file.rad_fullwidth;
      }
    }
  }

  // Poloidal.
  for(int i=0; i<input_file.num_rad_bins; i++){
    for(int k=0; k<input_file.num_par_bins; k++){
      for(int j=0; j<input_file.num_pol_bins; j++){
        grid[i][j][k].pol_mid = input_file.pol_halfwidth * (2 * double(j) / (double(input_file.num_pol_bins) - 1) - 1);
      }
    }
  }

  // Parallel.
  for(int i=0; i<input_file.num_rad_bins; i++){
    for(int j=0; j<input_file.num_pol_bins; j++){
      for(int k=0; k<input_file.num_par_bins; k++){
        grid[i][j][k].par_mid = input_file.par_halfwidth * (2 * double(k) / (double(input_file.num_par_bins) - 1) - 1);
      }
    }
  }

  return grid;
}

void PlasmaGrid::assign_background(vector<vector<vector<Cube>>> &grid){

  // Te and ne are entered in distance from omp, subtract start_window from them
  // so that zero is the origin (or at least very close to it).
  //cout << "Te size: " << input_file.te_rad.size() << endl;
  for(int m=0; m < int(input_file.te_rad.size()); m++){
    input_file.te_rad[m] -= input_file.start_window;
  }
  for(int m=0; m < int(input_file.ne_rad.size()); m++){
    input_file.ne_rad[m] -= input_file.start_window;
  }
  vector <double>::iterator index_closest_to_zero;
  index_closest_to_zero = lower_bound(input_file.te_rad.begin(), input_file.te_rad.end(), 0);

  // Finally, before setting up the entire SOL, assign Te, ne to the correct
  // Cubes that represent the middle of the windowed region (i.e. the normal idea
  // where we know Teu, neu and then get everything from that). Currently assumes
  // Te, ne have the same R data. To separate the two just copy and paste and
  // swap te with ne.
  int mid_par_bin = input_file.num_par_bins / 2;  // Intentionally do integer division here.
  cout << "mid_par_bin = " << mid_par_bin << endl;

  cout << "Closest te_rad index to zero is " << (index_closest_to_zero-input_file.te_rad.begin()) << endl;
  cout << "Assigning upstream values... ";
  vector <double>::iterator closest_te_rad;
  for(int k=0; k<input_file.num_par_bins; k++){

    // If we're in the middle of the SOL/simulation box.
    if(k == mid_par_bin){
      for(int i=0; i<input_file.num_rad_bins; i++){
        for(int j=0; j<input_file.num_pol_bins; j++){

          // Find closest te_rad to this Cube.rad_mid, assign that te to the Cube.
          closest_te_rad = lower_bound(input_file.te_rad.begin(), input_file.te_rad.end(), grid[i][j][k].rad_mid);
          //cout << "Assigning " << input_file.te[(closest_te_rad -input_file.te_rad.begin())] << ", "
          //     << input_file.ne[(closest_te_rad -input_file.te_rad.begin())] << " to " << i << j << k << endl;
          grid[i][j][k].te = input_file.te[closest_te_rad - input_file.te_rad.begin()];
          grid[i][j][k].ne = input_file.ne[closest_te_rad - input_file.te_rad.begin()];
        }
      }
    }
  }
  cout << "done." << endl;

  // Now for each Cube in the grid, assign plasma parameters. Dependent on SOL option.
  switch(input_file.sol_option){

    // A sheath-limited, 1D isothermal fluid model, constant source. Very basic.
    case 0:
      cout << "SOL Option 0." << endl << "Assigning flux tube values... " << endl;
      // Equations defining this model:
      // M = L/y +/- 1/2 * sqrt((2L/y + 2)(2L/y - 2))
      // Choose - if y is positive, + if y is negative.
      //
      // n(M) = n0 / (1 + M^2)
      // n0 is density at y = 0 for this tube.
      //
      // V(M) = -kT/e * ln(1 + M^2)
      //
      // E(M) = kT/eL * M(1 + M^2)/(1 - M^2)
      double conn_length = input_file.par_halfwidth;
      double y;  // Parallel variable. -L < y < L.
      double n0, t0;
      for(int i=0; i<input_file.num_rad_bins; i++){
        n0 = grid[i][0][mid_par_bin].ne;
        t0 = grid[i][0][mid_par_bin].te;
        for(int k=0; k<input_file.num_par_bins; k++){
          y = grid[i][0][k].par_mid;
          for(int j=0; j<input_file.num_pol_bins; j++){

            // First find the Mach number.
            if(y > 0){
              grid[i][j][k].mach = conn_length / y - 0.5 * sqrt((2*conn_length / y + 2)
                                  * (2*conn_length / y - 2));
            }
            else if(y < 0){
              grid[i][j][k].mach = conn_length / y + 0.5 * sqrt((2*conn_length / y + 2)
                                  * (2*conn_length / y- 2));
            }
            else if(is_zero(y)){
              grid[i][j][k].mach = 0.0;
            }

            // Now find the density, n(M).
            grid[i][j][k].ne = n0 / (1 + pow(grid[i][j][k].mach, 2.0));

            // Apply a constant temperature (isothermal) along each flux tube.
            grid[i][j][k].te = t0;
          }
        }
      }

      // If collector probe needs it own plasma, reconstruct it here.
      if(input_file.cp_option){

        double cp_pol_min, cp_pol_max, cp_rad_min;

        // Assume probe is centered at pol = 0. This is okay since it's kinda arbitrary.
        cp_pol_min = -input_file.probe_width / 2.0;
        cp_pol_max =  input_file.probe_width / 2.0;

        // Assume probe goes from its tip all the way out (to the wall).
        cp_rad_min = input_file.probe_tip;

        // These are the halfway points between each probe face and its corresponding
        // wall. If the probe is at Y=0, then these will be the same.
        double cp_yzero_right, cp_yzero_left, cp_y, n0_right, n0_left, conn_length_right, conn_length_left;

        // Maybe these three lines are redundant? It hurts to think about for some reason.
        conn_length_right = abs(input_file.par_halfwidth - input_file.probe_par_loc) / 2.0;
        conn_length_left  = abs(input_file.probe_par_loc + input_file.par_halfwidth) / 2.0;
        cp_yzero_right = input_file.probe_par_loc + conn_length_right;
        cp_yzero_left = -input_file.par_halfwidth + conn_length_left;

        //cp_yzero_right = (input_file.par_halfwidth - input_file.probe_par_loc) / 2.0;
        //cp_yzero_left = (-input_file.par_halfwidth - input_file.probe_par_loc) / 2.0;
        cout << "probe_par_loc        = " << input_file.probe_par_loc << endl;
        cout << "conn_length          = " << conn_length << endl;
        cout << "conn_length_right    = " << conn_length_right << endl;
        cout << "conn_length_left     = " << conn_length_left << endl;
        cout << "cp_yzero_right = " << cp_yzero_right << endl;
        cout << "cp_yzero_left  = " << cp_yzero_left  << endl;

        // Before we begin we need to define new mid_par_bin so that we can set
        // the appropriate n0's for the left/right. I.e. what parallel bin number
        // is closest to cp_yzero_right and left?
        int mid_par_bin_right, mid_par_bin_left, probe_bin;
        double tmp_dist;
        double old_tmp_left  = 9999999;
        double old_tmp_right = 9999999;
        double old_tmp_probe = 9999999;
        for(int k=0; k<input_file.num_par_bins; k++){
          tmp_dist = abs(grid[0][0][k].par_mid - cp_yzero_right);
          //cout << "tmp_dist right = " << tmp_dist << endl;
          if(tmp_dist < old_tmp_right){
            mid_par_bin_right = k;
          }
          old_tmp_right = tmp_dist;
          tmp_dist = abs(grid[0][0][k].par_mid - cp_yzero_left);
          if(tmp_dist < old_tmp_left){
            mid_par_bin_left = k;
          }
          old_tmp_left = tmp_dist;
          tmp_dist = abs(grid[0][0][k].par_mid - input_file.probe_par_loc);
          if(tmp_dist < old_tmp_probe){
            probe_bin = k;
          }
          old_tmp_probe = tmp_dist;
        }

        cout << "mid_par_bin_right = " << mid_par_bin_right << endl;
        cout << "mid_par_bin_left  = " << mid_par_bin_left  << endl;
        cout << "probe_bin         = " << probe_bin         << endl;

        for(int i=0; i<input_file.num_rad_bins; i++){

          // If we're in the radial range of the probe.
          if(grid[i][0][0].rad_mid > cp_rad_min){
            n0_right = grid[i][0][mid_par_bin_right].ne;
            n0_left  = grid[i][0][mid_par_bin_left].ne;
            t0 = grid[i][0][mid_par_bin].te;
            for(int j=0; j<input_file.num_pol_bins; j++){

              // If we're in the poloidal range of the probe.
              if(grid[0][j][0].pol_mid > cp_pol_min && grid[0][j][0].pol_mid < cp_pol_max){
                for(int k=0; k<input_file.num_par_bins; k++){

                  // Define a new y so that instead of the parallel location in the entire
                  // simulation volume, it's just the parallel location between the
                  // probe face and the wall.
                  y = grid[i][0][k].par_mid;
                  if(y < input_file.probe_par_loc){
                    //cp_y = conn_length_left * ((y + conn_length) / (conn_length_left + conn_length) - 1.0);
                    cp_y = conn_length_left * ((2 * (y + conn_length)) / (input_file.probe_par_loc + conn_length) - 1);
                  }
                  else if (y > input_file.probe_par_loc){
                    //cp_y = -conn_length_right * ((y - conn_length) / (conn_length_right - conn_length) - 1.0);
                    cp_y = conn_length_right * ((2 * (y - input_file.probe_par_loc)) / (conn_length - input_file.probe_par_loc) - 1);
                  }
                  //cout << "cp_y = " << cp_y << endl;

                  // Mach number. Middle of left or right probe tube.
                  if(is_zero(cp_y)){
                    //cout << "y    = " << y << endl;
                    //cout << "cp_y = " << cp_y << endl;
                    //cout << "Middle = ";
                    grid[i][j][k].mach = 0.0;
                    //cout << grid[i][j][k].mach << endl;
                  }

                  // Mach number. Right half of right side of probe tube.
                  else if(cp_y > 0.0 && y > input_file.probe_par_loc){
                    //cout << "y    = " << y << endl;
                    //cout << "cp_y = " << cp_y << endl;
                    //cout << "Right side right half = ";
                    grid[i][j][k].mach = conn_length_right / cp_y - 0.5 * sqrt((2*conn_length_right / cp_y + 2)
                                        * (2*conn_length_right / cp_y - 2));
                    //cout << grid[i][j][k].mach << endl;
                  }

                  // Mach number. Left half of right side of probe tube.
                  else if(cp_y < 0.0 && y > input_file.probe_par_loc){
                    //cout << "y    = " << y << endl;
                    //cout << "cp_y = " << cp_y << endl;
                    //cout << "Right side left half = ";
                    grid[i][j][k].mach = conn_length_right / cp_y + 0.5 * sqrt((2*conn_length_right / cp_y + 2)
                                        * (2*conn_length_right / cp_y - 2));
                    //cout << grid[i][j][k].mach << endl;
                  }

                  // Mach number. Right half of left probe tube.
                  else if(cp_y > 0.0 && y < input_file.probe_par_loc){
                    //cout << "Left side right half" << endl;
                    grid[i][j][k].mach = conn_length_left / cp_y - 0.5 * sqrt((2*conn_length_left / cp_y + 2)
                                        * (2*conn_length_left / cp_y - 2));
                  }

                  // Mach number. Left half of left probe tube.
                  else if(cp_y < 0.0 && y < input_file.probe_par_loc){
                    //cout << "Left side left half" << endl;
                    grid[i][j][k].mach = conn_length_left / cp_y + 0.5 * sqrt((2*conn_length_left / cp_y + 2)
                                        * (2*conn_length_left / cp_y- 2));
                  }

                  // Now find the density, n(M), using the right n0 value.
                  if(y > input_file.probe_par_loc){
                    //cout << "Right density" << endl;
                    //cout << "Old = " << grid[i][j][k].ne << ". New = ";
                    grid[i][j][k].ne = n0_right / (1 + pow(grid[i][j][k].mach, 2.0));
                    //cout << grid[i][j][k].ne << endl;
                  }
                  else if(y < input_file.probe_par_loc){
                    grid[i][j][k].ne = n0_left / (1 + pow(grid[i][j][k].mach, 2.0));
                  }
                  // Apply a constant temperature (isothermal) along each flux tube.
                  grid[i][j][k].te = t0;

                  // Finally get rid of the discontinuity at Y=0 (probably only
                  // necessary to make plots pretty, but not sure. Can't hurt to
                  // do this though) by using the average between the two points
                  // closest to each probe face.
                  grid[i][j][probe_bin].ne = 0.5 * (grid[i][j][probe_bin-1].ne
                                                    + grid[i][j][probe_bin+1].ne);
                }
              }
            }
          }
        }
      }

    cout << "done." << endl;

    // Print out the background near P=0 to a file for plotting in Python.
    ofstream grid_file;
    grid_file.open("grid.txt");
    for(int k=0; k<input_file.num_par_bins; k++){
      grid_file << grid[0][input_file.num_pol_bins/2][k].par_mid << " ";
    }
    grid_file << "\n";
    for(int i=0; i<input_file.num_rad_bins; i++){
      grid_file << grid[i][input_file.num_pol_bins/2][0].rad_mid << " ";
      for(int k=0; k<input_file.num_par_bins; k++){
        grid_file << grid[i][input_file.num_pol_bins/2][k].ne << " ";
      }
      grid_file << "\n";
    }
    grid_file.close();

    break;
  }
}
