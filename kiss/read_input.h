#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
using namespace std;

#ifndef INPUT_FILE
#define INPUT_FILE
// Function to read in values from the input file, as doubles, ints or vectors.
//void read_input_file(string input_filename);

class InputFile{
  public:
    int num_pol_bins, num_rad_bins, num_par_bins, sol_option, cp_option, inj_opt,
        num_imps;
    double probe_width, probe_tip, time_step, par_halfwidth, rad_fullwidth,
           pol_halfwidth, start_window, probe_par_loc, launch_rmin, launch_rmax,
           launch_pmin, launch_pmax, launch_ymin, launch_ymax, imp_amu;
    vector <double> te_rad, te, ne_rad, ne;

    InputFile();
    InputFile(string input_filename);
    double read_double(ifstream&, const string);
    void read_2d_vectors(ifstream&, const string, vector <double>&, vector <double>&);
    void read_input_file(string);
};

#endif
