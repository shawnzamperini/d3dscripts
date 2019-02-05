#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
using namespace std;


class InputFile{
  public:
    int num_pol_bins, num_rad_bins, num_par_bins, sol_option, cp_option;
    double probe_width, probe_tip, time_step, par_halfwidth, rad_fullwidth,
           pol_halfwidth, start_window, probe_par_loc;
    vector <double> te_rad, te, ne_rad, ne;

    InputFile();
    InputFile(string input_filename);
    double read_double(ifstream&, const string);
    void read_2d_vectors(ifstream&, const string, vector <double>&, vector <double>&);
    void read_input_file(string );
};

InputFile::InputFile(){
  //cout << "No input filename." << endl;
}
InputFile::InputFile(string input_filename){
  read_input_file(input_filename);
}
double InputFile::read_double(ifstream& input_file, const string tag){
  string line;
  double tag_value = 0;
  int bar_count;

  // First go to start of file...
  input_file.clear();
  input_file.seekg(0, ios::beg);

  // ... scan through each line until you come across the tag...
  while(getline(input_file, line)){

    // line.find will return string::npos if no match, so != indicates a match
    if (line.find(tag, 0) != string::npos){
      //cout << line << endl;

      // Get the substring at the end, which will be everything after the third bar.
      bar_count = 0;
      for (int i=0; i<int(line.size()); i++){
        if (line[i] == '|'){
          if (bar_count == 2){
            line = line.substr(i, string::npos);
            //cout << line << endl;
            break;  // We got the number (still surrounded by bars), so break out.
          }
          else{
            bar_count++;  // Move to the next bar.
          }
        }
      }
      break;  // We got the number (still surrounded by bars), so break out.
    }
  }

  // Currently we have, for example, "| 0.015 |", so let's trim off the two bars
  // and cast to a double.
  line.erase(remove(line.begin(), line.end(), '|'), line.end());
  //cout << line << endl;
  tag_value = stod(line);
  return tag_value;
}
void InputFile::read_2d_vectors(ifstream& input_file, const string tag,
                                vector <double> &vec_x, vector <double> &vec_y){

  int lines_to_read;
  double x, y;

  // Get number of lines to read.
  lines_to_read = read_double(input_file, tag);

  // At this point, the "cursor"? will be at the start of the array we want,
  // so go ahead and read the next lines_to_read lines into two vectors.
  for (int i=0; i<lines_to_read; i++){
    input_file >> x >> y;
    //cout << x << " " << y << endl;
    vec_x.push_back(x);
    vec_y.push_back(y);
  }

//  for (int i=0; i<vec_x.size(); i++){
//    cout << vec_x[i] << ", " << vec_y[i] << endl;
//  }
}
void InputFile::read_input_file(string input_filename){

  // Define input file as input only since we aren't writing to it.
  ifstream input_file (input_filename);

  // Only perform operations if the file is successfully opened.
  if (input_file.is_open()){
    cout << "Opened input file: " << input_filename << endl;

    // Find desired tag names in the input file, starting from the beginning.
    probe_width   =     read_double(input_file, "*T000");
    probe_tip     =     read_double(input_file, "*T001");
    time_step     =     read_double(input_file, "*T002");
    par_halfwidth =     read_double(input_file, "*T003");
    rad_fullwidth =     read_double(input_file, "*T004");
    pol_halfwidth =     read_double(input_file, "*T005");
    num_pol_bins  = int(read_double(input_file, "*T006"));
    num_rad_bins  = int(read_double(input_file, "*T007"));
    num_par_bins  = int(read_double(input_file, "*T008"));
    sol_option    = int(read_double(input_file, "*T009"));
    start_window  =     read_double(input_file, "*T010");
    cp_option     = int(read_double(input_file, "*T011"));
    probe_par_loc =     read_double(input_file, "*T012");
    read_2d_vectors(input_file, "*A000", te_rad, te);
    read_2d_vectors(input_file, "*A001", ne_rad, ne);

    // Make sure the pol and par bins are odd numbers so that there can be a
    // middle bin for P=0 and Y=0.
    if (num_pol_bins % 2 == 0){
      cout << "Using " << num_pol_bins + 1 << " poloidal bins instead of "
           << num_pol_bins << "." << endl;
      num_pol_bins++;
    }
    if (num_par_bins % 2 == 0){
      cout << "Using " << num_par_bins + 1 << " parallel bins instead of "
           << num_par_bins << "." << endl;
      num_par_bins++;
    }
  }

  // Print out error if error opening file.
  else{
    cout << "Error opening input file." << endl;
  }
}

/*
// Function to read a double from input file. It will return to the beginning
// of the file, and search for the tag, then get the value from the line.
double read_double(ifstream& input_file, const string tag){

  string line;
  double tag_value = 0;
  int bar_count;

  // First go to start of file...
  input_file.clear();
  input_file.seekg(0, ios::beg);

  // ... scan through each line until you come across the tag...
  while(getline(input_file, line)){

    // line.find will return string::npos if no match, so != indicates a match
    if (line.find(tag, 0) != string::npos){
      //cout << line << endl;

      // Get the substring at the end, which will be everything after the third bar.
      bar_count = 0;
      for (int i=0; i<line.size(); i++){
        if (line[i] == '|'){
          if (bar_count == 2){
            line = line.substr(i, string::npos);
            //cout << line << endl;
            break;  // We got the number (still surrounded by bars), so break out.
          }
          else{
            bar_count++;  // Move to the next bar.
          }
        }
      }
      break;  // We got the number (still surrounded by bars), so break out.
    }
  }

  // Currently we have, for example, "| 0.015 |", so let's trim off the two bars
  // and cast to a double.
  line.erase(remove(line.begin(), line.end(), '|'), line.end());
  //cout << line << endl;
  tag_value = stod(line);
  return tag_value;
}

// This function will read in a 2D array from the input file. It actually uses
// the above read_double to get the number of lines to be read for the array.
void read_2d_vectors(ifstream& input_file, const string tag,
                      vector <double> vec_x, vector <double> vec_y){

  int lines_to_read;
  double x, y;

  // Get number of lines to read.
  lines_to_read = read_double(input_file, tag);

  // At this point, the "cursor"? will be at the start of the array we want,
  // so go ahead and read the next lines_to_read lines into two vectors.
  for (int i=0; i<lines_to_read; i++){
    input_file >> x >> y;
    vec_x.push_back(x);
    vec_y.push_back(y);
  }

//  for (int i=0; i<vec_x.size(); i++){
//    cout << vec_x[i] << ", " << vec_y[i] << endl;
//  }
}

// Function to read in the input file. The input file is a format of columns
// of Tag, Name, Value, each separated by bars (|). Each input variable is
// uniquely located by its tag, and depending on its type, its value is loaded
// using one of the above functions.
void read_input_file(string input_filename){

  // Define variables to be used.
  int num_pol_bins, num_rad_bins, num_par_bins;
  double probe_width, probe_tip, time_step, par_halfwidth, rad_halfwidth,
         pol_halfwidth;
  vector <double> te_rad, te, ne_rad, ne;

  // Define input file as input only since we aren't writing to it.
  ifstream input_file (input_filename);

  // Only perform operations if the file is successfully opened.
  if (input_file.is_open()){
    cout << "File opened." << endl;

    // Find desired tag names in the input file, starting from the beginning.
    probe_width   = read_double(input_file, "*T000");
    probe_tip     = read_double(input_file, "*T001");
    time_step     = read_double(input_file, "*T002");
    par_halfwidth = read_double(input_file, "*T003");
    rad_halfwidth = read_double(input_file, "*T004");
    pol_halfwidth = read_double(input_file, "*T005");
    num_pol_bins  = read_double(input_file, "*T006");
    num_rad_bins  = read_double(input_file, "*T007");
    num_par_bins  = read_double(input_file, "*T008");
    read_2d_vectors(input_file, "*A000", te_rad, te);
    read_2d_vectors(input_file, "*A001", ne_rad, ne);
  }

  // Print out error if error opening file.
  else{
    cout << "Error opening input file." << endl;
  }
}
*/
