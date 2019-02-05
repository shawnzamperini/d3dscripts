#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


int main(){
  vector <int> my_vec;
  vector<int>::iterator low;

  for(int i=1; i<=5; i++){
    my_vec.push_back(i);
  }

  for(int i=0; i<5; i++){
    my_vec[i] -= 3;
    cout << my_vec[i] << endl;
  }

  cout << my_vec.size() << endl;

  low = lower_bound(my_vec.begin(), my_vec.end(), 0);
  cout << "Low is " << (low-my_vec.begin()) << endl;
}
