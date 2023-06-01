#include <iostream>
#include "Graph.h"
using namespace std;

int main() {
  Graph obj(4);

  //obj.print(cout);

  // Create an empty vector
  vector<int> vect;
  
  vect.push_back(0);
  vect.push_back(1);
  vect.push_back(2);
  vect.push_back(3);
  vect.push_back(4);
  vect.push_back(5);
  int x = 5;

  obj.fill(vect, true);

  obj.SIR(0, 0.3, vect, x);

  obj.print(cout);

  return 0;
}