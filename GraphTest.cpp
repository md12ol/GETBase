#include <iostream>
#include "Graph.h"
using namespace std;

int main() {
  Graph obj(5);

  vector<int> test1;

  for(int x = 0; x < 2; x++){
    test1.push_back(x);
  }
  
  for (auto i = test1.begin(); i != test1.end(); ++i){
    cout << *i << " ";
  }
  
  return 0;
}