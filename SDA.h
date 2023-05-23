#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;

class SDA {
public:
    SDA();           //creates an unallocated bitspray
    explicit SDA(int states, int numChars, int len);      //create a bitspray with buffer S states
    SDA(SDA &other);  //copy constructor
    ~SDA();                //destructor

    int create();
    int randomize();
    int copy(SDA &other);
    int print();
    int print(ostream &aus);
    static int destroy();
    int twoPtCrossover(SDA &other);
    int oneStateCrossover(SDA &other);
    int mutate(int numMuts);
    int getBitsVec(int len, vector<int> &rtn);
    int printBitsVec(int len, ostream &aus);

private:
    int initInput;
    int numStates;
    int initState;
    int curState;
    int numChars;
    double zeroProb = 0.95;
    vector<int> buf;
    vector<vector<int> > transitions;
    vector<vector<vector<int> > > responses;
};
