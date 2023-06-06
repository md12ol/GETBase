#include <iostream>
#include "SDA.h"
using namespace std;

int checker(int numVals, SDA &check){
    vector<int> vec(numVals);
    check.getBitsVec(numVals, vec);

    return 0;
}

int main() {
    SDA obj(3,2,25);// Initiallize SDA
    obj.create();// Create SDA
    cout << "Here is the SDA generated: " << endl;
    obj.print(cout);// print the SDA that was created

    cout << endl << "Print resulting Bit Vector: " << endl;
    obj.printBitsVec(25, cout);// Call bit vector printing method

    checker(25, obj);// call cjecker method and pass the desired length of SDA's outpur and pass the SDA

    return 0;
}
