#include <iostream>
#include "SDA.h"
using namespace std;

int checker(int numVals, SDA &check){
    vector<int> vec(numVals);
    check.getBitsVec(numVals, vec);

    if(vec.size() <= numVals){// if output form SDA is less than or equal to desired length
        for(int x : vec){// go through vector
            switch(vec[x]){// create switch for different cases
                case 0:// if reading a  zero
                cout << "A " << endl;// print an "A"
                break;
                case 1:// if reading a one
                cout << "B " << endl;// Print a "B"
                break;
            }
        }
    }else cout << "Error output from SDA is to great" << endl;// if ouput from SDA is greater than desired length print error
    cout << endl;
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
