#include <iostream>
#include "SDA.h"
using namespace std;

int checker(int numVals, SDA &check){
    vector<int> vec(numVals);
    check.getBitsVec(numVals, vec);

    cout << endl << "Print resulting SDA bitvector in alphabet form where 0 = A and 1 = B: " << endl;

    if(vec.size() <= numVals){// if output form SDA is less than or equal to desired length
        for(int x : vec){// go through vector
            switch(x){// create switch for different cases
                case 0:// if reading a  zero
                cout << "A ";// print an "A"
                break;
                case 1:// if reading a one
                cout << "B ";// Print a "B"
                break;
            }
        }
    }else cout << "Error output from SDA is too large" << endl;// if ouput from SDA is greater than desired length print error
    cout << endl;
    return 0;
}

int main() {
    SDA obj(3,2,25);// Initiallize SDA
    obj.create();// Create SDA
    cout << "Here is the SDA generated: " << endl;
    obj.print();// print the SDA that was created

    cout << endl << "Print resulting Bit Vector: " << endl;
    obj.printBitsVec(25, cout);// Call bit vector printing method

    checker(25, obj);// call cjecker method and pass the desired length of SDA's outpur and pass the SDA


    SDA obj2(3, 2, 13);// Create the two SDAs for the test
    obj2.create();
    SDA obj3(3, 2, 24);
    obj3.create();

    cout << "Here is the first SDA created: " << endl;// print both of the SDAs to the terminal
    obj2.print();
    cout << endl << "Here is the second SDA created: " << endl;
    obj3.print();

    obj2.mutate(5);// call the mutation function for both SDAs
    obj3.mutate(3);


    cout << endl << "Here is the first SDA again but mutated: " <<endl;// print both of the mutated SDAs
    obj2.print();
    cout << endl << "Here is the second SDA again but mutated: " << endl;
    obj3.print();

    obj2.twoPtCrossover(obj2);// perform two-point crossover using the mutated SDAs

    cout << endl << "Here is the first SDA after the crossover: " << endl;
    obj2.print();
    cout << endl << "Here is the second SDA after the crossover: " << endl;
    obj3.print();

    return 0;
}
