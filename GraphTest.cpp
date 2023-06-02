#include <iostream>
#include "Graph.h"
#include <cmath>
using namespace std;


int printStats(vector<int> &epiProfile, int &epiLength){
    
    int totI = 0;// initialize int to count total infected
  
    cout << "Profile: ";// printout what is the numbers are displaying
    for (int x = 0; x < epiLength; ++x) {// go through elements of vector
        totI += epiProfile[x];// sum elements of vector
        cout << epiProfile[x] << " ";// print out contents of vector
    }
    cout << endl;// end line
    cout <<"Total Infected: " << totI << endl;// printout the total number of infected

    return 0;
}

int runEpis(int numEpis, float alpha, int numNodes, vector<int> &weights) { //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    


    Graph obj(numNodes);// Initiallize "Graph" object
    obj.fill(weights, true);
    int totI = 4;
    int avgLen = 0;
    int avgTotI = 0;
    vector<int> avgPro = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    vector<int> vect2 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    cout << "The program will run " << numEpis << " on the following graph: "<< endl;
    // Find out how to print graph!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for(int x = 0; x < numEpis; ++x){// initialize for loop and have it run for the number of epidemics
      cout << "Epidemic: " << x+1 << endl;// print what epidemic is being performed
      int epiLen = obj.SIR(0, alpha, vect2, totI);// run the SIR program for the epidemic

      avgTotI += totI;// add total infected during epidemic to total infected
      avgLen += epiLen;// add epidemic length to total epidemic length
      for(int y = 0; y < epiLen; y++){// add epidemic profile to total epidemic profile
        avgPro[y] += vect2[y];
      }

      cout << "Length: " << epiLen << endl;// print out epidemic length
      printStats(vect2, epiLen);// can pass vect2 as epiProfile since passed as reference
    }// if performed all required epiodes break while loop

    cout << "Summary:" << endl;// print out summary information for the 30 runs performed
    cout << "Average Length: " << avgLen/numEpis << endl;// print average epidemic length
    cout << "Average Total Infections: " << avgTotI/numEpis << endl;
    cout << "Average Profile: ";
    for(int x = 0; x < avgPro.size(); ++x){
      cout << round(avgPro[x]/numEpis) << " ";
    }
    cout << endl;
    return 0;
}

int main() {
  //Graph obj(4);

  //obj.print(cout);

  // Create an empty vector
  vector<int> vect;
  // Fill vector
  vect.push_back(1);
  vect.push_back(1);
  vect.push_back(2);
  vect.push_back(3);
  vect.push_back(4);
  vect.push_back(4);
  vect.push_back(1);
  vect.push_back(2);
  vect.push_back(3);
  vect.push_back(4);
  vect.push_back(1);
  vect.push_back(2);
  vect.push_back(3);
  vect.push_back(4);
  vect.push_back(0);
  //int x = 5;
  //obj.fill(vect, true);
  //obj.SIR(0, 0.3, vect, x);
  //obj.print(cout);
  runEpis(30, 0.25, 6, vect);
  return 0;
}
