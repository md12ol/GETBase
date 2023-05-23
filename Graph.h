#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

class Graph {
public:
    Graph();
    explicit Graph(int nn);

    vector<int> fill(string fname);
    int fill(vector<int> &weights, bool diag);
    int SIR(double alpha, int p0, vector<int> &profile, int &ttl);
    void print(ostream &out);
    vector<int> weightHist();
    int hammy_distance(Graph &other, int penalty);
    static int quadForm(int A, int B, int C);

protected:
    static int infect(int nin, double alpha);
    int numNodes;
    int numEdges;
    int totWeight;
    int maxWeight;
    vector<int> state;
    vector<int> nin;
    vector<vector<int> > adjM;
};


#endif
