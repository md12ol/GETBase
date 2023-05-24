#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <bitset>
#include <random>
#include <algorithm>

// Variant Controls
const int DNALen = 128;

using namespace std;

class Graph {
public:
    Graph();
    explicit Graph(int nn);

    vector<int> fill(const string filename);
    int fill(const vector<int> &weights, bool diag);
    int SIR(int p0, double alpha, vector<int> &epiProfile, int &totInf);
    int print(ostream &out);
    vector<int> weightHist();
    int SIRwithVariants(int p0, double alpha, double varProb, int &varCnt, int maxVars, int maxLen, vector<int> varProfs[],
                        vector<bitset<DNALen>> &variants, int varParents[], int varStart[], int initBits,
                        int minEdits, int maxEdits);

protected:
    int numNodes;
    int numEdges;
    int totWeight;
    int maxWeight;

    vector<int> state;
    vector<int> numInfNeighs;
    vector<vector<int>> potStrains;
    vector<vector<int>> adjM;
    vector<bitset<DNALen>> immunity;

    static int quadForm(int A, int B, int C);
    static bool infect(int numInfNeighs, double alpha);
    static int variantInfect(bitset<DNALen> &immunity, double alpha, vector<int> &potStrains,
                             vector<bitset<DNALen>> &allStrains, int maxVars);
    static bool compareSeverity(pair<double, int> severity1, pair<double, int> severity2);
    static int newVariant(bitset<DNALen> &origVar, bitset<DNALen> &newVar, vector<int> &rndIdxVec, int minEdits,
                          int maxEdits);
};


#endif
