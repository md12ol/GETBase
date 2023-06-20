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
    int SIRwithVariants(int p0, double varAlphas[], bool coupled, double newVarProb, int &varCnt, int maxVars, int maxLen,
                               vector<int> varProfs[], vector<bitset<DNALen>> &varDNAs, int varParents[],
                               int varStarts[], int varInfSeverity[], int initBits, int minEdits, int maxEdits,
                               double alphaDelta, int &totInf);

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
    int variantInfect(bitset<DNALen> &immStr, double varAlphas[], vector<int> &potVars,
                             vector<bitset<DNALen>> &varStrs, int varInfSeverity[], int maxVars, bool coupled) const;
    static bool compareSeverity(pair<double, int> severity1, pair<double, int> severity2);
    static int newVariant(bitset<DNALen> &origVar, const double &origVarAlpha, bitset<DNALen> &newVar, double &newVarAlpha,
               vector<int> &rndIdxVec, int minEdits, int maxEdits, double alphaDelta, bool coupled);
};

#endif // GRAPH_H
