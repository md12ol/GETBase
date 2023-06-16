#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <bitset>
#include <numeric>
#include <algorithm>
#include <cassert>

#include "Graph.h"
#include "SDA.h"

// Generic Controls
int seed = 9107819;         // Random Number Seed
const bool verbose = true;  // Print Information to User?
const int numNodes = 256;

// Genetic Algorithm Controls 
int popsize = 50;
int initRunNum = 1;
int runs = 30;
int generations = 10000;                                    // Number of Mating Events
const int numReports = 100;                                 // Reporting Interval
const int reportEvery = (int) (generations / numReports);   // Report Every
int tournSize = 7;                                          // Tournament Size
double crossoverRate = 1.0;
double mutationRate = 1.0;
int maxMuts = 2;                                            // Maximum Number of Mutations

// Epidemic Controls
const int profLen = 16;         // Profile Length
char *pathToProfile;
int profileNum = 1;
int numSampEpis = 30;           // Number of Sample Epidemics
const double alpha = 0.5;       // Probability of Infection
const int minEpiLen = 3;        // Minimum Epidemic Length
const int shortEpiRetrys = 5;   // Retries for Short Epidemics

// Self-Driving Automata (SDA) Controls
int SDANumStates = 12;
const int SDAOutLen = numNodes * (numNodes - 1) / 2;
const int SDANumChars = 2;
const int SDAMaxRespLen = 2;

// Variant Controls
int initOneBits = 32;
const int maxNumVars = numNodes;
const int maxEpiLen = 128;

// Genetic Algorithm Variables
SDA *SDAPop;                // Stores the Population of SDAs
vector<double> fits;        // Fitness of each SDA
bool *dead;                 // Which SDAs Failed the Necrotic Filter (true)
int ctrlFitnessFctn;        // Program Control for Which Fitness Function to Use
double maxFit = numNodes;
double minFit = 0.0;

// Epidemic Variables
int *bestEpiLen;                    // Best Epidemic for each SDA
double targetProfile[profLen + 1];  // Profile dictionary

// Epidemic Variants Variables
int bestVarCount;                       // The Number of Variants in the Best Epidemic
int bestVarParents[maxNumVars];         // The Parents of Each Variant in the Best Epidemic
int bestVarStart[maxNumVars];           // The Variant Lengths (start, end) for Variants in the Best Epidemic
vector<int> bestVarProfs[maxNumVars];   // The Variant Profiles for the Variants in the Best Epidemic
bitset<DNALen> bestVarDNA[maxNumVars];  // The Variant DNA Strings for the Variants in the Best Epidemic
double variantProb;                     // Probability of Generating a new Variant
int minEdits;                           // Minimum Number of Edits to New Variant String
int maxEdits;                           // Maximum Number of Edits to New Variant String

// Other Variables
char outRoot[50] = "./Output/";     // The Root Output Directory
char pathToOut[200];                // Variable to Store Path to Output for Opening File Streams
const static bool diagFill = true;  // Should the SDA Fill the UTAM Diagonally
Graph network(numNodes);
vector<int> SDAOutput;
vector<int> epiProfile;

// Method Declarations
int getArgs(char *args[]);                              // Get Command Line Arguments
void createReadMe(ostream &outStrm);                    // Make readme.dat
void initAlg();                                         // Initialize the Algorithm
void cmdLineIntro(ostream &outStrm);                    // Print Intro to User
void cmdLineRun(int run, ostream &outStrm);             // Print Column Headers to User
void initPop();                                         // Initialize the Population
bool necroticFilter();
double fitnessPrep(int idx, bool final);                // Prepare to Calculate Fitness
double fitness(int idx, bool final);                    // Calculate the Fitness
void report(ostream &outStrm);
vector<double> calcStats(const vector<int> &goodIdxs,
                         bool biggerBetter);
void matingEvent();
vector<int> tournSelect(int size, bool decreasing);
bool compareFitness(int popIdx1, int popIdx2);
template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec,
                      const vector<int> &idxs,
                      const string &msg,
                      const string &sep,
                      bool newline);
void reportBest(ostream &outStrm);

int getArgs(char *args[]) {
    seed = stoi(args[1]);
    ctrlFitnessFctn = stoi(args[2]);
    variantProb = stod(args[3]);
    // Ensure that variants are not used in combination with profile matching fitness.
    assert(ctrlFitnessFctn != 1 || variantProb == 0.0);
    initRunNum = stoi(args[4]);
    runs = stoi(args[5]);
    popsize = stoi(args[6]);
    generations = stoi(args[7]);
    tournSize = stoi(args[8]);
    crossoverRate = stod(args[9]);
    mutationRate = stod(args[10]);
    maxMuts = stoi(args[11]);
    numSampEpis = stoi(args[12]);
    SDANumStates = stoi(args[13]);
    if (ctrlFitnessFctn == 1) { // Profile Matching
        pathToProfile = args[14];
        profileNum = stoi(args[15]);
    } else if (variantProb > 0.0) { // Variants
        initOneBits = stoi(args[14]);
        minEdits = stoi(args[15]);
        maxEdits = stoi(args[16]);
    }
    cout << "Arguments Captured!" << endl;
    return 0;
}

void createReadMe(ostream &outStrm) {
    outStrm << "This file contains the info about the files in this folder and the parameter settings for "
               "this experiment." << endl;
    outStrm << "Graph Evolution Tool." << endl;
    outStrm << "Fitness Function: ";
    if (ctrlFitnessFctn == 0) {
        outStrm << "Epidemic Length" << endl;
    } else if (ctrlFitnessFctn == 1) { // TODO: Not Yet Implemented.
        outStrm << "Profile Matching" << endl;
        outStrm << "Profile: ";
        for (double i: targetProfile) {
            outStrm << i << " ";
        }
        outStrm << endl;
    } else if (ctrlFitnessFctn == 2) {
        outStrm << "Epidemic Spread" << endl;
    }

    outStrm << "Epidemic Model Used: ";
    if (variantProb > 0.0) {
        outStrm << "SIR with Variants" << endl;
    } else {
        outStrm << "SIR" << endl;
    }

    outStrm << "Representation: SDAs Filling an Adjacency Matrix" << endl;
    outStrm << endl;
    outStrm << "The parameter settings are as follows: " << endl;
    outStrm << "Number of sample epidemics: " << numSampEpis << endl;
    outStrm << "Alpha: " << alpha << endl;
    outStrm << "Minimum epidemic length: " << minEpiLen << endl;
    outStrm << "Re-tries for short epidemics: " << shortEpiRetrys << endl;
    outStrm << "Runs: " << runs << endl;
    outStrm << "Mating events: " << generations << endl;
    outStrm << "Population size: " << popsize << endl;
    outStrm << "Number of vertices: " << numNodes << endl;
    outStrm << "Maximum number of mutations: " << maxMuts << endl;
    outStrm << "Tournament size: " << tournSize << endl;
    outStrm << "Number of States: " << SDANumStates << endl;
    outStrm << "Variant Probability: " << variantProb << endl;
    outStrm << "Edit Range: " << minEdits << " to " << maxEdits << endl;
    outStrm << endl;
    outStrm << "The file descriptions are as follows: " << endl;
    outStrm << "best##.dat -> the best fitness and it's associated data for one or more runs" << endl;
    outStrm << "run##.dat -> population statistics for each run" << endl;
}

void initAlg() {
    srand48(seed);

    SDAPop = new SDA[popsize];
    fits.reserve(popsize);
    bestEpiLen = new int[popsize];
    dead = new bool[popsize];
    for (int idx = 0; idx < popsize; ++idx) {
        SDAPop[idx] = SDA(SDANumStates, SDANumChars, SDAMaxRespLen, SDAOutLen);
        fits.push_back(-1);
        dead[idx] = false;
        bestEpiLen[idx] = -1;
    }
    epiProfile.reserve(maxEpiLen);
    for (int idx = 0; idx < maxEpiLen; ++idx) {
        epiProfile.push_back(-1);
    }
    SDAOutput.reserve(SDAOutLen);
    for (int i = 0; i < SDAOutLen; ++i) {
        SDAOutput.push_back(-1);
    }

    if (ctrlFitnessFctn == 1) { // TODO: Not Yet Implemented.
        fstream inp;
        char buf[20];
        inp.open(pathToProfile, ios::in);      //open input file
        for (int i = 0; i < profLen; i++) {
            targetProfile[i] = 0;  //pre-fill missing values
        }
        targetProfile[0] = 1;  //put in patient zero
        for (int i = 0; i < profLen; i++) {      //loop over input values
            inp.getline(buf, 19);  //read in the number
            targetProfile[i + 1] = strtod(buf, nullptr);    //translate the number
        }
        inp.close();
    }
}

void cmdLineIntro(ostream &outStrm) {
    outStrm << "Graph Evolution Tool." << endl;
    outStrm << "Fitness Function: ";

    if (ctrlFitnessFctn == 0) {
        outStrm << "Epidemic Length" << endl;
    } else if (ctrlFitnessFctn == 1) { // TODO: Not Yet Implemented.
        outStrm << "Profile Matching" << endl;
        outStrm << "Profile: ";
        for (double i: targetProfile) {
            outStrm << i << " ";
        }
        outStrm << endl;
    } else if (ctrlFitnessFctn == 2) {
        outStrm << "Epidemic Spread" << endl;
    }

    outStrm << "Epidemic Model Used: ";
    if (variantProb > 0.0) {
        outStrm << "SIR with Variants" << endl;
    } else {
        outStrm << "SIR" << endl;
    }
    outStrm << "Representation: SDAs Filling an Adjacency Matrix" << endl;
    outStrm << "Check readme.dat for more information about parameters and output.";
    outStrm << endl;
}

void cmdLineRun(int run, ostream &outStrm) {
    outStrm << endl << "Beginning Run " << run << "." << endl;
    outStrm << left << setw(5) << "Run";
    outStrm << left << setw(4) << "RI";
    outStrm << left << setw(10) << "Mean";
    outStrm << left << setw(12) << "95% CI";
    outStrm << left << setw(10) << "SD";
    outStrm << left << setw(12) << "Best Fit";
    outStrm << left << setw(10) << "Best Epi";
    outStrm << left << setw(8) << "# Dead";
    outStrm << endl;
    outStrm << left << setw(5) << run;
    outStrm << left << setw(4) << "0";
}

void initPop() {
    // Generate the Initial Population
    for (int idx = 0; idx < popsize; idx++) {
        dead[idx] = false;
        SDAPop[idx].fillOutput(SDAOutput);
        while (necroticFilter()) {
            SDAPop[idx].randomize();
            SDAPop[idx].fillOutput(SDAOutput);
        }
        fits[idx] = fitnessPrep(idx, false);
    }

    // Determine the Min or Max fits
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing Fitness
        minFit = fits[0];
        for (double i: fits) {
            if (i < minFit) {
                minFit = i;
            }
        }
    } else { // Minimizing Fitness
        maxFit = fits[0];
        for (double i: fits) {
            if (i > maxFit) {
                maxFit = i;
            }
        }
    }
}

bool necroticFilter() {
    int len = SDAOutLen;
    int count[2] = {0, 0};
    int bounds[2] = {1 * numNodes, 6 * numNodes};

    for (int val: SDAOutput) {
        count[val]++;
    }
    if (count[1] < bounds[0] || count[1] > bounds[1]) {
        return true; // DEAD
    } else {
        return false;
    }
}

double fitnessPrep(int idx, bool final) {
    SDAPop[idx].fillOutput(SDAOutput);
    if (necroticFilter()) {
        dead[idx] = true;
        if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) {
            return minFit;
        } else {
            return maxFit;
        }
    }

    return fitness(idx, final);
}

double fitness(int idx, bool final) {//compute the fitness
    network.fill(SDAOutput, diagFill);
    int len, totInf;   //maximum, length, and total removed
    int cnt;             //counter for tries
    double trials[numSampEpis];  //stores squared error for each trial
    double accu = 0.0;         //accumulator
    double mean = 0.0;
    int best_epi = 0;

    int numVars = 0;
    vector<int> varProfs[maxNumVars];
    vector<bitset<DNALen>> variants(maxNumVars);
    bitset<DNALen> emptyBS(0);
    for (int i = 0; i < maxNumVars; ++i) {
        variants.push_back(emptyBS);
    }
    int varParents[maxNumVars];
    int varStart[maxNumVars];

    if (ctrlFitnessFctn != 0) {
        // TODO: Not Yet Implemented.
        cout << "ERROR!  NOT IMPLEMENTED!!" << endl;
    } else { // Epidemic Length Fitness
        if (variantProb == 0.0) { // ... without variants
            for (int epi = 0; epi < (final ? 10 * numSampEpis : numSampEpis); ++epi) {
                cnt = 0;
                do {
                    totInf = 0;
                    len = network.SIR(0, alpha, epiProfile, totInf);
                    cnt += 1;
                } while (len < minEpiLen && cnt < shortEpiRetrys);
                if (final) {
                    if (len > best_epi) {
                        best_epi = len;
                        bestVarCount = 0;
                        bestVarStart[0] = 0;
                        bestVarProfs[0] = epiProfile;
                        bestVarParents[0] = -1;
                        bestVarDNA[0] = bitset<DNALen>();
                    }
                } else {
                    trials[epi] = len;
                }
            }
            if (!final) {
                int longest = 0;
                for (double trial: trials) {//loop over trials
                    if (trial > longest) {
                        longest = (int) trial;
                    }
                    accu += trial;
                }
                mean = accu / numSampEpis;
                bestEpiLen[idx] = longest;
            }
        } else { // ... with variants
            for (int epi = 0; epi < (final ? 10 * numSampEpis : numSampEpis); ++epi) {
                cnt = 0;
                do {
                    len = network.SIRwithVariants(0, alpha, variantProb, numVars, maxNumVars, maxEpiLen, varProfs,
                                                  variants, varParents, varStart, initOneBits, minEdits, maxEdits);
                    cnt += 1;
                } while (len < minEpiLen && cnt < shortEpiRetrys);
                if (final) {
                    if (len > best_epi) {
                        best_epi = len;
                        bestVarCount = numVars;
                        for (int var = 0; var <= bestVarCount; ++var) {
                            bestVarStart[var] = varStart[var];
                            bestVarProfs[var] = varProfs[var];
                            bestVarParents[var] = varParents[var];
                            bestVarDNA[var] = variants[var];
                        }
                    }
                } else {
                    trials[epi] = len;
                }
            }
            if (!final) {
                int longest = 0;
                for (double trial: trials) {//loop over trials
                    if (trial > longest) {
                        longest = (int) trial;
                    }
                    accu += trial;
                }
                mean = accu / numSampEpis;
                bestEpiLen[idx] = longest;
            }
        }
    }
    return mean;
}

void report(ostream &outStrm) {//make a statistical report
    vector<int> goodIdxs;
    int deaths = 0;

    for (int i = 0; i < popsize; i++) {
        if (!dead[i]) {
            goodIdxs.push_back(i);
        } else {
            deaths++;
        }
    }

    vector<double> stats = calcStats(goodIdxs, ctrlFitnessFctn != 1);
    double mean = stats[0];
    double stdDev = stats[1];
    double CI95 = stats[2];
    double bestVal = stats[3];
    double worstVal = stats[4];

    // Find the Best and Update the Fitness of Dead Members of the Population
    int bestIdx = 0;
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing
        for (int i = 1; i < popsize; i++) {
            if (fits[i] > fits[bestIdx]) {
                bestIdx = i;
            }
        }
        if (worstVal != minFit) {
            for (int i = 0; i < popsize; i++) {
                if (dead[i]) { // Update to Minimum Fitness
                    fits[i] = worstVal;
                }
            }
            minFit = worstVal;
        }
    } else { // Minimizing
        for (int i = 1; i < popsize; i++) {
            if (fits[i] < fits[bestIdx]) {
                bestIdx = i; //find best fitness
            }
        }
        if (bestVal != maxFit) {
            for (int i = 0; i < popsize; i++) {
                if (dead[i]) { // Update to Maximum Fitness
                    fits[i] = bestVal;
                }
            }
            maxFit = bestVal;
        }
    }

    // Print Report
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing
        outStrm << left << setw(10) << mean;
        outStrm << left << setw(12) << CI95;
        outStrm << left << setw(10) << stdDev;
        outStrm << left << setw(12) << bestVal;
        outStrm << left << setw(10) << bestEpiLen[bestIdx];
        outStrm << left << setw(8) << deaths << endl;
        if (verbose) {
            cout << left << setw(10) << mean;
            cout << left << setw(12) << CI95;
            cout << left << setw(10) << stdDev;
            cout << left << setw(12) << bestVal;
            cout << left << setw(10) << bestEpiLen[bestIdx];
            cout << left << setw(8) << deaths << endl;
        }
    } else { // Minimizing
        // TODO: Not Yet Implemented.
        outStrm << left << setw(10) << mean;
        outStrm << left << setw(12) << CI95;
        outStrm << left << setw(10) << stdDev;
        outStrm << left << setw(12) << bestVal << "\t";
        outStrm << "Dead: " << deaths << endl;
        if (verbose) {
            cout << left << setw(10) << mean;
            cout << left << setw(12) << CI95;
            cout << left << setw(10) << stdDev;
            cout << left << setw(12) << bestVal << "\t";
            cout << "Dead: " << deaths << endl;
        }
    }
}

vector<double> calcStats(const vector<int> &goodIdxs, bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);
    double worstVal = (biggerBetter ? MAXFLOAT : 0.0);

    for (int idx: goodIdxs) {
        sum += fits[idx];
        if ((biggerBetter && fits[idx] > bestVal) || (!biggerBetter && fits[idx] < bestVal)) {
            bestVal = fits[idx];
        }
        if ((biggerBetter && fits[idx] < worstVal) || (!biggerBetter && fits[idx] > worstVal)) {
            worstVal = fits[idx];
        }
    }

    double mean = sum / (double) goodIdxs.size();
    double stdDevSum = 0.0;
    for (int idx: goodIdxs) {
        stdDevSum += pow(fits[idx] - mean, 2);
    }
    double stdDev = sqrt(stdDevSum / ((double) goodIdxs.size() - 1.0));
    double CI95 = 1.96 * (stdDev / sqrt((double) goodIdxs.size()));

    return {mean, stdDev, CI95, bestVal, worstVal}; // {mean, stdDev, 95CI, best, worst}
}

void matingEvent() {
    int numMuts;
    vector<int> tournIdxs;
    // Selection
    tournIdxs = tournSelect(tournSize, ctrlFitnessFctn == 1); // TODO: Will not work for profile matching.

    // Copy the Parents -> Children
    SDAPop[tournIdxs[0]].copy(SDAPop[tournIdxs[tournSize - 2]]);
    SDAPop[tournIdxs[1]].copy(SDAPop[tournIdxs[tournSize - 1]]);

    // Crossover
    if (drand48() < crossoverRate) SDAPop[tournIdxs[0]].twoPointCrossover(SDAPop[tournIdxs[1]]);

    // Mutation
    if (drand48() < mutationRate) {
        numMuts = (int) lrand48() % (maxMuts + 1) + 1;
        SDAPop[tournIdxs[0]].mutate(numMuts);
        numMuts = (int) lrand48() % (maxMuts + 1) + 1;
        SDAPop[tournIdxs[1]].mutate(numMuts);
    }

    // Reset dead SDAs
    SDAPop[tournIdxs[0]].fillOutput(SDAOutput);
    dead[tournIdxs[0]] = necroticFilter();
    SDAPop[tournIdxs[1]].fillOutput(SDAOutput);
    dead[tournIdxs[1]] = necroticFilter();

    if (!dead[tournIdxs[0]]) {
        fits[tournIdxs[0]] = fitnessPrep(tournIdxs[0], false);
    } else {
        fits[tournIdxs[0]] = minFit;
    }

    if (!dead[tournIdxs[1]]) {
        fits[tournIdxs[1]] = fitnessPrep(tournIdxs[1], false);
    } else {
        fits[tournIdxs[1]] = minFit;
    }
}

vector<int> tournSelect(int size, bool decreasing) {
    vector<int> tournIdxs;
    int idxToAdd;

    tournIdxs.reserve(size);
    if (size == popsize) {
        for (int idx = 0; idx < size; idx++) {
            tournIdxs.push_back(idx);
        }
    } else {
        do {
            idxToAdd = (int) lrand48() % popsize;
            if (count(tournIdxs.begin(), tournIdxs.end(), idxToAdd) == 0) {
                tournIdxs.push_back(idxToAdd);
            }
        } while (tournIdxs.size() < tournSize);
    }

    sort(tournIdxs.begin(), tournIdxs.end(), compareFitness);
    if (decreasing) {
        reverse(tournIdxs.begin(), tournIdxs.end());
    }
    return tournIdxs;
}

bool compareFitness(int popIdx1, int popIdx2) {
    return fits[popIdx1] < fits[popIdx2];
}

template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs, const string &msg, const string &sep,
                      bool newline) {
    outp << msg;
    for (auto idx: idxs) {
        outp << vec[idx] << sep;
    }
    if (newline) outp << "\n";
    return 0;
}

void reportBest(ostream &outStrm) {//report the best graph
    int bestIdx;
    Graph G(numNodes);

    bestIdx = 0;
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing
        for (int i = 1; i < popsize; i++) {
            if (fits[i] > fits[bestIdx]) {
                bestIdx = i;
            }
        }
    } else { // Minimizing
        for (int i = 1; i < popsize; i++) {
            if (fits[i] < fits[bestIdx]) {
                bestIdx = i;
            }
        }
    }

    outStrm << "Best Fitness: " << fits[bestIdx] << endl;

    // Print the Profiles
    fitnessPrep(bestIdx, true);
    outStrm << "Epidemic Profile" << endl;
    if (variantProb > 0.0) {
        for (int i = 0; i <= bestVarCount; i++) {
            if (i == 0) {
                outStrm << "NA-->" << "V" << left << setw(2) << i << "\t";
            } else {
                outStrm << "V" << left << setw(2) << bestVarParents[i] << "->V" << left << setw(2) << i << "\t";
            }
            outStrm << "[" << left << setw(3) << bestVarStart[i] << "-";
            outStrm << left << setw(3) << bestVarStart[i] + bestVarProfs[i].size() << "]:\t";
            for (int j = 0; j < bestVarStart[i]; j++) {
                outStrm << "\t";
            }
            for (int j: bestVarProfs[i]) {
                outStrm << j << "\t";
            }
            outStrm << endl;
        }
        for (int i = 0; i <= bestVarCount; i++) {
            outStrm << "V" << i << "\t";
            outStrm << bestVarDNA[i] << endl;
        }
    } else {
        outStrm << left << setw(4) << "V0" << " ";
        outStrm << "[" << left << setw(2) << bestVarStart[0] << " ";
        outStrm << bestVarProfs[0].size() << "]: ";
        for (int j: bestVarProfs[0]) {
            outStrm << j << " ";
        }
        outStrm << endl;
    }

    // Write the SDA
    outStrm << "Self-Driving Automata" << endl;
    SDAPop[bestIdx].printSDA(outStrm);

    SDAPop[bestIdx].fillOutput(SDAOutput);
    G.fill(SDAOutput, diagFill);
    outStrm << "Graph" << endl;
    G.print(outStrm);
    outStrm << endl;
}