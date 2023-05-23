#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <bitset>
#include <numeric>
#include <algorithm>

#include "Graph.h"
#include "SDA.h"

// Generic Controls
int seed = 9107819;                                // Random Number Seed
#define verbose true                                // Print Information ot User?
#define RIs 100                                     // Reporting Interval
#define RE ((long)mevs/RIs)                         // Report Every
#define verts 256

// Genetic Algorithm Controls
#define popsize 50
#define runs 30
#define mevs 10000                                  // Number of Mating Events
#define tsize 7                                     // Tournament Size
#define maxMuts 2                                       // Maximum Number of Mutations

// Epidemic Controls
#define PL 16                                       // Profile Length
#define NSE 30                                      // Number of Sample Epidemics
#define alpha 0.5                                   // Probability of Infection
#define mepl 3                                      // Minimum Epidemic Length
#define rse 5                                       // Retries for Short Epidemics

// Self-Driving Automata (SDA) Controls
#define states 12
#define stringLen verts*(verts-1)/2 + 2

#define DNALen (int)128
#define initOneBits (int)32
#define maxNumVars (int)64
#define maxEpiLen (int)128

// Genetic Algorithm Variables
SDA *SDAPop;                               // Stores the Population of SDAs
vector<double> fit(popsize);                                // Fitness of each SDA
bool dead[popsize];                                 // Which SDAs Failed the Necrotic Filter (true)
int ctrlFitnessFctn;                                // Program Control for Which Fitness Function to Use
double maxFit = verts;
double minFit = 0.0;

// SDA Variables
int SDAOutput[stringLen];                           // Stores the Output from an SDA

// Epidemic Variables
int bestEpiLen[popsize];                               // Best Epidemic for each SDA
double PD[PL + 1];                                  // Profile dictionary
bool ctrlProfileMatching;                           // Program Control for Doing Profile Matching
bool ctrlVariants;                                  // Program Control for Having Epidemic Variants
bool ctrlEpiSpread;                                 // Program Control for Using Epidemic Spread Fitness

// Epidemic Variants Variables
int bestVarCount;                                   // The Number of Variants in the Best Epidemic
int bestVarParents[maxNumVars];                     // The Parents of Each Variant in the Best Epidemic
pair<int, int> bestVarLens[maxNumVars];             // The Variant Lengths (start, end) for Variants in the Best Epidemic
vector<int> bestVarProfs[maxNumVars];               // The Variant Profiles for the Variants in the Best Epidemic
bitset<DNALen> bestVarDNA[maxNumVars];              // The Variant DNA Strings for the Variants in the Best Epidemic
double variantProb;                                 // Probability of Generating a new Variant
int minEdits;                                       // Minimum Number of Edits to New Variant String
int maxEdits;                                       // Maximum Number of Edits to New Variant String

// Other Variables
char outRoot[20] = "./Output/";                // The Root Output Directory
char pathToOut[150];                                // Variable to Store Path to Output for Opening File Streams
static bool diagFill = true;                        // Should the SDA Fill the UTAM Diagonally
Graph G(verts);
vector<int> vals;
vector<int> profile;


// Method Declarations
int getArgs(char *args[]);                          // Get Command Line Arguments
void cmdLineIntro(ostream &aus);                    // Print Intro to User
void createReadMe(ostream &aus);                    // Make readme.dat
void cmdLineRun(int run, ostream &aus);             // Print Column Headers to User
void initAlg(const char *pLoc);                     // Initialize the Algorithm
void initPop();                                     // Initialize the Population
void matingEvent();
double fitnessPrep(int idx, SDA &A, bool final);    // Prepare to Calculate Fitness
double fitness(int idx, bool final);    // Calculate the Fitness
void report(ostream &aus);
void reportBest(ostream &aus);

void createReadMe(ostream &aus) {
    aus << "This file contains the info about the files in this folder and the parameter settings for this experiment."
        << endl;
    aus << "Graph Evolution Tool." << endl;
    aus << "Fitness Function: ";
    if (ctrlProfileMatching) { // TODO: Not Yet Implemented.
        aus << "Profile Matching" << endl;
        aus << "Profile: ";
        for (double i: PD) {
            aus << i << " ";
        }
        aus << endl;
    } else if (ctrlEpiSpread) {
        aus << "Epidemic Spread" << endl;
    } else {
        aus << "Epidemic Length" << endl;
    }

    aus << "Epidemic Model Used: ";
    if (ctrlVariants) {
        aus << "SIR with Variants" << endl;
    } else {
        aus << "SIR" << endl;
    }

    aus << "Representation: SDAs Filling an Adjacency Matrix" << endl;
    aus << endl;
    aus << "The parameter settings are as follows: " << endl;
    aus << "Number of sample epidemics: " << NSE << endl;
    aus << "Alpha: " << alpha << endl;
    aus << "Minimum epidemic length: " << mepl << endl;
    aus << "Re-tries for short epidemics: " << rse << endl;
    aus << "Runs: " << runs << endl;
    aus << "Mating events: " << mevs << endl;
    aus << "Population size: " << popsize << endl;
    aus << "Number of vertices: " << verts << endl;
    aus << "Maximum number of mutations: " << maxMuts << endl;
    aus << "Tournament size: " << tsize << endl;
    aus << "Number of States: " << states << endl;
    aus << "Variant Probability: " << variantProb << endl;
    aus << "Edit Range: " << minEdits << " to " << maxEdits << endl;
    aus << endl;
    aus << "The file descriptions are as follows: " << endl;
    aus << "best.dat -> the best fitness and it's associated data for each run" << endl;
    aus << "run##.dat -> population statistics for each run" << endl;
}

void cmdLineIntro(ostream &aus) {
    aus << "Graph Evolution Tool." << endl;
    aus << "Fitness Function: ";
    if (ctrlProfileMatching) { // TODO: Not Yet Implemented.
        aus << "Profile Matching" << endl;
        aus << "Profile: ";
        for (double i: PD) {
            aus << i << " ";
        }
        aus << endl;
    } else if (ctrlEpiSpread) {
        aus << "Epidemic Spread" << endl;
    } else {
        aus << "Epidemic Length" << endl;
    }

    aus << "Epidemic Model Used: ";
    if (ctrlVariants) {
        aus << "SIR with Variants" << endl;
    } else {
        aus << "SIR" << endl;
    }
    aus << "Representation: SDAs Filling an Adjacency Matrix" << endl;
    aus << "Check readme.dat for more information about parameters and output.";
    aus << endl;
}

void cmdLineRun(int run, ostream &aus) {
    aus << endl << "Beginning Run " << run << " of " << runs << endl;
    aus << left << setw(5) << "Run";
    aus << left << setw(4) << "RI";
    aus << left << setw(10) << "Mean";
    aus << left << setw(12) << "95% CI";
    aus << left << setw(10) << "SD";
    aus << left << setw(12) << "Best Fit";
    aus << left << setw(10) << "Best Epi";
    aus << left << setw(8) << "# Dead";
    aus << endl;
    aus << left << setw(5) << run;
    aus << left << setw(4) << "0";
}

void initAlg(const char *pLoc) {
    srand48(seed);

    SDAPop = new SDA[popsize];
    profile.reserve(verts);
    for (int idx = 0; idx < popsize; ++idx) {
        SDAPop[idx] = SDA(states, 2, stringLen - 2); // TODO: Make a parameter.
        profile.push_back(-1);
    }
    vals.reserve(stringLen - 2);
    for (int i = 0; i < stringLen - 2; ++i) {
        vals.push_back(-1);
    }

    if (ctrlProfileMatching) { // TODO: Not Yet Implemented.
        fstream inp;
        char buf[20];
        inp.open(pLoc, ios::in);      //open input file
        for (int i = 0; i < PL; i++) {
            PD[i] = 0;  //pre-fill missing values
        }
        PD[0] = 1;  //put in patient zero
        for (int i = 0; i < PL; i++) {      //loop over input values
            inp.getline(buf, 19);  //read in the number
            PD[i + 1] = strtod(buf, nullptr);    //translate the number
        }
        inp.close();
    }
}

int getArgs(char *args[]) {
    string arg = args[1]; // ctrlVariants
    try {
        size_t pos;
        ctrlVariants = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[2]; // variantProb
    try {
        variantProb = std::stod(arg);
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[3]; // minEdits
    try {
        size_t pos;
        minEdits = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[4]; // maxEdits
    try {
        size_t pos;
        maxEdits = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[5]; // ctrlEpiSpread
    try {
        size_t pos;
        ctrlEpiSpread = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[6]; // ctrlProfileMatching
    try {
        size_t pos;
        ctrlProfileMatching = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    cout << "Arguments Captured!" << endl;
    return 0;
}

bool necroticFilter() {
    int len = stringLen - 2;
    int count[2] = {0, 0};
    int bounds[2] = {1 * verts, 6 * verts};

    for (int i = 0; i < len; i++) {
        count[vals[i]]++;
    }
    if (count[1] < bounds[0] || count[1] > bounds[1]) {
        return true; // DEAD
    } else {
        return false;
    }
}

double fitness(int idx, bool final) {//compute the fitness
    G.fill(vals, diagFill);
    int len, ttl;   //maximum, length, and total removed
    int cnt;             //counter for tries
    double trials[NSE];  //stores squared error for each trial
    double accu = 0.0;         //accumulator
    int best_epi = 0;

    if (ctrlFitnessFctn != 0 || ctrlVariants) {
        // TODO: Not Yet Implemented.
        cout << "ERROR!  NOT IMPLEMENTED!!" << endl;
    } else {
        for (int epi = 0; epi < (final ? 10 * NSE : NSE); ++epi) {
            cnt = 0;
            do {
                ttl = 0;
                len = G.SIR(alpha, 0, profile, ttl);
                cnt += 1;
            } while (len < mepl && cnt < rse);
            if (final) {
                if (len > best_epi) {
                    best_epi = len;
                    bestVarCount = 0;
                    bestVarLens[0].first = 0;
                    bestVarLens[0].second = len;
                    bestVarProfs[0] = profile;
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
            accu = accu / NSE;
            bestEpiLen[idx] = longest;
        }
    }
    return accu;
}

double fitnessPrep(int idx, SDA &A, bool final) {
    A.getBitsVec(stringLen - 2, vals);
    if (necroticFilter()) {
        dead[idx] = true;
        if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) {
            return minFit;
        } else {
            return maxFit;
        }
    }

    double fi = fitness(idx, final);
    return fi;
}

void initPop() {
    // Generate the Initial Population
    int size = stringLen - 2;

    for (int i = 0; i < popsize; i++) {
        dead[i] = false;
        SDAPop[i].create();
        SDAPop[i].getBitsVec(size, vals);
        while (necroticFilter()){
            SDAPop[i].randomize();
            SDAPop[i].getBitsVec(size, vals);
        }
        fit[i] = fitnessPrep(i, SDAPop[i], false);
    }

    // Determine the Min or Max fit
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing Fitness
        minFit = fit[0];
        for (double i: fit) {
            if (i < minFit) {
                minFit = i;
            }
        }
    } else { // Minimizing Fitness
        maxFit = fit[0];
        for (double i: fit) {
            if (i > maxFit) {
                maxFit = i;
            }
        }
    }
}

bool compareFitness(int popIdx1, int popIdx2) {
    return fit[popIdx1] < fit[popIdx2];
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
        } while (tournIdxs.size() < tsize);
    }

    sort(tournIdxs.begin(), tournIdxs.end(), compareFitness);
    if (decreasing) {
        reverse(tournIdxs.begin(), tournIdxs.end());
    }
    return tournIdxs;
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

void matingEvent() {//run a mating event
    int rp;   //loop index, random position, swap variable
    vector<int> tournIdxs;
    tournIdxs = tournSelect(tsize, ctrlFitnessFctn == 1); // TODO: Will not work for profile matching.

//    printIdxsOfVector(cout, fit, tournIdxs, "Fit Vals in T: ", ", ", true);

    SDAPop[tournIdxs[0]].copy(SDAPop[tournIdxs[tsize - 2]]);
    SDAPop[tournIdxs[1]].copy(SDAPop[tournIdxs[tsize - 1]]);
    SDAPop[tournIdxs[0]].twoPtCrossover(SDAPop[tournIdxs[1]]);
    rp = (int) lrand48() % maxMuts + 1;
    SDAPop[tournIdxs[0]].mutate(rp);
    rp = (int) lrand48() % maxMuts + 1;
    SDAPop[tournIdxs[1]].mutate(rp);

    // reset dead SDAs
    int size = stringLen - 2;
    SDAPop[tournIdxs[1]].getBitsVec(size, vals);
    dead[tournIdxs[0]] = necroticFilter();
    SDAPop[tournIdxs[1]].getBitsVec(size, vals);
    dead[tournIdxs[1]] = necroticFilter();

    if (!dead[tournIdxs[0]]) {
        fit[tournIdxs[0]] = fitnessPrep(tournIdxs[0], SDAPop[tournIdxs[0]], false);
    } else {
        fit[tournIdxs[0]] = minFit;
    }

    if (!dead[tournIdxs[1]]) {
        fit[tournIdxs[1]] = fitnessPrep(tournIdxs[1], SDAPop[tournIdxs[1]], false);
    } else {
        fit[tournIdxs[1]] = minFit;
    }
}

vector<double> calcStats(const vector<int> &goodIdxs, bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);
    double worstVal = (biggerBetter ? MAXFLOAT : 0.0);

    for (int idx: goodIdxs) {
        sum += fit[idx];
        if ((biggerBetter && fit[idx] > bestVal) || (!biggerBetter && fit[idx] < bestVal)) {
            bestVal = fit[idx];
        }
        if ((biggerBetter && fit[idx] < worstVal) || (!biggerBetter && fit[idx] > worstVal)) {
            worstVal = fit[idx];
        }
    }

    double mean = sum / (double) goodIdxs.size();
    double stdDevSum = 0.0;
    for (int idx: goodIdxs) {
        stdDevSum += pow(fit[idx] - mean, 2);
    }
    double stdDev = sqrt(stdDevSum / ((double) goodIdxs.size() - 1.0));
    double CI95 = 1.96 * (stdDev / sqrt((double) goodIdxs.size()));

    return {mean, stdDev, CI95, bestVal, worstVal}; // {mean, stdDev, 95CI, best, worst}
}

void report(ostream &aus) {//make a statistical report
    vector<int> goodIdxs;
    int deaths = 0;

    for (int i = 0; i < popsize; i++) {
        if (!dead[i]) {
            goodIdxs.push_back(i);
        } else {
            deaths++;
        }
    }

    vector<double> stats = calcStats(goodIdxs, ctrlProfileMatching != 1);
    double mean = stats[0];
    double stdDev = stats[1];
    double CI95 = stats[2];
    double bestVal = stats[3];
    double worstVal = stats[4];

    // Find the Best and Update the Fitness of Dead Members of the Population
    int bestIdx = 0;
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing
        for (int i = 1; i < popsize; i++) {
            if (fit[i] > fit[bestIdx]) {
                bestIdx = i;
            }
        }
        if (worstVal != minFit) {
            for (int i = 0; i < popsize; i++) {
                if (dead[i]) { // Update to Minimum Fitness
                    fit[i] = worstVal;
                }
            }
            minFit = worstVal;
        }
    } else { // Minimizing
        for (int i = 1; i < popsize; i++) {
            if (fit[i] < fit[bestIdx]) {
                bestIdx = i; //find best fitness
            }
        }
        if (bestVal != maxFit) {
            for (int i = 0; i < popsize; i++) {
                if (dead[i]) { // Update to Maximum Fitness
                    fit[i] = bestVal;
                }
            }
            maxFit = bestVal;
        }
    }

    // Print Report
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing
        aus << left << setw(10) << mean;
        aus << left << setw(12) << CI95;
        aus << left << setw(10) << stdDev;
        aus << left << setw(12) << bestVal;
        aus << left << setw(10) << bestEpiLen[bestIdx];
        aus << left << setw(8) << deaths << endl;
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
        aus << left << setw(10) << mean;
        aus << left << setw(12) << CI95;
        aus << left << setw(10) << stdDev;
        aus << left << setw(8) << worstVal << "\t";
        aus << "Dead: " << deaths << endl;
        if (verbose) {
            cout << left << setw(10) << mean;
            cout << left << setw(12) << CI95;
            cout << left << setw(10) << stdDev;
            cout << left << setw(12) << worstVal << "\t";
            cout << "Dead: " << deaths << endl;
        }
    }
}

void reportBest(ostream &aus) {//report the best graph
    int bestIdx;
    Graph G(verts);

    bestIdx = 0;
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing
        for (int i = 1; i < popsize; i++) {
            if (fit[i] > fit[bestIdx]) {
                bestIdx = i;
            }
        }
    } else { // Minimizing
        for (int i = 1; i < popsize; i++) {
            if (fit[i] < fit[bestIdx]) {
                bestIdx = i;
            }
        }
    }

    aus << "Best Fitness: " << fit[bestIdx] << endl;

    // Print the Profiles
    fitnessPrep(bestIdx, SDAPop[bestIdx], true);
    aus << "Epidemic Profile" << endl;
    if (ctrlVariants) { // TODO: Not Yet Implemented.
        for (int i = 0; i <= bestVarCount; i++) {
            if (i == 0) {
                aus << "NA-->" << "V" << left << setw(2) << i << "\t";
            } else {
                aus << "V" << left << setw(2) << bestVarParents[i] << "->V" << left << setw(2) << i << "\t";
            }
            aus << "[" << left << setw(3) << bestVarLens[i].first << "-";
            aus << left << setw(3) << bestVarLens[i].second << "]:\t";
            for (int j = 0; j < bestVarLens[i].first; j++) {
                aus << "\t";
            }
            for (int j = 0; j <= bestVarLens[i].second - bestVarLens[i].first; j++) {
                aus << bestVarProfs[i].at(j) << "\t";
            }
            aus << endl;
        }
        for (int i = 0; i <= bestVarCount; i++) {
            aus << "V" << i << "\t";
            aus << bestVarDNA[i] << endl;
        }
    } else {
        aus << left << setw(4) << "V0" << " ";
        aus << "[" << left << setw(2) << bestVarLens[0].first << " ";
        aus << bestVarLens[0].second << "]: ";
        for (int j = 0; j < bestVarLens[0].second - bestVarLens[0].first; j++) {
            aus << bestVarProfs[0].at(j) << " ";
        }
        aus << endl;
    }

    // Write the SDA
    aus << "Self-Driving Automata" << endl;
    SDAPop[bestIdx].print(aus);

    SDAPop[bestIdx].getBitsVec(stringLen-2, vals);
    G.fill(vals, diagFill);
    aus << "Graph" << endl;
    G.print(aus);
    aus << endl;
}