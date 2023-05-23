#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <bitset>
#include <numeric>

#include "setu.h"
#include "stat.h"
#include "SDA.h"

// Generic Controls
#define RNS 91207819                                // Random Number Seed
#define verbose true                                // Print Information ot User?
#define RIs 100                                     // Reporting Interval
#define RE ((long)mevs/RIs)                         // Report Every
#define verts 256

// Genetic Algorithm Controls
#define popsize 50
#define runs 30
#define mevs 10000                                  // Number of Mating Events
#define tsize 7                                     // Tournament Size
#define MNM 2                                       // Maximum Number of Mutations

// Epidemic Controls
#define PL 16                                       // Profile Length
#define NSE 30                                      // Number of Sample Epidemics
#define alpha 0.5                                   // Probability of Infection
#define mepl 3                                      // Minimum Epidemic Length
#define rse 5                                       // Retries for Short Epidemics

// Self-Driving Automata (SDA) Controls
#define states 12
#define stringLen verts*(verts-1)/2 + 2

// Genetic Algorithm Variables
SDA *SDAPop[popsize];                               // Stores the Population of SDAs
int netPop[popsize][stringLen - 2];                 // Stores the Upper Triangular Adjacency Matrices (UTAMs)
double fit[popsize];                                // Fitness of each SDA
bool dead[popsize];                                 // Which SDAs Failed the Necrotic Filter (true)
int ctrlFitnessFctn;                                // Program Control for Which Fitness Function to Use
double maxFit = verts;
double minFit = 0.0;

// SDA Variables
int SDAOutput[stringLen];                           // Stores the Output from an SDA

// Epidemic Variables
int bestEpi[popsize];                               // Best Epidemic for each SDA
double PD[PL + 1];                                  // Profile dictionary
bool ctrlProfileMatching;                           // Program Control for Doing Profile Matching
bool ctrlVariants;                                  // Program Control for Having Epidemic Variants
bool ctrlEpiSpread;                                 // Program Control for Using Epidemic Spread Fitness

// Epidemic Variants Variables
int numVars;                                        // Number of Variants
int bestVarCount;                                   // The Number of Variants in the Best Epidemic
int bestVarParents[maxNumVars];                     // The Parents of Each Variant in the Best Epidemic
pair<int, int> bestVarLens[maxNumVars];             // The Variant Lengths (start, end) for Variants in the Best Epidemic
vector<int> bestVarProfs[maxNumVars];               // The Variant Profiles for the Variants in the Best Epidemic
bitset<DNALen> bestVarDNA[maxNumVars];              // The Variant DNA Strings for the Variants in the Best Epidemic
double variantProb;                                 // Probability of Generating a new Variant
int minEdits;                                       // Minimum Number of Edits to New Variant String
int maxEdits;                                       // Maximum Number of Edits to New Variant String

// Other Variables
int dx[popsize];                                    // Used to Sort the Population

// Method Declarations
void cmdLineIntro(ostream &aus);                    // Print Intro to User
void createReadMe(ostream &aus);                    // Make readme.dat
void cmdLineRun(int run, ostream &aus);             // Print Column Headers to User
void initAlg(const char *pLoc);                     // Initialize the Algorithm
void initPop();                                     // Initialize the Population
void matingEvent();
void createUpTri(SDA &A, int *upTri);               // Pull Output from the SDA and Generate the UTAM
double fitnessPrep(int idx, SDA &A, bool final);    // Prepare to Calculate Fitness
double fitness(int *upTri, int idx, bool final);    // Calculate the Fitness
void report(ostream &aus);
void reportBest(ostream &aus);

/**
 * Run the program: initialize a population, evolve the population for the specified number of mating events,
 * and output results.  This will be repeated runs times.  The command line arguments are explained below:
 * 1. 0 -> No Variants (any fitness function); 1 -> Variants (only epidemic length or spread fitness)
 * 2. The Probability of Generating a New Variant
 * 3. Minimum Number of Edits to Variant String for a New Variant
 * 4. Maximum Number of Edits to Variant String for a New Variant
 * 5. 0 -> Profile Matching or Epidemic Length Fitness; 1 -> Epidemic Spread Fitness
 * 6. 0 -> Epidemic Spread or Epidemic Length Fitness; 1 -> Profile Matching Fitness (w/o variants)
 *
 * Thus, for arguments 5 and 6 the possibilities are:
 * Mode 0: 0, 0 -> Epidemic Length Fitness with or without Variants (see argument 1)
 * Mode 1: 0, 1 -> Profile Matching Fitness without Variants
 * Mode 2: 1, 0 -> Epidemic Spread Fitness with ot without Variants (see argument 1)
 * Mode -1: 1, 1 -> Not Possible/Implemented TODO: Not Yet Implemented.
 *
 * @param argc Number of Command Line Arguments
 * @param argv Those Arguments, Explained Above
 * @return Hopefully a Zero :)
 */
int main(int argc, char *argv[]) {
    char filename[200];
    fstream runStats, expStats, readMe; // For File Output
    char *pLoc;                         // Location of the Profile
    int pNum;                           // Profile Number

    getArgs(argv);

    // Determine Fitness Function
    if (!ctrlVariants) {
        if (ctrlProfileMatching) { // Mode 1: Profile Matching Fitness w/o Variants TODO: Not Yet Implemented.
            ctrlFitnessFctn = 1;
        } else if (ctrlEpiSpread) {
            ctrlFitnessFctn = 2; // Mode 2: Epidemic Spread Fitness w/o Variants
        } else {
            ctrlFitnessFctn = 0; // Mode 0: Epidemic Length Fitness w/o Variants
        }
    } else {
        if (ctrlEpiSpread) {
            ctrlFitnessFctn = 2; // Mode 2: Epidemic Spread Fitness with Variants
        } else {
            ctrlFitnessFctn = 0; // Mode 0: Epidemic Length Fitness with Variants
        }
    }

    if (ctrlProfileMatching) { // TODO: Not Yet Implemented.
        pLoc = argv[2];
        pNum = stoi(argv[3]);
    }

    // Create Directory for Output
    if (ctrlProfileMatching) { // TODO: Not Yet Implemented.
        sprintf(pathToOut, "%sOutput - Profile%d %dS, %02dP, %dM/", outRoot, pNum, states, popsize, maxMuts);
    } else if (ctrlVariants) {
        if (ctrlEpiSpread) {
            sprintf(pathToOut, "%sOutput - ES %03dP, %02dI, %.4f%%, %02d-%02dE/",
                    outRoot, popsize, initOneBits, variantProb, minEdits, maxEdits);
        } else {
            sprintf(pathToOut, "%sOutput - EL %03dP, %02dI, %.4f%%, %02d-%02dE/",
                    outRoot, popsize, initOneBits, variantProb, minEdits, maxEdits);
        }
    } else {
        if (ctrlEpiSpread) {
            sprintf(pathToOut, "%sOutput - ES Base/", outRoot);
        } else {
            sprintf(pathToOut, "%sOutput - EL Base/", outRoot);
        }
    }
    mkdir(pathToOut, 0777);

    // Okay, Let's Get Started!
    initAlg(pLoc);

    // Determine the Location of the Different Output Files
    sprintf(filename, "%sbest.dat", pathToOut);
    expStats.open(filename, ios::out);
    sprintf(filename, "%sreadme.dat", pathToOut);
    readMe.open(filename, ios::out);

    // Generate the Readme File
    createReadMe(readMe);
    readMe.close();
    cmdLineIntro(cout);

    for (int run = 1; run <= runs; run++) {
        sprintf(filename, "%srun%02d.dat", pathToOut, run);
        runStats.open(filename, ios::out);
        if (verbose) cmdLineRun(run, cout);
        initPop();
        report(runStats); // Initial Report
        for (int mev = 0; mev < mevs; mev++) { // Evolve
            matingEvent();
            if ((mev + 1) % RE == 0) { // Is it time to report?
                if (verbose) {
                    cout << left << setw(5) << run;
                    cout << left << setw(4) << (mev + 1) / RE;
                }
                report(runStats);
            }
        }
        runStats.close();
        reportBest(expStats);
        cout << "Done run " << run << " of " << runs << endl;
    }
    expStats.close();
    return (0);
}