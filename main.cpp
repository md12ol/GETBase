#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <bitset>
#include <numeric>
#include <algorithm>

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
 * 1. The Output Root (usually "./Output/")
 * 2. 0 -> No Variants (any fitness function); 1 -> Variants (only epidemic length or spread fitness)
 * 3. The Probability of Generating a New Variant
 * 4. Minimum Number of Edits to Variant String for a New Variant
 * 5. Maximum Number of Edits to Variant String for a New Variant
 * 6. 0 -> Profile Matching or Epidemic Length Fitness; 1 -> Epidemic Spread Fitness
 * 7. 0 -> Epidemic Spread or Epidemic Length Fitness; 1 -> Profile Matching Fitness (w/o variants)
 *
 * Thus, for arguments 6 and 7 the possibilities are:
 * Mode 0: 0, 0 -> Epidemic Length Fitness with or without Variants (see argument 2)
 * Mode 1: 0, 1 -> Profile Matching Fitness without Variants
 * Mode 2: 1, 0 -> Epidemic Spread Fitness with ot without Variants (see argument 2)
 * Mode -1: 1, 1 -> Not Possible/Implemented TODO: Not Yet Implemented.
 *
 * @param argc Number of Command Line Variables
 * @param argv Those Variables, Explained Above
 * @return Hopefully a Zero :)
 */
int main(int argc, char *argv[]) {
    fstream stat, best, readme, iGOut;  // For File Output
    char filename[60];                  // I Love C++!
    char *outLoc = new char[45];
    char *outRoot = argv[1];            // Root Directory for Output (usually ./Output/
    char *pLoc;                         // Location of the Profile
    int pNum;                           // Profile Number

    // Get Command Line Arguments
    ctrlVariants = (int) strtol(argv[2], nullptr, 10) == 1;
    variantProb = strtod(argv[3], nullptr);
    minEdits = (int) strtol(argv[4], nullptr, 10);
    maxEdits = (int) strtol(argv[5], nullptr, 10);
    ctrlEpiSpread = (int) strtol(argv[6], nullptr, 10) == 1;
    ctrlProfileMatching = (int) strtol(argv[7], nullptr, 10) == 1;

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
        sprintf(outLoc, "%sOutput - Profile%d %dS, %02dP, %dM/", outRoot, pNum, states, popsize, MNM);
    } else if (ctrlVariants) {
        if (ctrlEpiSpread) {
            sprintf(outLoc, "%sOutput - ES %03dP, %02dI, %.4f%%, %02d-%02dE/",
                    outRoot, popsize, initOneBits, variantProb, minEdits, maxEdits);
        } else {
            sprintf(outLoc, "%sOutput - EL %03dP, %02dI, %.4f%%, %02d-%02dE/",
                    outRoot, popsize, initOneBits, variantProb, minEdits, maxEdits);
        }
    } else {
        if (ctrlEpiSpread) {
            sprintf(outLoc, "%sOutput - ES Base/", outRoot);
        } else {
            sprintf(outLoc, "%sOutput - EL Base/", outRoot);
        }
    }
    mkdir(outLoc, 0777);

    // Okay, Let's Get Started!
    initAlg(pLoc);

    // Determine the Location of the Different Output Files
    sprintf(filename, "%sbest.dat", outLoc);
    best.open(filename, ios::out);
    sprintf(filename, "%sreadme.dat", outLoc);
    readme.open(filename, ios::out);

    // Generate the Readme File
    createReadMe(readme);
    readme.close();

    if (verbose) {
        cmdLineIntro(cout);
    }
    if (!verbose) {
        cout << "Started" << endl;
    }
    for (int run = 1; run <= runs; run++) {
        sprintf(filename, "%srun%02d.dat", outLoc, run);
        stat.open(filename, ios::out);
        if (verbose) cmdLineRun(run, cout);
        initPop();
        report(stat); // Initial Report
        for (int mev = 0; mev < mevs; mev++) { // Evolve
            matingEvent();
            if ((mev + 1) % RE == 0) { // Is it time to report?
                if (verbose) {
                    cout << left << setw(5) << run;
                    cout << left << setw(4) << (mev + 1) / RE;
                }
                report(stat);
            }
        }
        stat.close();
        reportBest(best);
        cout << "Done run " << run << " of " << runs << endl;
    }
    best.close();
    delete[] outLoc;
    return (0);
}

void createReadMe(ostream &aus) {
    aus << "This file contains the info about the files in this folder." << endl;
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
    aus << "Maximum number of mutations: " << MNM << endl;
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
    srand48(RNS);

    for (auto &i: SDAPop) {
        i = new SDA(states);
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

bool necroticFilter(const int *upTri) {
    int len = stringLen - 2;
    int count[2] = {0, 0};
    int bounds[2] = {1 * verts, 6 * verts};
    if (ctrlFitnessFctn == 2) {
        bounds[0] = 1 * verts;
        bounds[1] = 2.5 * verts;
    }

    for (int i = 0; i < len; i++) {
        count[upTri[i]]++;
    }
    if (count[1] < bounds[0] || count[1] > bounds[1]) {
        return true; // DEAD
    } else {
        return false;
    }
}

double fitness(int *upTri, int idx, bool final) {//compute the fitness
    graph G(verts);      //scratch graph
    G.UTAM(upTri);
    int max, len, ttl;   //maximum, length, and total removed
    int cnt;             //counter for tries
    double prof[verts];  //profile variable
    double trials[NSE];  //stores squared error for each trial
    int en;              //epidemic number
    double delta;        //difference between profile and trial
    double accu = 0.0;         //accumulator
    int best_epi = 0;
    vector<int> varBaseProf;
    numVars = 0;
    int var_parents[maxNumVars];
    pair<int, int> var_lens[maxNumVars];
    vector<int> var_profs[maxNumVars];
    bitset<DNALen> variants[verts];
    vector<double> v_lengths;
    v_lengths.clear();
    v_lengths.reserve(NSE);
    vector<double> v_spreads;
    v_spreads.clear();
    v_spreads.reserve(NSE);

    if (ctrlFitnessFctn == 0) { //  Epidemic length
        for (en = 0; en < (final ? 10 * NSE : NSE); en++) {
            cnt = 0;
            do {
                if (!ctrlVariants) {
                    G.SIR(0, max, len, ttl, alpha, varBaseProf);
                } else {
                    G.varSIR(0, numVars, var_profs, variants,
                             var_parents, var_lens, 0.5,
                             minEdits, maxEdits, variantProb);
                    v_lengths.clear();
                    for (int i = 0; i <= numVars; i++) {
                        v_lengths.push_back((double) (var_lens[i].second));
                    }
                    sort(v_lengths.begin(), v_lengths.end());
                    len = (int) v_lengths.at(numVars);
                }
                cnt++;
            } while (len < mepl && cnt < rse);
            if (final) {
                if (!ctrlVariants) {
                    if (len > best_epi) {
                        best_epi = len;
                        bestVarCount = 0;
                        bestVarLens[0].first = 0;
                        bestVarLens[0].second = len;
                        bestVarProfs[0] = varBaseProf;
                        bestVarParents[0] = -1;
                        bestVarDNA[0] = bitset<DNALen>();
                    }
                } else {
                    if (len > best_epi) {
                        best_epi = len;
                        for (int i = 0; i <= numVars; i++) {
                            bestVarCount = numVars;
                            bestVarLens[i] = var_lens[i];
                            bestVarProfs[i] = var_profs[i];
                            bestVarParents[i] = var_parents[i];
                            bestVarDNA[i] = variants[i];
                        }
                    }
                }
                accu = -1;
            } else {
                trials[en] = len;
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
            bestEpi[idx] = longest;
        }
    } else if (ctrlFitnessFctn == 1) { // Profile Matching
        for (en = 0; en < NSE; en++) {//loop over epidemics
            cnt = 0;
            do {
                G.SIRProfile(0, max, len, ttl, alpha, prof);
                cnt++;
            } while (len < mepl && cnt < rse);
            trials[en] = 0;  //zero the current squared error
            if (len < PL + 1) {
                len = PL + 1;  //find length of epi/prof (longer)
            }
            for (int i = 0; i < len; i++) {//loop over time periods
                delta = prof[i] - PD[i];
                trials[en] += delta * delta;
            }
            trials[en] = sqrt(trials[en] / len); //convert to RMS error
        }
    } else { //  Epidemic Spread
        for (en = 0; en < (final ? 10 * NSE : NSE); en++) {
            cnt = 0;
            do {
                if (!ctrlVariants) {
                    G.SIR(0, max, len, ttl, alpha, varBaseProf);
                } else {
                    G.varSIR(0, numVars, var_profs, variants,
                             var_parents, var_lens, 0.5,
                             minEdits, maxEdits, variantProb);
                    v_lengths.clear();
                    ttl = 0;
                    for (int i = 0; i <= numVars; i++) {
                        v_lengths.push_back((double) (var_lens[i].second));
                        ttl += accumulate(var_profs[i].begin(), var_profs[i].end(), 0);
                    }
                    sort(v_lengths.begin(), v_lengths.end());
                    len = (int) v_lengths.at(numVars);
                }
                cnt++;
            } while (len < mepl && cnt < rse);
            if (final) {
                if (!ctrlVariants) {
                    if (ttl > best_epi) {
                        best_epi = ttl;
                        bestVarCount = 0;
                        bestVarLens[0].first = 0;
                        bestVarLens[0].second = len;
                        bestVarProfs[0] = varBaseProf;
                        bestVarParents[0] = -1;
                        bestVarDNA[0] = bitset<DNALen>();
                    }
                } else {
                    if (ttl > best_epi) {
                        best_epi = ttl;
                        for (int i = 0; i < maxNumVars; i++) {
                            bestVarCount = numVars;
                            bestVarLens[i] = var_lens[i];
                            bestVarProfs[i] = var_profs[i];
                            bestVarParents[i] = var_parents[i];
                            bestVarDNA[i] = variants[i];
                        }
                    }
                }
                accu = -1;
            } else {
                trials[en] = ttl;
            }
        }
        if (!final) {
            int furthest = 0;
            for (double trial: trials) {//loop over trials
                if (trial > furthest) {
                    furthest = (int) trial;
                }
                accu += trial;
            }
            accu = accu / NSE;
            bestEpi[idx] = furthest;
        }
    }
    return accu;  //return the fitness value
}

void createUpTri(SDA &A, int *upTri) {//unpack the queue
    int h, t;  //head and tail of queue

    for (t = 0; t < stringLen; t++) {
        SDAOutput[t] = 0;        //clear the queue
    }
    A.reset(SDAOutput, h, t);     //reset the self driving automata
    while (t < stringLen - 2) {
        A.next(SDAOutput, h, t, stringLen);  //run the automata
    }

    for (int i = 0; i < stringLen - 2; i++) {
        upTri[i] = SDAOutput[i];
    }
}

double fitnessPrep(int idx, SDA &A, bool final) {
    createUpTri(A, netPop[idx]);
    if (necroticFilter(netPop[idx])) {
        dead[idx] = true;
        if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) {
            return minFit;
        } else {
            return maxFit;
        }
    }
    double fi = fitness(netPop[idx], idx, final);
    return fi;
}

void printUpTri(int *upTri) {
    for (int i = 0; i < stringLen - 2; i++) {
        cout << upTri[i] << " ";
    }
    cout << endl;
}

void initPop() {
    // Generate the Initial Population
    for (int i = 0; i < popsize; i++) {
        dead[i] = false;
        do {
            SDAPop[i]->randomize();
            createUpTri(*SDAPop[i], netPop[i]); // Stored in SDAOutput
        } while (necroticFilter(netPop[i])); // Until the SDA Passes the Necrotic Filter
        fit[i] = fitnessPrep(i, *SDAPop[i], false);
        dx[i] = i;
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

void matingEvent() {//run a mating event
    int rp;   //loop index, random position, swap variable

    //perform tournament selection
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // duration, ctrlEpiSpread
        tselect(fit, dx, tsize, popsize); // lowest first
    } else { // profile matching
        Tselect(fit, dx, tsize, popsize); // highest first
    }

    SDAPop[dx[0]]->copy(*SDAPop[dx[tsize - 2]]);
    SDAPop[dx[1]]->copy(*SDAPop[dx[tsize - 1]]);
    SDAPop[dx[0]]->tpc(*SDAPop[dx[1]]);
    rp = (int) lrand48() % MNM + 1;
    SDAPop[dx[0]]->mutate(rp);
    rp = (int) lrand48() % MNM + 1;
    SDAPop[dx[1]]->mutate(rp);
    // reset dead SDAs
    dead[dx[0]] = false;
    dead[dx[1]] = false;
    fit[dx[0]] = fitnessPrep(dx[0], *SDAPop[dx[0]], false);
    fit[dx[1]] = fitnessPrep(dx[1], *SDAPop[dx[1]], false);
}

void culling() {
    for (int i = 0; i < popsize; i++) {
        if (dead[i]) {
            do {
                SDAPop[i]->randomize();
                createUpTri(*SDAPop[i], netPop[i]);
            } while (necroticFilter(netPop[i]));
            fit[i] = fitnessPrep(i, *SDAPop[i], false);
            minFit = fit[0];
            dead[i] = false;
        }
    }
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) {
        for (double i: fit) {
            if (i < minFit) {
                minFit = i;
            }
        }
    } else {
        for (double i: fit) {
            if (i > maxFit) {
                maxFit = i;
            }
        }
    }
}

void report(ostream &aus) {//make a statistical report
    dset D;
    int deaths = 0;
    for (bool i: dead) {
        if (i) {
            deaths++;
        }
    }

    // Gather the Fitness from the Alive Members of the Population
    double good_fit[popsize - deaths];
    for (double &i: good_fit) {
        i = 0.0;
    }
    int cnt = 0;
    for (int i = 0; i < popsize; i++) {
        if (!dead[i]) {
            good_fit[cnt++] = fit[i];
        }
    }

    // Do the Stats!
    D.add(good_fit, popsize - deaths);

    // Find the Best and Update the Fitness of Dead Members of the Population
    int bestIdx = 0;
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing
        for (int i = 1; i < popsize; i++) {
            if (fit[i] > fit[bestIdx]) {
                bestIdx = i;
            }
        }
        if (D.Rmin() != minFit) {
            for (int i = 0; i < popsize; i++) {
                if (dead[i]) { // Update to Minimum Fitness
                    fit[i] = D.Rmin();
                }
            }
            minFit = D.Rmin();
        }
    } else { // Minimizing
        for (int i = 1; i < popsize; i++) {
            if (fit[i] < fit[bestIdx]) {
                bestIdx = i; //find best fitness
            }
        }
        if (D.Rmax() != maxFit) {
            for (int i = 0; i < popsize; i++) {
                if (dead[i]) { // Update to Maximum Fitness
                    fit[i] = D.Rmax();
                }
            }
            maxFit = D.Rmax();
        }
    }

    // Print Report
    if (ctrlFitnessFctn == 0 || ctrlFitnessFctn == 2) { // Maximizing
        aus << left << setw(10) << D.Rmu();
        aus << left << setw(12) << D.RCI95();
        aus << left << setw(10) << D.Rsg();
        aus << left << setw(12) << D.Rmax();
        aus << left << setw(10) << bestEpi[bestIdx];
        aus << left<< setw(8) << deaths << endl;
        if (verbose) {
            cout << left << setw(10) << D.Rmu();
            cout << left << setw(12) << D.RCI95();
            cout << left << setw(10) << D.Rsg();
            cout << left << setw(12) << D.Rmax();
            cout << left << setw(10) << bestEpi[bestIdx];
            cout << left << setw(8) << deaths << endl;
        }
    } else { // Minimizing
        // TODO: Not Yet Implemented.
        aus << left << setw(10) << D.Rmu();
        aus << left << setw(12) << D.RCI95();
        aus << left << setw(10) << D.Rsg();
        aus << left << setw(8) << D.Rmin() << "\t";
        aus << "Dead: " << deaths << endl;
        if (verbose) {
            cout << left << setw(10) << D.Rmu();
            cout << left << setw(12) << D.RCI95();
            cout << left << setw(10) << D.Rsg();
            cout << left << setw(12) << D.Rmin() << "\t";
            cout << "Dead: " << deaths << endl;
        }
    }
}

void reportBest(ostream &aus) {//report the best graph
    int bestIdx;
    graph G(verts);
    double En;

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
    fitnessPrep(bestIdx, *SDAPop[bestIdx], true);
    aus << "Epidemic Profile" << endl;
    if (ctrlVariants) {
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
        for (int j = 0; j <= bestVarLens[0].second - bestVarLens[0].first; j++) {
            aus << bestVarProfs[0].at(j) << " ";
        }
        aus << endl;
    }

    // Write the SDA
    aus << "Self-Driving Automata" << endl;
    SDAPop[bestIdx]->print(aus);

    G.empty(verts);
    G.UTAM(netPop[bestIdx]);
    aus << "Graph" << endl;
    G.write(aus);
    aus << endl;
}

