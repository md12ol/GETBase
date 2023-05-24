#include "main.h"

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
        sprintf(pathToOut, "%sOutput - Profile%d %dS, %02dP, %dM/", outRoot, pNum, SDANumStates, popsize, maxMuts);
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
        for (int mev = 0; mev < generations; mev++) { // Evolve
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