#include "SDA.h"

#define VERBOSEBIT false

SDA::SDA() {
    initInput = -1;
    numStates = -1;
    initState = -1;
    curState = -1;
    if (VERBOSEBIT) cout << "SDA Made" << endl;
}

SDA::SDA(int states, int numChars, int len) {
    initInput = -1;
    numStates = states;
    initState = -1;
    curState = -1;
    this->numChars = numChars;
    buf.reserve(len);
    transitions.reserve(states);

    for (vector<int> v: transitions) {
        v.reserve(numChars);
    }
    responses.reserve(states);
    for (vector<vector<int> > v: responses) {
        v.reserve(numChars);
    }

    if (VERBOSEBIT) cout << "SDA Made w " << states << " SDANumStates" << endl;
}

SDA::SDA(SDA &other) {
    copy(other);
}

SDA::~SDA() {//destructor
    destroy();    //call the deallocation routines.
}

int SDA::create() {
    initInput = (int) lrand48() % numChars;
    initState = 0;
    curState = -1;

    transitions.clear();
    responses.clear();
    buf.clear();

    for (int i = 0; i < buf.capacity(); ++i) {
        buf.push_back(-1);
    }

    vector<int> oneState;
    for (int s = 0; s < numStates; ++s) {
        oneState.clear();
        for (int i = 0; i < numChars; ++i) {
            oneState.push_back((int) lrand48() % numStates);
        }
        transitions.push_back(oneState);
    }

    vector<int> oneResponse;
    vector<vector<int>> oneStateResps;
    int respSize;
    for (int s = 0; s < numStates; ++s) {
        oneStateResps.clear();
        for (int t = 0; t < numChars; ++t) {
            oneResponse.clear();
            respSize = (int) lrand48() % 2 + 1;
            for (int i = 0; i < respSize; ++i) {
                oneResponse.push_back((int) lrand48() % numChars);
                //TODO: Check!
            }
            oneStateResps.push_back(oneResponse);
        }
        responses.push_back(oneStateResps);
    }
    return 0;
}

int SDA::randomize() {
    initInput = (int) lrand48() % numChars;
    initState = 0;
    curState = -1;

    int respSize;
    vector<int> oneResponse;
    for (int s = 0; s < numStates; ++s) {
        for (int i = 0; i < numChars; ++i) {
            transitions[s][i] = (int) lrand48() % numStates;
            oneResponse.clear();
            respSize = (int) lrand48() % 2 + 1;
            for (int j = 0; j < respSize; ++j) {
                oneResponse.push_back((int) lrand48() % numChars);
            }
            responses[s][i] = oneResponse;
        }
    }
    return 0;
}

int SDA::copy(SDA &other) {// Creates a copy of an SDA
    initInput = other.initInput;
    numStates = other.numStates;
    initState = other.initState;
    curState = other.curState;
    this->numChars = other.numChars;

    transitions = other.transitions;
    responses = other.responses;
    if (VERBOSEBIT) cout << "SDA Copied" << endl;
    return 0;
}

// int SDA::print() {// THIS IS POINTLESS AS FAR AS I CAN TELL
//     print(cout);
//     return 0;
// }

int SDA::print(ostream &aus) {
    aus << initState << " <- " << initInput << endl;
    for (int s = 0; s < numStates; ++s) {
        if (transitions[s].size() > numChars) {
            aus << "ERROR!  More transitions than MAXVAL!" << endl;
        }
        if (responses[s].size() > numChars) {
            aus << "ERROR!  More responses than MAXVAL!" << endl;
        }
        for (int t = 0; t < numChars; ++t) {
            aus << s << " + " << t << " -> " << transitions.at(s).at(t) << " [";
            if (responses[s][t].size() > 2) {
                aus << "ERROR!  Response length more than 2!" << endl;
            }
            for (int v: responses.at(s).at(t)) {
                aus << " " << v;
            }
            aus << " ]" << endl;
        }
    }
    if (transitions.size() > numStates) {
        aus << "ERROR!  More transitions than the number of SDANumStates!" << endl;
    }
    if (responses.size() > numStates) {
        aus << "ERROR!  More responses than the number of SDANumStates!" << endl;
    }
    if (VERBOSEBIT) cout << "SDA Printed" << endl;
    return 0;
}

int SDA::destroy() {
    return 0;
}

int SDA::twoPtCrossover(SDA &other) {
    int cp1, cp2;
    int swapInt;
    vector<int> swapVec;

    if (numStates != other.numStates) {
        return 1;
    }

    do {
        cp1 = (int) lrand48() % numStates;
        cp2 = (int) lrand48() % numStates;
        if (cp1 > cp2) {
            swapInt = cp1;
            cp1 = cp2;
            cp2 = swapInt;
        }
    } while (cp1 == cp2);

    if (cp1 == 0) {
        swapInt = initInput;
        initInput = other.initInput;
        other.initInput = swapInt;
    }

    for (int s = cp1; s < cp2; s++) {
        swapVec = transitions.at(s);
        transitions.at(s) = other.transitions.at(s);
        other.transitions.at(s) = swapVec;
        swapVec = responses.at(s).at(0);
        responses.at(s).at(0) = other.responses.at(s).at(0);
        other.responses.at(s).at(0) = swapVec;
        swapVec = responses.at(s).at(1);
        responses.at(s).at(1) = other.responses.at(s).at(1);
        other.responses.at(s).at(1) = swapVec;
    }
    return 0;
}

int SDA::oneStateCrossover(SDA &other) {
    int state, swapInt;
    vector<int> swapVec;

    if (numStates != other.numStates) {
        return 1;
    }

    state = (int) lrand48() % numStates;
    if (state == 0) {
        swapInt = initInput;
        initInput = other.initInput;
        other.initInput = swapInt;
    }

    swapVec = transitions.at(state);
    transitions.at(state) = other.transitions.at(state);
    other.transitions.at(state) = swapVec;

    return 0;
}

int SDA::mutate(int numMuts) {
    int mutPt;
    vector<int> oneResponse;
    int respSize;

    for (int m = 0; m < numMuts; ++m) {
        mutPt = (int) lrand48() % (2 * numStates + 1);

        if (mutPt == 0) {
            initInput = (int) lrand48() % numChars;
            return 0;
        }
        mutPt = (mutPt - 1) / 2;
        int transNum = (int) lrand48() % numChars;
        if ((int) lrand48() % 2 == 0) { // Mutate transition
            transitions.at(mutPt).at(transNum) = (int) lrand48() % numStates;
        } else { // Mutate response
            oneResponse.clear();
            respSize = (int) lrand48() % 2 + 1;
            for (int i = 0; i < respSize; ++i) {
                oneResponse.push_back((int) lrand48() % numChars);
            }
            responses.at(mutPt).at(transNum) = oneResponse;
        }
    }
    return 0;
}

int SDA::getBitsVec(int len, vector<int> &rtn) {
    int nextBit;
    int head = 0;
    int tail = 0;
    curState = initState;
    rtn[head] = initInput;
    head += 1;

    while (head < len) {
        nextBit = rtn[tail];
        tail += 1;
        for (int i: responses[curState][nextBit]) {
            if (head < len) {
                rtn[head] = i;
                head += 1;
            }
        }
        curState = transitions[curState][nextBit];
    }
    return 0;
}

int SDA::printBitsVec(int len, ostream &aus) {
    vector<int> vec(len);
    getBitsVec(len, vec);
    for (int i: vec) {
        aus << i << " ";
    }
    aus << endl;
    return 0;
}