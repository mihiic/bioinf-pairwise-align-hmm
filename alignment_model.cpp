#include "alignment_model.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

AlignmentModel::AlignmentModel(unsigned long sample) {
    sampleSize = sample;
    transitions = {
            {0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0},   // t(0, 0), t(0, 1), t(0, 2), t(0, 3), t(0, 4)
            {0, 0.7,       0.1,       0.1,       0.1}, // t(1, 0), t(1, 1), t(1, 2), t(1, 3), t(1, 4)
            {0, 0.25,      0.65,      0,         0.1}, // t(2, 0), t(2, 1), t(2, 2), t(2, 3), t(2, 4)
            {0, 0.25,      0,         0.65,      0.1}, // t(3, 0), t(3, 1), t(3, 2), t(3, 3), t(3, 4)
            {0, 0,         0,         0,         0}    // t(4, 0), t(4, 1), t(4, 2), t(4, 3), t(4, 4)
    };

    emissions = {
            {0,    0.25,       0.25,       0.25,       0.25 / 1.0},  // {"-, -"}, {"-, A"}, {"-, T"}, {"-, C"}, {"-, G"},
            {0.25, 1.0 / 7.0,  1.0 / 28.0, 1.0 / 28.0, 1.00 / 28.0}, // {"A, -"}, {"A, A"}, {"A, T"}, {"A, C"}, {"A, G"},
            {0.25, 1.0 / 28.0, 1.0 / 7.0,  1.0 / 28.0, 1.0 / 28.0},  // {"T, -"}, {"T, A"}, {"T, T"}, {"T, C"}, {"T, G"},
            {0.25, 1.0 / 28.0, 1.0 / 28.0, 1.0 / 7.0,  1.0 / 28.0},  // {"C, -"}, {"C, A"}, {"C, T"}, {"C, C"}, {"C, G"},
            {0.25, 1.0 / 28.0, 1.0 / 28.0, 1.0 / 28.0, 1.0 / 7.0}    // {"G, -"}, {"G, A"}, {"G, T"}, {"G, C"}, {"G, G"},
    };
}

AlignmentModel::~AlignmentModel() {

}

double AlignmentModel::getEmissionProb(char i, char j) {
    return emissions[charToIndex(i)][charToIndex(j)];
}

int AlignmentModel::charToIndex(char x) {
    switch (x) {
        case '-':
            return 0;
        case 'A':
            return 1;
        case 'T':
            return 2;
        case 'C':
            return 3;
        case 'G':
            return 4;
        default:
            break;
    }
    return 0;
}

void AlignmentModel::run(char *filenameA, char *filenameB) {
    ifstream fileA(filenameA);
    ifstream fileB(filenameB);

    L1 = 0;
    L2 = 1;

    calculateDimensions(fileA, fileB);
    for (int i = 1; i <= totalLines; i++) {
        cout << "current line: " << i << endl;
        if (i == totalLines) {
            reserveMemory(leftover1, leftover2);
        } else {
            reserveMemory(maxL1, maxL2);
        }
    }

    fileA.close();
    fileB.close();
}

void AlignmentModel::reserveMemory(unsigned long l1, unsigned long l2) {
    viterbi.resize(0);
    point.resize(0);

    viterbi.resize(l1 + 2, {});
    point.resize(l1 + 2, {});
    for (long i = 0; i < l1 + 2; i++) {
        viterbi[i].resize(l2 + 2, {});
        point[i].resize(l2 + 2, {});
        for (long j = 0; j < l2 + 2L; j++) {
            viterbi[i][j].resize(5, 0.0);
            point[i][j].resize(5, -1);
        }
    }

    viterbi[0][0][0] = 1;
}

void AlignmentModel::calculateDimensions(ifstream &fileA, ifstream &fileB) {
    char c;
    bool firstLine = false;
    while (fileA.get(c)) {
        if (c != '\n') {
            if (firstLine) {
                L1 += 1;
            }
        } else {
            firstLine = true;
        }
    }

    firstLine = false;
    while (fileB.get(c)) {
        if (c != '\n') {
            if (firstLine) {
                L2 += 1;
            }
        } else {
            firstLine = true;
        }
    }

    if (L1 > L2) {
        maxL1 = sampleSize;
        maxL2 = (unsigned long)ceil(L2 / (L1 / (double)maxL1));
    } else {
        maxL2 = sampleSize;
        maxL1 = (unsigned long)ceil(L1 / (L2 / (double)maxL2));
    }

    leftover1 = L1 % maxL1;
    leftover2 = L2 % maxL2;

    totalLines = (L1 / maxL1) + 1;

    fileA.clear();
    fileB.clear();
    fileA.seekg(0, ios::beg);
    fileB.seekg(0, ios::beg);

    // skip comment lines in both files
    while (fileA.get(c)) {
        if (c == '\n') {
            break;
        }
    }
    fileA.get(c);

    while (fileB.get(c)) {
        if (c == '\n') {
            break;
        }
    }
}
